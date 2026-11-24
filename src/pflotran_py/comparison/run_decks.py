"""Run the closed-batch decks through PFLOTRAN in the container.

Each deck gets its own working directory, because PFLOTRAN writes its output
files next to the input using a fixed prefix, and running two decks in one
directory would have them overwrite each other's results.

Two edits are made to every deck on the way in, neither of which changes the
science:

* The ``DATABASE`` line is repointed at the thermodynamic database as the
  container sees it. The decks carry a path valid on the host, and the
  repository is mounted somewhere else inside the image.
* The chemistry-level ``OUTPUT`` block containing only
  ``WATER_ACTIVITY_COEFFICIENT`` is removed. PFLOTRAN version 6 rejects it as
  an unknown species. The water-activity sandboxes are unaffected: they receive
  water activity internally through ``ACTIVITY_WATER`` and
  ``ACTIVITY_COEFFICIENTS TIMESTEP``, not through that output request. The same
  edit is made by the repository's existing integration tests.
"""

import logging
import os
import re
import shutil
import subprocess

logger = logging.getLogger(__name__)

CONTAINER_IMAGE = "pflotran-py-test"
CONTAINER_WORKDIR = "/work"
CONTAINER_DATABASE_DIR = f"{CONTAINER_WORKDIR}/sandbox"
CONTAINER_PFLOTRAN = "/opt/pflotran-py/pflotran"


def prepare_deck(source_path, destination_path):
    """Copy a deck, repointing its database and removing the v6-hostile block."""
    with open(source_path) as handle:
        text = handle.read()

    # Repoint the database into the container, keeping whichever file the deck
    # asked for. Rewriting the whole path to a fixed filename would silently
    # substitute a different database: the closed-batch decks use a patched copy
    # carrying gas-phase methane, and forcing them back to the stock file makes
    # PFLOTRAN reject the deck for a reason that looks like a chemistry error.
    def _repoint(match):
        requested = match.group(1).strip()
        return f"  DATABASE {CONTAINER_DATABASE_DIR}/{os.path.basename(requested)}"

    text = re.sub(r"^\s*DATABASE\s+(.*)$", _repoint, text, flags=re.MULTILINE)
    text = re.sub(
        r"\n[ \t]*OUTPUT\n[ \t]*WATER_ACTIVITY_COEFFICIENT\n[ \t]*/\n",
        "\n",
        text,
    )

    os.makedirs(os.path.dirname(destination_path), exist_ok=True)
    with open(destination_path, "w") as handle:
        handle.write(text)
    return destination_path


def run_deck(
    deck_path, run_root, repo_root, timeout_seconds=3600, container_runtime="docker"
):
    """Run one deck and return its working directory.

    Parameters
    ----------
    deck_path : str
        Path to the deck on the host.
    run_root : str
        Directory under which each deck gets its own subdirectory.
    repo_root : str
        Repository root, mounted into the container so the deck can reach the
        thermodynamic database.

    Returns
    -------
    dict
        ``{"name", "workdir", "returncode", "log"}``. A non-zero return code is
        reported rather than raised, so that one failing condition does not stop
        the rest of the set from running.
    """
    name = os.path.splitext(os.path.basename(deck_path))[0]
    workdir = os.path.join(run_root, name)
    os.makedirs(workdir, exist_ok=True)

    prepare_deck(deck_path, os.path.join(workdir, "sim.in"))

    relative_workdir = os.path.relpath(workdir, repo_root)
    command = [
        container_runtime,
        "run",
        "--rm",
        "-v",
        f"{repo_root}:{CONTAINER_WORKDIR}",
        "-w",
        f"{CONTAINER_WORKDIR}/{relative_workdir}",
        "--entrypoint",
        CONTAINER_PFLOTRAN,
        CONTAINER_IMAGE,
        "-input_prefix",
        "sim",
    ]

    logger.info("Running %s", name)
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout_seconds,
    )

    log_path = os.path.join(workdir, "pflotran.log")
    with open(log_path, "w") as handle:
        handle.write(completed.stdout)

    if completed.returncode != 0:
        logger.warning(
            "%s exited with code %s; see %s", name, completed.returncode, log_path
        )

    return {
        "name": name,
        "workdir": workdir,
        "returncode": completed.returncode,
        "log": completed.stdout,
    }


def run_all(deck_dir, run_root, repo_root, **kwargs):
    """Run every deck in a directory, one after another.

    Sequential rather than parallel. Each run is chemistry-heavy and memory
    hungry, and the repository's existing batch runner makes the same choice
    for the same reason.
    """
    decks = sorted(
        os.path.join(deck_dir, name)
        for name in os.listdir(deck_dir)
        if name.endswith(".in")
    )

    results = []
    for deck in decks:
        results.append(run_deck(deck, run_root, repo_root, **kwargs))
    return results


def main():
    import argparse
    import time

    parser = argparse.ArgumentParser(
        description="Run the closed-batch decks through PFLOTRAN in the container."
    )
    parser.add_argument("--deck-dir", default="decks")
    parser.add_argument("--run-root", default="runs")
    parser.add_argument("--repo-root", default=os.getcwd())
    parser.add_argument("--only", nargs="+", help="Run only these deck basenames.")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument(
        "--clean", action="store_true", help="Delete the run root before starting."
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.clean and os.path.isdir(args.run_root):
        shutil.rmtree(args.run_root)

    decks = sorted(
        os.path.join(args.deck_dir, name)
        for name in os.listdir(args.deck_dir)
        if name.endswith(".in")
    )
    if args.only:
        decks = [
            d for d in decks if os.path.splitext(os.path.basename(d))[0] in args.only
        ]

    results = []
    for deck in decks:
        started = time.monotonic()
        result = run_deck(
            deck, args.run_root, args.repo_root, timeout_seconds=args.timeout
        )
        result["seconds"] = time.monotonic() - started
        results.append(result)
        status = (
            "ok" if result["returncode"] == 0 else f"FAILED ({result['returncode']})"
        )
        print(f"  {result['name']:32s} {status:20s} {result['seconds']:6.1f} s")

    failed = [r for r in results if r["returncode"] != 0]
    print()
    print(f"{len(results) - len(failed)} of {len(results)} runs completed.")
    if failed:
        print("Failed:")
        for result in failed:
            print(f"  {result['name']}  -> {result['workdir']}/pflotran.log")


if __name__ == "__main__":
    main()
