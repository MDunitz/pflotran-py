"""Make repo packages importable in every environment that runs the tests.

CI installs the package (``pip install -e .``), but the integration
container only installs ``requirements.txt`` and mounts the source, so it
has no editable install. Adding both the repo root (for top-level
``exploratory.*`` modules) and ``src`` (for the ``pflotran_py`` package)
keeps imports resolving in CI, the container, and local runs alike.
"""

import os
import sys

_HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(_HERE, ".."))
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
