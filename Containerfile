# Containerfile for pflotran-py integration testing
#
# Provides PFLOTRAN v6.0 with custom AWINHIBIT sandboxes + Python test deps.
# Works with both Podman and Docker.
#
# Build:
#   docker build -t pflotran-py-test -f Containerfile .
#
# Run integration tests:
#   docker run --rm -v $(pwd):/work -w /work pflotran-py-test \
#       pytest tests/ -v --tb=short -m integration

FROM pshuai/jupyter-pflotran-multiplatform:base_v6

USER root

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Build PFLOTRAN with pflotran-py custom AWINHIBIT reaction sandboxes.
COPY sandbox/ /tmp/sandbox/
COPY scripts/patch_pflotran_sandboxes.py /tmp/scripts/patch_pflotran_sandboxes.py
RUN python3 /tmp/scripts/patch_pflotran_sandboxes.py \
        --pflotran-src /scratch/pflotran/src/pflotran \
        --sandbox-dir /tmp/sandbox && \
    cd /scratch/pflotran/src/pflotran && \
    export PETSC_DIR=/scratch/petsc PETSC_ARCH=petsc-arch && \
    make clean >/dev/null 2>&1 || true && \
    make -j$(nproc) pflotran && \
    mkdir -p /opt/pflotran-py && \
    cp pflotran /opt/pflotran-py/pflotran

WORKDIR /work

ENTRYPOINT ["python3", "-m"]
CMD ["pytest", "tests/", "-v", "--tb=short"]
