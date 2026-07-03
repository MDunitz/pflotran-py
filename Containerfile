# Containerfile for pflotran-py integration testing
#
# Provides PFLOTRAN v6.0 binary + Python test dependencies.
# Works with both Podman and Docker.
#
# Build:
#   podman build -t pflotran-py-test -f Containerfile .
#
# Run tests:
#   podman run --rm -v $(pwd):/work:Z -w /work pflotran-py-test \
#       pytest tests/ -v --tb=short
#
# Run integration tests only:
#   podman run --rm -v $(pwd):/work:Z -w /work pflotran-py-test \
#       pytest tests/test_integration.py -v --tb=short -m integration

FROM pshuai/jupyter-pflotran-multiplatform:base_v6

USER root

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /work

ENTRYPOINT ["python3", "-m"]
CMD ["pytest", "tests/", "-v", "--tb=short"]
