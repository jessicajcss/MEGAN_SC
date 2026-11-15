# Reproducible production environment for MEGAN-SC.
# Build (from the repository root):  docker build -t megan-sc .
# Run:    docker run --rm -it -v "$PWD":/work megan-sc bash
#         then, inside:  python reproduce.py all
FROM mambaorg/micromamba:1.5-jammy

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /work
ENV PYTHONUNBUFFERED=1 \
    PROJ_NETWORK=ON
# Earth Engine auth is interactive on first use:  earthengine authenticate
CMD ["bash"]
