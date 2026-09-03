FROM python:3.12-slim

LABEL org.opencontainers.image.title="SangerFlow"
LABEL org.opencontainers.image.description="Auditable Sanger chromatogram analysis"

WORKDIR /opt/sangerflow
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN python -m pip install --no-cache-dir .

RUN useradd --create-home --uid 10001 sangerflow
USER sangerflow
WORKDIR /data
ENTRYPOINT ["sangerflow"]

