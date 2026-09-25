# Image reproductible d'orfeval : versions exactes (requirements.txt), utilisateur non root.
#   docker build -t orfeval .
#   docker run --rm -v "$PWD/results/docker:/app/results" orfeval            # démo
#   docker run --rm -p 8501:8501 --entrypoint streamlit orfeval \
#       run app/streamlit_app.py --server.address 0.0.0.0                   # dashboard
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

WORKDIR /app

# Dépendances d'abord (couche mise en cache tant que requirements.txt ne change pas).
COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-deps .

COPY config ./config
COPY data/demo ./data/demo
COPY app ./app

RUN useradd --create-home --uid 1000 orfeval \
    && mkdir -p /app/results \
    && chown orfeval /app/results
USER orfeval

EXPOSE 8501
ENTRYPOINT ["orfeval"]
CMD ["demo", "--output", "/app/results/demo"]
