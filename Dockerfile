FROM python:3.12.9-slim-bookworm@sha256:48a11b7ba705fd53bf15248d1f94d36c39549903c5d59edcfa2f3f84126e7b44
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/matplotlib HOME=/tmp
WORKDIR /workspace
COPY requirements.lock ./
RUN python -m pip install --no-cache-dir --only-binary=:all: -r requirements.lock && python -m pip check
COPY benchmark benchmark
COPY data data
COPY tests tests
CMD ["python", "-m", "unittest", "discover", "-s", "tests", "-v"]
