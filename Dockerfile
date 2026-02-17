FROM ubuntu:20.04 as builder

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
    && apt-get install -y wget \
    && apt-get autoremove

WORKDIR /tools

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && tar zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz \
    && rm ncbi-blast-2.13.0+-x64-linux.tar.gz

#####

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy NCBI BLAST+ binaries from builder
WORKDIR /tools/ncbi-blast-2.13.0+
COPY --from=builder /tools/ncbi-blast-2.13.0+ .

WORKDIR /tools
RUN ln -s ncbi-blast-2.13.0+ blast

# Create necessary directories
RUN install -d -m 755 /var/www/app && \
    install -d -m 1777 /var/tmp && \
    install -d -m 755 /var/www/conf

WORKDIR /var/www/app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY www/app/*.py ./
COPY www/conf/*.json /var/www/conf/

# Environment variables
ENV DATA_DIR=/data/blast/ \
    BIN_PATH=/tools/blast/bin/ \
    TMP_DIR=/var/tmp/ \
    CONF_DIR=/var/www/conf/

# Health check
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
    CMD curl -fsS http://localhost:8000/ || exit 1

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
