# Use Python 3.9 slim image as base
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install system dependencies
# CRITICAL: openbabel and autoscan-vina are essential for the docking engine
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    openbabel \
    autoscan-vina \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY . /app/

# Install Python dependencies
# Phase 2: Added meeko and rdkit for enhanced molecular preparation
RUN pip install --no-cache-dir -e .

# Create data directories
RUN mkdir -p /app/data/receptors /app/data/ligands

# AutoScan runtime defaults
ENV AUTOSCAN_DATA_DIR=/app/data \
    AUTOSCAN_CONFIG_DIR=/app/config

# Set the entrypoint to the CLI
ENTRYPOINT ["autoscan"]
CMD ["--help"]


