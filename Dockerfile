# CRISPR Analysis Suite - Docker Container
# Multi-platform support for Mac, Windows, Linux

FROM python:3.10-slim

# Set metadata
LABEL maintainer="your.email@example.com"
LABEL description="CRISPR Analysis Suite - gRNA Design, Primer Design, Indel Analysis"
LABEL version="1.0.0"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0

# Create app directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p /app/data /app/cache /app/output

# Set proper permissions
RUN chmod -R 755 /app

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Default command - run GUI
CMD ["streamlit", "run", "crispr_gui.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.headless=true"]
