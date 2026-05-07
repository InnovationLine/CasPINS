# CasPINS - Deployment Guide

This guide covers all deployment options for CasPINS (Cas-Primer-Indel Suite).

## Quick Start Options

### Option 1: Local Installation (Recommended for most users)

```bash
# Clone the repository
git clone https://github.com/InnovationLine/CasPINS.git
cd CasPINS

# Install dependencies
pip install -r requirements.txt

# Run the GUI
python run.py gui
# Or: streamlit run src/gui/app.py
```

### Option 2: Docker (Recommended for servers/HPC)

```bash
# Using Docker Compose (easiest)
docker-compose up -d

# Access at http://localhost:8501

# Or build and run manually
docker build -t caspins.
docker run -p 8501:8501 -v $(pwd)/data:/app/data caspins
```

### Option 3: Streamlit Cloud (Zero installation - web demo)

1. Fork this repository to your GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Sign in with GitHub
4. Click "New app"
5. Select your forked repository
6. Set main file: `streamlit_app.py`
7. Click "Deploy"

Your app will be available at: `https://[your-app-name].streamlit.app`

### Option 4: PyPI Installation

```bash
pip install caspins

# Run GUI
caspins-gui

# Or use CLI tools
caspins-grna GENE_NAME --species human
caspins-primers GENE_NAME
caspins-analyze
```

---

## Deployment Configurations

### Streamlit Cloud Configuration

The repository includes `.streamlit/config.toml` for cloud deployment:

```toml
[server]
maxUploadSize = 200  # MB - for AB1 files

[theme]
primaryColor = "#1E88E5"
backgroundColor = "#FFFFFF"
```

### Docker Configuration

The `Dockerfile` and `docker-compose.yml` are pre-configured:

```yaml
# docker-compose.yml
services:
  crispr-gui:
    build:.
    ports:
      - "8501:8501"
    volumes:
      -./data:/app/data    # Persistent data
      -./output:/app/output
```

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `STREAMLIT_SERVER_PORT` | 8501 | Port for web interface |
| `STREAMLIT_SERVER_ADDRESS` | 0.0.0.0 | Bind address |
| `DATA_DIR` |./data | Data directory path |

---

## Platform-Specific Notes

### Windows

```batch
# Run using batch file
run_gui.bat

# Or via Python
python run.py gui
```

### macOS / Linux

```bash
# Make scripts executable
chmod +x run.py

# Run
./run.py gui
# Or
python run.py gui
```

### HPC Clusters (Slurm)

```bash
#!/bin/bash
#SBATCH --job-name=crispr-gui
#SBATCH --output=crispr_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=4G

module load python/3.10
source venv/bin/activate

streamlit run src/gui/app.py --server.port $((8500 + SLURM_ARRAY_TASK_ID))
```

---

## Binder (For Paper Reviewers)

Click the Binder badge to launch without installation:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/InnovationLine/CasPINS/main)

Configuration is in `.binder/`:
- `requirements.txt` - Python dependencies
- `postBuild` - Setup script
- `start` - Launch command

---

## Production Deployment

### Using Nginx as Reverse Proxy

```nginx
server {
    listen 80;
    server_name crispr.yourdomain.com;

    location / {
        proxy_pass http://localhost:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
    }
}
```

### Using HTTPS with Let's Encrypt

```bash
# Install certbot
sudo apt install certbot python3-certbot-nginx

# Get certificate
sudo certbot --nginx -d crispr.yourdomain.com
```

### Systemd Service

```ini
# /etc/systemd/system/caspins.service
[Unit]
Description=CasPINS GUI
After=network.target

[Service]
User=www-data
WorkingDirectory=/opt/CasPINS
ExecStart=/opt/CasPINS/venv/bin/streamlit run src/gui/app.py
Restart=always

[Install]
WantedBy=multi-user.target
```

---

## Troubleshooting

### Port Already in Use

```bash
# Find process using port 8501
lsof -i :8501
# Or on Windows
netstat -ano | findstr :8501

# Use different port
streamlit run src/gui/app.py --server.port 8502
```

### Docker Build Fails

```bash
# Clear Docker cache
docker system prune -a

# Rebuild without cache
docker-compose build --no-cache
```

### Memory Issues

```bash
# Increase Streamlit memory limit
streamlit run src/gui/app.py --server.maxMessageSize 500
```

---

## Support

- **Issues**: https://github.com/InnovationLine/CasPINS/issues
- **Documentation**: https://caspins.readthedocs.io
- **Email**: rinkidsgpt@gmail.com
