# CasPINS – Deployment Guide

CasPINS (Cas-Primer-Indel Suite) is a Streamlit-based web application.
It runs locally in any modern browser and does **not** require cloud
infrastructure for standard use.

---

## Option 1: Local installation (recommended for most users)

**Requirements**: Python 3.8 or later, pip, ~500 MB disk space.

```bash
# Clone the repository
git clone https://github.com/InnovationLine/CasPINS.git
cd CasPINS

# Install dependencies
pip install -r requirements.txt

# Launch the GUI
python run.py gui
```

The application opens automatically at **http://localhost:8501**.
On Windows you can also double-click `run_gui.bat`.

### Virtual environment (optional but recommended)

```bash
python -m venv crispr_env

# Activate — Windows
.\crispr_env\Scripts\activate

# Activate — macOS / Linux
source crispr_env/bin/activate

pip install -r requirements.txt
python run.py gui
```

---

## Option 2: Docker (recommended for servers and HPC)

**Requirements**: Docker and Docker Compose installed.

```bash
# Clone the repository
git clone https://github.com/InnovationLine/CasPINS.git
cd CasPINS

# Start with Docker Compose (includes persistent data volume)
docker-compose up -d

# Access at http://localhost:8501
```

To stop:

```bash
docker-compose down
```

Data is persisted in the `./data`, `./output`, `./cache`, and `./logs`
directories, which are mounted into the container via the volumes defined
in `docker-compose.yml`.

### Manual Docker build (alternative)

```bash
docker build -t caspins .
docker run -p 8501:8501 \
  -v "$(pwd)/data":/app/data \
  -v "$(pwd)/logs":/app/logs \
  caspins
```

---

## Option 3: Institutional server (Nginx reverse proxy)

Expose CasPINS over HTTPS at a custom domain by placing Nginx in front
of the Streamlit server.

**Nginx configuration** (`/etc/nginx/sites-available/caspins`):

```nginx
server {
    listen 80;
    server_name caspins.yourdomain.com;

    location / {
        proxy_pass http://localhost:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_read_timeout 86400;
    }
}
```

**Enable HTTPS with Let's Encrypt**:

```bash
sudo apt install certbot python3-certbot-nginx
sudo certbot --nginx -d caspins.yourdomain.com
```

**Run as a systemd service** (`/etc/systemd/system/caspins.service`):

```ini
[Unit]
Description=CasPINS web application
After=network.target

[Service]
User=www-data
WorkingDirectory=/opt/CasPINS
ExecStart=/opt/CasPINS/crispr_env/bin/streamlit run src/gui/app.py \
          --server.port 8501 --server.address 0.0.0.0
Restart=always

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl enable caspins
sudo systemctl start caspins
```

---

## Option 4: HPC clusters (Slurm)

Submit an interactive Streamlit session with port forwarding from your
local machine.

**Slurm submission script**:

```bash
#!/bin/bash
#SBATCH --job-name=caspins
#SBATCH --output=caspins_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2

module load python/3.10
source /path/to/CasPINS/crispr_env/bin/activate

streamlit run src/gui/app.py \
  --server.port 8501 \
  --server.address 0.0.0.0
```

Then SSH tunnel from your laptop:

```bash
ssh -L 8501:compute-node-name:8501 youruser@cluster.example.edu
```

Open `http://localhost:8501` in your local browser.

---

## Environment variables

| Variable | Default | Description |
|---|---|---|
| `STREAMLIT_SERVER_PORT` | `8501` | Web interface port |
| `STREAMLIT_SERVER_ADDRESS` | `0.0.0.0` | Bind address |
| `CRISPR_DATA_DIR` | `./data` | Data directory (optional) |

---

## Troubleshooting

### Port already in use

```bash
# Find the process using port 8501
# Linux / macOS
lsof -i :8501

# Windows
netstat -ano | findstr :8501

# Run on a different port
streamlit run src/gui/app.py --server.port 8502
```

### Docker build fails

```bash
docker system prune -a           # clear cached layers
docker-compose build --no-cache  # full rebuild
```

### High memory usage with large gene sequences

```bash
# Increase Streamlit's message size limit
streamlit run src/gui/app.py --server.maxMessageSize 500
```

### Windows: "python not found"

Ensure Python 3.8+ is on `PATH`, or use the full path:

```batch
C:\Users\YourName\AppData\Local\Programs\Python\Python311\python.exe run.py gui
```

---

## Support

- **Issues and feature requests**: https://github.com/InnovationLine/CasPINS/issues
- **Documentation**: https://github.com/InnovationLine/CasPINS/tree/main/docs
- **Email**: rinkidsgpt@gmail.com
