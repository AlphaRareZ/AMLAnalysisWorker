FROM python:3.11-slim

WORKDIR /app

# 1. Install system dependencies FIRST (Better Docker cache utilization)
RUN apt-get update && apt-get install -y \
    xvfb \
    libgl1-mesa-glx \
    libgl1 \
    libxrender1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# 2. Copy and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 3. Copy the rest of your application code
COPY . .

# 4. Run with xvfb-run using the -a flag
CMD ["xvfb-run", "-a", "python", "main.py"]