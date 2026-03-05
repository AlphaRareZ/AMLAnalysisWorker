FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
# Install system dependencies for headless 3D rendering
RUN apt-get update && apt-get install -y \
    xvfb \
    libgl1-mesa-glx \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*
COPY . .
CMD ["xvfb-run","python", "main.py"]