#!/bin/bash
# Startup script for Gene Extractor application

echo "================================================"
echo "   NCBI RefSeq Gene Extractor - Enterprise     "
echo "================================================"
echo ""

# Check if .env exists
if [ ! -f .env ]; then
    echo "⚠️  Warning: .env file not found. Creating from template..."
    cp .env.example .env
    echo "✅ Created .env file. Please edit it with your credentials."
    echo ""
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "🐍 Python version: $PYTHON_VERSION"

# Create necessary directories
mkdir -p data/cache data/logs
echo "✅ Directories created"

# Check if dependencies are installed
if ! python3 -c "import streamlit" 2>/dev/null; then
    echo "📦 Installing dependencies..."
    pip install -r requirements.txt
else
    echo "✅ Dependencies already installed"
fi

echo ""
echo "🚀 Starting application..."
echo "🌐 Access the app at: http://localhost:8501"
echo ""
echo "Press Ctrl+C to stop"
echo ""

# Start Streamlit
streamlit run app.py --server.port 8501 --server.address 0.0.0.0
