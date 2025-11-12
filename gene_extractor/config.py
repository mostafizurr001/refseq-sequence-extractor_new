"""Configuration management for Gene Extractor application."""
import os
from pathlib import Path
from typing import Dict, Any

class Config:
    """Application configuration."""
    
    # Base paths
    BASE_DIR = Path(__file__).parent
    DATA_DIR = BASE_DIR / "data"
    CACHE_DIR = DATA_DIR / "cache"
    LOGS_DIR = DATA_DIR / "logs"
    
    # Database
    DATABASE_PATH = DATA_DIR / "gene_extractor.db"
    
    # NCBI Configuration
    NCBI_EMAIL = os.getenv("NCBI_EMAIL", "your_email@example.com")
    NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")
    NCBI_RATE_LIMIT = int(os.getenv("NCBI_RATE_LIMIT", "3"))  # requests per second
    NCBI_MAX_RETRIES = int(os.getenv("NCBI_MAX_RETRIES", "3"))
    
    # Authentication
    APP_USERNAME = os.getenv("APP_USERNAME", "admin")
    APP_PASSWORD = os.getenv("APP_PASSWORD", "admin123")
    SESSION_TIMEOUT_MINUTES = int(os.getenv("SESSION_TIMEOUT_MINUTES", "30"))
    
    # Email Configuration
    SMTP_SENDER = os.getenv("SMTP_SENDER", "")
    SMTP_PASSWORD = os.getenv("SMTP_PASSWORD", "")
    SMTP_HOST = os.getenv("SMTP_HOST", "smtp.gmail.com")
    SMTP_PORT = int(os.getenv("SMTP_PORT", "465"))
    
    # Cache settings
    CACHE_ENABLED = os.getenv("CACHE_ENABLED", "true").lower() == "true"
    CACHE_EXPIRY_DAYS = int(os.getenv("CACHE_EXPIRY_DAYS", "7"))
    
    # Application settings
    MAX_BATCH_SIZE = int(os.getenv("MAX_BATCH_SIZE", "100"))
    ENABLE_ANALYTICS = os.getenv("ENABLE_ANALYTICS", "true").lower() == "true"
    
    # Supported organisms
    COMMON_ORGANISMS = [
        "Homo sapiens (Human)",
        "Mus musculus (Mouse)",
        "Rattus norvegicus (Rat)",
        "Danio rerio (Zebrafish)",
        "Drosophila melanogaster (Fruit fly)",
        "Caenorhabditis elegans (Nematode)",
        "Saccharomyces cerevisiae (Yeast)",
        "Escherichia coli (E. coli)",
        "Arabidopsis thaliana (Thale cress)",
        "Bos taurus (Cow)",
        "Oryza sativa (Rice)",
        "Zea mays (Maize)",
        "Schizosaccharomyces pombe (Fission yeast)",
        "Plasmodium falciparum (Malaria parasite)",
        "Chlamydomonas reinhardtii (Green algae)",
        "Other (type manually)"
    ]
    
    @classmethod
    def ensure_directories(cls):
        """Ensure all required directories exist."""
        cls.DATA_DIR.mkdir(exist_ok=True)
        cls.CACHE_DIR.mkdir(exist_ok=True)
        cls.LOGS_DIR.mkdir(exist_ok=True)
    
    @classmethod
    def to_dict(cls) -> Dict[str, Any]:
        """Export configuration as dictionary."""
        return {
            "ncbi_email": cls.NCBI_EMAIL,
            "ncbi_rate_limit": cls.NCBI_RATE_LIMIT,
            "cache_enabled": cls.CACHE_ENABLED,
            "cache_expiry_days": cls.CACHE_EXPIRY_DAYS,
            "max_batch_size": cls.MAX_BATCH_SIZE,
            "session_timeout": cls.SESSION_TIMEOUT_MINUTES,
            "analytics_enabled": cls.ENABLE_ANALYTICS
        }

# Initialize directories on import
Config.ensure_directories()
