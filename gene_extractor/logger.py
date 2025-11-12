"""Comprehensive logging configuration."""
import logging
import sys
from pathlib import Path
from logging.handlers import RotatingFileHandler
from datetime import datetime
from config import Config

class AppLogger:
    """Application logger with file and console handlers."""
    
    _loggers = {}
    
    @classmethod
    def get_logger(cls, name: str) -> logging.Logger:
        """Get or create a logger instance."""
        if name in cls._loggers:
            return cls._loggers[name]
        
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        
        # Remove existing handlers
        logger.handlers.clear()
        
        # Console handler (INFO level)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
        # File handler (DEBUG level) - rotating
        log_file = Config.LOGS_DIR / f"{name}.log"
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10*1024*1024,  # 10MB
            backupCount=5
        )
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        
        # Error file handler (ERROR level only)
        error_log_file = Config.LOGS_DIR / f"{name}_errors.log"
        error_handler = RotatingFileHandler(
            error_log_file,
            maxBytes=10*1024*1024,
            backupCount=3
        )
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(file_formatter)
        logger.addHandler(error_handler)
        
        cls._loggers[name] = logger
        return logger

# Create default loggers
app_logger = AppLogger.get_logger("app")
ncbi_logger = AppLogger.get_logger("ncbi")
auth_logger = AppLogger.get_logger("auth")
db_logger = AppLogger.get_logger("database")
analytics_logger = AppLogger.get_logger("analytics")
