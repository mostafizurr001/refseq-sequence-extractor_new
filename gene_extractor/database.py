"""Database operations for job history and caching."""
import sqlite3
import json
from datetime import datetime, timedelta
from typing import List, Dict, Optional, Any
from pathlib import Path
from config import Config
from logger import db_logger

class Database:
    """SQLite database manager for job history and caching."""
    
    def __init__(self, db_path: Path = Config.DATABASE_PATH):
        self.db_path = db_path
        self.init_database()
    
    def get_connection(self) -> sqlite3.Connection:
        """Get database connection."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn
    
    def init_database(self):
        """Initialize database tables."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                # Job history table
                cursor.execute("""
                    CREATE TABLE IF NOT EXISTS job_history (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        job_id TEXT UNIQUE NOT NULL,
                        username TEXT NOT NULL,
                        genes TEXT NOT NULL,
                        organism TEXT NOT NULL,
                        sequence_type TEXT NOT NULL,
                        output_format TEXT NOT NULL,
                        status TEXT NOT NULL,
                        total_genes INTEGER DEFAULT 0,
                        successful_genes INTEGER DEFAULT 0,
                        failed_genes INTEGER DEFAULT 0,
                        zip_path TEXT,
                        metadata_path TEXT,
                        error_message TEXT,
                        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                        completed_at TIMESTAMP
                    )
                """)
                
                # Results cache table
                cursor.execute("""
                    CREATE TABLE IF NOT EXISTS results_cache (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        cache_key TEXT UNIQUE NOT NULL,
                        gene TEXT NOT NULL,
                        organism TEXT NOT NULL,
                        sequence_type TEXT NOT NULL,
                        output_format TEXT NOT NULL,
                        accession TEXT,
                        sequence_data TEXT,
                        metadata TEXT,
                        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                        last_accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                        access_count INTEGER DEFAULT 1
                    )
                """)
                
                # Audit log table
                cursor.execute("""
                    CREATE TABLE IF NOT EXISTS audit_log (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        username TEXT NOT NULL,
                        action TEXT NOT NULL,
                        details TEXT,
                        ip_address TEXT,
                        timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                    )
                """)
                
                # Analytics table
                cursor.execute("""
                    CREATE TABLE IF NOT EXISTS analytics (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        event_type TEXT NOT NULL,
                        event_data TEXT,
                        timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                    )
                """)
                
                # Create indexes
                cursor.execute("CREATE INDEX IF NOT EXISTS idx_job_username ON job_history(username)")
                cursor.execute("CREATE INDEX IF NOT EXISTS idx_cache_key ON results_cache(cache_key)")
                cursor.execute("CREATE INDEX IF NOT EXISTS idx_audit_username ON audit_log(username)")
                cursor.execute("CREATE INDEX IF NOT EXISTS idx_analytics_type ON analytics(event_type)")
                
                conn.commit()
                db_logger.info("Database initialized successfully")
        except Exception as e:
            db_logger.error(f"Database initialization error: {e}")
            raise
    
    # Job History Methods
    def create_job(self, job_id: str, username: str, genes: List[str], 
                   organism: str, sequence_type: str, output_format: str) -> bool:
        """Create a new job entry."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    INSERT INTO job_history 
                    (job_id, username, genes, organism, sequence_type, output_format, status, total_genes)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (job_id, username, json.dumps(genes), organism, sequence_type, output_format, "running", len(genes)))
                conn.commit()
                db_logger.info(f"Job created: {job_id} for user {username}")
                return True
        except Exception as e:
            db_logger.error(f"Error creating job: {e}")
            return False
    
    def update_job_status(self, job_id: str, status: str, successful: int = 0, 
                          failed: int = 0, zip_path: str = None, 
                          metadata_path: str = None, error_message: str = None):
        """Update job status."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    UPDATE job_history 
                    SET status = ?, successful_genes = ?, failed_genes = ?,
                        zip_path = ?, metadata_path = ?, error_message = ?,
                        completed_at = CURRENT_TIMESTAMP
                    WHERE job_id = ?
                """, (status, successful, failed, zip_path, metadata_path, error_message, job_id))
                conn.commit()
                db_logger.info(f"Job updated: {job_id} - status: {status}")
        except Exception as e:
            db_logger.error(f"Error updating job: {e}")
    
    def get_user_jobs(self, username: str, limit: int = 50) -> List[Dict]:
        """Get job history for a user."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    SELECT * FROM job_history 
                    WHERE username = ? 
                    ORDER BY created_at DESC 
                    LIMIT ?
                """, (username, limit))
                rows = cursor.fetchall()
                return [dict(row) for row in rows]
        except Exception as e:
            db_logger.error(f"Error fetching jobs: {e}")
            return []
    
    def get_job_by_id(self, job_id: str) -> Optional[Dict]:
        """Get job details by ID."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("SELECT * FROM job_history WHERE job_id = ?", (job_id,))
                row = cursor.fetchone()
                return dict(row) if row else None
        except Exception as e:
            db_logger.error(f"Error fetching job: {e}")
            return None
    
    # Cache Methods
    def get_cached_result(self, cache_key: str) -> Optional[Dict]:
        """Get cached result if available and not expired."""
        if not Config.CACHE_ENABLED:
            return None
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                expiry_date = datetime.now() - timedelta(days=Config.CACHE_EXPIRY_DAYS)
                cursor.execute("""
                    SELECT * FROM results_cache 
                    WHERE cache_key = ? AND created_at > ?
                """, (cache_key, expiry_date.strftime('%Y-%m-%d %H:%M:%S')))
                row = cursor.fetchone()
                
                if row:
                    # Update access count and timestamp
                    cursor.execute("""
                        UPDATE results_cache 
                        SET last_accessed = CURRENT_TIMESTAMP, 
                            access_count = access_count + 1
                        WHERE cache_key = ?
                    """, (cache_key,))
                    conn.commit()
                    db_logger.info(f"Cache hit for key: {cache_key}")
                    return dict(row)
                return None
        except Exception as e:
            db_logger.error(f"Error fetching cache: {e}")
            return None
    
    def save_to_cache(self, cache_key: str, gene: str, organism: str,
                      sequence_type: str, output_format: str, accession: str,
                      sequence_data: str, metadata: Dict):
        """Save result to cache."""
        if not Config.CACHE_ENABLED:
            return
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    INSERT OR REPLACE INTO results_cache 
                    (cache_key, gene, organism, sequence_type, output_format, 
                     accession, sequence_data, metadata)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (cache_key, gene, organism, sequence_type, output_format,
                      accession, sequence_data, json.dumps(metadata)))
                conn.commit()
                db_logger.info(f"Cached result for: {cache_key}")
        except Exception as e:
            db_logger.error(f"Error saving to cache: {e}")
    
    def clear_expired_cache(self):
        """Clear expired cache entries."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                expiry_date = datetime.now() - timedelta(days=Config.CACHE_EXPIRY_DAYS)
                cursor.execute("""
                    DELETE FROM results_cache 
                    WHERE created_at < ?
                """, (expiry_date.strftime('%Y-%m-%d %H:%M:%S'),))
                deleted = cursor.rowcount
                conn.commit()
                db_logger.info(f"Cleared {deleted} expired cache entries")
                return deleted
        except Exception as e:
            db_logger.error(f"Error clearing cache: {e}")
            return 0
    
    # Audit Methods
    def log_action(self, username: str, action: str, details: str = None, ip_address: str = None):
        """Log user action to audit trail."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    INSERT INTO audit_log (username, action, details, ip_address)
                    VALUES (?, ?, ?, ?)
                """, (username, action, details, ip_address))
                conn.commit()
        except Exception as e:
            db_logger.error(f"Error logging action: {e}")
    
    def get_audit_logs(self, username: str = None, limit: int = 100) -> List[Dict]:
        """Get audit logs."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                if username:
                    cursor.execute("""
                        SELECT * FROM audit_log 
                        WHERE username = ? 
                        ORDER BY timestamp DESC 
                        LIMIT ?
                    """, (username, limit))
                else:
                    cursor.execute("""
                        SELECT * FROM audit_log 
                        ORDER BY timestamp DESC 
                        LIMIT ?
                    """, (limit,))
                rows = cursor.fetchall()
                return [dict(row) for row in rows]
        except Exception as e:
            db_logger.error(f"Error fetching audit logs: {e}")
            return []
    
    # Analytics Methods
    def track_event(self, event_type: str, event_data: Dict = None):
        """Track analytics event."""
        if not Config.ENABLE_ANALYTICS:
            return
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("""
                    INSERT INTO analytics (event_type, event_data)
                    VALUES (?, ?)
                """, (event_type, json.dumps(event_data) if event_data else None))
                conn.commit()
        except Exception as e:
            db_logger.error(f"Error tracking event: {e}")
    
    def get_analytics_summary(self, days: int = 30) -> Dict:
        """Get analytics summary."""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                since_date = datetime.now() - timedelta(days=days)
                
                # Total queries
                cursor.execute("""
                    SELECT COUNT(*) as total FROM analytics 
                    WHERE event_type = 'query_executed' 
                    AND timestamp > ?
                """, (since_date.strftime('%Y-%m-%d %H:%M:%S'),))
                total_queries = cursor.fetchone()[0]
                
                # Total genes searched
                cursor.execute("""
                    SELECT COUNT(*) as total FROM job_history 
                    WHERE created_at > ?
                """, (since_date.strftime('%Y-%m-%d %H:%M:%S'),))
                total_jobs = cursor.fetchone()[0]
                
                # Cache hit rate
                cursor.execute("""
                    SELECT COUNT(*) as hits FROM analytics 
                    WHERE event_type = 'cache_hit' 
                    AND timestamp > ?
                """, (since_date.strftime('%Y-%m-%d %H:%M:%S'),))
                cache_hits = cursor.fetchone()[0]
                
                return {
                    "total_queries": total_queries,
                    "total_jobs": total_jobs,
                    "cache_hits": cache_hits,
                    "cache_hit_rate": (cache_hits / total_queries * 100) if total_queries > 0 else 0
                }
        except Exception as e:
            db_logger.error(f"Error fetching analytics: {e}")
            return {}

# Global database instance
db = Database()
