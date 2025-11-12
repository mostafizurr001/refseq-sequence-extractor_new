"""Analytics and usage tracking."""
import streamlit as st
from typing import Dict, List
from datetime import datetime, timedelta
from database import db
from logger import analytics_logger
import pandas as pd

class Analytics:
    """Analytics dashboard and tracking."""
    
    @staticmethod
    def display_dashboard():
        """Display analytics dashboard."""
        st.markdown("### 📊 Usage Analytics")
        
        # Time period selector
        col1, col2 = st.columns([3, 1])
        with col2:
            period = st.selectbox("Period", ["7 days", "30 days", "90 days", "All time"])
        
        days = {
            "7 days": 7,
            "30 days": 30,
            "90 days": 90,
            "All time": 36500
        }[period]
        
        # Get analytics data
        stats = db.get_analytics_summary(days=days)
        
        # Display metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(
                label="Total Queries",
                value=stats.get("total_queries", 0)
            )
        
        with col2:
            st.metric(
                label="Total Jobs",
                value=stats.get("total_jobs", 0)
            )
        
        with col3:
            st.metric(
                label="Cache Hits",
                value=stats.get("cache_hits", 0)
            )
        
        with col4:
            cache_rate = stats.get("cache_hit_rate", 0)
            st.metric(
                label="Cache Hit Rate",
                value=f"{cache_rate:.1f}%"
            )
        
        # Recent activity
        st.markdown("#### Recent Activity")
        username = st.session_state.get("username")
        recent_jobs = db.get_user_jobs(username, limit=10)
        
        if recent_jobs:
            df = pd.DataFrame(recent_jobs)
            display_cols = ['job_id', 'organism', 'status', 'total_genes', 'successful_genes', 'created_at']
            available_cols = [col for col in display_cols if col in df.columns]
            st.dataframe(df[available_cols], use_container_width=True)
        else:
            st.info("No recent activity")
        
        # Audit logs
        with st.expander("📝 View Audit Log"):
            logs = db.get_audit_logs(username, limit=50)
            if logs:
                df_logs = pd.DataFrame(logs)
                st.dataframe(df_logs, use_container_width=True)
            else:
                st.info("No audit logs available")
    
    @staticmethod
    def track_query(gene_count: int, organism: str):
        """Track a query execution."""
        try:
            db.track_event("query_executed", {
                "gene_count": gene_count,
                "organism": organism,
                "timestamp": datetime.now().isoformat()
            })
            analytics_logger.info(f"Query tracked: {gene_count} genes, {organism}")
        except Exception as e:
            analytics_logger.error(f"Error tracking query: {e}")
    
    @staticmethod
    def track_download(file_type: str, file_size: int):
        """Track a file download."""
        try:
            db.track_event("file_download", {
                "file_type": file_type,
                "file_size": file_size,
                "timestamp": datetime.now().isoformat()
            })
            analytics_logger.info(f"Download tracked: {file_type}, {file_size} bytes")
        except Exception as e:
            analytics_logger.error(f"Error tracking download: {e}")
    
    @staticmethod
    def track_error(error_type: str, error_message: str):
        """Track an error occurrence."""
        try:
            db.track_event("error", {
                "error_type": error_type,
                "error_message": error_message,
                "timestamp": datetime.now().isoformat()
            })
            analytics_logger.warning(f"Error tracked: {error_type}")
        except Exception as e:
            analytics_logger.error(f"Error tracking error: {e}")
