"""Enhanced authentication and session management."""
import streamlit as st
from datetime import datetime, timedelta
from typing import Optional
from config import Config
from logger import auth_logger
from database import db

class AuthManager:
    """Manage user authentication and sessions."""
    
    @staticmethod
    def init_session():
        """Initialize session state variables."""
        if "auth_ok" not in st.session_state:
            st.session_state.auth_ok = False
        if "username" not in st.session_state:
            st.session_state.username = None
        if "last_activity" not in st.session_state:
            st.session_state.last_activity = None
        if "login_attempts" not in st.session_state:
            st.session_state.login_attempts = 0
        if "locked_until" not in st.session_state:
            st.session_state.locked_until = None
    
    @staticmethod
    def check_session_timeout() -> bool:
        """Check if session has timed out."""
        if not st.session_state.auth_ok:
            return False
        
        if st.session_state.last_activity:
            timeout = timedelta(minutes=Config.SESSION_TIMEOUT_MINUTES)
            if datetime.now() - st.session_state.last_activity > timeout:
                auth_logger.warning(f"Session timeout for user: {st.session_state.username}")
                AuthManager.logout()
                return True
        
        # Update last activity
        st.session_state.last_activity = datetime.now()
        return False
    
    @staticmethod
    def is_account_locked() -> tuple[bool, Optional[str]]:
        """Check if account is locked due to failed attempts."""
        if st.session_state.locked_until:
            if datetime.now() < st.session_state.locked_until:
                remaining = (st.session_state.locked_until - datetime.now()).seconds
                return True, f"Account locked. Try again in {remaining} seconds."
            else:
                st.session_state.locked_until = None
                st.session_state.login_attempts = 0
        return False, None
    
    @staticmethod
    def authenticate(username: str, password: str) -> tuple[bool, str]:
        """Authenticate user credentials."""
        # Check if account is locked
        locked, message = AuthManager.is_account_locked()
        if locked:
            return False, message
        
        # Validate credentials
        if username == Config.APP_USERNAME and password == Config.APP_PASSWORD:
            st.session_state.auth_ok = True
            st.session_state.username = username
            st.session_state.last_activity = datetime.now()
            st.session_state.login_attempts = 0
            
            # Log successful login
            auth_logger.info(f"Successful login: {username}")
            db.log_action(username, "login", "Successful authentication")
            db.track_event("user_login", {"username": username})
            
            return True, "Login successful!"
        else:
            st.session_state.login_attempts += 1
            
            # Lock account after 5 failed attempts
            if st.session_state.login_attempts >= 5:
                st.session_state.locked_until = datetime.now() + timedelta(minutes=15)
                auth_logger.warning(f"Account locked due to failed attempts: {username}")
                db.log_action(username, "login_failed", f"Account locked after {st.session_state.login_attempts} attempts")
                return False, "Too many failed attempts. Account locked for 15 minutes."
            
            auth_logger.warning(f"Failed login attempt for: {username}")
            db.log_action(username, "login_failed", f"Attempt {st.session_state.login_attempts}")
            
            remaining = 5 - st.session_state.login_attempts
            return False, f"Invalid credentials. {remaining} attempts remaining."
    
    @staticmethod
    def logout():
        """Logout user and clear session."""
        username = st.session_state.get("username")
        if username:
            auth_logger.info(f"User logged out: {username}")
            db.log_action(username, "logout", "User logged out")
        
        st.session_state.auth_ok = False
        st.session_state.username = None
        st.session_state.last_activity = None
    
    @staticmethod
    def require_auth():
        """Require authentication before accessing the app."""
        AuthManager.init_session()
        
        # Check session timeout
        if AuthManager.check_session_timeout():
            st.warning("⏰ Session expired. Please login again.")
            st.rerun()
        
        # Show login form if not authenticated
        if not st.session_state.auth_ok:
            st.markdown("""
            <div style='text-align: center; padding: 2rem;'>
                <h2>🔐 Login Required</h2>
                <p>Please authenticate to access the Gene Extractor application.</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Check if account is locked
            locked, lock_message = AuthManager.is_account_locked()
            if locked:
                st.error(lock_message)
                st.stop()
            
            with st.form("login_form", clear_on_submit=False):
                col1, col2, col3 = st.columns([1, 2, 1])
                with col2:
                    username = st.text_input("👤 Username", key="login_username")
                    password = st.text_input("🔑 Password", type="password", key="login_password")
                    
                    col_a, col_b, col_c = st.columns([1, 1, 1])
                    with col_b:
                        submitted = st.form_submit_button("🚀 Login", use_container_width=True)
            
            if submitted:
                success, message = AuthManager.authenticate(username, password)
                if success:
                    st.success(message)
                    st.rerun()
                else:
                    st.error(message)
            
            st.stop()
        
        # Show logout button in sidebar for authenticated users
        with st.sidebar:
            st.markdown(f"**👤 User:** {st.session_state.username}")
            if st.session_state.last_activity:
                st.caption(f"Last activity: {st.session_state.last_activity.strftime('%H:%M:%S')}")
            
            if st.button("🚪 Logout", use_container_width=True):
                AuthManager.logout()
                st.rerun()
