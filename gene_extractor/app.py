

import streamlit as st
import pandas as pd
import io
import zipfile
import uuid
from datetime import datetime
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

# Import custom modules
from config import Config
from auth import AuthManager
from database import db
from ncbi_service import ncbi_service
from email_service import EmailService
from analytics import Analytics
from logger import app_logger
from utils import (
    parse_gene_input, is_refseq_accession, nucleotide_counts,
    format_timestamp, validate_email, parse_batch_file,
    create_fasta_format, get_organism_short_name, calculate_sequence_stats,
    sanitize_filename
)

# Page configuration
st.set_page_config(
    page_title="Gene Extractor Enterprise",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better UI
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background: #f8f9fa;
        padding: 1.5rem;
        border-radius: 8px;
        border-left: 4px solid #667eea;
        margin: 1rem 0;
    }
    .success-box {
        background: #d4edda;
        border: 1px solid #c3e6cb;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .warning-box {
        background: #fff3cd;
        border: 1px solid #ffeaa7;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .info-box {
        background: #d1ecf1;
        border: 1px solid #bee5eb;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .stButton>button {
        width: 100%;
        border-radius: 5px;
        height: 3em;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# Ensure authentication
AuthManager.require_auth()

# Initialize session state
if 'current_page' not in st.session_state:
    st.session_state.current_page = 'extractor'
if 'job_results' not in st.session_state:
    st.session_state.job_results = None

# Header
st.markdown("""
<div class="main-header">
    <h1>🧬 NCBI RefSeq Multi-Gene Sequence Extractor</h1>
    <p>Enterprise-Grade Bioinformatics Tool with Advanced Features</p>
</div>
""", unsafe_allow_html=True)

# Sidebar navigation
with st.sidebar:
    st.markdown("## 📋 Navigation")
    page = st.radio(
        "Select Module",
        ["🔬 Gene Extractor", "📊 Job History", "📈 Analytics", "⚙️ Settings"],
        label_visibility="collapsed"
    )
    
    st.markdown("---")
    
    # Quick stats
    st.markdown("### 📌 Quick Stats")
    username = st.session_state.get("username")
    recent_jobs = db.get_user_jobs(username, limit=5)
    st.metric("Recent Jobs", len(recent_jobs))
    
    stats = db.get_analytics_summary(days=7)
    st.metric("Queries (7d)", stats.get("total_queries", 0))
    st.metric("Cache Hit Rate", f"{stats.get('cache_hit_rate', 0):.1f}%")

# ====================
# GENE EXTRACTOR PAGE
# ====================
if page == "🔬 Gene Extractor":
    st.session_state.current_page = 'extractor'
    
    # Input method selection
    st.markdown("### 📝 Input Method")
    input_method = st.radio(
        "Choose input method",
        ["✍️ Manual Entry", "📁 Batch Upload (CSV/Excel)"],
        horizontal=True
    )
    
    genes = []
    
    if input_method == "✍️ Manual Entry":
        col1, col2 = st.columns([2, 1])
        
        with col1:
            gene_input = st.text_area(
                "🔍 Enter Gene Symbols or RefSeq Accessions",
                height=150,
                help="Enter one gene per line or comma-separated. Examples: BRCA1, TP53, NM_000546",
                placeholder="BRCA1\nTP53\nEGFR"
            )
            if gene_input:
                genes = parse_gene_input(gene_input)
                if genes:
                    st.success(f"✅ {len(genes)} unique genes identified")
                    if len(genes) > Config.MAX_BATCH_SIZE:
                        st.warning(f"⚠️ Maximum batch size is {Config.MAX_BATCH_SIZE}. Only first {Config.MAX_BATCH_SIZE} will be processed.")
                        genes = genes[:Config.MAX_BATCH_SIZE]
        
        with col2:
            st.markdown("#### 💡 Tips")
            st.info("""
            - One gene per line or comma-separated
            - Mix gene symbols and accessions
            - Duplicates are automatically removed
            - Supports RefSeq accessions (NM_, NP_, etc.)
            """)
    
    else:  # Batch Upload
        st.markdown("#### 📁 Upload Batch File")
        uploaded_file = st.file_uploader(
            "Upload CSV or Excel file with gene list",
            type=['csv', 'xlsx', 'xls'],
            help="File should contain a column named 'gene', 'genes', or 'gene_symbol'"
        )
        
        if uploaded_file:
            file_type = "csv" if uploaded_file.name.endswith('.csv') else "excel"
            genes, error = parse_batch_file(uploaded_file.read(), file_type)
            
            if error:
                st.error(f"❌ {error}")
            else:
                st.success(f"✅ {len(genes)} genes loaded from {uploaded_file.name}")
                if len(genes) > Config.MAX_BATCH_SIZE:
                    st.warning(f"⚠️ Maximum batch size is {Config.MAX_BATCH_SIZE}. Only first {Config.MAX_BATCH_SIZE} will be processed.")
                    genes = genes[:Config.MAX_BATCH_SIZE]
                
                # Show preview
                with st.expander("👁️ Preview Genes"):
                    st.write(genes[:20])
                    if len(genes) > 20:
                        st.caption(f"... and {len(genes) - 20} more")
    
    # Configuration section
    st.markdown("### ⚙️ Configuration")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        organism_choice = st.selectbox(
            "🌱 Organism",
            Config.COMMON_ORGANISMS,
            help="Select organism or choose 'Other' to enter manually"
        )
        if organism_choice == "Other (type manually)":
            organism = st.text_input("Type Organism Name", "Homo sapiens")
        else:
            organism = get_organism_short_name(organism_choice)
    
    with col2:
        sequence_type = st.selectbox(
            "📦 Sequence Type",
            ["Nucleotide", "Protein"]
        )
    
    with col3:
        output_format = st.selectbox(
            "📄 Output Format",
            ["GenBank", "FASTA"],
            help="GenBank required for 3'UTR extraction"
        )
    
    # Advanced options
    with st.expander("🔧 Advanced Options"):
        col1, col2 = st.columns(2)
        
        with col1:
            prefer_variant1 = st.checkbox(
                "Prefer transcript variant 1",
                value=True,
                help="Prioritize transcript variant 1 in search results"
            )
            
            exclude_variant_x = st.checkbox(
                "Exclude 'transcript variant X*' records",
                value=True,
                help="Filter out uncertain variant designations"
            )
        
        with col2:
            if output_format == "GenBank":
                selected_types = st.multiselect(
                    "🔖 Filter Feature Types",
                    ["CDS", "gene", "exon", "mRNA", "3'UTR", "5'UTR"],
                    default=["CDS", "gene"],
                    help="Select which features to display in annotations"
                )
            else:
                selected_types = []
    
    # Email notification
    col1, col2 = st.columns([1, 2])
    with col1:
        send_email = st.checkbox("✉️ Email results", value=False)
    with col2:
        if send_email:
            user_email = st.text_input("📧 Email Address", placeholder="your.email@example.com")
            if user_email and not validate_email(user_email):
                st.error("❌ Invalid email format")
        else:
            user_email = ""
    
    # Execute button
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        execute_button = st.button("🚀 Extract Sequences", type="primary", use_container_width=True)
    
    # Main execution logic
    if execute_button:
        if not genes:
            st.error("❌ Please enter at least one gene symbol or RefSeq accession.")
            st.stop()
        
        if send_email and (not user_email or not validate_email(user_email)):
            st.error("❌ Please provide a valid email address.")
            st.stop()
        
        # Create job
        job_id = f"job_{uuid.uuid4().hex[:8]}_{int(datetime.now().timestamp())}"
        username = st.session_state.get("username")
        
        # Initialize job in database
        db.create_job(job_id, username, genes, organism, sequence_type, output_format)
        db.log_action(username, "job_started", f"Job {job_id}: {len(genes)} genes")
        Analytics.track_query(len(genes), organism)
        
        app_logger.info(f"Starting job {job_id} for user {username}: {len(genes)} genes")
        
        # Progress tracking
        st.markdown("### 🔄 Processing")
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        metadata_rows = []
        zip_buffer = io.BytesIO()
        successful_count = 0
        failed_count = 0
        
        with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zipf:
            for idx, gene in enumerate(genes):
                progress = (idx + 1) / len(genes)
                progress_bar.progress(progress)
                status_text.markdown(f"🔎 Processing **{gene}** ({idx + 1}/{len(genes)})...")
                
                try:
                    # Check if it's an accession
                    is_accession = is_refseq_accession(gene)
                    
                    # Get sequence (with caching)
                    result = ncbi_service.get_sequence_with_cache(
                        gene, organism, sequence_type, output_format,
                        prefer_variant1, exclude_variant_x, is_accession
                    )
                    
                    if not result["success"]:
                        st.warning(f"⚠️ No sequence found for `{gene}`: {result.get('error', 'Unknown error')}")
                        failed_count += 1
                        continue
                    
                    sequence_data = result["sequence_data"]
                    accession = result["accession"]
                    from_cache = result.get("from_cache", False)
                    
                    if from_cache:
                        status_text.markdown(f"💾 Using cached data for **{gene}**")
                    
                    # Save to ZIP
                    filename = sanitize_filename(f"{gene}_{organism.replace(' ', '_')}.{output_format.lower()}")
                    zipf.writestr(filename, sequence_data)
                    
                    # Process based on format
                    if sequence_type == "Nucleotide" and output_format == "GenBank":
                        utr3, cds_start, cds_end, full_seq, gb = ncbi_service.extract_utr3_from_genbank(sequence_data)
                        
                        # Full sequence stats
                        full_stats = nucleotide_counts(full_seq)
                        metadata_rows.append({
                            "Gene": gene,
                            "Accession": accession,
                            "Region": "FULL",
                            "Organism": organism,
                            **full_stats
                        })
                        
                        # 3'UTR stats
                        if utr3:
                            utr_stats = nucleotide_counts(utr3)
                            metadata_rows.append({
                                "Gene": gene,
                                "Accession": accession,
                                "Region": "3'UTR",
                                "Organism": organism,
                                **utr_stats
                            })
                            
                            # Save 3'UTR as FASTA
                            utr_filename = sanitize_filename(f"{gene}_{organism.replace(' ', '_')}_3UTR.fasta")
                            utr_fasta = create_fasta_format(
                                f"{gene}|{organism}|3'UTR_from_{cds_end+1 if cds_end else 'NA'}",
                                utr3
                            )
                            zipf.writestr(utr_filename, utr_fasta)
                    
                    elif output_format == "FASTA":
                        lines = [ln.strip() for ln in sequence_data.splitlines() if ln.strip()]
                        seq = "".join(ln for ln in lines if not ln.startswith(">"))
                        stats = nucleotide_counts(seq) if sequence_type == "Nucleotide" else {}
                        
                        metadata_rows.append({
                            "Gene": gene,
                            "Accession": accession,
                            "Region": "FULL",
                            "Organism": organism,
                            **stats
                        })
                    
                    successful_count += 1
                    
                except Exception as e:
                    app_logger.error(f"Error processing {gene}: {e}")
                    st.error(f"❌ Error processing `{gene}`: {str(e)}")
                    Analytics.track_error("processing_error", str(e))
                    failed_count += 1
        
        progress_bar.progress(1.0)
        status_text.markdown("✅ **Processing complete!**")
        
        # Update job status
        zip_buffer.seek(0)
        zip_path = f"job_{job_id}.zip"
        db.update_job_status(job_id, "completed", successful_count, failed_count, zip_path, None)
        db.log_action(username, "job_completed", f"Job {job_id}: {successful_count} successful, {failed_count} failed")
        
        # Display results
        st.markdown("---")
        st.markdown("### 📊 Results Summary")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Genes", len(genes))
        with col2:
            st.metric("✅ Successful", successful_count)
        with col3:
            st.metric("❌ Failed", failed_count)
        with col4:
            success_rate = (successful_count / len(genes) * 100) if genes else 0
            st.metric("Success Rate", f"{success_rate:.1f}%")
        
        # Download section
        st.markdown("### 💾 Download Results")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.download_button(
                "📦 Download Sequences (ZIP)",
                zip_buffer.getvalue(),
                f"sequences_{job_id}.zip",
                "application/zip",
                use_container_width=True
            )
            Analytics.track_download("zip", len(zip_buffer.getvalue()))
        
        # Metadata section
        if metadata_rows:
            df = pd.DataFrame(metadata_rows)
            
            st.markdown("### 📋 Metadata & Analysis")
            
            # Display dataframe
            st.dataframe(df, use_container_width=True)
            
            # Statistics if nucleotide sequences
            if sequence_type == "Nucleotide" and 'GC%' in df.columns:
                st.markdown("#### 📈 Sequence Statistics")
                
                # Visualizations
                col1, col2 = st.columns(2)
                
                with col1:
                    # GC content distribution
                    fig_gc = px.histogram(
                        df[df['Region'] == 'FULL'],
                        x='GC%',
                        title='GC Content Distribution',
                        labels={'GC%': 'GC Content (%)'},
                        color_discrete_sequence=['#667eea']
                    )
                    st.plotly_chart(fig_gc, use_container_width=True)
                
                with col2:
                    # Sequence length distribution
                    fig_len = px.box(
                        df,
                        x='Region',
                        y='Length',
                        title='Sequence Length by Region',
                        color='Region',
                        color_discrete_sequence=['#667eea', '#764ba2']
                    )
                    st.plotly_chart(fig_len, use_container_width=True)
                
                # Nucleotide composition
                full_seqs = df[df['Region'] == 'FULL']
                if not full_seqs.empty:
                    total_nucleotides = {
                        'A': full_seqs['A'].sum(),
                        'T': full_seqs['T'].sum(),
                        'G': full_seqs['G'].sum(),
                        'C': full_seqs['C'].sum()
                    }
                    
                    fig_composition = go.Figure(data=[go.Pie(
                        labels=list(total_nucleotides.keys()),
                        values=list(total_nucleotides.values()),
                        hole=0.4,
                        marker_colors=['#667eea', '#764ba2', '#f093fb', '#4facfe']
                    )])
                    fig_composition.update_layout(title='Overall Nucleotide Composition')
                    st.plotly_chart(fig_composition, use_container_width=True)
            
            # Export options
            col1, col2 = st.columns(2)
            
            with col1:
                csv_bytes = df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    "📄 Download Metadata (CSV)",
                    csv_bytes,
                    f"metadata_{job_id}.csv",
                    "text/csv",
                    use_container_width=True
                )
            
            with col2:
                xlsx_buffer = io.BytesIO()
                with pd.ExcelWriter(xlsx_buffer, engine='xlsxwriter') as writer:
                    df.to_excel(writer, index=False, sheet_name='Metadata')
                xlsx_bytes = xlsx_buffer.getvalue()
                
                st.download_button(
                    "📊 Download Metadata (Excel)",
                    xlsx_bytes,
                    f"metadata_{job_id}.xlsx",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True
                )
        
        # Email results
        if send_email and user_email:
            with st.spinner("📧 Sending email..."):
                success, message = EmailService.send_results_email(
                    user_email, username,
                    zip_buffer.getvalue(),
                    csv_bytes if metadata_rows else None,
                    xlsx_bytes if metadata_rows else None,
                    len(genes), successful_count, organism
                )
                if success:
                    st.success(f"✅ {message} Check your inbox at {user_email}")
                    db.log_action(username, "email_sent", f"Results sent to {user_email}")
                else:
                    st.error(f"❌ {message}")

# ====================
# JOB HISTORY PAGE
# ====================
elif page == "📊 Job History":
    st.markdown("## 📊 Job History")
    
    username = st.session_state.get("username")
    jobs = db.get_user_jobs(username, limit=100)
    
    if not jobs:
        st.info("📭 No job history yet. Start by extracting some sequences!")
    else:
        st.success(f"Found {len(jobs)} jobs in your history")
        
        # Filters
        col1, col2, col3 = st.columns(3)
        with col1:
            status_filter = st.selectbox("Filter by Status", ["All", "completed", "running", "failed"])
        with col2:
            organism_filter = st.selectbox("Filter by Organism", ["All"] + list(set(j['organism'] for j in jobs)))
        with col3:
            limit = st.number_input("Show records", min_value=10, max_value=100, value=20)
        
        # Apply filters
        filtered_jobs = jobs[:limit]
        if status_filter != "All":
            filtered_jobs = [j for j in filtered_jobs if j['status'] == status_filter]
        if organism_filter != "All":
            filtered_jobs = [j for j in filtered_jobs if j['organism'] == organism_filter]
        
        # Display jobs
        for job in filtered_jobs:
            with st.expander(f"🔬 Job {job['job_id']} - {job['organism']} ({format_timestamp(job['created_at'])})"):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Status", job['status'].upper())
                    st.metric("Total Genes", job['total_genes'])
                
                with col2:
                    st.metric("Successful", job['successful_genes'])
                    st.metric("Failed", job['failed_genes'])
                
                with col3:
                    st.metric("Sequence Type", job['sequence_type'])
                    st.metric("Format", job['output_format'])
                
                # Show genes
                if job['genes']:
                    import json
                    genes_list = json.loads(job['genes'])
                    st.markdown("**Genes:**")
                    st.code(", ".join(genes_list[:10]) + ("..." if len(genes_list) > 10 else ""))
                
                # Re-run option
                if job['status'] == 'completed':
                    if st.button(f"🔄 Re-run this job", key=f"rerun_{job['id']}"):
                        st.info("Feature coming soon: Re-run job with same parameters")

# ====================
# ANALYTICS PAGE
# ====================
elif page == "📈 Analytics":
    Analytics.display_dashboard()

# ====================
# SETTINGS PAGE
# ====================
elif page == "⚙️ Settings":
    st.markdown("## ⚙️ Settings & Configuration")
    
    # Display current configuration
    st.markdown("### 📋 Current Configuration")
    
    config_dict = Config.to_dict()
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### NCBI Settings")
        st.info(f"**Email:** {config_dict['ncbi_email']}")
        st.info(f"**Rate Limit:** {config_dict['ncbi_rate_limit']} req/sec")
    
    with col2:
        st.markdown("#### Application Settings")
        st.info(f"**Session Timeout:** {config_dict['session_timeout']} minutes")
        st.info(f"**Max Batch Size:** {config_dict['max_batch_size']} genes")
    
    st.markdown("#### Cache Settings")
    col1, col2 = st.columns(2)
    with col1:
        st.info(f"**Cache Enabled:** {'✅ Yes' if config_dict['cache_enabled'] else '❌ No'}")
    with col2:
        st.info(f"**Cache Expiry:** {config_dict['cache_expiry_days']} days")
    
    # Cache management
    st.markdown("### 🗄️ Cache Management")
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("🧹 Clear Expired Cache", use_container_width=True):
            deleted = db.clear_expired_cache()
            st.success(f"✅ Cleared {deleted} expired cache entries")
            db.log_action(st.session_state.username, "cache_cleared", f"{deleted} entries")
    
    # Help section
    st.markdown("### 📚 Help & Documentation")
    
    with st.expander("🔍 How to use this tool"):
        st.markdown("""
        **Basic Usage:**
        1. Enter gene symbols or RefSeq accessions (one per line or comma-separated)
        2. Select organism, sequence type, and output format
        3. Configure advanced options if needed
        4. Click "Extract Sequences" to start processing
        5. Download results as ZIP file with metadata
        
        **Batch Upload:**
        - Upload CSV or Excel file with a column containing gene names
        - Supported column names: 'gene', 'genes', 'gene_symbol', 'gene_name', 'symbol'
        - Maximum batch size: {Config.MAX_BATCH_SIZE} genes
        
        **Features:**
        - Automatic caching of results for faster repeated queries
        - 3'UTR extraction for nucleotide sequences (GenBank format)
        - Nucleotide composition analysis and visualization
        - Email delivery of results with professional reports
        - Job history tracking and audit logging
        """)
    
    with st.expander("⚡ Performance Tips"):
        st.markdown("""
        - Use cache to avoid re-fetching same sequences
        - Batch upload for large gene lists (faster than manual entry)
        - Enable "Prefer transcript variant 1" for more consistent results
        - Results are cached for {Config.CACHE_EXPIRY_DAYS} days
        - NCBI API rate limit: {Config.NCBI_RATE_LIMIT} requests per second
        """)
    
    with st.expander("🔐 Security & Privacy"):
        st.markdown("""
        - All user actions are logged for audit purposes
        - Session timeout: {Config.SESSION_TIMEOUT_MINUTES} minutes of inactivity
        - Account lockout after 5 failed login attempts (15 minutes)
        - Email results are sent securely via encrypted SMTP
        - No sequences or metadata are permanently stored on servers
        """)

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 1rem;'>
    <p>🧬 <strong>Gene Extractor Enterprise Edition</strong> | Powered by NCBI RefSeq Database</p>
    <p style='font-size: 0.8rem;'>For support or feature requests, contact your system administrator</p>
</div>
""", unsafe_allow_html=True)
