# 🧬 NCBI RefSeq Multi-Gene Sequence Extractor - Enterprise Edition

## Overview

An enterprise-grade bioinformatics tool for extracting gene sequences from NCBI RefSeq database with advanced features including batch processing, caching, analytics, and comprehensive audit logging.

## ✨ Features

### Core Functionality
- **Multi-gene sequence extraction** from NCBI RefSeq
- **Support for 16+ organisms** with manual entry option
- **Nucleotide and protein sequences**
- **GenBank and FASTA formats**
- **3'UTR extraction** from mRNA sequences
- **Nucleotide composition analysis**

### Enterprise Features
- ✅ **Batch upload** via CSV/Excel files
- ✅ **Result caching** for faster repeated queries
- ✅ **Job history** tracking and management
- ✅ **Session management** with timeout
- ✅ **Audit logging** of all user actions
- ✅ **Usage analytics** dashboard
- ✅ **Rate limiting** for NCBI API
- ✅ **Professional email** reports with HTML templates
- ✅ **Interactive visualizations** (GC content, length distributions)
- ✅ **Account security** with lockout protection

## 🚀 Quick Start

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Installation

```bash
# Navigate to the application directory
cd /app/gene_extractor

# Install dependencies
pip install -r requirements.txt

# Copy environment template
cp .env.example .env

# Edit .env file with your credentials
nano .env
```

### Configuration

Edit `.env` file and set the following:

1. **NCBI Email** (required by NCBI):
   ```
   NCBI_EMAIL=your_email@example.com
   ```

2. **Authentication** (change defaults):
   ```
   APP_USERNAME=admin
   APP_PASSWORD=your_secure_password
   ```

3. **SMTP Settings** (for email delivery):
   ```
   SMTP_SENDER=your_email@gmail.com
   SMTP_PASSWORD=your_app_password
   ```

4. **Optional: NCBI API Key** (for higher rate limits):
   ```
   NCBI_API_KEY=your_ncbi_api_key
   ```
   Get your API key from: https://www.ncbi.nlm.nih.gov/account/settings/

### Running the Application

```bash
streamlit run app.py
```

The application will start on `http://localhost:8501`

## 📖 User Guide

### Basic Usage

1. **Login** with your credentials
2. **Enter genes** manually or upload a batch file
3. **Configure** organism, sequence type, and format
4. **Extract** sequences and download results
5. **View** metadata with composition analysis

### Batch Upload

1. Prepare CSV or Excel file with gene list
2. Ensure file has a column named: `gene`, `genes`, `gene_symbol`, or `gene_name`
3. Upload file in the "Batch Upload" section
4. Review preview and proceed with extraction

**Example CSV format:**
```csv
gene
BRCA1
TP53
EGFR
KRAS
```

### Job History

- View all previous extraction jobs
- Filter by status, organism, or date
- Re-run previous jobs with same parameters
- Download results from completed jobs

### Analytics Dashboard

- View usage statistics (queries, cache hit rate)
- Monitor recent activity
- Access audit logs
- Track performance metrics

## 🔧 Advanced Configuration

### Cache Settings

```env
CACHE_ENABLED=true          # Enable/disable caching
CACHE_EXPIRY_DAYS=7         # Cache expiration period
```

### Rate Limiting

```env
NCBI_RATE_LIMIT=3           # Requests per second (3 without API key, 10 with key)
NCBI_MAX_RETRIES=3          # Retry failed requests
```

### Session Management

```env
SESSION_TIMEOUT_MINUTES=30  # Inactive session timeout
```

### Batch Processing

```env
MAX_BATCH_SIZE=100          # Maximum genes per batch
```

## 📊 Database Schema

The application uses SQLite with the following tables:

- **job_history**: Track all extraction jobs
- **results_cache**: Cache NCBI query results
- **audit_log**: Log all user actions
- **analytics**: Track usage metrics

## 🔐 Security Features

### Authentication
- Username/password authentication
- Session timeout after inactivity
- Account lockout after failed attempts
- Logout functionality

### Audit Trail
- All actions are logged with timestamp
- User activity tracking
- IP address logging (when available)
- Full audit history available in Analytics page

### Data Privacy
- Sequences are cached temporarily
- Automatic cache expiration
- No permanent storage of sensitive data
- Secure email delivery via encrypted SMTP

## 📧 Email Configuration

### Gmail Setup

1. Enable 2-factor authentication in Gmail
2. Generate an App Password:
   - Go to Google Account → Security → App Passwords
   - Create password for "Mail"
3. Use the generated password in `.env`:
   ```
   SMTP_SENDER=your_email@gmail.com
   SMTP_PASSWORD=generated_app_password
   ```

### Other Email Providers

Update SMTP settings accordingly:
```env
SMTP_HOST=smtp.yourprovider.com
SMTP_PORT=465  # or 587 for TLS
```

## 🐛 Troubleshooting

### NCBI Connection Issues

**Error:** "NCBI API rate limit exceeded"
- **Solution**: Get an API key from NCBI or reduce rate limit

**Error:** "No sequences found"
- **Solution**: Check gene symbol spelling, try with RefSeq accession

### Email Delivery Issues

**Error:** "SMTP authentication failed"
- **Solution**: Use App Password instead of regular password (Gmail)

**Error:** "Connection timeout"
- **Solution**: Check firewall settings, verify SMTP port

### Performance Issues

**Slow queries:**
- Enable caching in settings
- Get NCBI API key for faster rate limits
- Process smaller batches

### Database Issues

**Database locked:**
- Close other connections to the database
- Restart the application

## 📝 Logs

Logs are stored in `/app/gene_extractor/data/logs/`:

- `app.log`: General application logs
- `ncbi.log`: NCBI API interactions
- `auth.log`: Authentication events
- `database.log`: Database operations
- `analytics.log`: Analytics tracking
- `*_errors.log`: Error-only logs

## 🔄 Maintenance

### Clear Expired Cache

1. Go to Settings page
2. Click "Clear Expired Cache"
3. Confirm deletion

### Database Backup

```bash
cp data/gene_extractor.db data/gene_extractor.db.backup
```

### View Logs

```bash
# View recent errors
tail -f data/logs/app_errors.log

# View all application logs
tail -f data/logs/app.log
```

## 🎯 Best Practices

1. **Use batch upload** for large gene lists (>10 genes)
2. **Enable caching** to reduce NCBI API calls
3. **Get NCBI API key** for production use
4. **Regular backup** of database file
5. **Monitor logs** for errors and issues
6. **Clear cache** periodically to save disk space
7. **Use strong passwords** for authentication

## 📚 Resources

- [NCBI RefSeq Database](https://www.ncbi.nlm.nih.gov/refseq/)
- [Biopython Documentation](https://biopython.org/)
- [Streamlit Documentation](https://docs.streamlit.io/)
- [NCBI API Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)

## 🤝 Support

For issues, questions, or feature requests:
- Check the troubleshooting section
- Review application logs
- Contact your system administrator

## 📄 License

This is an enterprise application for internal use.

## 🔖 Version

**Version:** 2.0.0 (Enterprise Edition)
**Last Updated:** 2025

---

**Note:** This application requires proper NCBI credentials and SMTP configuration for full functionality.
