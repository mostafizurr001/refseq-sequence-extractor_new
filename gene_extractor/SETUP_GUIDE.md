# 🚀 Setup & Testing Guide

## Quick Setup (5 minutes)

### Step 1: Configure Environment

```bash
cd /app/gene_extractor

# Create .env file from template
cp .env.example .env

# Edit the .env file (minimum required changes)
nano .env
```

**Minimum Required Configuration:**

```env
# Required: Your email (NCBI requirement)
NCBI_EMAIL=your.email@example.com

# Required: Change default credentials
APP_USERNAME=admin
APP_PASSWORD=YourSecurePassword123

# Optional: For email functionality
SMTP_SENDER=your.email@gmail.com
SMTP_PASSWORD=your_app_password
```

### Step 2: Install Dependencies

```bash
# Install Python packages
pip install -r requirements.txt
```

### Step 3: Test the Application

```bash
# Quick test of core modules
python3 -c "from config import Config; print('✅ Config loaded')"
python3 -c "from database import db; print('✅ Database initialized')"
python3 -c "from ncbi_service import ncbi_service; print('✅ NCBI service ready')"
```

### Step 4: Run the Application

```bash
# Option 1: Using the startup script
./run.sh

# Option 2: Direct Streamlit command
streamlit run app.py --server.port 8501 --server.address 0.0.0.0
```

The application will be available at: **http://localhost:8501**

## 🧪 Testing Guide

### Test 1: Basic Authentication

1. Open http://localhost:8501
2. Login with credentials from `.env` file:
   - Username: `admin` (or your configured username)
   - Password: `YourSecurePassword123` (or your configured password)
3. Verify successful login

### Test 2: Single Gene Extraction

1. After login, go to "🔬 Gene Extractor"
2. Enter a test gene: `BRCA1`
3. Select organism: `Homo sapiens (Human)`
4. Keep defaults: Nucleotide, GenBank format
5. Click "🚀 Extract Sequences"
6. Verify:
   - Progress bar shows
   - Results summary displays
   - Download buttons appear
   - Metadata table shows statistics

### Test 3: Multiple Genes

Enter multiple genes (one per line):
```
BRCA1
TP53
EGFR
KRAS
```

Expected: All 4 genes should be processed successfully

### Test 4: Batch Upload

1. Create a test CSV file (`test_genes.csv`):
```csv
gene
BRCA1
TP53
EGFR
MYC
PTEN
```

2. Select "📁 Batch Upload (CSV/Excel)"
3. Upload the file
4. Verify gene list is loaded
5. Process and download results

### Test 5: Cache Functionality

1. Search for `BRCA1` first time (fetches from NCBI)
2. Search for `BRCA1` again immediately
3. Look for "💾 Using cached data" message
4. Verify faster processing time

### Test 6: Job History

1. Go to "📊 Job History" page
2. Verify previous jobs are listed
3. Check job details and statistics
4. Test filters (status, organism)

### Test 7: Analytics Dashboard

1. Go to "📈 Analytics" page
2. Verify metrics display:
   - Total Queries
   - Total Jobs
   - Cache Hit Rate
3. Check recent activity table
4. View audit logs

### Test 8: 3'UTR Extraction

1. Search for gene: `TP53`
2. Organism: `Homo sapiens`
3. Format: **GenBank** (required for UTR)
4. After processing, check:
   - Both FULL and 3'UTR rows in metadata
   - 3'UTR FASTA file in ZIP download
   - Nucleotide composition for both regions

### Test 9: Visualizations

After processing multiple genes:
1. Scroll to "📈 Sequence Statistics" section
2. Verify charts appear:
   - GC Content Distribution (histogram)
   - Sequence Length by Region (box plot)
   - Overall Nucleotide Composition (pie chart)

### Test 10: Session Timeout

1. Login to application
2. Wait for 30+ minutes (or set SESSION_TIMEOUT_MINUTES=1 for testing)
3. Try to perform an action
4. Verify session timeout message appears
5. Verify automatic redirect to login

## 🎯 Sample Test Data

### Quick Test Genes (Homo sapiens)
```
BRCA1
TP53
EGFR
KRAS
MYC
```

### RefSeq Accessions Test
```
NM_000546
NM_007294
NM_005228
```

### Multi-Organism Test
- `INS` - Homo sapiens (Human)
- `Ins2` - Mus musculus (Mouse)
- `insulin` - Danio rerio (Zebrafish)

## 📊 Expected Results

### Successful Extraction
- **Status**: ✅ Successful
- **Files in ZIP**:
  - `GENENAME_Homo_sapiens.genbank`
  - `GENENAME_Homo_sapiens_3UTR.fasta` (if applicable)
- **Metadata CSV**:
  - Gene, Accession, Region, Length, GC%, A, T, G, C counts

### Performance Benchmarks
- **Single gene**: 2-5 seconds (first time), <1 second (cached)
- **10 genes**: 20-50 seconds (first time), 5-10 seconds (with cache)
- **Cache hit rate**: Should reach 70%+ with repeated queries

## 🔧 Troubleshooting Tests

### Test NCBI Connection
```python
python3 -c "
from Bio import Entrez
from config import Config
Entrez.email = Config.NCBI_EMAIL
handle = Entrez.esearch(db='nuccore', term='BRCA1[Gene] AND Homo sapiens[Organism]', retmax=1)
result = Entrez.read(handle)
print('✅ NCBI connection working')
print(f'Found {len(result[\"IdList\"])} results')
"
```

### Test Database
```python
python3 -c "
from database import db
jobs = db.get_user_jobs('admin', limit=10)
print(f'✅ Database working: {len(jobs)} jobs found')
"
```

### Test Email (if configured)
```python
python3 -c "
from email_service import EmailService
from config import Config
if Config.SMTP_SENDER and Config.SMTP_PASSWORD:
    print('✅ Email configuration found')
else:
    print('⚠️ Email not configured (optional)')
"
```

## 🐛 Common Issues & Solutions

### Issue: "NCBI email not configured"
**Solution**: Set `NCBI_EMAIL` in `.env` file

### Issue: "No sequences found"
**Solution**: 
- Check gene symbol spelling
- Try with RefSeq accession (e.g., NM_000546)
- Verify organism name is correct

### Issue: "Rate limit exceeded"
**Solution**:
- Get NCBI API key: https://www.ncbi.nlm.nih.gov/account/
- Add to `.env`: `NCBI_API_KEY=your_key_here`
- This increases rate limit from 3 to 10 requests/second

### Issue: Email not sending
**Solution**:
- For Gmail: Use App Password, not regular password
- Enable 2FA and generate app password
- Update SMTP_PASSWORD in `.env`

### Issue: Database locked
**Solution**:
```bash
# Stop application
# Remove lock (if exists)
rm data/gene_extractor.db-wal data/gene_extractor.db-shm
# Restart application
```

## 📈 Performance Optimization

### For Production Use:

1. **Get NCBI API Key**
   - Increases rate limit to 10 req/sec
   - Register at: https://www.ncbi.nlm.nih.gov/account/

2. **Enable Caching**
   ```env
   CACHE_ENABLED=true
   CACHE_EXPIRY_DAYS=7
   ```

3. **Optimize Batch Size**
   ```env
   MAX_BATCH_SIZE=100  # Adjust based on your needs
   ```

4. **Monitor Logs**
   ```bash
   tail -f data/logs/app.log
   tail -f data/logs/ncbi.log
   ```

## ✅ Verification Checklist

- [ ] Application starts without errors
- [ ] Login works with configured credentials
- [ ] Single gene extraction succeeds
- [ ] Multiple genes process correctly
- [ ] Batch upload accepts CSV/Excel
- [ ] Results download as ZIP
- [ ] Metadata displays with statistics
- [ ] Visualizations render (charts)
- [ ] Job history shows previous jobs
- [ ] Analytics dashboard displays metrics
- [ ] Cache functionality works
- [ ] Session timeout occurs after inactivity
- [ ] Audit logs capture user actions
- [ ] Settings page shows configuration

## 🎓 Advanced Testing

### Load Testing (10+ genes)
```python
# Create large test file
genes = ["BRCA1", "TP53", "EGFR", "KRAS", "MYC", "PTEN", "AKT1", "BRAF", "PIK3CA", "NRAS"]
```

### Multi-Organism Testing
Test same gene across different organisms to verify organism filtering works correctly.

### Edge Cases
- Empty input
- Invalid gene symbols
- Non-existent accessions
- Special characters in input
- Very long gene lists (>50)

## 📞 Need Help?

Check the logs for detailed error information:
```bash
# Application logs
cat data/logs/app.log

# Error logs only
cat data/logs/app_errors.log

# NCBI API logs
cat data/logs/ncbi.log
```

## 🎉 Success Criteria

Your setup is complete when:
1. ✅ All 10 basic tests pass
2. ✅ No errors in application logs
3. ✅ Cache hit rate improves over time
4. ✅ Job history tracks all operations
5. ✅ Email delivery works (if configured)
6. ✅ Visualizations render correctly
7. ✅ Session management functions properly
8. ✅ Audit logs capture all actions

---

**Ready to use!** Your enterprise-grade gene extractor is now fully operational. 🚀
