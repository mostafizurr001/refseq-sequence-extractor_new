#!/usr/bin/env python3
"""
Verification script to test core modules and configuration
Run this after setup to ensure everything is working correctly
"""

import sys
import os

def test_imports():
    """Test all module imports."""
    print("🔍 Testing module imports...")
    try:
        from config import Config
        print("  ✅ config.py imported successfully")
        
        from logger import app_logger
        print("  ✅ logger.py imported successfully")
        
        from database import db
        print("  ✅ database.py imported successfully")
        
        from auth import AuthManager
        print("  ✅ auth.py imported successfully")
        
        from ncbi_service import ncbi_service
        print("  ✅ ncbi_service.py imported successfully")
        
        from email_service import EmailService
        print("  ✅ email_service.py imported successfully")
        
        from analytics import Analytics
        print("  ✅ analytics.py imported successfully")
        
        from utils import parse_gene_input, is_refseq_accession
        print("  ✅ utils.py imported successfully")
        
        return True
    except Exception as e:
        print(f"  ❌ Import error: {e}")
        return False

def test_configuration():
    """Test configuration loading."""
    print("\n⚙️  Testing configuration...")
    try:
        from config import Config
        
        config_dict = Config.to_dict()
        print(f"  ✅ NCBI Email: {config_dict['ncbi_email']}")
        print(f"  ✅ Rate Limit: {config_dict['ncbi_rate_limit']} req/sec")
        print(f"  ✅ Cache Enabled: {config_dict['cache_enabled']}")
        print(f"  ✅ Max Batch Size: {config_dict['max_batch_size']}")
        
        # Check if directories exist
        if Config.DATA_DIR.exists():
            print(f"  ✅ Data directory exists: {Config.DATA_DIR}")
        else:
            print(f"  ⚠️  Data directory missing: {Config.DATA_DIR}")
        
        return True
    except Exception as e:
        print(f"  ❌ Configuration error: {e}")
        return False

def test_database():
    """Test database initialization."""
    print("\n🗄️  Testing database...")
    try:
        from database import db
        from config import Config
        
        # Check if database file exists
        if Config.DATABASE_PATH.exists():
            print(f"  ✅ Database file exists: {Config.DATABASE_PATH}")
        else:
            print(f"  ⚠️  Database file will be created on first run")
        
        # Test database operations
        test_username = "test_user"
        jobs = db.get_user_jobs(test_username, limit=1)
        print(f"  ✅ Database query successful (found {len(jobs)} jobs)")
        
        return True
    except Exception as e:
        print(f"  ❌ Database error: {e}")
        return False

def test_logger():
    """Test logging system."""
    print("\n📝 Testing logger...")
    try:
        from logger import app_logger
        from config import Config
        
        # Test log write
        app_logger.info("Verification test log entry")
        
        # Check if log directory exists
        if Config.LOGS_DIR.exists():
            log_files = list(Config.LOGS_DIR.glob("*.log"))
            print(f"  ✅ Log directory exists: {Config.LOGS_DIR}")
            print(f"  ✅ Found {len(log_files)} log files")
        else:
            print(f"  ⚠️  Log directory will be created on first run")
        
        return True
    except Exception as e:
        print(f"  ❌ Logger error: {e}")
        return False

def test_dependencies():
    """Test required dependencies."""
    print("\n📦 Testing dependencies...")
    
    dependencies = [
        ('streamlit', 'Streamlit'),
        ('pandas', 'Pandas'),
        ('Bio', 'Biopython'),
        ('plotly', 'Plotly'),
        ('xlsxwriter', 'XlsxWriter'),
        ('openpyxl', 'OpenPyXL'),
    ]
    
    all_ok = True
    for module, name in dependencies:
        try:
            __import__(module)
            print(f"  ✅ {name} installed")
        except ImportError:
            print(f"  ❌ {name} NOT installed")
            all_ok = False
    
    return all_ok

def test_ncbi_connection():
    """Test NCBI connection (optional)."""
    print("\n🧬 Testing NCBI connection (this may take a few seconds)...")
    try:
        from Bio import Entrez
        from config import Config
        
        Entrez.email = Config.NCBI_EMAIL
        if Config.NCBI_API_KEY:
            Entrez.api_key = Config.NCBI_API_KEY
        
        # Test simple search
        handle = Entrez.esearch(db="nuccore", term="BRCA1[Gene] AND Homo sapiens[Organism]", retmax=1)
        result = Entrez.read(handle)
        handle.close()
        
        if result.get("IdList"):
            print(f"  ✅ NCBI connection working (found {len(result['IdList'])} results)")
            return True
        else:
            print("  ⚠️  NCBI search returned no results (may be rate limited)")
            return True
    except Exception as e:
        print(f"  ⚠️  NCBI connection test failed: {e}")
        print("     This is OK if you haven't configured NCBI_EMAIL yet")
        return True  # Don't fail on this

def test_utils():
    """Test utility functions."""
    print("\n🔧 Testing utility functions...")
    try:
        from utils import parse_gene_input, is_refseq_accession, nucleotide_counts
        
        # Test gene parsing
        test_input = "BRCA1, TP53\nEGFR"
        genes = parse_gene_input(test_input)
        assert len(genes) == 3, f"Gene parsing failed: got {len(genes)} genes instead of 3"
        print("  ✅ Gene parsing works")
        
        # Test accession detection
        assert is_refseq_accession("NM_000546"), "Accession detection failed"
        assert not is_refseq_accession("BRCA1"), "Accession detection failed"
        print("  ✅ Accession detection works")
        
        # Test nucleotide counting
        stats = nucleotide_counts("ATGCATGC")
        assert stats["A"] == 2, "Nucleotide counting failed"
        assert stats["Length"] == 8, "Nucleotide counting failed"
        print("  ✅ Nucleotide counting works")
        
        return True
    except Exception as e:
        print(f"  ❌ Utils error: {e}")
        return False

def main():
    """Run all verification tests."""
    print("=" * 60)
    print("  🧬 Gene Extractor - Setup Verification")
    print("=" * 60)
    
    results = []
    
    results.append(("Imports", test_imports()))
    results.append(("Configuration", test_configuration()))
    results.append(("Database", test_database()))
    results.append(("Logger", test_logger()))
    results.append(("Dependencies", test_dependencies()))
    results.append(("Utils", test_utils()))
    results.append(("NCBI Connection", test_ncbi_connection()))
    
    # Summary
    print("\n" + "=" * 60)
    print("  📊 Verification Summary")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {status} - {name}")
    
    print("\n" + "=" * 60)
    print(f"  Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("  ✅ All tests passed! System is ready.")
        print("\n  🚀 Start the application with: ./run.sh")
        print("     or: streamlit run app.py")
        return 0
    else:
        print("  ⚠️  Some tests failed. Check the output above.")
        print("     You may need to install dependencies or configure .env")
        return 1

if __name__ == "__main__":
    sys.exit(main())
