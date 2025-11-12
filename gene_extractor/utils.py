"""Utility functions for the application."""
import re
import hashlib
from typing import List, Dict, Tuple
from datetime import datetime
import pandas as pd

def is_refseq_accession(token: str) -> bool:
    """Check if token is a RefSeq accession."""
    return bool(re.match(r"^[A-Z]{2,3}_\d+(\.\d+)?$", token))

def parse_gene_input(gene_input: str) -> List[str]:
    """Parse gene input from text area."""
    genes = [g.strip() for g in gene_input.replace(",", "\n").splitlines() if g.strip()]
    return list(dict.fromkeys(genes))  # Remove duplicates while preserving order

def nucleotide_counts(seq: str) -> Dict[str, float]:
    """Calculate nucleotide composition."""
    seq = seq.upper()
    length = len(seq)
    
    if length == 0:
        return {
            "A": 0, "T": 0, "G": 0, "C": 0,
            "Length": 0, "GC%": 0.0,
            "AT%": 0.0, "N_count": 0
        }
    
    a_count = seq.count("A")
    t_count = seq.count("T")
    g_count = seq.count("G")
    c_count = seq.count("C")
    n_count = seq.count("N")
    
    gc_percent = ((g_count + c_count) / length * 100) if length > 0 else 0.0
    at_percent = ((a_count + t_count) / length * 100) if length > 0 else 0.0
    
    return {
        "A": a_count,
        "T": t_count,
        "G": g_count,
        "C": c_count,
        "Length": length,
        "GC%": round(gc_percent, 2),
        "AT%": round(at_percent, 2),
        "N_count": n_count
    }

def generate_cache_key(gene: str, organism: str, sequence_type: str, output_format: str) -> str:
    """Generate cache key for results."""
    key_string = f"{gene}_{organism}_{sequence_type}_{output_format}".lower()
    return hashlib.md5(key_string.encode()).hexdigest()

def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} TB"

def format_timestamp(timestamp: str) -> str:
    """Format timestamp for display."""
    try:
        dt = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        return dt.strftime('%b %d, %Y %I:%M %p')
    except:
        return timestamp

def validate_email(email: str) -> bool:
    """Validate email format."""
    pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return bool(re.match(pattern, email))

def parse_batch_file(file_content: bytes, file_type: str) -> Tuple[List[str], str]:
    """Parse uploaded batch file (CSV or Excel)."""
    try:
        if file_type == "csv":
            df = pd.read_csv(pd.io.common.BytesIO(file_content))
        else:  # Excel
            df = pd.read_excel(pd.io.common.BytesIO(file_content))
        
        # Look for gene column (case-insensitive)
        gene_col = None
        for col in df.columns:
            if col.lower() in ['gene', 'genes', 'gene_symbol', 'gene_name', 'symbol']:
                gene_col = col
                break
        
        if gene_col is None:
            # If no gene column found, use first column
            gene_col = df.columns[0]
        
        genes = df[gene_col].dropna().astype(str).tolist()
        genes = [g.strip() for g in genes if g.strip()]
        
        return genes, None
    except Exception as e:
        return [], f"Error parsing file: {str(e)}"

def create_fasta_format(header: str, sequence: str, line_width: int = 70) -> str:
    """Format sequence in FASTA format."""
    wrapped = "\n".join([sequence[i:i+line_width] for i in range(0, len(sequence), line_width)])
    return f">{header}\n{wrapped}\n"

def get_organism_short_name(organism_full: str) -> str:
    """Extract short organism name."""
    if "(" in organism_full:
        return organism_full.split("(")[0].strip()
    return organism_full

def calculate_sequence_stats(sequences: List[Dict]) -> Dict:
    """Calculate aggregate statistics for multiple sequences."""
    if not sequences:
        return {}
    
    total_length = sum(s.get('Length', 0) for s in sequences)
    avg_gc = sum(s.get('GC%', 0) for s in sequences) / len(sequences) if sequences else 0
    
    return {
        "total_sequences": len(sequences),
        "total_length": total_length,
        "avg_length": total_length // len(sequences) if sequences else 0,
        "avg_gc_percent": round(avg_gc, 2),
        "min_length": min(s.get('Length', 0) for s in sequences),
        "max_length": max(s.get('Length', 0) for s in sequences)
    }

def sanitize_filename(filename: str) -> str:
    """Sanitize filename for safe file system usage."""
    # Remove or replace invalid characters
    filename = re.sub(r'[<>:"/\\|?*]', '_', filename)
    # Remove leading/trailing spaces and dots
    filename = filename.strip('. ')
    # Limit length
    if len(filename) > 255:
        filename = filename[:255]
    return filename
