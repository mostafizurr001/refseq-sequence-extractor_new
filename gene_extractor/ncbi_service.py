"""NCBI API service with rate limiting and caching."""
import time
import io
from typing import List, Dict, Optional, Tuple
from Bio import Entrez, SeqIO
from config import Config
from logger import ncbi_logger
from database import db
from utils import generate_cache_key, nucleotide_counts

class NCBIService:
    """Handle NCBI Entrez API interactions with rate limiting."""
    
    def __init__(self):
        Entrez.email = Config.NCBI_EMAIL
        if Config.NCBI_API_KEY:
            Entrez.api_key = Config.NCBI_API_KEY
        self.last_request_time = 0
        self.request_interval = 1.0 / Config.NCBI_RATE_LIMIT
    
    def _rate_limit(self):
        """Enforce rate limiting."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.request_interval:
            time.sleep(self.request_interval - elapsed)
        self.last_request_time = time.time()
    
    def _retry_request(self, func, *args, **kwargs):
        """Retry failed requests."""
        for attempt in range(Config.NCBI_MAX_RETRIES):
            try:
                self._rate_limit()
                return func(*args, **kwargs)
            except Exception as e:
                ncbi_logger.warning(f"Request failed (attempt {attempt + 1}/{Config.NCBI_MAX_RETRIES}): {e}")
                if attempt == Config.NCBI_MAX_RETRIES - 1:
                    raise
                time.sleep(2 ** attempt)  # Exponential backoff
    
    def search_ids_for_gene(self, gene: str, organism: str, db: str, 
                           prefer_variant1: bool = True, 
                           exclude_variant_x: bool = True) -> List[str]:
        """Search for NCBI IDs matching gene and organism."""
        ncbi_logger.info(f"Searching for gene: {gene}, organism: {organism}, db: {db}")
        
        try:
            # Build base query
            base = f"{gene}[Gene] AND {organism}[Organism] AND RefSeq[Filter]"
            
            if db == "nuccore":
                base += " AND biomol_mrna[PROP]"
                if prefer_variant1:
                    base += ' AND "transcript variant 1"[Title]'
                if exclude_variant_x:
                    base += ' NOT "transcript variant X"[Title]'
            
            # Primary search
            def search():
                handle = Entrez.esearch(db=db, term=base, retmax=20)
                rec = Entrez.read(handle)
                handle.close()
                return rec.get("IdList", [])
            
            ids = self._retry_request(search)
            
            # Fallback searches if no results
            if not ids and prefer_variant1 and db == "nuccore":
                ncbi_logger.info(f"No variant 1 found, searching without variant filter")
                base_fallback = f"{gene}[Gene] AND {organism}[Organism] AND RefSeq[Filter] AND biomol_mrna[PROP]"
                
                def search_fallback():
                    handle = Entrez.esearch(db=db, term=base_fallback, retmax=20)
                    rec = Entrez.read(handle)
                    handle.close()
                    return rec.get("IdList", [])
                
                ids = self._retry_request(search_fallback)
            
            if not ids:
                ncbi_logger.info(f"No results with gene filter, trying broader search")
                broader_query = f"{gene}[All Fields] AND {organism}[Organism] AND RefSeq[Filter]"
                
                def search_broad():
                    handle = Entrez.esearch(db=db, term=broader_query, retmax=20)
                    rec = Entrez.read(handle)
                    handle.close()
                    return rec.get("IdList", [])
                
                ids = self._retry_request(search_broad)
            
            # Choose best variant if nuccore
            if ids and db == "nuccore":
                chosen = self._choose_variant1_id(ids, db, exclude_variant_x)
                return [chosen] if chosen else ids[:1]
            
            return ids[:1] if ids else []
            
        except Exception as e:
            ncbi_logger.error(f"Error searching for gene {gene}: {e}")
            raise
    
    def _choose_variant1_id(self, id_list: List[str], db: str, exclude_variant_x: bool) -> Optional[str]:
        """Choose the best variant ID from a list."""
        if not id_list:
            return None
        
        try:
            def get_summaries():
                handle = Entrez.esummary(db=db, id=",".join(id_list))
                summaries = Entrez.read(handle)
                handle.close()
                return summaries
            
            summaries = self._retry_request(get_summaries)
            
            # First, look for transcript variant 1
            for s in summaries:
                title = s.get("Title", "") or s.get("Caption", "")
                uid = s.get("Id") or s.get("Uid")
                title_l = title.lower()
                
                if "transcript variant 1" in title_l:
                    if exclude_variant_x and "transcript variant x" in title_l:
                        continue
                    ncbi_logger.info(f"Selected variant 1: {uid}")
                    return uid
            
            # If no variant 1, exclude variant X if requested
            if exclude_variant_x:
                for s in summaries:
                    title = s.get("Title", "") or s.get("Caption", "")
                    uid = s.get("Id") or s.get("Uid")
                    if "transcript variant x" not in (title or "").lower():
                        ncbi_logger.info(f"Selected non-X variant: {uid}")
                        return uid
            
            # Default to first
            default_id = summaries[0].get("Id") or summaries[0].get("Uid")
            ncbi_logger.info(f"Using default (first) ID: {default_id}")
            return default_id
            
        except Exception as e:
            ncbi_logger.warning(f"Error choosing variant: {e}")
            return id_list[0]
    
    def search_by_accession(self, accession: str, db: str) -> List[str]:
        """Search for sequence by RefSeq accession."""
        ncbi_logger.info(f"Searching by accession: {accession}, db: {db}")
        
        try:
            search_terms = [
                f"{accession}[Accession]",
                f"{accession}[All Fields] AND RefSeq[Filter]",
                accession
            ]
            
            for query in search_terms:
                def search():
                    handle = Entrez.esearch(db=db, term=query, retmax=1)
                    rec = Entrez.read(handle)
                    handle.close()
                    return rec.get("IdList", [])
                
                ids = self._retry_request(search)
                if ids:
                    ncbi_logger.info(f"Found accession with query: {query}")
                    return ids
            
            return []
            
        except Exception as e:
            ncbi_logger.error(f"Error searching accession {accession}: {e}")
            raise
    
    def fetch_sequence(self, rec_id: str, db: str, rettype: str) -> str:
        """Fetch sequence data from NCBI."""
        ncbi_logger.info(f"Fetching sequence: {rec_id}, db: {db}, type: {rettype}")
        
        try:
            def fetch():
                handle = Entrez.efetch(db=db, id=rec_id, rettype=rettype.lower(), retmode="text")
                data = handle.read()
                handle.close()
                return data
            
            return self._retry_request(fetch)
            
        except Exception as e:
            ncbi_logger.error(f"Error fetching sequence {rec_id}: {e}")
            raise
    
    def extract_utr3_from_genbank(self, gb_text: str) -> Tuple[str, Optional[int], Optional[int], str, object]:
        """Extract 3'UTR from GenBank format."""
        try:
            record_io = io.StringIO(gb_text)
            gb = SeqIO.read(record_io, "genbank")
            full_seq = str(gb.seq)
            
            cds_start = None
            cds_end = None
            
            # Find CDS feature
            for feat in gb.features:
                if feat.type == "CDS":
                    cds_start = int(feat.location.start)
                    cds_end = int(feat.location.end)
                    break
            
            utr3 = ""
            if cds_end is not None and cds_end < len(gb.seq):
                utr3 = str(gb.seq[cds_end:])
                ncbi_logger.info(f"Extracted 3'UTR: {len(utr3)} bp")
            
            return utr3, cds_start, cds_end, full_seq, gb
            
        except Exception as e:
            ncbi_logger.error(f"Error extracting 3'UTR: {e}")
            raise
    
    def get_sequence_with_cache(self, gene: str, organism: str, sequence_type: str,
                                output_format: str, prefer_variant1: bool = True,
                                exclude_variant_x: bool = True, is_accession: bool = False) -> Dict:
        """Get sequence with caching support."""
        # Generate cache key
        cache_key = generate_cache_key(gene, organism, sequence_type, output_format)
        
        # Try to get from cache
        cached = db.get_cached_result(cache_key)
        if cached:
            ncbi_logger.info(f"Using cached result for {gene}")
            db.track_event("cache_hit", {"gene": gene, "organism": organism})
            return {
                "success": True,
                "gene": gene,
                "accession": cached["accession"],
                "sequence_data": cached["sequence_data"],
                "from_cache": True
            }
        
        # Fetch from NCBI
        db.track_event("ncbi_fetch", {"gene": gene, "organism": organism})
        
        try:
            db_type = "nuccore" if sequence_type == "Nucleotide" else "protein"
            
            # Search for IDs
            if is_accession:
                ids = self.search_by_accession(gene, db_type)
            else:
                ids = self.search_ids_for_gene(gene, organism, db_type, prefer_variant1, exclude_variant_x)
            
            if not ids:
                return {
                    "success": False,
                    "gene": gene,
                    "error": "No sequences found"
                }
            
            rec_id = ids[0]
            sequence_data = self.fetch_sequence(rec_id, db_type, output_format)
            
            # Save to cache
            metadata = {"db_type": db_type, "rec_id": rec_id}
            db.save_to_cache(cache_key, gene, organism, sequence_type, output_format, rec_id, sequence_data, metadata)
            
            return {
                "success": True,
                "gene": gene,
                "accession": rec_id,
                "sequence_data": sequence_data,
                "from_cache": False
            }
            
        except Exception as e:
            ncbi_logger.error(f"Error fetching sequence for {gene}: {e}")
            return {
                "success": False,
                "gene": gene,
                "error": str(e)
            }

# Global NCBI service instance
ncbi_service = NCBIService()
