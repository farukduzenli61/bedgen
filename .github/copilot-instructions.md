# BED File Generator Web Application

## Project Overview
FastAPI web application for generating BED files from gene position data. Processes UCSC Table Browser or similar gene coordinate files by filtering contigs, merging gene positions, and extending genomic regions.

## Architecture
- **Backend**: FastAPI application ([app.py](app.py)) with single-page HTML frontend embedded
- **Core Logic**: `process_bed_file()` function handles all gene position processing
- **Data Flow**: File upload → pandas processing → BED file download
- No external database; all processing is in-memory with temporary files

## Running the Application

**Start server**:
```bash
python app.py
# Or with auto-reload for development:
uvicorn app:app --reload --host 0.0.0.0 --port 8000
```

**Access**: Navigate to http://localhost:8000

**Dependencies**:
```bash
pip install -r requirements.txt
```

## Key Components

### Input Format
Tab-separated file with columns: `chr`, `start`, `end`, `gene`
- Example source: UCSC Table Browser gene annotations
- Handles reverse strand genes (where start > end)

### Processing Pipeline ([app.py](app.py#L21-L60))
1. Filter alt/fix chromosome contigs
2. Fix reverse strand positions (ensure start < end)
3. Merge multiple entries per gene (min start, max end)
4. Extend regions by specified base pairs
5. Natural sort chromosomes (chr1, chr2, ..., chr10, not chr1, chr10, chr2)

### API Endpoints
- `GET /` - Web interface
- `POST /process` - Main processing endpoint (multipart form-data)
  - Parameters: `file`, `extend_bp`, `output_filename`, `remove_alt`, `remove_fix`
- `GET /health` - Health check

## Code Conventions

**Data Processing Pattern**:
- All genomic data uses pandas DataFrames
- Chromosome sorting with `natsort.index_natsorted()` for natural ordering
- Temporary files cleaned up via FileResponse background tasks

**Error Handling**:
- Input validation at API level with HTTPException
- Try-except blocks around file operations with proper cleanup

**File Handling**:
```python
# Always use temporary files for uploads/downloads
with tempfile.NamedTemporaryFile(...) as tmp:
    # Process
output_path = tempfile.mktemp(suffix='.bed')
```

## Related Files
- [create_bed_file_merge_gene_position.ipynb](create_bed_file_merge_gene_position.ipynb) - Original Jupyter notebook with analysis workflow
- [README.md](README.md) - User documentation and usage examples

## Notes for AI Agents
- Gene position data can have start > end for reverse strand genes - always normalize
- Use `natsort` for chromosome sorting to avoid chr1, chr10, chr2 ordering issues
- BED format is 0-based start, but input files may vary - document assumptions
- When extending regions, ensure start doesn't go negative (use `.clip(lower=0)`)
- FastAPI File uploads require `python-multipart` package
- Frontend is embedded in app.py - no separate static files needed for deployment
