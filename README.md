# BED File Generator - Web Application

A FastAPI-based web application for generating BED files from gene position data. This tool processes gene coordinates from UCSC Table Browser or similar sources, merges overlapping gene regions, and optionally extends genomic regions.

## Features

- 📤 **Upload gene position files** (TSV/TXT format with chr, start, end, gene columns)
- 🧹 **Automatic filtering** of alternative and fix chromosomes
- 🔄 **Gene position merging** for multiple entries of the same gene
- ↔️ **Region extension** with customizable base pair length
- 📊 **Natural chromosome sorting** (chr1, chr2, ..., chr10, chr11, etc.)
- 💾 **Direct BED file download** with custom filename

## Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Setup

1. **Clone or navigate to the project directory**:
   ```bash
   cd /Users/farukduzenli/Documents/Projeler/rnaseq-batcheffect
   ```

2. **Create a virtual environment** (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On macOS/Linux
   # or
   venv\Scripts\activate  # On Windows
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Starting the Server

Run the FastAPI application:
```bash
python app.py
```

Or using uvicorn directly:
```bash
uvicorn app:app --reload --host 0.0.0.0 --port 8000
```

The application will start on: **http://localhost:8000**

### Using the Web Interface

1. Open your browser and navigate to `http://localhost:8000`
2. Upload your gene position file (TSV/TXT format)
3. Configure options:
   - **Extend regions**: Number of base pairs to add to each end (default: 1000)
   - **Output filename**: Name for the generated BED file
   - **Remove alt/fix chromosomes**: Filter out alternative assemblies
4. Click "Generate BED File" to process and download

### Input File Format

Your input file should be tab-separated with the following columns:
```
chr     start      end        gene
chr1    1000000    1005000    GENE1
chr1    2000000    2010000    GENE2
chr2    500000     505000     GENE3
```

Example sources:
- [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
- GENCODE/Ensembl gene annotations
- Custom gene coordinate files

### API Endpoints

#### `GET /`
Returns the web interface HTML page

#### `POST /process`
Process uploaded gene file and return BED file

**Form parameters:**
- `file`: Gene position file (multipart/form-data)
- `extend_bp`: Base pairs to extend (integer, default: 1000)
- `output_filename`: Output file name (string, default: "merged_genes.bed")
- `remove_alt`: Remove alt chromosomes (boolean, default: true)
- `remove_fix`: Remove fix chromosomes (boolean, default: true)

**Response:** BED file download

#### `GET /health`
Health check endpoint - returns `{"status": "healthy"}`

## How It Works

The application performs the following steps:

1. **Read input file** with gene coordinates (chr, start, end, gene)
2. **Filter contigs** - Remove chromosomes containing 'alt' or 'fix'
3. **Fix reverse strand genes** - Ensure start < end for all genes
4. **Merge gene positions** - Combine multiple entries for the same gene (min start, max end)
5. **Extend regions** - Add specified base pairs to both ends
6. **Sort naturally** - Order by chromosome (chr1, chr2, ..., chr10) and start position
7. **Output BED format** - chr, start, end, gene (tab-separated)

## Development

### Project Structure
```
rnaseq-batcheffect/
├── app.py                                   # FastAPI application
├── requirements.txt                          # Python dependencies
├── create_bed_file_merge_gene_position.ipynb # Original Jupyter notebook
└── README.md                                 # This file
```

### Running in Development Mode
```bash
uvicorn app:app --reload --host 127.0.0.1 --port 8000
```

The `--reload` flag enables auto-restart on code changes.

## Example Use Cases

### Parkinson's Disease Gene Analysis
```python
# Upload gene list from UCSC with Parkinson's disease genes
# Extend regions by 1kb for regulatory element analysis
# Output: parkinson_genes_plus1kb.bed
```

### Custom Genomic Region Analysis
```python
# Input: GENCODE gene annotations
# Extend: 5000 bp for promoter regions
# Filter: Remove alt/fix chromosomes for cleaner analysis
```

## Troubleshooting

### File Upload Issues
- Ensure file is tab-separated
- Check that column order is: chr, start, end, gene
- Verify file encoding is UTF-8

### Processing Errors
- Start positions should be less than end positions (or will be auto-corrected)
- Gene names should not contain special characters causing parsing issues
- Chromosome names should be standard (chr1, chr2, etc.)

## License

This project is open source and available for academic and research purposes.

## Credits

Based on the original Jupyter notebook for RNA-seq batch effect analysis and gene position processing.
