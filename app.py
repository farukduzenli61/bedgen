from fastapi import FastAPI, Form, HTTPException
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import pandas as pd
from natsort import index_natsorted
import numpy as np
import os
import tempfile
import requests
from typing import List, Dict

app = FastAPI(title="BED Dosyası Oluşturucu", description="BED dosyalarını gen pozisyon verilerinden oluşturun")

# Enable CORS for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Standard chromosomes (chr1-22, chrX, chrY)
STANDARD_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def normalize_chromosome(seq_region: str) -> str:
    """Convert Ensembl format (17) to BED format (chr17)"""
    if seq_region.startswith("chr"):
        return seq_region
    return f"chr{seq_region}"


def is_standard_chromosome(chrom: str) -> bool:
    """Check if chromosome is standard (chr1-22, chrX, chrY)"""
    return chrom in STANDARD_CHROMOSOMES


def fetch_gene_positions(gene_symbols: List[str]) -> Dict:
    """
    Fetch gene positions from Ensembl REST API.

    Uses POST /lookup/symbol/homo_sapiens for batch queries.
    Returns: {"found": [...], "not_found": [...]}
    """
    if not gene_symbols:
        return {"found": [], "not_found": []}

    # Remove duplicates and empty strings
    unique_symbols = list(set(s.strip().upper() for s in gene_symbols if s.strip()))

    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    found = []
    not_found = []

    try:
        response = requests.post(
            url,
            headers=headers,
            json={"symbols": unique_symbols}
        )

        if response.status_code == 200:
            data = response.json()

            for symbol in unique_symbols:
                if symbol in data and data[symbol]:
                    gene_data = data[symbol]
                    chrom = normalize_chromosome(gene_data.get("seq_region_name", ""))
                    start = gene_data.get("start", 0)
                    end = gene_data.get("end", 0)
                    strand = gene_data.get("strand", 1)

                    # Fix reverse strand (ensure start < end)
                    if start > end:
                        start, end = end, start

                    found.append({
                        "gene": gene_data.get("display_name", symbol),
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "is_standard": is_standard_chromosome(chrom)
                    })
                else:
                    not_found.append(symbol)
        else:
            # If batch request fails, try individual lookups
            for symbol in unique_symbols:
                try:
                    single_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
                    single_response = requests.get(single_url, headers=headers)

                    if single_response.status_code == 200:
                        gene_data = single_response.json()
                        chrom = normalize_chromosome(gene_data.get("seq_region_name", ""))
                        start = gene_data.get("start", 0)
                        end = gene_data.get("end", 0)
                        strand = gene_data.get("strand", 1)

                        if start > end:
                            start, end = end, start

                        found.append({
                            "gene": gene_data.get("display_name", symbol),
                            "chr": chrom,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "is_standard": is_standard_chromosome(chrom)
                        })
                    else:
                        not_found.append(symbol)
                except Exception:
                    not_found.append(symbol)

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Ensembl API error: {str(e)}")

    return {"found": found, "not_found": not_found}


def process_bed_file(df: pd.DataFrame, extend_bp: int = 0, remove_alt: bool = True, remove_fix: bool = True) -> pd.DataFrame:
    """
    Process gene position data to create BED file
    
    Args:
        df: DataFrame with columns chr, start, end, gene
        extend_bp: Number of base pairs to extend on each side
        remove_alt: Remove chromosomes containing 'alt'
        remove_fix: Remove chromosomes containing 'fix'
    """
    # Remove unwanted contigs
    if remove_alt:
        df = df[~df["chr"].str.contains("alt", na=False)]
    if remove_fix:
        df = df[~df["chr"].str.contains("fix", na=False)]
    
    # Fix reverse strand genes (where start > end)
    for idx in df.index:
        if int(df.loc[idx, 'start']) > int(df.loc[idx, 'end']):
            start = df.loc[idx, 'start']
            end = df.loc[idx, 'end']
            df.at[idx, 'start'] = end
            df.at[idx, 'end'] = start
    
    # Merge multiple entries for the same gene
    merged = df.groupby("gene").agg({
        "chr": "first",
        "start": "min",
        "end": "max"
    }).reset_index()
    
    # Extend start and end regions
    if extend_bp > 0:
        merged["start"] = merged["start"].astype(int) - extend_bp
        merged["end"] = merged["end"].astype(int) + extend_bp
        
        # Ensure start is not negative
        merged["start"] = merged["start"].clip(lower=0)
    
    # Reorder columns for BED format
    merged = merged[["chr", "start", "end", "gene"]]
    
    # Sort by chromosome and start position
    merged = merged.sort_values(by="start", ascending=True)
    merged = merged.sort_values(
        by=["chr", "start"],
        key=lambda col: np.argsort(index_natsorted(col)) if col.name == "chr" else col
    )
    
    return merged


@app.get("/", response_class=HTMLResponse)
async def read_root():
    """Serve the main HTML page"""
    return """
    <!DOCTYPE html>
    <html>
    <head>
        <title>BED Dosyası Oluşturucu</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                max-width: 900px;
                margin: 50px auto;
                padding: 20px;
                background-color: #f5f5f5;
            }
            .container {
                background-color: white;
                padding: 30px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            h1 { color: #333; text-align: center; }
            h3 { color: #555; margin-top: 25px; margin-bottom: 10px; }
            .form-group { margin-bottom: 20px; }
            label {
                display: block;
                margin-bottom: 5px;
                font-weight: bold;
                color: #555;
            }
            textarea, input[type="number"], input[type="text"] {
                width: 100%;
                padding: 10px;
                border: 1px solid #ddd;
                border-radius: 4px;
                box-sizing: border-box;
                font-family: monospace;
            }
            textarea { min-height: 120px; resize: vertical; }
            button {
                background-color: #4CAF50;
                color: white;
                padding: 12px 30px;
                border: none;
                border-radius: 4px;
                cursor: pointer;
                font-size: 16px;
                width: 100%;
                margin-top: 10px;
            }
            button:hover { background-color: #45a049; }
            button:disabled { background-color: #ccc; cursor: not-allowed; }
            button.secondary {
                background-color: #2196F3;
            }
            button.secondary:hover { background-color: #1976D2; }
            .info {
                background-color: #e7f3ff;
                padding: 15px;
                border-left: 4px solid #2196F3;
                margin-bottom: 20px;
            }
            .warning {
                background-color: #fff3cd;
                padding: 15px;
                border-left: 4px solid #ffc107;
                margin-top: 15px;
            }
            .error {
                background-color: #f8d7da;
                padding: 15px;
                border-left: 4px solid #dc3545;
                margin-top: 15px;
            }
            .success {
                background-color: #d4edda;
                padding: 15px;
                border-left: 4px solid #28a745;
                margin-top: 15px;
            }
            .hidden { display: none; }
            #loading {
                text-align: center;
                margin-top: 20px;
            }
            .spinner {
                border: 4px solid #f3f3f3;
                border-top: 4px solid #3498db;
                border-radius: 50%;
                width: 40px;
                height: 40px;
                animation: spin 1s linear infinite;
                margin: 0 auto;
            }
            @keyframes spin {
                0% { transform: rotate(0deg); }
                100% { transform: rotate(360deg); }
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin-top: 15px;
                font-size: 14px;
            }
            th, td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            th { background-color: #4CAF50; color: white; }
            tr:nth-child(even) { background-color: #f9f9f9; }
            .checkbox-item {
                display: flex;
                align-items: center;
                padding: 8px;
                border-bottom: 1px solid #eee;
            }
            .checkbox-item input { margin-right: 10px; }
            .checkbox-item label { font-weight: normal; margin: 0; }
            .step {
                border: 2px solid #ddd;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 20px;
            }
            .step.active { border-color: #4CAF50; }
            .step-header {
                display: flex;
                align-items: center;
                margin-bottom: 15px;
            }
            .step-number {
                background-color: #4CAF50;
                color: white;
                width: 30px;
                height: 30px;
                border-radius: 50%;
                display: flex;
                align-items: center;
                justify-content: center;
                font-weight: bold;
                margin-right: 10px;
            }
            .step.disabled .step-number { background-color: #ccc; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>BED Dosyası Oluşturucu</h1>

            <div class="info">
                <strong>Gen sembollerinden BED dosyası oluşturun:</strong>
                <ul>
                    <li>Gen sembollerini virgül veya satır satır girin (örn: BRCA1, TP53, EGFR)</li>
                    <li>Pozisyonlar Ensembl API'den otomatik çekilir (GRCh38)</li>
                    <li>Standart olmayan kromozomları dahil edip etmeyeceğinizi seçebilirsiniz</li>
                </ul>
            </div>

            <!-- Step 1: Gene Input -->
            <div class="step active" id="step1">
                <div class="step-header">
                    <div class="step-number">1</div>
                    <h3 style="margin: 0;">Gen Sembolleri</h3>
                </div>
                <div class="form-group">
                    <label for="genes">Gen sembollerini girin (virgül veya satır satır):</label>
                    <textarea id="genes" placeholder="BRCA1, TP53, EGFR, ACTB&#10;veya&#10;BRCA1&#10;TP53&#10;EGFR&#10;ACTB"></textarea>
                </div>
                <button type="button" class="secondary" onclick="fetchPositions()">Pozisyonları Getir</button>
            </div}
            <div id="loading" class="hidden">
                <div class="spinner"></div>
                <p>Ensembl API'dan pozisyonlar alınıyor...</p>
            </div>

            <div id="errorBox" class="error hidden"></div>
            <div id="warningBox" class="warning hidden"></div>

            <!-- Step 2: Results & Non-standard chromosomes -->
            <div class="step disabled hidden" id="step2">
                <div class="step-header">
                    <div class="step-number">2</div>
                    <h3 style="margin: 0;">Bulunan Genler</h3>
                </div>

                <div id="resultsTable"></div>

                <div id="nonStandardSection" class="hidden">
                    <h3>Standart Olmayan Kromozomlar</h3>
                    <p style="color: #666;">Aşağıdaki genler standart kromozomlar (chr1-22, X, Y) dışındadır. Dahil etmek istediklerinizi seçin:</p>
                    <div id="nonStandardList"></div>
                </div>
            </div>

            <!-- Step 3: Extend Regions -->
            <div class="step disabled hidden" id="step3">
                <div class="step-header">
                    <div class="step-number">3</div>
                    <h3 style="margin: 0;">Bölgeleri Genişlet</h3>
                </div>

                <div class="form-group">
                    <label for="extend_bp">Genişletme miktarı (bp):</label>
                    <input type="number" id="extend_bp" value="1000" min="0" step="100">
                    <small style="color: #666;">Start'tan çıkarılacak ve End'e eklenecek baz çifti sayısı</small>
                </div>

                <button type="button" class="secondary" onclick="previewExtension()">Önizleme</button>

                <div id="previewSection" class="hidden" style="margin-top: 20px;">
                    <h4>Genişletme Karşılaştırması (Rastgele 3-5 gen):</h4>
                    <table id="comparisonTable">
                        <thead>
                            <tr>
                                <th>Gen</th>
                                <th>Kromozom</th>
                                <th>Orijinal Start</th>
                                <th>Yeni Start</th>
                                <th>Orijinal End</th>
                                <th>Yeni End</th>
                            </tr>
                        </thead>
                        <tbody id="comparisonBody"></tbody>
                    </table>
                </div>
            </div>

            <!-- Step 4: Download BED -->
            <div class="step disabled hidden" id="step4">
                <div class="step-header">
                    <div class="step-number">4</div>
                    <h3 style="margin: 0;">BED Dosyasını İndir</h3>
                </div>

                <div class="form-group">
                    <label for="output_filename">Çıktı dosya adı:</label>
                    <input type="text" id="output_filename" value="genes.bed">
                </div>

                <button type="button" onclick="generateBed()">BED Dosyasını İndir</button>
            </div>

            <div id="successBox" class="success hidden"></div>
        </div>

        <script>
            let foundGenes = [];
            let nonStandardGenes = [];

            async function fetchPositions() {
                const genes = document.getElementById('genes').value.trim();
                if (!genes) {
                    showError('Lütfen en az bir gen sembolü girin.');
                    return;
                }

                showLoading(true);
                hideMessages();

                try {
                    const formData = new FormData();
                    formData.append('genes', genes);

                    const response = await fetch('/fetch-positions', {
                        method: 'POST',
                        body: formData
                    });

                    if (!response.ok) {
                        const error = await response.json();
                        throw new Error(error.detail || 'API hatasi');
                    }

                    const data = await response.json();
                    foundGenes = data.genes || [];
                    nonStandardGenes = data.non_standard || [];
                    const notFound = data.not_found || [];

                    // Show warnings for not found genes
                    if (notFound.length > 0) {
                        showWarning('Bulunamayan genler: ' + notFound.join(', '));
                    }

                    if (foundGenes.length === 0 && nonStandardGenes.length === 0) {
                        showError('Hiçbir gen bulunamadı. Lütfen gen sembollerini kontrol edin.');
                        return;
                    }

                    // Display results
                    displayResults();

                    // Show step 2 and 3
                    document.getElementById('step2').classList.remove('hidden', 'disabled');
                    document.getElementById('step2').classList.add('active');
                    document.getElementById('step3').classList.remove('hidden', 'disabled');
                    document.getElementById('step3').classList.add('active');

                } catch (error) {
                    showError('Hata: ' + error.message);
                } finally {
                    showLoading(false);
                }
            }

            function displayResults() {
                // Display found genes table
                let html = '<table><thead><tr><th>Gen</th><th>Kromozom</th><th>Başlangıç</th><th>Bitiş</th><th>Strand</th></tr></thead><tbody>';

                foundGenes.forEach(gene => {
                    html += `<tr>
                        <td>${gene.gene}</td>
                        <td>${gene.chr}</td>
                        <td>${gene.start.toLocaleString()}</td>
                        <td>${gene.end.toLocaleString()}</td>
                        <td>${gene.strand === 1 ? '+' : '-'}</td>
                    </tr>`;
                });

                html += '</tbody></table>';
                document.getElementById('resultsTable').innerHTML = html;

                // Display non-standard chromosomes if any
                if (nonStandardGenes.length > 0) {
                    document.getElementById('nonStandardSection').classList.remove('hidden');

                    let nsHtml = '';
                    nonStandardGenes.forEach((gene, idx) => {
                        nsHtml += `<div class="checkbox-item">
                            <input type="checkbox" id="ns_${idx}" value="${idx}">
                            <label for="ns_${idx}">${gene.gene} (${gene.chr}: ${gene.start.toLocaleString()} - ${gene.end.toLocaleString()})</label>
                        </div>`;
                    });
                    document.getElementById('nonStandardList').innerHTML = nsHtml;
                } else {
                    document.getElementById('nonStandardSection').classList.add('hidden');
                }
            }

            function getSelectedGenes() {
                let selectedGenes = [...foundGenes];
                nonStandardGenes.forEach((gene, idx) => {
                    const checkbox = document.getElementById('ns_' + idx);
                    if (checkbox && checkbox.checked) {
                        selectedGenes.push(gene);
                    }
                });
                return selectedGenes;
            }

            function previewExtension() {
                const selectedGenes = getSelectedGenes();
                if (selectedGenes.length === 0) {
                    showError('Onizleme icin gen secilmedi.');
                    return;
                }

                const extendBp = parseInt(document.getElementById('extend_bp').value) || 0;

                // Randomly select 3-5 genes for preview
                const sampleSize = Math.min(Math.floor(Math.random() * 3) + 3, selectedGenes.length);
                const shuffled = [...selectedGenes].sort(() => 0.5 - Math.random());
                const sampleGenes = shuffled.slice(0, sampleSize);

                // Build comparison table
                let html = '';
                sampleGenes.forEach(gene => {
                    const newStart = Math.max(0, gene.start - extendBp);
                    const newEnd = gene.end + extendBp;
                    const startDiff = gene.start - newStart;
                    const endDiff = newEnd - gene.end;

                    html += `<tr>
                        <td><strong>${gene.gene}</strong></td>
                        <td>${gene.chr}</td>
                        <td>${gene.start.toLocaleString()}</td>
                        <td style="color: #28a745;">${newStart.toLocaleString()} <small>(-${startDiff.toLocaleString()})</small></td>
                        <td>${gene.end.toLocaleString()}</td>
                        <td style="color: #28a745;">${newEnd.toLocaleString()} <small>(+${endDiff.toLocaleString()})</small></td>
                    </tr>`;
                });

                document.getElementById('comparisonBody').innerHTML = html;
                document.getElementById('previewSection').classList.remove('hidden');

                // Show step 4
                document.getElementById('step4').classList.remove('hidden', 'disabled');
                document.getElementById('step4').classList.add('active');

                showSuccess('Önizleme hazırlandı. ' + sampleSize + ' gen gösteriliyor. İndirmek için aşağıdaki butonu kullanın.');
            }

            async function generateBed() {
                const selectedGenes = getSelectedGenes();
                if (selectedGenes.length === 0) {
                    showError('Indirilecek gen yok.');
                    return;
                }

                const extendBp = document.getElementById('extend_bp').value;
                const outputFilename = document.getElementById('output_filename').value || 'genes.bed';

                showLoading(true);
                hideMessages();

                try {
                    const formData = new FormData();
                    formData.append('genes_json', JSON.stringify(selectedGenes));
                    formData.append('extend_bp', extendBp);
                    formData.append('output_filename', outputFilename);

                    const response = await fetch('/generate-bed', {
                        method: 'POST',
                        body: formData
                    });

                    if (!response.ok) {
                        const error = await response.json();
                        throw new Error(error.detail || 'Dosya oluşturulamadı');
                    }

                    const blob = await response.blob();
                    const url = window.URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = outputFilename;
                    document.body.appendChild(a);
                    a.click();
                    window.URL.revokeObjectURL(url);
                    document.body.removeChild(a);

                    showSuccess('BED dosyası başarıyla oluşturuldu! (' + selectedGenes.length + ' gen)');

                } catch (error) {
                    showError('Hata: ' + error.message);
                } finally {
                    showLoading(false);
                }
            }

            function showLoading(show) {
                document.getElementById('loading').classList.toggle('hidden', !show);
            }

            function hideMessages() {
                document.getElementById('errorBox').classList.add('hidden');
                document.getElementById('warningBox').classList.add('hidden');
                document.getElementById('successBox').classList.add('hidden');
            }

            function showError(msg) {
                const el = document.getElementById('errorBox');
                el.textContent = msg;
                el.classList.remove('hidden');
            }

            function showWarning(msg) {
                const el = document.getElementById('warningBox');
                el.textContent = msg;
                el.classList.remove('hidden');
            }

            function showSuccess(msg) {
                const el = document.getElementById('successBox');
                el.textContent = msg;
                el.classList.remove('hidden');
            }
        </script>
    </body>
    </html>
    """


@app.post("/fetch-positions")
async def fetch_positions(genes: str = Form(...)):
    """
    Fetch gene positions from Ensembl API.
    Returns found genes, not found genes, and non-standard chromosome genes.
    """
    # Parse gene symbols (comma or newline separated)
    gene_list = []
    for line in genes.replace(",", "\n").split("\n"):
        symbol = line.strip()
        if symbol:
            gene_list.append(symbol)

    if not gene_list:
        raise HTTPException(status_code=400, detail="No gene symbols provided")

    result = fetch_gene_positions(gene_list)

    # Separate standard and non-standard chromosome genes
    standard_genes = [g for g in result["found"] if g["is_standard"]]
    non_standard_genes = [g for g in result["found"] if not g["is_standard"]]

    return JSONResponse({
        "genes": standard_genes,
        "non_standard": non_standard_genes,
        "not_found": result["not_found"]
    })


@app.post("/generate-bed")
async def generate_bed(
    genes_json: str = Form(...),
    extend_bp: int = Form(1000),
    output_filename: str = Form("genes.bed")
):
    """
    Generate BED file from gene position data.
    genes_json: JSON array of gene objects with chr, start, end, gene fields
    """
    import json

    try:
        genes = json.loads(genes_json)
    except json.JSONDecodeError:
        raise HTTPException(status_code=400, detail="Invalid JSON format")

    if not genes:
        raise HTTPException(status_code=400, detail="No genes provided")

    # Create DataFrame from gene list
    df = pd.DataFrame(genes)

    # Ensure required columns exist
    required_cols = ["chr", "start", "end", "gene"]
    for col in required_cols:
        if col not in df.columns:
            raise HTTPException(status_code=400, detail=f"Missing column: {col}")

    # Process the data (extend regions, sort)
    # Skip alt/fix filtering since we already filtered by standard chromosomes
    merged = process_bed_file(df, extend_bp, remove_alt=False, remove_fix=False)

    # Save to temporary output file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as tmp:
        merged.to_csv(tmp.name, sep="\t", index=False, header=False)
        output_path = tmp.name

    return FileResponse(
        output_path,
        media_type="text/plain",
        filename=output_filename,
        background=lambda: os.unlink(output_path)
    )


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy"}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
