#!/usr/bin/env bash
set -euo pipefail

# Setup VFDB and REBASE databases for mgraft_v2
# Usage: bash setup_databases.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="${SCRIPT_DIR}/databases"

echo "=== Setting up databases in ${DB_DIR} ==="
mkdir -p "${DB_DIR}/vfdb" "${DB_DIR}/rebase"

# ============ VFDB (setB) ============
cd "${DB_DIR}/vfdb"

VFDB_URL="http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz"
VFDB_FASTA="VFDB_setB_pro.fas"

if [[ -f "${VFDB_FASTA}" ]]; then
    echo "VFDB FASTA already exists, skipping download"
else
    echo "Downloading VFDB setB..."
    wget -q "${VFDB_URL}" -O "${VFDB_FASTA}.gz"
    gunzip -f "${VFDB_FASTA}.gz"
    echo "Downloaded: ${VFDB_FASTA}"
fi

# Build DIAMOND database
if [[ -f "vfdb.dmnd" ]]; then
    echo "VFDB DIAMOND database already exists"
else
    echo "Building DIAMOND database for VFDB..."
    diamond makedb --in "${VFDB_FASTA}" -d vfdb --quiet
    echo "Built: vfdb.dmnd"
fi

echo "VFDB setup complete: $(grep -c '^>' "${VFDB_FASTA}") proteins"

# ============ REBASE ============
echo ""
cd "${DB_DIR}/rebase"

# REBASE protein sequences with organism info and recseq in headers
REBASE_URL="ftp://ftp.neb.com/pub/rebase/protein_org_seqs.txt"
REBASE_RAW="protein_org_seqs.txt"
REBASE_FILTERED="rebase_with_recseq.fasta"

if [[ -f "${REBASE_RAW}" ]]; then
    echo "REBASE file already exists, skipping download"
else
    echo "Downloading REBASE protein sequences..."
    wget -q "${REBASE_URL}" -O "${REBASE_RAW}"
    echo "Downloaded: ${REBASE_RAW}"
fi

# Filter REBASE to only include enzymes with annotated recognition sequences
# protein_org_seqs.txt has key:value pairs in header, e.g. RecSeq:GATC
echo "Filtering REBASE to enzymes with recognition sequences..."
python3 << 'PYTHON_SCRIPT'
import re
from pathlib import Path

raw_file = Path("protein_org_seqs.txt")
filtered_fasta = Path("rebase_with_recseq.fasta")

kept = 0
total = 0
with open(raw_file) as fin, open(filtered_fasta, 'w') as fout:
    current_header = None
    current_seq = []
    has_recseq = False
    
    for line in fin:
        line = line.rstrip('\n\r')
        if not line:
            continue
        if line.startswith('>'):
            # Write previous record if it had recseq
            if current_header and has_recseq and current_seq:
                fout.write(current_header + '\n')
                # Join sequence lines without spaces
                seq = ''.join(s.replace(' ', '') for s in current_seq)
                fout.write(seq + '\n')
                kept += 1
            
            # Parse header - check for RecSeq: field with actual value
            total += 1
            # Look for RecSeq:\s*(\S+) pattern
            match = re.search(r'RecSeq:\s*([A-Za-z]+)', line)
            has_recseq = match is not None
            current_header = line
            current_seq = []
        else:
            # Sequence line - remove record separator markers
            if line.strip() == '<>':
                continue
            # Remove <> from end of sequence lines
            line = line.rstrip('<>')
            if line:  # Only add non-empty lines
                current_seq.append(line)
    
    # Write last record
    if current_header and has_recseq and current_seq:
        fout.write(current_header + '\n')
        seq = ''.join(s.replace(' ', '') for s in current_seq)
        fout.write(seq + '\n')
        kept += 1

print(f"Filtered: {kept}/{total} sequences have recognition sequences")
PYTHON_SCRIPT

# Build DIAMOND database
if [[ -f "rebase.dmnd" ]]; then
    echo "REBASE DIAMOND database already exists"
else
    echo "Building DIAMOND database for REBASE..."
    diamond makedb --in "${REBASE_FILTERED}" -d rebase --quiet
    echo "Built: rebase.dmnd"
fi

echo "REBASE setup complete: $(grep -c '^>' "${REBASE_FILTERED}") proteins with recognition sequences"

# ============ Summary ============
echo ""
echo "=== Database setup complete ==="
echo "VFDB:   ${DB_DIR}/vfdb/vfdb.dmnd"
echo "REBASE: ${DB_DIR}/rebase/rebase.dmnd"
echo ""
echo "Update your config.yaml to use these paths:"
echo "  vfdb_db: \"workflow/resources/databases/vfdb/vfdb.dmnd\""
echo "  rebase_db: \"workflow/resources/databases/rebase/rebase.dmnd\""
