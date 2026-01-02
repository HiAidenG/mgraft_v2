#!/usr/bin/env bash
set -euo pipefail

# Install CRISPRDetect 3.0 for mgraft_v2
# Usage: bash setup_crisprdetect.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="${SCRIPT_DIR}/crisprdetect"

echo "=== Installing CRISPRDetect 3.0 to ${INSTALL_DIR} ==="

if [[ -d "${INSTALL_DIR}" ]]; then
    echo "CRISPRDetect directory exists, checking installation..."
    if [[ -f "${INSTALL_DIR}/CRISPRDetect3" ]]; then
        echo "CRISPRDetect already installed at ${INSTALL_DIR}/CRISPRDetect3"
        exit 0
    fi
fi

mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

# Clone CRISPRDetect 3
echo "Cloning CRISPRDetect 3..."
if [[ -d "CRISPRDetect3" ]]; then
    echo "Repository already exists, pulling latest..."
    cd CRISPRDetect3
    git pull
    cd ..
else
    git clone https://github.com/ambarishbiswas/CRISPRDetect_3.0
fi

# Create wrapper script in the install dir
echo "Creating wrapper script..."
cat > "${INSTALL_DIR}/CRISPRDetect3" << 'WRAPPER'
#!/usr/bin/env bash
# CRISPRDetect3 wrapper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
perl "${SCRIPT_DIR}/CRISPRDetect3/CRISPRDetect.pl" "$@"
WRAPPER
chmod +x "${INSTALL_DIR}/CRISPRDetect3"

# Create symlink for Perl modules
cd "${INSTALL_DIR}/CRISPRDetect_3.0"

# Check dependencies
echo ""
echo "=== Checking dependencies ==="
echo -n "Perl: "
perl --version | head -1

echo -n "BioPerl: "
if perl -MBio::Seq -e 'print "installed\n"' 2>/dev/null; then
    echo "OK"
else
    echo "MISSING - install with: conda install -c bioconda perl-bioperl"
fi

echo -n "cd-hit: "
if command -v cd-hit &> /dev/null; then
    echo "$(cd-hit 2>&1 | head -1)"
else
    echo "MISSING - install with: conda install -c bioconda cd-hit"
fi

echo ""
echo "=== CRISPRDetect installation complete ==="
echo "Path: ${INSTALL_DIR}/CRISPRDetect3"
echo ""
echo "Test with:"
echo "  ${INSTALL_DIR}/CRISPRDetect3 -f test.fasta -o test_output"
