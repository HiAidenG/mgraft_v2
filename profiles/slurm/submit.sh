#!/usr/bin/env bash
# Snakemake SLURM submit script

set -euo pipefail

# Parse Snakemake job properties
JOBSCRIPT="$1"
shift

# Extract job properties
RULE=$(cat "$JOBSCRIPT" | grep -oP '(?<=# properties = ).*' | python3 -c "import sys, json; print(json.load(sys.stdin).get('rule', 'unknown'))")
THREADS=$(cat "$JOBSCRIPT" | grep -oP '(?<=# properties = ).*' | python3 -c "import sys, json; print(json.load(sys.stdin).get('threads', 1))")
RESOURCES=$(cat "$JOBSCRIPT" | grep -oP '(?<=# properties = ).*' | python3 -c "import sys, json; print(json.dumps(json.load(sys.stdin).get('resources', {})))")

# Parse resources
MEM_MB=$(echo "$RESOURCES" | python3 -c "import sys, json; print(json.load(sys.stdin).get('mem_mb', 4000))")
RUNTIME=$(echo "$RESOURCES" | python3 -c "import sys, json; print(json.load(sys.stdin).get('runtime', 60))")

# Convert memory to GB for SLURM
MEM_GB=$(python3 -c "import math; print(math.ceil($MEM_MB / 1024))")

# Submit job
sbatch --job-name="mgraft_v2.$RULE" \
       --cpus-per-task="$THREADS" \
       --mem="${MEM_GB}G" \
       --time="${RUNTIME}" \
       --output="logs/slurm-%x.%j.out" \
       --error="logs/slurm-%x.%j.err" \
       "$JOBSCRIPT"
