# Benchmark Tool Installation Guide

This document describes how to install all tools required for the LLmap comparative benchmark campaign (Phase 11).

## Quick Start

```bash
# Check if tools are already installed at correct versions
./check_versions.sh

# Install missing tools
./install_tools.sh
```

## Tools Required

| Tool | Version | Purpose |
|------|---------|---------|
| LLmap | 1.0.0 | This work |
| minimap2 | 2.28 | Long + short read alignment (gold standard) |
| BWA-MEM2 | 2.2.1 | Short-read alignment (gold standard) |
| Winnowmap2 | 2.03 | Long-read repeat-aware alignment |
| STAR | 2.7.11b | RNA-seq spliced alignment |
| Bowtie2 | 2.5.4 | Short-read BWT-based alignment |
| samtools | 1.21 | BAM/SAM processing |
| seqkit | 2.8.2 | FASTQ manipulation |
| htslib | 1.21 | BAM library (dependency) |

All version pins are defined in [`datasets/tools.yaml`](../datasets/tools.yaml).

## Installation Methods

### Method 1: HPC Module System (Hummel-2)

On Hummel-2 and similar HPC clusters with a module system:

```bash
# Load all modules
module load minimap2/2.28
module load bwa-mem2/2.2.1
module load winnowmap/2.03
module load star/2.7.11b
module load bowtie2/2.5.4
module load samtools/1.21
module load htslib/1.21

# Or use the install script
./install_tools.sh --module-only
```

Add these to your `~/.bashrc` or SLURM job script for persistence.

### Method 2: Conda/Mamba

For local machines or when modules are unavailable:

```bash
# Create and populate conda environment
./install_tools.sh --conda-only --env llmap-bench

# Activate
conda activate llmap-bench
```

Manual installation:

```bash
# Create environment
conda create -n llmap-bench python=3.11 -y
conda activate llmap-bench

# Add bioconda channel
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install tools (mamba is faster if available)
mamba install -y \
  minimap2=2.28 \
  bwa-mem2=2.2.1 \
  winnowmap=2.03 \
  star=2.7.11b \
  bowtie2=2.5.4 \
  samtools=1.21 \
  htslib=1.21 \
  seqkit=2.8.2
```

### Method 3: Manual Installation

If neither module nor conda works:

#### minimap2

```bash
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -xjf -
export PATH="$PWD/minimap2-2.28_x64-linux:$PATH"
```

#### BWA-MEM2

```bash
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar -xjf -
export PATH="$PWD/bwa-mem2-2.2.1_x64-linux:$PATH"
```

#### Winnowmap2

```bash
git clone https://github.com/marbl/Winnowmap.git
cd Winnowmap && git checkout v2.03
make -j$(nproc)
export PATH="$PWD/bin:$PATH"
```

#### STAR

```bash
curl -L https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip -o STAR.zip
unzip STAR.zip
export PATH="$PWD/STAR-2.7.11b/Linux_x86_64_static:$PATH"
```

#### Bowtie2

```bash
curl -L https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip -o bowtie2.zip
unzip bowtie2.zip
export PATH="$PWD/bowtie2-2.5.4-linux-x86_64:$PATH"
```

#### samtools / htslib

```bash
# htslib
curl -L https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 | tar -xjf -
cd htslib-1.21 && ./configure --prefix=$HOME/.local && make -j$(nproc) && make install
cd ..

# samtools
curl -L https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 | tar -xjf -
cd samtools-1.21 && ./configure --prefix=$HOME/.local && make -j$(nproc) && make install
cd ..
```

## LLmap Installation

LLmap is built from source:

```bash
# Clone (if not already)
git clone https://github.com/schlein-lab/LLmap.git
cd LLmap

# Build
mkdir -p build && cd build
cmake -DLLMAP_ENABLE_CUDA=ON ..   # or CUDA=OFF for CPU-only
cmake --build . -j$(nproc)

# Test
ctest --output-on-failure

# Add to PATH
export LLMAP_ROOT="$PWD/.."
export LLMAP_BIN="$PWD/llmap"
export PATH="$PWD:$PATH"
```

For benchmark runners, set `LLMAP_BIN` to point to your built binary:

```bash
export LLMAP_BIN=/path/to/llmap/build/llmap
```

## Verification

After installation, verify all tools:

```bash
./check_versions.sh
```

Expected output when all tools are correct:

```
Checking benchmark tool versions...

  [OK] minimap2     2.28
  [OK] bwa_mem2     2.2.1
  [OK] winnowmap    2.03
  [OK] star         2.7.11b
  [OK] bowtie2      2.5.4
  [OK] samtools     1.21
  [OK] seqkit       2.8.2
  [OK] llmap        1.0.0

Results written to: ../reports/tool_versions.json

PASS: all tools present at expected versions
```

The script exits with code 0 on success, 1 on any mismatch or missing tool.

## SLURM Job Setup

For SLURM jobs on Hummel-2, add to your job script:

```bash
#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00

# Load modules
module load minimap2/2.28
module load bwa-mem2/2.2.1
module load winnowmap/2.03
module load star/2.7.11b
module load bowtie2/2.5.4
module load samtools/1.21
module load htslib/1.21

# Set LLmap path
export LLMAP_BIN=/beegfs/u/$USER/llmap/build/llmap
export PATH="$(dirname $LLMAP_BIN):$PATH"

# Verify before run
benchmarks/runners/check_versions.sh --json-only

# ... rest of benchmark script
```

## Troubleshooting

### Module not found

Check available versions:
```bash
module avail minimap2
```

Use closest available version and note in results.

### Conda solve fails

Try mamba (faster solver):
```bash
conda install -n base mamba
mamba install -n llmap-bench minimap2=2.28
```

Or install tools one at a time.

### Version mismatch

If exact version unavailable, document the actual version in your run manifest. Minor version differences (e.g., 2.28 vs 2.28-r1199) are acceptable.

### Permission denied

Ensure scripts are executable:
```bash
chmod +x check_versions.sh install_tools.sh
```
