# GRAFT - Genomic Read Analysis and Filtering Tool

[![Python Version](https://img.shields.io/badge/python-3.12%2B-blue.svg)](https://www.python.org/)
[![GitHub last commit](https://img.shields.io/github/last-commit/michi-sxc/GRAFT.svg)](https://github.com/michi-sxc/graft2/commits/main)
[![GitHub issues](https://img.shields.io/github/issues/michi-sxc/GRAFT.svg)](https://github.com/michi-sxc/GRAFT/issues)
[![License](https://img.shields.io/github/license/michi-sxc/GRAFT.svg)](https://github.com/michi-sxc/GRAFT/blob/main/LICENSE)
[![Maintenance](https://img.shields.io/badge/maintenance-active-green.svg)](https://github.com/michi-sxc/GRAFT)


## Overview

GRAFT is a Python-based toolkit for analyzing ancient DNA data through an interactive web interface. The application integrates sequence data visualization, quality assessment, and filtering capabilities specific to degraded DNA analysis. Built on Dash and Plotly, GRAFT processes BAM, SAM, FASTA, and FASTQ files.
## Core Functions

* Interactive sequence data visualization
* Read quality and length distribution analysis
* C>T transition pattern assessment
* MAPQ-based filtering system
* Cross-format file conversion
* Post-mortem damage modeling

## Technical Requirements

* Python 3.9+
* Redis server
* Core dependencies:
  ```
  dash==2.17.1
  dash-bootstrap-components==1.6.0
  dash_dangerously_set_inner_html==0.0.2
  dash-iconify==0.1.2
  pysam==0.22.1
  biopython==1.81
  numpy==1.26.4
  plotly==5.22.0
  pandas==2.2.1
  celery==5.4.0
  PyYAML==6.0.2
  ansi2html==1.9.2
  ```

## Installation

### Unix-Based Systems (Linux/MacOS)

1. Clone repository:
   ```bash
   git clone https://github.com/michi-sxc/GRAFT.git
   cd GRAFT
   ```

2. Set directory permissions:
   ```bash
   chmod 755 .
   chmod -R u+w .
   ```

3. Install requirements:
   ```bash
   pip install -r requirements.txt
   ```

4. Run index.py:
   ```bash
   python index.py
   ```

### Windows (WSL)

1. Install WSL:
   ```powershell
   wsl --install
   ```

2. Execute WSL setup:
   ```bash
   git clone https://github.com/michi-sxc/GRAFT.git
   cd GRAFT
   chmod 755 .
   chmod -R u+w .
   python setup.py
   ```

Note: If running under WSL on a Windows filesystem (NTFS), the setup script will automatically relocate the project to the WSL filesystem to ensure proper permissions handling.

### File Permissions Structure

The setup creates the following permission structure:
* Directory permissions: 755 (drwxr-xr-x)
* Python files: 644 (-rw-r--r--)
* Executable scripts: 755 (-rwxr-xr-x)
* Config files: 644 (-rw-r--r--)

## Usage Instructions

### Data Analysis

1. Upload sequence files through the web interface
2. Select analysis parameters:
   * MAPQ thresholds (0-255)
   * Read length constraints
   * Mismatch tolerances
   * C>T transition filters

3. Export filtered data in BAM/SAM/FASTA/FASTQ format

### Sequence Visualization

* Read length distribution plots
* GC content analysis
* Damage pattern visualization
* Quality score distribution
* Mismatch frequency assessment

### File Management

* Format conversion between BAM/SAM/FASTA/FASTQ
* Batch processing capabilities
* Automatic quality control

## Temporary Files

GRAFT creates temporary files during operation in:
* Linux/WSL: `/tmp/GRAFT_temp/`
* MacOS: `/private/tmp/GRAFT_temp/`

Ensure appropriate write permissions for these directories.

## Performance Considerations

* Memory usage scales with input file size
* Recommended minimum: 8GB RAM
* Processing speed dependent on CPU cores

## Issue Reporting

Submit issues via GitHub with:
* Input file characteristics
* Error messages
* System specifications
* Steps to reproduce

## License

MIT License

## Contact

[michael.schneider6@mailbox.tu-dresden.de](mailto:michael.schneider6@mailbox.tu-dresden.de)
