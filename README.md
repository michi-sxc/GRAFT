# GRAFT - Genomic Read Analysis and Filtering Tool

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation and Setup](#installation-and-setup)
  - [Mac/Linux Users](#maclinux-users)
  - [Windows Users (Using WSL)](#windows-users-using-wsl)
  - [Opening WSL in VSCode](#opening-wsl-in-vscode)
- [Usage](#usage)
  - [Uploading Files](#uploading-files)
  - [Data Filtering](#data-filtering)
  - [Data Visualization](#data-visualization)
  - [File Conversion](#file-conversion)
  - [Centrifuge Analysis](#centrifuge-analysis)
- [Configuration](#configuration)
- [Folder Structure](#folder-structure)
- [Known Issues](#known-issues)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Overview



GRAFT is an interactive dashboard written in Python for genomic data analysis especially for working with aDNA data, focusing on BAM, SAM, FASTA, and FASTQ files. Built with Dash and Plotly, GRAFT provides functionalities for data upload, visualization, filtering, and batch conversion of genomic read files. It supports exploring read length distributions, CG content, damage patterns, and mismatch frequencies.
## Features

- **Interactive Data Visualization**: Visualize read length distributions, CG content, damage patterns, and mismatch frequencies.
- **Advanced Filtering**: Apply filters on Mapping Quality (MAPQ), read mismatch counts, and C>T changes.
- **File Upload and Conversion**: Support for BAM, SAM, FASTA, and FASTQ files. Batch conversion between these formats.
- **Centrifuge Integration**: Perform taxonomic classification on FASTQ files.

## Requirements

- Python 3.9+
- Dependencies (see `requirements.txt`):
  - dash
  - dash-bootstrap-components
  - dash-iconify
  - pysam
  - biopython
  - numpy
  - plotly
  - pandas
  - celery
  - redis
  - yaml
- **Redis**: Required for Celery task queue.
- **Centrifuge**: For taxonomic classification features only. Dashboard works perfectly fine without Centrifuge being installed.

## Installation and Setup

### Mac/Linux Users

1. **Clone the repository**:

   ```bash
   git clone https://github.com/michi-sxc/GRAFT.git
   cd GRAFT
   ```

2. **Install dependencies**:

   ```bash
   pip install -r requirements.txt
   ```

3. **Install and configure Redis (required for Celery)**:

   - On Ubuntu/Debian:
     ```bash
     sudo apt-get install redis-server
     ```
   - On macOS using Homebrew:
     ```bash
     brew install redis
     ```

4. **Start Redis**:

   ```bash
   redis-server
   ```

5. **Install Centrifuge (required for taxonomic classification)**:

   - **Install Bioconda** (if not already installed):
     ```bash
     conda config --add channels defaults
     conda config --add channels bioconda
     conda config --add channels conda-forge
     ```
   - **Install Centrifuge**:
     ```bash
     conda install -c bioconda centrifuge
     ```
   - Alternatively, you can install Centrifuge from source following the instructions on the [Centrifuge GitHub repository](https://github.com/DaehwanKimLab/centrifuge).

6. **Configure Celery**:

   - Start a Celery worker:
     ```bash
     celery -A app.celery_app worker --loglevel=info
     ```

7. **Running the Application**:

   - To start the dashboard, run:
     ```bash
     python app.py
     ```
     The dashboard will be accessible at [http://127.0.0.1:8050/](http://127.0.0.1:8050/).

### Windows Users (Using WSL)

1. **Install WSL (Windows Subsystem for Linux)**:

   - Open PowerShell as Administrator and run:
     ```powershell
     wsl --install
     ```
   - Restart your computer if prompted.

2. **Set up WSL**:

   - Open a WSL terminal (e.g., Ubuntu).
   - **Clone the repository**:
     ```bash
     git clone https://github.com/michi-sxc/GRAFT.git
     cd GRAFT
     ```

3. **Install dependencies**:

   ```bash
   pip install -r requirements.txt
   ```

4. **Install and configure Redis**:

   - In the WSL terminal:
     ```bash
     sudo apt-get install redis-server
     ```

5. **Start Redis**:

   ```bash
   redis-server
   ```

6. **Install Centrifuge (required for taxonomic classification)**:

   - **Install Bioconda**:
     ```bash
     conda config --add channels defaults
     conda config --add channels bioconda
     conda config --add channels conda-forge
     ```
   - **Install Centrifuge**:
     ```bash
     conda install -c bioconda centrifuge
     ```

7. **Configure Celery**:

   - Start a Celery worker:
     ```bash
     celery -A app.celery_app worker --loglevel=info
     ```

8. **Running the Application**:

   - To start the dashboard, run:
     ```bash
     python app.py
     ```
     The dashboard will be accessible at [http://127.0.0.1:8050/](http://127.0.0.1:8050/).

### Opening WSL in VSCode (optional):

1. **Install VSCode (if not already installed)**:
   - Download and install [Visual Studio Code](https://code.visualstudio.com/).

2. **Install the WSL Extension**:
   - Open VSCode and go to the Extensions view by clicking on the Extensions icon in the Activity Bar on the side of the window.
   - Search for "Remote - WSL" and click Install.

3. **Open the Project in WSL**:
   - Open a WSL terminal and navigate to your project directory:
     ```bash
     cd /path/to/GRAFT
     ```
   - Launch VSCode from the terminal:
     ```bash
     code .
     ```
     This will open the project in VSCode with WSL integration, allowing you to edit and manage your project files seamlessly.

## Usage

### Uploading Files

- Click "Drag and Drop or Click to Select Files" to upload BAM, SAM, FASTA, or FASTQ files.
- Uploaded files are displayed in the sidebar under "File Selection".

### Data Filtering

- **MAPQ Range**: Adjust the slider to select mapping quality for filtering reads.
- **Mismatch Filter**: Filter reads based on the number of mismatches (NM tag).
- **C>T Change Filter**:
  - Show only C>T changes: Display reads with C>T changes.
  - Exclude C>T changed reads: Exclude reads with C>T changes.
  - Subtract C>T changes from NM: Adjust the mismatch count by subtracting C>T changes.

### Data Visualization

- Navigate through different tabs to view:
  - Read Length Distribution
  - CG Content Distribution
  - Damage Patterns
  - Mismatch Frequencies
  - Alignment Statistics

- **Interactive Plots**: Click and drag on the histograms to select specific ranges or data points. The selections will update the other visualizations accordingly (not working with read length distribution).

### File Conversion

- Navigate to the "Convert" tab to batch convert between BAM, SAM, FASTA, and FASTQ formats.
- Upload files you wish to convert.
- Select the desired output format.
- Click "Convert Files" to start the conversion.
- Converted files can be downloaded directly or as a ZIP archive if multiple files are converted.
- Please note that fasta/q to BAM/SAM conversions will generate dummy headers and empty metadata as fasta/q files do not contain BAM/SAM specific data.

### Centrifuge Analysis

- Select a FASTQ file from the "File Selection" dropdown.
- Ensure the Centrifuge database path is correctly configured in the settings.
- Click "Run Centrifuge Analysis" to perform taxonomic classification.
- The results will be displayed once the analysis is complete.

## Configuration

- **Settings Panel**: Click "Settings" in the navigation bar to configure:
  - **Theme**: Choose between Dark or Light mode.
  - **Centrifuge Database Path**: Specify the path to your Centrifuge database.
  - **Maximum Read Length**: Set the maximum read length for processing.
  - **Minimum Quality Score**: Define the minimum base quality score for filtering.

- **Config File**: Settings are saved in `config.yaml`. You can edit this file directly or through the Settings panel.

## Folder Structure

- `app.py`: Main application script.
- `config.yaml`: Configuration file for application settings.
- `requirements.txt`: Python dependencies.
- `assets/`: Folder for custom CSS and assets used in the dashboard.
- `templates/`: HTML templates for the Dash application.
- `README.md`: Project documentation.

## Known Issues

- **Performance**: Processing large BAM files can be slow and may consume significant memory. Consider downsampling or using a machine with higher specifications.
- **File Compatibility**: Ensure uploaded files conform to the appropriate format specifications. Corrupted or improperly formatted files may cause errors.
- **Centrifuge Database**: The database for Centrifuge can be large. Ensure you have sufficient disk space and that the path is correctly set in the settings.
- **Data Selecting**: Selecting data from the read length histogram might not always work. This is a limitation of the plotly library.
- **Exporting data**: The export file function is still under development and not fully implemented.
- **Centrifuge implementation**: Proper centrifuge taxonomic classification is not working yet. This is a placeholder module with basic functionality.

## Contributing

Contributions are welcome! Please fork the repository and create a pull request with your improvements.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For questions or suggestions, please contact [michael.schneider6@mailbox.tu-dresden.de](mailto:michael.schneider6@mailbox.tu-dresden.de).
