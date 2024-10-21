#!/usr/bin/env python3

import subprocess
import sys
import os
import shutil
from datetime import datetime

# Colors for terminal output
class TerminalColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

# Utility function to print progress messages
def log(message, level='INFO'):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if level == 'INFO':
        print(f"{TerminalColors.OKBLUE}[{timestamp}] [INFO] {message}{TerminalColors.ENDC}")
    elif level == 'SUCCESS':
        print(f"{TerminalColors.OKGREEN}[{timestamp}] [SUCCESS] {message}{TerminalColors.ENDC}")
    elif level == 'WARNING':
        print(f"{TerminalColors.WARNING}[{timestamp}] [WARNING] {message}{TerminalColors.ENDC}")
    elif level == 'ERROR':
        print(f"{TerminalColors.FAIL}[{timestamp}] [ERROR] {message}{TerminalColors.ENDC}")

# Function to execute shell commands and handle errors
def run_command(command):
    log(f"Executing: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    if process.returncode != 0:
        log(f"Error executing command: {command}\n{error.decode('utf-8')}", level='ERROR')
        sys.exit(1)
    log(f"Successfully executed: {command}", level='SUCCESS')
    return output.decode('utf-8')

def install_python_packages():
    log("Installing Python packages...")
    packages = [
        "pysam==0.22.1",
        "biopython==1.81",
        "numpy==1.26.4",
        "dash==2.17.1",
        "dash-bootstrap-components==1.6.0",
        "dash-iconify==0.1.2",
        "dash_dangerously_set_inner_html==0.0.2",
        "plotly==5.22.0",
        "pandas==2.2.1",
        "celery==5.4.0",
        "PyYAML==6.0.2",
        "ansi2html==1.9.2"
    ]
    try:
        run_command(f"{sys.executable} -m pip install {' '.join(packages)}")
        log("Python packages installed successfully.", level='SUCCESS')
    except Exception as e:
        log(str(e), level='ERROR')

def install_dust_module():
    log("Installing DUST module...")
    dust_dir = 'DUST'
    if not os.path.exists(dust_dir):
        log(f"Error: {dust_dir} directory not found.", level='ERROR')
        sys.exit(1)

    build_dir = 'dust_build'
    os.makedirs(build_dir, exist_ok=True)

    files_to_copy = ['dust_module.c', 'sdust.c', 'sdust.h', 'kdq.h', 'ketopt.h', 'kalloc.h', 'kseq.h', 'kvec.h']
    for file in files_to_copy:
        src_path = os.path.join(dust_dir, file)
        if os.path.exists(src_path):
            shutil.copy(src_path, build_dir)
            log(f"Copied {file} to build directory.", level='INFO')
        else:
            log(f"Warning: {file} not found in {dust_dir}. Build might fail.", level='WARNING')

    with open(os.path.join(build_dir, 'setup.py'), 'w') as f:
        f.write('''
from setuptools import setup, Extension

module = Extension('dust_module',
                   sources=['dust_module.c', 'sdust.c'],
                   include_dirs=['.', './DUST'],  
                   libraries=['z'])

setup(name='dust_module',
      version='1.0',
      description='DUST module',
      ext_modules=[module])
''')
    
    os.chdir(build_dir)
    try:
        run_command(f"{sys.executable} -m pip install --no-deps .")
        log("DUST module installed successfully.", level='SUCCESS')
    except Exception as e:
        log(str(e), level='ERROR')
    finally:
        os.chdir('..')
        shutil.rmtree(build_dir)

def install_centrifuge():
    log("Checking if Centrifuge is installed...")
    if not os.path.exists('centrifuge'):
        log("Centrifuge is not installed. User interaction required.", level='WARNING')
        if input("Do you want to install Centrifuge? (y/n): ").lower() == 'y':
            try:
                run_command("conda install -c conda-forge -c bioconda centrifuge")
                log("Centrifuge installed successfully.", level='SUCCESS')
            except Exception as e:
                log(str(e), level='ERROR')
        else:
            log("Skipping Centrifuge installation.", level='WARNING')
    else:
        log("Centrifuge is already installed.", level='SUCCESS')

def install_mapdamage():
    log("Installing mapDamage...")
    if not os.path.exists('mapDamage'):
        try:
            run_command("git clone https://github.com/ginolhac/mapDamage.git")
            log("Cloned mapDamage repository.", level='SUCCESS')
        except Exception as e:
            log(str(e), level='ERROR')
            return
    os.chdir('mapDamage')
    try:
        run_command(f"{sys.executable} -m pip install .")
        log("mapDamage installed successfully.", level='SUCCESS')
    except Exception as e:
        log(str(e), level='ERROR')
    finally:
        os.chdir('..')

def check_r_dependencies():
    log("Checking R dependencies...")
    if run_command("R --version").startswith("R version"):
        log("R is installed. Installing R packages...", level='INFO')
        r_packages = ["ggplot2", "Rcpp", "gam"]
        for package in r_packages:
            try:
                run_command(f"R -e \"if (!requireNamespace('{package}', quietly = TRUE)) install.packages('{package}', repos='http://cran.rstudio.com/')\"")
                log(f"R package '{package}' installed successfully.", level='SUCCESS')
            except Exception as e:
                log(str(e), level='ERROR')
        try:
            run_command("sudo apt-get install -y r-cran-rcppgsl libgsl-dev")
            log("R dependencies installed successfully.", level='SUCCESS')
        except Exception as e:
            log(str(e), level='ERROR')
    else:
        log("R is not installed. User interaction required.", level='WARNING')
        if input("Do you want to install R? (y/n): ").lower() == 'y':
            try:
                run_command("sudo apt-get update && sudo apt-get install -y r-base")
                check_r_dependencies()
            except Exception as e:
                log(str(e), level='ERROR')
        else:
            log("Skipping R and its dependencies installation.", level='WARNING')

def main():
    log("Starting setup...")
    install_python_packages()
    install_dust_module()
    install_centrifuge()
    install_mapdamage()
    check_r_dependencies()
    log("Setup completed successfully!", level='SUCCESS')

if __name__ == "__main__":
    main()
