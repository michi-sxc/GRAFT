#!/usr/bin/env python3

import subprocess
import sys
import os

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    if process.returncode != 0:
        print(f"Error executing command: {command}")
        print(error.decode('utf-8'))
        sys.exit(1)
    return output.decode('utf-8')

def install_python_packages():
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
    run_command(f"{sys.executable} -m pip install {' '.join(packages)}")

def install_dust_module():
    with open('setup.py', 'w') as f:
        f.write('''
from setuptools import setup, Extension
module = Extension('dust_module',
                   sources=['dust_module.c', 'sdust.c'],
                   libraries=['z'])  
setup(name='dust_module',
      version='1.0',
      description='DUST module',
      ext_modules=[module])
''')
    run_command(f"{sys.executable} setup.py install")

def install_centrifuge():
    if not os.path.exists('centrifuge'):
        print("Centrifuge is not installed. Do you want to install it? (y/n)")
        if input().lower() == 'y':
            if not run_command("conda --version").startswith("conda"):
                print("Bioconda is not installed. Do you want to install it? (y/n)")
                if input().lower() == 'y':
                    run_command("conda install -c conda-forge -c bioconda centrifuge")
                else:
                    print("Skipping Centrifuge installation.")
            else:
                run_command("conda install -c conda-forge -c bioconda centrifuge")
        else:
            print("Skipping Centrifuge installation.")

def install_mapdamage():
    if not os.path.exists('mapDamage'):
        run_command("git clone https://github.com/ginolhac/mapDamage.git")
    os.chdir('mapDamage')
    run_command(f"{sys.executable} setup.py install")
    os.chdir('..')

def check_r_dependencies():
    if run_command("R --version").startswith("R version"):
        r_packages = ["ggplot2", "Rcpp", "gam"]
        for package in r_packages:
            run_command(f"R -e \"if (!requireNamespace('{package}', quietly = TRUE)) install.packages('{package}', repos='http://cran.rstudio.com/')\"")
        
        run_command("sudo apt-get install -y r-cran-rcppgsl libgsl-dev")
    else:
        print("R is not installed. Do you want to install it? (y/n)")
        if input().lower() == 'y':
            run_command("sudo apt-get update && sudo apt-get install -y r-base")
            check_r_dependencies()
        else:
            print("Skipping R and its dependencies installation.")

def main():
    install_python_packages()
    install_dust_module()
    install_centrifuge()
    install_mapdamage()
    check_r_dependencies()
    print("Setup completed successfully!")

if __name__ == "__main__":
    main()