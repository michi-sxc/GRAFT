#!/usr/bin/env python3

import subprocess
import sys
import os
import shutil
from datetime import datetime
from pathlib import Path
import platform

class SetupLogger:
    """Simple logging with timestamps"""
    @staticmethod
    def log(message: str, level: str = 'INFO') -> None:
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"[{timestamp}] [{level}] {message}")

class CommandExecutor:
    """Handles command execution"""
    @staticmethod
    def run_command(command, check=True):
        SetupLogger.log(f"Executing: {' '.join(command)}")
        try:
            result = subprocess.run(
                command,
                check=check,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if result.stdout:
                SetupLogger.log(result.stdout.strip())
            if result.stderr:
                SetupLogger.log(result.stderr.strip(), 'WARNING')
            return result
        except subprocess.CalledProcessError as e:
            SetupLogger.log(f"Command failed: {str(e)}", 'ERROR')
            if check:
                sys.exit(1)
            return None

class SystemDependencyInstaller:
    """Handles system-specific dependency installation"""
    
    @staticmethod
    def install_homebrew():
        """Install Homebrew on MacOS if not present"""
        try:
            CommandExecutor.run_command(['brew', '--version'], check=False)
        except FileNotFoundError:
            SetupLogger.log("Installing Homebrew...")
            install_cmd = '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
            subprocess.run(install_cmd, shell=True, check=True)

    @staticmethod
    def install_linux_deps():
        """Install dependencies on Linux"""
        deps = [
            'build-essential',
            'gcc',
            'python3-dev',
            'zlib1g-dev',
            'liblzma-dev',
            'libzstd-dev'
        ]
        
        CommandExecutor.run_command(['sudo', 'apt-get', 'update'])
        CommandExecutor.run_command(['sudo', 'apt-get', 'install', '-y'] + deps)

    @staticmethod
    def install_macos_deps():
        """Install dependencies on MacOS"""
        deps = [
            'gcc',
            'python',
            'zlib',
            'xz',
            'zstd'
        ]
        
        # Install Xcode Command Line Tools if needed
        try:
            CommandExecutor.run_command(['xcode-select', '--version'], check=False)
        except FileNotFoundError:
            SetupLogger.log("Installing Xcode Command Line Tools...")
            CommandExecutor.run_command(['xcode-select', '--install'])

        # Install Homebrew if needed
        SystemDependencyInstaller.install_homebrew()
        
        # Update Homebrew and install dependencies
        CommandExecutor.run_command(['brew', 'update'])
        for dep in deps:
            CommandExecutor.run_command(['brew', 'install', dep])

class DependencyInstaller:
    """Handles installation of Python packages and DUST module"""
    def __init__(self):
        self.python_packages = [
            "wheel>=0.40.0",
            "setuptools>=69.0.0",
            "numpy>=1.26.4",
            "pysam==0.22.1",
            "biopython==1.81",
            "dash==2.17.1",
            "dash-bootstrap-components==1.6.0",
            "dash-iconify==0.1.2",
            "dash_dangerously_set_inner_html==0.0.2",
            "plotly==5.22.0",
            "pandas==2.2.1",
            "psutil==6.0.0",
            "PyYAML==6.0.2",
            "scipy==1.15.1",
        ]

    def install_system_deps(self):
        """Install system dependencies based on platform"""
        system = platform.system()
        if system == 'Linux':
            SystemDependencyInstaller.install_linux_deps()
        elif system == 'Darwin':  # MacOS
            SystemDependencyInstaller.install_macos_deps()
        else:
            SetupLogger.log(f"Unsupported operating system: {system}", 'WARNING')

    def install_python_packages(self):
        """Install required Python packages"""
        SetupLogger.log("Installing Python packages...")
        
        # Set environment variables for MacOS
        if platform.system() == 'Darwin':
            # Use Homebrew OpenSSL for package installation
            brew_prefix = subprocess.check_output(['brew', '--prefix']).decode().strip()
            os.environ['CFLAGS'] = f'-I{brew_prefix}/opt/openssl/include'
            os.environ['LDFLAGS'] = f'-L{brew_prefix}/opt/openssl/lib'

        for package in self.python_packages:
            CommandExecutor.run_command([
                sys.executable, '-m', 'pip', 'install',
                '--no-cache-dir', package
            ])

    def install_dust_module(self):
        """Install DUST module from source"""
        SetupLogger.log("Installing DUST module...")
        dust_dir = Path('DUST')
        if not dust_dir.exists():
            SetupLogger.log("Error: DUST directory not found", 'ERROR')
            sys.exit(1)

        build_dir = Path('dust_build')
        if build_dir.exists():
            shutil.rmtree(build_dir)
        build_dir.mkdir()

        # Copy required files
        for file in ['dust_module.c', 'sdust.c', 'sdust.h', 'kdq.h', 
                    'ketopt.h', 'kalloc.h', 'kseq.h', 'kvec.h']:
            src = dust_dir / file
            if src.exists():
                shutil.copy2(src, build_dir / file)

        # Create setup.py with platform-specific settings
        include_dirs = ['.']
        if platform.system() == 'Darwin':
            # Add Homebrew include directories for MacOS
            brew_prefix = subprocess.check_output(['brew', '--prefix']).decode().strip()
            include_dirs.extend([
                f'{brew_prefix}/include',
                f'{brew_prefix}/opt/zlib/include',
            ])

        setup_content = f'''
from setuptools import setup, Extension

module = Extension(
    'dust_module',
    sources=['dust_module.c', 'sdust.c'],
    include_dirs={include_dirs},
    libraries=['z'],
    extra_compile_args=['-O3']
)

setup(
    name='dust_module',
    version='1.0',
    description='DUST module for sequence complexity filtering',
    ext_modules=[module]
)
'''
        (build_dir / 'setup.py').write_text(setup_content)

        # Install module
        current_dir = os.getcwd()
        try:
            os.chdir(build_dir)
            CommandExecutor.run_command([
                sys.executable, '-m', 'pip', 'install', '--no-deps', '.'
            ])
        finally:
            os.chdir(current_dir)
            shutil.rmtree(build_dir)

def main():
    """Main setup function"""
    system = platform.system()
    if system not in ['Linux', 'Darwin']:
        SetupLogger.log(f"Unsupported operating system: {system}", 'ERROR')
        sys.exit(1)
        
    SetupLogger.log(f"Starting installation on {system}...")
    
    installer = DependencyInstaller()
    
    try:
        installer.install_system_deps()
        installer.install_python_packages()
        installer.install_dust_module()
        SetupLogger.log("Installation completed successfully!", 'SUCCESS')
    except KeyboardInterrupt:
        SetupLogger.log("Installation interrupted by user", 'WARNING')
        sys.exit(1)
    except Exception as e:
        SetupLogger.log(f"Installation failed: {str(e)}", 'ERROR')
        sys.exit(1)

if __name__ == "__main__":
    main()
