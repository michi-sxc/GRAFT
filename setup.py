#!/usr/bin/env python3

import subprocess
import sys
import os
import shutil
from datetime import datetime
from typing import Optional, List
import argparse
import platform
from pathlib import Path
import threading

class TerminalColors:
    """ANSI color codes for terminal output."""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

class SetupLogger:
    """Handles logging with timestamps and colored output."""
    
    @staticmethod
    def log(message: str, level: str = 'INFO') -> None:
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        color_map = {
            'INFO': TerminalColors.OKBLUE,
            'SUCCESS': TerminalColors.OKGREEN,
            'WARNING': TerminalColors.WARNING,
            'ERROR': TerminalColors.FAIL
        }
        color = color_map.get(level, TerminalColors.OKBLUE)
        print(f"{color}[{timestamp}] [{level}] {message}{TerminalColors.ENDC}")

class PathPermissionError(Exception):
    """Custom exception for permission-related errors."""
    pass

class FilePermissionHandler:
    """Handles file permission operations and checks."""
    
    @staticmethod
    def ensure_directory_permissions(path: Path) -> None:
        """
        Ensure directory has correct permissions (755).
        """
        try:
            if not path.exists():
                path.mkdir(parents=True, mode=0o755)
            else:
                path.chmod(0o755)
            SetupLogger.log(f"Set permissions for directory: {path}", 'SUCCESS')
        except Exception as e:
            raise PathPermissionError(f"Failed to set permissions for {path}: {str(e)}")

    @staticmethod
    def ensure_file_permissions(path: Path) -> None:
        """
        Ensure file has correct permissions (644).
        """
        try:
            path.chmod(0o644)
            SetupLogger.log(f"Set permissions for file: {path}", 'SUCCESS')
        except Exception as e:
            raise PathPermissionError(f"Failed to set permissions for {path}: {str(e)}")

    @staticmethod
    def check_write_permissions(path: Path) -> bool:
        """
        Check if we have write permissions for the given path.
        """
        if path.exists():
            return os.access(path, os.W_OK)
        return os.access(path.parent, os.W_OK)

    @staticmethod
    def make_executable(path: Path) -> None:
        """
        Make a file executable (755).
        """
        try:
            path.chmod(0o755)
            SetupLogger.log(f"Made file executable: {path}", 'SUCCESS')
        except Exception as e:
            raise PathPermissionError(f"Failed to make {path} executable: {str(e)}")

class CommandExecutor:
    """Handles command execution and error handling."""

    @staticmethod
    def run_command(command: List[str], check: bool = True, sudo: bool = False, timeout: int = None) -> Optional[str]:
        """
        Execute a command with optional sudo privileges.
        """
        if sudo and platform.system() == 'Linux' and os.geteuid() != 0:
            command = ['sudo'] + command

        SetupLogger.log(f"Executing: {' '.join(command)}")
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Function to read from a stream and log
            def reader(stream, log_level):
                for line in iter(stream.readline, ''):
                    if line:
                        SetupLogger.log(line.rstrip(), log_level)

            # Start threads to read stdout and stderr
            stdout_thread = threading.Thread(target=reader, args=(process.stdout, 'INFO'))
            stderr_thread = threading.Thread(target=reader, args=(process.stderr, 'ERROR'))

            stdout_thread.start()
            stderr_thread.start()

            # Wait for the process to complete
            process.wait(timeout=timeout)

            # Wait for threads to finish reading
            stdout_thread.join()
            stderr_thread.join()

            if process.returncode == 0:
                SetupLogger.log(f"Successfully executed: {' '.join(command)}", 'SUCCESS')
                return None
            else:
                SetupLogger.log(f"Error executing command: {' '.join(command)}", 'ERROR')
                if check:
                    sys.exit(1)
                return None

        except subprocess.TimeoutExpired as e:
            SetupLogger.log(f"Command timed out: {str(e)}", 'ERROR')
            process.kill()
            stdout_thread.join()
            stderr_thread.join()
            if check:
                sys.exit(1)
            return None
        except Exception as e:
            SetupLogger.log(f"Failed to execute command: {str(e)}", 'ERROR')
            if check:
                sys.exit(1)
            return None

class DependencyInstaller:
    """Handles installation of various dependencies."""
    
    def __init__(self):
        self.python_packages = [
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
        self.r_packages = ["ggplot2", "Rcpp", "gam"]

    def install_python_packages(self) -> None:
        """Install required Python packages."""
        SetupLogger.log("Installing Python packages...")
        try:
            for package in self.python_packages:
                SetupLogger.log(f"Installing Python package: {package}", 'INFO')
                CommandExecutor.run_command(
                    [sys.executable, '-m', 'pip', 'install', package]
                )
            SetupLogger.log("Python packages installed successfully.", 'SUCCESS')
        except Exception as e:
            SetupLogger.log(f"Failed to install Python packages: {str(e)}", 'ERROR')
            sys.exit(1)
    
    def install_dust_module(self) -> None:
        """Install DUST module from source with proper permission handling."""
        SetupLogger.log("Installing DUST module...")
        dust_dir = Path('DUST')
        if not dust_dir.exists():
            SetupLogger.log(f"Error: {dust_dir} directory not found.", 'ERROR')
            sys.exit(1)

        build_dir = Path('dust_build')
        
        # Check permissions before proceeding
        if not FilePermissionHandler.check_write_permissions(Path.cwd()):
            SetupLogger.log("Insufficient permissions in current directory.", 'ERROR')
            raise PathPermissionError("Cannot write to current directory")

        # Create build directory with proper permissions
        try:
            FilePermissionHandler.ensure_directory_permissions(build_dir)
        except PathPermissionError as e:
            SetupLogger.log(str(e), 'ERROR')
            sys.exit(1)

        files_to_copy = ['dust_module.c', 'sdust.c', 'sdust.h', 'kdq.h', 
                        'ketopt.h', 'kalloc.h', 'kseq.h', 'kvec.h']
        
        # Copy files with proper permissions
        for file in files_to_copy:
            src_path = dust_dir / file
            dst_path = build_dir / file
            if src_path.exists():
                try:
                    shutil.copy2(src_path, dst_path)
                    # FilePermissionHandler.ensure_file_permissions(dst_path)
                    SetupLogger.log(f"Copied {file}", 'INFO')
                except Exception as e:
                    SetupLogger.log(f"Error copying {file}: {str(e)}", 'ERROR')
                    shutil.rmtree(build_dir, ignore_errors=True)
                    sys.exit(1)
            else:
                SetupLogger.log(f"Warning: {file} not found in {dust_dir}.", 'WARNING')

        setup_content = '''
from setuptools import setup, Extension

module = Extension(
    'dust_module',
    sources=['dust_module.c', 'sdust.c'],
    include_dirs=['.', './DUST'],
    libraries=['z'],
    extra_compile_args=['-O3']
)

setup(
    name='dust_module',
    version='1.0',
    description='DUST module for sequence complexity filtering',
    ext_modules=[module],
    python_requires='>=3.6'
)
'''
        setup_path = build_dir / 'setup.py'
        setup_path.write_text(setup_content)
        # FilePermissionHandler.ensure_file_permissions(setup_path)

        current_dir = os.getcwd()
        try:
            os.chdir(build_dir)
            # Try installation without sudo first
            try:
                CommandExecutor.run_command([sys.executable, '-m', 'pip', 'install', '--no-deps', '.'])
            except Exception:
                # If regular installation fails, try with sudo
                SetupLogger.log("Regular installation failed, attempting with sudo...", 'WARNING')
                CommandExecutor.run_command([sys.executable, '-m', 'pip', 'install', '--no-deps', '.'], sudo=True)
            SetupLogger.log("DUST module installed successfully.", 'SUCCESS')
        except Exception as e:
            SetupLogger.log(f"Failed to install DUST module: {str(e)}", 'ERROR')
            sys.exit(1)
        finally:
            os.chdir(current_dir)
            shutil.rmtree(build_dir, ignore_errors=True)

    def check_and_install_r_dependencies(self) -> None:
        """Check and install R and required R packages."""
        SetupLogger.log("Checking R dependencies...")
        r_version = CommandExecutor.run_command(['R', '--version'], check=False)
        if not r_version:
            SetupLogger.log("R is not installed.", 'WARNING')
            if platform.system() == 'Linux':
                SetupLogger.log("Installing R base...", 'INFO')
                CommandExecutor.run_command(['sudo', 'apt-get', 'update'])
                CommandExecutor.run_command(['sudo', 'apt-get', 'install', '-y', 'r-base'])
            else:
                SetupLogger.log("Please install R manually for your operating system.", 'ERROR')
                return

        # Install system dependencies required by R packages
        if platform.system() == 'Linux':
            SetupLogger.log("Installing system dependencies for R packages...", 'INFO')
            system_packages = [
                'libxml2-dev',
                'libcurl4-openssl-dev',
                'libssl-dev',
                'libblas-dev',
                'liblapack-dev',
                'gfortran',
                'r-base-dev',
                'r-recommended',
            ]
            CommandExecutor.run_command(['sudo', 'apt-get', 'install', '-y'] + system_packages)

        # Ensure R user library path exists
        try:
            r_libs_user = subprocess.check_output(
                ['Rscript', '-e', 'cat(Sys.getenv("R_LIBS_USER"))'], text=True
            ).strip()
            if not r_libs_user:
                # Set a default user library path if R_LIBS_USER is empty
                r_libs_user = os.path.expanduser('~/R/library')
            os.makedirs(r_libs_user, exist_ok=True)
            SetupLogger.log(f"R user library path: {r_libs_user}", 'INFO')
        except Exception as e:
            SetupLogger.log(f"Failed to create R user library directory: {str(e)}", 'ERROR')
            sys.exit(1)

        for package in self.r_packages:
            SetupLogger.log(f"Installing R package: {package}", 'INFO')
            r_cmd = [
                'R',
                '--slave',
                '--no-restore',
                '-e',
                f"if (!requireNamespace('{package}', quietly = TRUE)) "
                f"{{ install.packages('{package}', repos='http://cran.rstudio.com/', "
                f"dependencies=TRUE, quiet=TRUE, lib='{r_libs_user}') }}"
            ]
            CommandExecutor.run_command(r_cmd)

        if platform.system() == 'Linux':
            SetupLogger.log("Installing additional R system packages...", 'INFO')
            CommandExecutor.run_command(
                ['sudo', 'apt-get', 'install', '-y', 'r-cran-rcppgsl', 'libgsl-dev']
            )

def is_wsl_on_ntfs() -> bool:
    """
    Check if the script is running under WSL on a Windows filesystem (NTFS).
    """
    if 'microsoft' in platform.uname().release.lower():
        current_path = Path.cwd()
        if current_path.is_mount() and str(current_path).startswith('/mnt/'):
            return True
    return False

def move_to_wsl_filesystem() -> None:
    """
    Move the project directory to the WSL filesystem and re-run the script.
    """
    home_dir = Path.home()
    projects_dir = home_dir / 'projects'
    projects_dir.mkdir(exist_ok=True)

    current_dir = Path.cwd()
    target_dir = projects_dir / current_dir.name

    SetupLogger.log(f"Copying project to WSL filesystem at {target_dir}...", 'INFO')

    if target_dir.exists():
        shutil.rmtree(target_dir)

    try:
        shutil.copytree(current_dir, target_dir, symlinks=True)
        SetupLogger.log("Project copied successfully.", 'SUCCESS')
    except Exception as e:
        SetupLogger.log(f"Failed to copy project: {str(e)}", 'ERROR')
        sys.exit(1)

    # Re-run the script from the new location
    new_script = target_dir / sys.argv[0]
    SetupLogger.log("Re-running setup script from the WSL filesystem...", 'INFO')

    os.chdir(target_dir)

    # Execute the script in the new location
    new_command = [sys.executable, str(new_script)] + sys.argv[1:]
    try:
        subprocess.check_call(new_command)
        sys.exit(0)
    except subprocess.CalledProcessError as e:
        SetupLogger.log(f"Failed to re-run setup script: {str(e)}", 'ERROR')
        sys.exit(1)

def main():
    """Main setup function with permission handling."""
    parser = argparse.ArgumentParser(description='Setup script for bioinformatics pipeline')
    parser.add_argument('--skip-r', action='store_true', help='Skip R dependencies installation')
    parser.add_argument('--skip-dust', action='store_true', help='Skip DUST module installation')
    parser.add_argument('--source-dir', help='Source directory to copy from')
    parser.add_argument('--force-sudo', action='store_true', help='Force using sudo for installations')
    args = parser.parse_args()

    SetupLogger.log("Starting setup...")

    # Check if running under WSL on Windows filesystem
    if is_wsl_on_ntfs():
        SetupLogger.log("Detected WSL running on Windows filesystem.", 'WARNING')
        SetupLogger.log("Moving project to WSL filesystem for proper permissions.", 'INFO')
        move_to_wsl_filesystem()
        return  # The script will exit after re-running from new location

    # Check if we have write permissions in the current directory
    if not FilePermissionHandler.check_write_permissions(Path.cwd()):
        SetupLogger.log("Warning: Insufficient permissions in current directory.", 'WARNING')
        if args.force_sudo:
            SetupLogger.log("Proceeding with sudo permissions...", 'INFO')
        else:
            SetupLogger.log("Please run the script from a directory where you have write permissions, or use --force-sudo", 'ERROR')
            sys.exit(1)

    try:
        installer = DependencyInstaller()

        # Install Python packages
        installer.install_python_packages()
        
        if not args.skip_dust:
            installer.install_dust_module()
            
        if not args.skip_r:
            installer.check_and_install_r_dependencies()
            
        SetupLogger.log("Setup completed successfully!", 'SUCCESS')
        
    except KeyboardInterrupt:
        SetupLogger.log("\nSetup interrupted by user.", 'WARNING')
        sys.exit(1)
    except PathPermissionError as e:
        SetupLogger.log(f"Permission error: {str(e)}", 'ERROR')
        SetupLogger.log("Try running with --force-sudo flag", 'INFO')
        sys.exit(1)
    except Exception as e:
        SetupLogger.log(f"Setup failed: {str(e)}", 'ERROR')
        sys.exit(1)

if __name__ == "__main__":
    main()
