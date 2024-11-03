import os
import shutil
import tempfile
from pathlib import Path
import atexit
import logging
from datetime import datetime, timedelta

class TempDirectoryManager:
    def __init__(self, base_dir=None, max_age_hours=24):
        self.base_dir = base_dir or tempfile.gettempdir()
        self.temp_dirs = set()
        self.max_age = timedelta(hours=max_age_hours)
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Register cleanup on exit
        atexit.register(self.cleanup_all)
    
    def create_temp_dir(self):
        """Create a new temporary directory and track it"""
        temp_dir = tempfile.mkdtemp(dir=self.base_dir)
        self.temp_dirs.add(temp_dir)
        return temp_dir
    
    def cleanup_dir(self, dir_path):
        """Clean up a specific temporary directory"""
        try:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)
                self.temp_dirs.discard(dir_path)
                self.logger.info(f"Cleaned up temporary directory: {dir_path}")
        except Exception as e:
            self.logger.error(f"Error cleaning up directory {dir_path}: {str(e)}")
    
    def cleanup_all(self):
        """Clean up all tracked temporary directories"""
        for dir_path in list(self.temp_dirs):
            self.cleanup_dir(dir_path)
    
    def cleanup_old_dirs(self):
        """Clean up temporary directories older than max_age"""
        current_time = datetime.now()
        
        for dir_path in list(self.temp_dirs):
            try:
                dir_stat = os.stat(dir_path)
                dir_time = datetime.fromtimestamp(dir_stat.st_mtime)
                
                if current_time - dir_time > self.max_age:
                    self.cleanup_dir(dir_path)
            except FileNotFoundError:
                self.temp_dirs.discard(dir_path)
            except Exception as e:
                self.logger.error(f"Error checking directory age {dir_path}: {str(e)}")
    
    def clean_orphaned_temps(self):
        """Clean up orphaned temporary directories matching our pattern"""
        temp_pattern = "tmp*"
        base_path = Path(self.base_dir)
        
        for temp_dir in base_path.glob(temp_pattern):
            if temp_dir.is_dir() and temp_dir.stem.startswith('tmp'):
                try:
                    dir_stat = os.stat(temp_dir)
                    dir_time = datetime.fromtimestamp(dir_stat.st_mtime)
                    
                    if datetime.now() - dir_time > self.max_age:
                        shutil.rmtree(temp_dir)
                        self.logger.info(f"Cleaned up orphaned temp directory: {temp_dir}")
                except Exception as e:
                    self.logger.error(f"Error cleaning orphaned directory {temp_dir}: {str(e)}")