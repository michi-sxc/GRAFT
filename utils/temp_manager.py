import os
import shutil
import tempfile
from pathlib import Path
import logging
from datetime import datetime, timedelta

class TempDirectoryManager:
    def __init__(self, base_dir=None, max_age_hours=24):
        self.base_dir = Path(base_dir if base_dir is not None else tempfile.gettempdir())
        self.temp_dirs = set()
        self.max_age = timedelta(hours=max_age_hours)

        # Setup logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

        # Register cleanup on exit - use cautiously
        # atexit.register(self.cleanup_all) # Can be problematic in some server environments 
        self.logger.info("TempDirectoryManager initialized.")
        # Consider periodic cleanup instead of relying solely on atexit
        # self.schedule_cleanup() # If using a scheduler library

    def create_temp_dir(self, prefix="graft_temp_"):
        """Create a new temporary directory and track it"""
        try:
            self.base_dir.mkdir(parents=True, exist_ok=True)
            # Pass self.base_dir explicitly to mkdtemp
            temp_dir_path = tempfile.mkdtemp(prefix=prefix, dir=self.base_dir)
            self.temp_dirs.add(temp_dir_path)
            self.logger.info(f"Created temporary directory: {temp_dir_path} within {self.base_dir}")
            return temp_dir_path
        except Exception as e:
            self.logger.error(f"Failed to create temporary directory in {self.base_dir}: {e}")
            # Fallback to system default if custom fails
            temp_dir_path = tempfile.mkdtemp(prefix=prefix)
            self.temp_dirs.add(temp_dir_path)
            self.logger.warning(f"Fell back to system default temp dir: {temp_dir_path}")
            return temp_dir_path

    def cleanup_dir(self, dir_path):
        """Clean up a specific temporary directory"""
        if dir_path in self.temp_dirs:
            try:
                if os.path.exists(dir_path):
                    shutil.rmtree(dir_path)
                    self.logger.info(f"Cleaned up tracked temporary directory: {dir_path}")
                self.temp_dirs.discard(dir_path)
            except Exception as e:
                self.logger.error(f"Error cleaning up tracked directory {dir_path}: {str(e)}")
        else:
             self.logger.warning(f"Attempted to clean up untracked directory: {dir_path}")


    def cleanup_all(self):
        """Clean up all tracked temporary directories"""
        self.logger.info(f"Starting cleanup of {len(self.temp_dirs)} tracked temporary directories.")
        # Iterate over a copy because self.temp_dirs might be modified during iteration
        for dir_path in list(self.temp_dirs):
            self.cleanup_dir(dir_path)
        self.logger.info("Finished cleanup of tracked directories.")

    def cleanup_old_dirs(self, base_path=None, prefix="graft_temp_"):
        """Clean up temporary directories older than max_age in a given base path."""
        current_time = datetime.now()
        search_path = Path(base_path or self.base_dir)
        cleaned_count = 0
        self.logger.info(f"Starting cleanup of old directories in {search_path} with prefix '{prefix}' older than {self.max_age}.")

        if not search_path.is_dir():
            self.logger.warning(f"Base path for cleanup does not exist or is not a directory: {search_path}")
            return

        for item in search_path.glob(f'{prefix}*'):
            if item.is_dir():
                try:
                    dir_stat = item.stat()
                    dir_time = datetime.fromtimestamp(dir_stat.st_mtime)

                    if current_time - dir_time > self.max_age:
                        self.logger.info(f"Found old directory: {item} (age: {current_time - dir_time}). Attempting removal.")
                        shutil.rmtree(item)
                        # If it was tracked, remove it from tracking
                        if str(item) in self.temp_dirs:
                            self.temp_dirs.discard(str(item))
                        cleaned_count += 1
                        self.logger.info(f"Successfully removed old directory: {item}")
                except FileNotFoundError:
                     self.logger.warning(f"Directory {item} not found during cleanup check (possibly removed concurrently).")
                     if str(item) in self.temp_dirs:
                         self.temp_dirs.discard(str(item))
                except Exception as e:
                    self.logger.error(f"Error checking or cleaning directory {item}: {str(e)}")

        self.logger.info(f"Finished cleanup of old directories. Removed {cleaned_count} directories.")

    # Consider adding a method to schedule cleanup_old_dirs periodically if the app runs long term
    # def schedule_cleanup(self):
    #     # Implementation using APScheduler or similar
    #     pass