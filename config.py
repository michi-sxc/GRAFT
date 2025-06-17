import os
import yaml
from utils.temp_manager import TempDirectoryManager
import logging 
import tempfile

logger = logging.getLogger(__name__)

# --- Configuration Loading ---
config_data = {} 
try:
    with open('config.yaml', 'r') as f:
        loaded_yaml = yaml.safe_load(f)
        if loaded_yaml is not None: # Check if loading resulted in something
            config_data = loaded_yaml
            logger.info("Successfully loaded 'config.yaml'.")
        else:
            logger.warning("'config.yaml' was found but loaded as None (e.g., empty or malformed). Using default empty config.")
            # config_data remains {}
except FileNotFoundError:
    logger.info("'config.yaml' not found. Using default empty config.")
    # config_data remains {}
except yaml.YAMLError as e:
    logger.error(f"Error parsing 'config.yaml': {e}. Using default empty config.")
    # config_data remains {}


# --- Constants ---
# Color Palette
colors = {
    'primary': '#7E57C2',    # Softer Deep Purple (Material Design 300)
    'secondary': '#26A69A',  # Softer Teal (Material Design 300)
    
    'background': '#212529', # Bootstrap's default dark bg, slightly darker than #343a40
    'surface': '#2c3034',    # A bit lighter than background for cards/surfaces
    'surface_alt': '#343a40', # For slightly different surfaces (not used)

    'plot_bg': '#2c3034',    # Match card surface for plots
    'grid': 'rgba(255, 255, 255, 0.1)', # Softer grid lines
    
    'on_primary': '#FFFFFF',
    'on_secondary': '#000000',
    'on_background': '#E0E0E0', # Off-white for text on main background
    'on_surface': '#F5F5F5',   # Brighter white for text on cards
    'muted': '#9E9E9E',      # Softer muted text (grey 500)
    
    'highlight': '#EC407A',  # Softer Pink (Material Design 300)
    'highlight2': '#FFB300', # Softer Amber (Material Design 600)
    
    'line': '#42A5F5',       # Softer Blue for plot lines (Material Design 300)
    'marker': '#FF7043',     # Softer Orange for markers (Material Design 400)

    # Colors for specific elements (for contrast)
    'input_bg': '#343a40',
    'input_border': '#555e67',
    'input_focus_border': '#7E57C2', # Primary color
    'input_focus_shadow': 'rgba(126, 87, 194, 0.25)',

    'stats_tab_color': '#2c3034', # Match surface
}

# --- Temporary Directory Management (Singleton-like approach) ---
_temp_manager_instance = None
_custom_temp_dir_instance = None

# Define a FIXED base directory for all temp operations of this app session type
APP_TEMP_BASE = os.path.join(tempfile.gettempdir(), "graft_app_sessions")
os.makedirs(APP_TEMP_BASE, exist_ok=True) # Ensure this base exists

def get_temp_manager():
    global _temp_manager_instance
    if _temp_manager_instance is None:
        # The TempDirectoryManager will create subdirectories within APP_TEMP_BASE
        _temp_manager_instance = TempDirectoryManager(base_dir=APP_TEMP_BASE, max_age_hours=12)
        logger.info(f"TempDirectoryManager initialized with base: {APP_TEMP_BASE}")
    return _temp_manager_instance

def get_custom_temp_dir():
    global _custom_temp_dir_instance
    if _custom_temp_dir_instance is None:
        manager = get_temp_manager()
        # Creates a specific directory for *this run* of the app (or worker)

        _custom_temp_dir_instance = manager.create_temp_dir(prefix="graft_main_run_")
        os.environ["TMPDIR"] = _custom_temp_dir_instance
        logger.info(f"Main application run temporary directory set to: {_custom_temp_dir_instance}")
    return _custom_temp_dir_instance

# Initialize them
temp_manager = get_temp_manager()
CUSTOM_TEMP_DIR = get_custom_temp_dir() # used for uploads etc.

# TMPDIR environment variable for dependencies like pysam
os.environ["TMPDIR"] = CUSTOM_TEMP_DIR

print(f"Using temporary directory: {CUSTOM_TEMP_DIR}") # for debugging

# --- Other Constants ---
MAX_MISMATCH_DISTANCE = 50
RECORDS_PER_PAGE_VIEWER = 20