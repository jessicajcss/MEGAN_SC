"""Make ``src/`` importable in tests without an editable install (reproducibility)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
