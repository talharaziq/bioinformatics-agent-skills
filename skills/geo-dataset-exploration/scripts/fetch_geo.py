#!/usr/bin/env python3
"""Run GEO metadata export without requiring an editable package install."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"src"))
from bioskills.cli import main
if __name__=="__main__": raise SystemExit(main(["geo-fetch",*sys.argv[1:]]))
