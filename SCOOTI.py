import argparse
from scooti import __version__

parser = argparse.ArgumentParser(description="SCOOTI CLI")
parser.add_argument("--version", action="version", version=f"SCOOTI {__version__}")

args = parser.parse_args()

