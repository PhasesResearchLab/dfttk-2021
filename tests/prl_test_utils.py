import argparse
import shutil
import os

MODULE_DIR = os.path.abspath(os.path.dirname(__file__))

parser = argparse.ArgumentParser(description='Utilities for testing with DFTTK')
parser.add_argument('-c', '--clean', action='store_true', help='Clean the temporary directories created in testing')
args = parser.parse_args()

if args.clean:
    for f in os.listdir(MODULE_DIR):
        if f.startswith('tmp'):
            shutil.rmtree(os.path.join(MODULE_DIR, f))
    scratch_dir = os.path.join(MODULE_DIR, 'scratch_dir')
    if os.path.exists(scratch_dir):
        shutil.rmtree(scratch_dir)
