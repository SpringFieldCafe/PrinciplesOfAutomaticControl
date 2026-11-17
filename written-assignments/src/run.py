"""Recompute chapter data and figures. Invoke with --chapter ch07, or --all."""
from pathlib import Path
import argparse
import os
import runpy
import sys

ROOT = Path(__file__).resolve().parents[1]
STEPS = {
    'ch02': ['compute.py'],
    'ch04': ['compute.py'],
    'ch05': ['compute.py', 'redraw_axes.py'],
    'ch06': ['calculate.py'],
    'ch07': ['enhance.py', 'damping_intersections.py'],
    'ch08': ['calculate.py', 'settling_metrics.py'],
    'ch09': ['compute.py', 'refine_figures.py'],
    'ch10': ['compute.py', 'refine_figures.py'],
}
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--chapter', choices=STEPS)
    group.add_argument('--all', action='store_true')
    args = parser.parse_args()
    for name in ('MPLCONFIGDIR', 'PYTHONPYCACHEPREFIX', 'TEMP', 'TMP'):
        folder = ROOT / '.local' / 'cache' / name.lower()
        folder.mkdir(parents=True, exist_ok=True)
        os.environ[name] = str(folder)
    sys.dont_write_bytecode = True
    for chapter in STEPS if args.all else [args.chapter]:
        for script in STEPS[chapter]:
            print(f'{chapter}: {script}', flush=True)
            runpy.run_path(str(ROOT / 'src' / chapter / script), run_name='__main__')
if __name__ == '__main__':
    main()
