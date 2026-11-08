#!/usr/bin/env python3
"""Run the packaged Python suite and record failures and skips explicitly."""
import argparse
import json
from pathlib import Path
import sys
import time
import unittest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parent), pattern='test_*.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = dict(tests=result.testsRun, failures=[str(t) for t, _ in result.failures],
                  errors=[str(t) for t, _ in result.errors],
                  skipped=[dict(test=str(t), reason=reason) for t, reason in result.skipped],
                  seconds=round(time.monotonic()-started, 3))
    report['status'] = 'PASS' if result.wasSuccessful() and not result.skipped else 'REVIEW'
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('[SUITE_REPORT] '+str(args.json))
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
