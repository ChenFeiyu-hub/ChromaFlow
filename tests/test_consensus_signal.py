import importlib.util
from pathlib import Path
import random
import unittest

spec = importlib.util.spec_from_file_location('consensus_signal', Path(__file__).resolve().parents[1]/'lib/steps/support/consensus_signal.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class OverlapTests(unittest.TestCase):
    def test_nested_long_peak_is_not_lost(self):
        index = mod.build_index([(0, 4000000, 9), (100, 200, 3)])
        self.assertEqual(mod.max_signal(*index, 3000000, 3000100), 9)
        self.assertIsNone(mod.max_signal(*index, 4000000, 4000100))

    def test_matches_brute_force_for_random_intervals(self):
        rng = random.Random(10)
        intervals = sorted((s, s+rng.randrange(1, 2000000), rng.random()*100)
                           for s in (rng.randrange(10000000) for _ in range(1000)))
        index = mod.build_index(intervals)
        for _ in range(100):
            start = rng.randrange(12000000)
            end = start+300
            expected = max((v for s,e,v in intervals if s<end and e>start), default=None)
            self.assertEqual(mod.max_signal(*index, start, end), expected)
