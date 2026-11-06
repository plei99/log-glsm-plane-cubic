"""Exercise imports, exact arithmetic, old entry points, and worker transport."""

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import multiprocessing
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ModuleInterfaceTests(unittest.TestCase):
    def test_fraction_interface_without_sage_initialization(self):
        result = subprocess.run(
            [sys.executable, "-c", """
import sys
from fractions import Fraction
from log_glsm.curves.elliptic import O3TwistedPlane
assert O3TwistedPlane.nonequivariant_i0_coefficient(3) == Fraction(1680)
assert 'sage.all' not in sys.modules
"""], capture_output=True, text=True, check=True,
        )
        self.assertEqual(result.stdout, "")

    def test_exact_hodge_request_survives_legacy_loads_and_pickle(self):
        from log_glsm.core.hodge import HodgeIntegralBackend, HodgeIntegralRequest
        from sage.repl.load import load
        from sage.all import QQ

        request = HodgeIntegralRequest(2, (2,), (2,))
        namespace = {"__name__": "not_a_command_line"}
        load(str(ROOT / "hodge_integrals.sage"), namespace)
        load(str(ROOT / "genus_four_marked.sage"), namespace)
        for item in (request, pickle.loads(pickle.dumps(request))):
            self.assertEqual(
                namespace["HodgeIntegralBackend"]().integral(item), QQ(7) / 5760
            )
            self.assertEqual(HodgeIntegralBackend().integral(item), QQ(7) / 5760)

    def test_vendored_chow_calculation_outside_repository(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, "-c", """
from log_glsm.core.hodge import AdmcyclesHodgeIntegralBackend
from sage.all import QQ
assert AdmcyclesHodgeIntegralBackend().integral(2, (2,), (2,)) == QQ(7)/5760
from log_glsm.targets.projective import ProjectiveCompleteIntersection
assert ProjectiveCompleteIntersection(3, (4,)).euler_characteristic() == 24
"""], cwd=directory, capture_output=True, text=True, check=True,
            )
        self.assertEqual(result.stdout, "")

    def test_spawned_workers_can_import_sage_functions(self):
        from log_glsm.cycles.marked import marked_multiplier

        with ProcessPoolExecutor(
            max_workers=1, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            result = pool.submit(marked_multiplier, 2, 1, ("psi",) * 5).result()
        self.assertEqual(result, 720)


if __name__ == "__main__":
    unittest.main()
