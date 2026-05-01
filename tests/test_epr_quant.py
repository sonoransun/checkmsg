import numpy as np
import pytest

from checkmsg.epr import (
    SpinSystem,
    absolute_spin_count,
    count_spins,
    double_integral,
)
from checkmsg.synthetic import generate_epr


def _build_pair(spin_count: float, ref_spins: float):
    sys_ = SpinSystem(name="probe", S=0.5, g=2.0023, linewidth_mT=0.3)
    fields = np.linspace(335, 345, 1001)
    sample = generate_epr(sys_, fields, 9.5, spin_count=spin_count, noise=0.0,
                          orientations=(7, 1), seed=0)
    reference = generate_epr(sys_, fields, 9.5, spin_count=ref_spins, noise=0.0,
                             orientations=(7, 1), seed=1)
    return sample, reference


def test_double_integral_scales_linearly_with_spin_count():
    a, _ = _build_pair(1.0, 1.0)
    b, _ = _build_pair(3.0, 1.0)
    da = double_integral(a)
    db = double_integral(b)
    assert abs(db / da - 3.0) / 3.0 < 0.10


def test_count_spins_recovers_known_ratio():
    sample, ref = _build_pair(2.5, 1.0)
    sc = count_spins(sample, ref, reference_spins=1.0e18)
    # Sample/reference spin ratio is 2.5; expect SpinCount.total ~ 2.5e18 within 10%.
    assert 2.0e18 < sc.total < 3.0e18
    assert 2.0 < sc.ratio_to_reference < 3.0


def test_count_spins_rejects_frequency_mismatch():
    sys_ = SpinSystem(name="p", S=0.5, g=2.0)
    fields = np.linspace(335, 345, 401)
    a = generate_epr(sys_, fields, 9.5, noise=0.0, orientations=(5, 1))
    fields_q = np.linspace(1245, 1255, 401)
    b = generate_epr(sys_, fields_q, 35.0, noise=0.0, orientations=(5, 1))
    with pytest.raises(ValueError, match="frequency mismatch"):
        count_spins(a, b, reference_spins=1.0)


def test_per_gram_reporting():
    sample, ref = _build_pair(1.0, 1.0)
    sc = count_spins(sample, ref, reference_spins=1.0e18, sample_mass_g=0.5)
    assert sc.per_gram is not None
    assert abs(sc.per_gram - sc.total / 0.5) < 1e-3


def test_absolute_count_formula_round_trip():
    sample, _ = _build_pair(1.0, 1.0)
    di = double_integral(sample)
    sc = absolute_spin_count(sample, cavity_Q=1.0, modulation_amplitude_mT=1.0,
                             microwave_power_mW=1.0, calibration_factor=1.0)
    # With Q=mod=sqrt(P)=cal=1, total should equal the bare double integral.
    assert abs(sc.total - di) / max(abs(di), 1e-12) < 1e-9


def test_absolute_count_rejects_zero_inputs():
    sample, _ = _build_pair(1.0, 1.0)
    with pytest.raises(ValueError):
        absolute_spin_count(sample, cavity_Q=0, modulation_amplitude_mT=1.0,
                            microwave_power_mW=1.0)
