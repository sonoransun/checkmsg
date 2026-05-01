import numpy as np
import pytest

from checkmsg import laicpms
from checkmsg.refdata.icpms_data import LAMBDA_235, LAMBDA_238, U238_OVER_U235
from checkmsg.synthetic import generate_laicpms_run


def _zircon_run(age_yr, *, u238_atoms=5e6, common_pb=1.0, sens=1000.0, seed=0):
    """Synthesize a zircon transient consistent with the given age (years)."""
    pb206_per_u238 = np.expm1(LAMBDA_238 * age_yr)
    pb207_per_u235 = np.expm1(LAMBDA_235 * age_yr)
    u235_atoms = u238_atoms / U238_OVER_U235
    signals = {
        "U238": u238_atoms * sens,
        "U235": u235_atoms * sens,
        "Pb206": u238_atoms * pb206_per_u238 * sens,
        "Pb207": u235_atoms * pb207_per_u235 * sens,
        "Pb204": common_pb * sens,
        "Pb208": common_pb * sens * 2,
    }
    return generate_laicpms_run(
        signals, blank_cps={k: 0.5 for k in signals},
        sample_duration_s=15, blank_duration_s=5,
        noise_factor=0.02, seed=seed,
    )


@pytest.mark.parametrize("age_Ma", [50.0, 100.0, 500.0, 2500.0])
def test_u_pb_age_round_trip(age_Ma):
    run = _zircon_run(age_yr=age_Ma * 1e6, seed=int(age_Ma))
    recovered = laicpms.u_pb_age(run)
    assert abs(recovered - age_Ma) < 5.0, f"got {recovered:.1f}, expected {age_Ma}"


def test_u_pb_age_requires_pb206_and_u238():
    run = generate_laicpms_run({"Pb208": 1000.0}, seed=10)
    with pytest.raises(ValueError, match="Pb206 and U238"):
        laicpms.u_pb_age(run)


def test_u_pb_age_rejects_discordant_data():
    """Synthesize a run with a 50 Ma 206/238 age but a 200 Ma 207/235 age — should reject."""
    age_206_yr = 50e6
    age_207_yr = 200e6
    pb206_per_u238 = np.expm1(LAMBDA_238 * age_206_yr)
    pb207_per_u235 = np.expm1(LAMBDA_235 * age_207_yr)
    u238_atoms = 5e6
    u235_atoms = u238_atoms / U238_OVER_U235
    sens = 1000.0
    signals = {
        "U238": u238_atoms * sens,
        "U235": u235_atoms * sens,
        "Pb206": u238_atoms * pb206_per_u238 * sens,
        "Pb207": u235_atoms * pb207_per_u235 * sens,
    }
    run = generate_laicpms_run(signals, blank_cps={k: 0.5 for k in signals},
                               sample_duration_s=15, blank_duration_s=5,
                               noise_factor=0.005, seed=99)
    with pytest.raises(ValueError, match="discordant"):
        laicpms.u_pb_age(run, discordance_tolerance=0.10)


def test_u_pb_age_uses_natural_u_ratio_when_u235_missing():
    """If U235 channel isn't measured, the function derives it from U238 via natural ratio."""
    age_yr = 100e6
    pb206_per_u238 = np.expm1(LAMBDA_238 * age_yr)
    pb207_per_u235 = np.expm1(LAMBDA_235 * age_yr)
    u238_atoms = 5e6
    u235_atoms = u238_atoms / U238_OVER_U235
    sens = 1000.0
    signals = {
        "U238": u238_atoms * sens,
        "Pb206": u238_atoms * pb206_per_u238 * sens,
        "Pb207": u235_atoms * pb207_per_u235 * sens,
    }
    run = generate_laicpms_run(signals, blank_cps={k: 0.5 for k in signals},
                               noise_factor=0.01, seed=11)
    age = laicpms.u_pb_age(run)
    assert abs(age - 100.0) < 5.0
