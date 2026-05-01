import numpy as np

from checkmsg.muon.source import MuonSource


def test_default_construction():
    src = MuonSource()
    assert src.mean_momentum_MeV > 0
    assert src.flux_per_s > 0


def test_sample_returns_correct_shapes():
    src = MuonSource(mean_momentum_MeV=300.0, momentum_FWHM_MeV=10.0)
    momenta, dirs = src.sample(1000, rng=np.random.default_rng(0))
    assert momenta.shape == (1000,)
    assert dirs.shape == (1000, 3)


def test_sampled_directions_are_unit_vectors():
    src = MuonSource()
    _, dirs = src.sample(500, rng=np.random.default_rng(0))
    norms = np.linalg.norm(dirs, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-10)


def test_sampled_momenta_match_mean_within_tolerance():
    src = MuonSource(mean_momentum_MeV=300.0, momentum_FWHM_MeV=10.0)
    momenta, _ = src.sample(5000, rng=np.random.default_rng(0))
    assert abs(momenta.mean() - 300.0) < 1.0


def test_angular_spread_within_tolerance():
    src = MuonSource(angular_spread_mrad=5.0,
                     direction=(0.0, 0.0, -1.0))
    _, dirs = src.sample(2000, rng=np.random.default_rng(0))
    # The deviation from the nominal direction (-z) should have RMS angle ≈ 5 mrad
    nominal = np.array([0.0, 0.0, -1.0])
    cos_th = np.clip(dirs @ nominal, -1.0, 1.0)
    theta_mrad = np.arccos(cos_th) * 1000.0
    rms = float(np.sqrt(np.mean(theta_mrad ** 2)))
    # Two perpendicular Gaussian deviates → expected RMS ≈ sqrt(2) * 5 mrad
    assert 5.0 < rms < 12.0


def test_expected_count_proportional_to_exposure():
    src = MuonSource(flux_per_s=1e9)
    assert src.expected_count(1.0) == 1_000_000_000
    assert src.expected_count(0.5) == 500_000_000
