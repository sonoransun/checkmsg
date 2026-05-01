import numpy as np

from checkmsg import uvvis
from checkmsg.refdata.chromophores import CHROMOPHORES, assign
from checkmsg.synthetic import PeakSpec, generate


def test_chromophore_table_has_known_entries():
    names = {c.name for c in CHROMOPHORES}
    assert any("Cr3+" in n for n in names)
    assert any("Fe2+/Ti4+" in n for n in names)


def test_assign_requires_all_bands_for_multi_band_chromophore():
    # Cr3+ ruby/spinel needs both 405 and 555 ± 20 nm.
    only_one = assign([405.0])
    assert not any("Cr3+ d-d (ruby" in c.name for _p, c in only_one)
    both = assign([405.0, 555.0])
    assert any("Cr3+ d-d (ruby" in c.name for _p, c in both)


def test_assign_finds_emerald_cr3_with_shifted_bands():
    # Emerald has bands at 430 + 605 (shifted vs ruby).
    out = assign([430.0, 605.0])
    assert any("emerald" in c.name for _p, c in out)


def test_uvvis_assign_bands_emerald_synthetic():
    axis = np.linspace(380, 800, 421)
    s = generate([PeakSpec(430, 0.85, 18, 8), PeakSpec(605, 0.95, 22, 10)],
                 axis, "uvvis", "nm", noise=0.005, seed=4)
    res = uvvis.assign_bands(s)
    chrom_names = [c.name for c in res.chromophores()]
    assert any("emerald" in n for n in chrom_names)


def test_uvvis_rejects_wrong_technique():
    import pytest

    from checkmsg.spectrum import Spectrum
    with pytest.raises(ValueError, match="expected uvvis"):
        uvvis.preprocess_uvvis(Spectrum(np.linspace(0, 1, 10), np.zeros(10), "raman", "cm-1"))
