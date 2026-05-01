import pytest

from checkmsg import minerals


def test_catalog_size():
    assert len(minerals.CATALOG) >= 40, "catalog must have ≥40 entries"


def test_get_canonical_name():
    profile = minerals.get("ruby")
    assert profile.name == "ruby"
    assert profile.species == "corundum"


def test_get_resolves_aliases():
    # 'imperial jade' is an alias of jadeite
    p = minerals.get("imperial jade")
    assert p.name == "jadeite"


def test_get_unknown_raises():
    with pytest.raises(KeyError):
        minerals.get("kryptonite")


def test_by_species_garnets():
    garnets = minerals.by_species("garnet")
    assert "pyrope" in garnets
    assert "almandine" in garnets
    assert "spessartine" in garnets


def test_by_color_blue_set():
    blues = minerals.by_color("blue")
    assert "sapphire_blue" in blues
    assert "tanzanite" in blues
    assert "aquamarine" in blues


def test_by_confusable_finds_listing_minerals():
    """Whoever lists 'ruby' as a confusable should appear here."""
    listing = minerals.by_confusable("ruby")
    assert "red_spinel" in listing


def test_resolve_confusables_returns_profiles():
    out = minerals.resolve_confusables("ruby")
    # Profiles, not names
    assert all(isinstance(v, minerals.MineralProfile) for v in out.values())
    assert "red_spinel" in out


def test_every_profile_has_at_least_one_diagnostic_signal():
    for name, p in minerals.CATALOG.items():
        signal_count = (len(p.raman_peaks_cm) + len(p.uvvis_bands_nm)
                        + len(p.xrf_signature) + len(p.libs_signature)
                        + len(p.epr_centers))
        assert signal_count > 0, f"{name} has no diagnostic data"


def test_every_profile_has_at_least_one_diagnostic_feature_string():
    for name, p in minerals.CATALOG.items():
        assert p.diagnostic_features, f"{name} missing diagnostic_features"


def test_every_confusable_resolves():
    """Every confusable listed in any profile must itself be in the catalog."""
    missing: list[tuple[str, str]] = []
    for name, p in minerals.CATALOG.items():
        for conf in p.confusables:
            try:
                minerals.get(conf)
            except KeyError:
                missing.append((name, conf))
    assert not missing, f"unresolvable confusables: {missing}"


def test_names_returns_sorted():
    nm = minerals.names()
    assert nm == sorted(nm)
