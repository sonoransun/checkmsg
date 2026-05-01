import numpy as np
import pytest

from checkmsg.epr import (
    Hyperfine,
    SpinSystem,
    _build_operators,
    _static_hamiltonian,
    spin_matrices,
)


def test_spin_half_matrices_are_pauli_over_two():
    Sx, Sy, Sz = spin_matrices(0.5)
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
    np.testing.assert_allclose(Sx, 0.5 * pauli_x, atol=1e-12)
    np.testing.assert_allclose(Sy, 0.5 * pauli_y, atol=1e-12)
    np.testing.assert_allclose(Sz, 0.5 * pauli_z, atol=1e-12)


@pytest.mark.parametrize("S", [0.5, 1.0, 1.5, 2.0, 2.5])
def test_spin_matrices_obey_commutator_algebra(S):
    Sx, Sy, Sz = spin_matrices(S)
    np.testing.assert_allclose(Sx @ Sy - Sy @ Sx, 1j * Sz, atol=1e-9)
    np.testing.assert_allclose(Sy @ Sz - Sz @ Sy, 1j * Sx, atol=1e-9)
    np.testing.assert_allclose(Sz @ Sx - Sx @ Sz, 1j * Sy, atol=1e-9)


@pytest.mark.parametrize("S", [0.5, 1.0, 1.5, 2.5])
def test_S_squared_eigenvalue(S):
    Sx, Sy, Sz = spin_matrices(S)
    S2 = Sx @ Sx + Sy @ Sy + Sz @ Sz
    expected = S * (S + 1)
    np.testing.assert_allclose(S2, expected * np.eye(int(2 * S + 1)), atol=1e-9)


def test_invalid_spin_raises():
    with pytest.raises(ValueError):
        spin_matrices(0.3)


def test_static_hamiltonian_is_hermitian():
    sys_ = SpinSystem(name="t", S=1.5, g=2.0, D_MHz=200.0, E_MHz=20.0,
                     hyperfine=(Hyperfine("14N", I=1.0, A_iso_MHz=50.0),))
    Sx, Sy, Sz, nuc_ops = _build_operators(sys_)
    H = _static_hamiltonian(sys_, Sx, Sy, Sz, nuc_ops)
    np.testing.assert_allclose(H, H.conj().T, atol=1e-12)


def test_dimension_matches_product_basis():
    sys_ = SpinSystem(name="t", S=2.5,
                     hyperfine=(Hyperfine("55Mn", I=2.5, A_iso_MHz=240.0),))
    Sx, _Sy, _Sz, _nuc = _build_operators(sys_)
    assert Sx.shape == (6 * 6, 6 * 6)
