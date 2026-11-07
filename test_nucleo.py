"""Testes do núcleo numérico:  python -m pytest -q"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modulos import (rep_analise, rep_canonical_ss, rep_canonical_tf, rep_discretize,  # noqa: E402
                     rep_feedback, rep_observer, rep_servo)
from nucleo import (ackermann_k, ackermann_l, discretize, matrix_exp, parse_complex_list,  # noqa: E402
                    safe_eval, servo_alt, servo_continuous_ufrn, servo_discrete_ufrn)

P = parse_complex_list


def poles(M):
    return np.sort_complex(np.round(np.linalg.eigvals(M), 4))


def test_parsers():
    assert P("-2+1i, -2-1i, -10") == [-2 + 1j, -2 - 1j, -10 + 0j]
    assert P("-1±0.1i, -5") == [-1 + 0.1j, -1 - 0.1j, -5 + 0j]
    assert safe_eval("1/3") == pytest.approx(1 / 3)
    assert safe_eval("0,5") == 0.5
    assert safe_eval("2^3") == 8
    with pytest.raises(ValueError):
        safe_eval("__import__('os')")


def test_matrix_exp_vs_scipy():
    scipy_linalg = pytest.importorskip("scipy.linalg")
    M = np.array([[0, 1, 0.3], [2, -1, 0.5], [0, 0, 0]]) * 1.7
    assert np.allclose(matrix_exp(M), scipy_linalg.expm(M), atol=1e-10)


def test_discretize_double_integrator():
    d = discretize(np.array([[0, 1], [0, 0.]]), np.array([[0], [1.]]), 1.0)
    assert np.allclose(d["G"], [[1, 1], [0, 1]])
    assert np.allclose(d["H"], [[0.5], [1]])


def test_ackermann_places_poles():
    A = np.array([[-4, 3, 0], [-1, -4, -1], [0, 0, -2.]])
    B = np.array([[0], [0], [3.]])
    C = np.array([[1, 0, 0.]])
    d = P("-2+1i, -2-1i, -10")
    assert np.allclose(poles(A + B @ ackermann_k(A, B, d)["K"]), np.sort_complex(np.array(d)))
    assert np.allclose(poles(A - ackermann_l(A, C, d)["L"] @ C), np.sort_complex(np.array(d)))


def test_uncontrollable_reports_error():
    assert "error" in ackermann_k(np.eye(2), np.array([[1], [0.]]), P("-1,-2"))


@pytest.mark.parametrize("G,H,C,desired,k1,k2", [
    ([[1, 0.6321], [0, 0.3679]], [[0.3679], [0.6321]], [[-2, 3]], "0.5+0.2i, 0.5-0.2i, -0.1",
     [-2.1576, -1.5298], 0.2523),
    ([[1, 0, 2], [0, 1, 0], [2, 3, -1]], [[-1], [2], [0]], [[1, 0, 1]],
     "0.5+0.5i, 0.5-0.5i, -0.5+0.5i, -0.5-0.5i", [-1.3698, -1.7370, -1.4219], -0.1042),
])
def test_discrete_servo_ufrn_matches_formulario(G, H, C, desired, k1, k2):
    d = P(desired)
    r = servo_discrete_ufrn(G, H, C, d)
    assert np.allclose(r["k1"][0], k1, atol=2e-4)
    assert r["k2"] == pytest.approx(k2, abs=2e-4)
    assert np.allclose(poles(r["Acl"]), np.sort_complex(np.array(d)), atol=1e-3)


def test_continuous_servo_both_conventions_reach_desired_poles():
    A = np.array([[0, 1, 0], [0, 0, 1], [-5, -7, -3.]])
    B = np.array([[0], [0], [1.]])
    C = np.array([[-2, 4, 3.]])
    d = P("-1+0.1i, -1-0.1i, -5, -12")
    for r in (servo_alt(A, B, C, d), servo_continuous_ufrn(A, B, C, d)):
        assert np.allclose(poles(r["Acl"]), np.sort_complex(np.array(d)), atol=1e-3)


def test_discrete_classic_servo_has_true_integrator():
    G = np.array([[1, 0.6321], [0, 0.3679]])
    H = np.array([[0.3679], [0.6321]])
    C = np.array([[-2, 3.]])
    r = servo_alt(G, H, C, P("0.5+0.2i, 0.5-0.2i, -0.1"), discrete=True)
    assert r["Aa"][0, 0] == 1.0           # integrador discreto: ξ(k+1) = ξ(k) + ...
    assert np.allclose(poles(r["Acl"]), np.sort_complex(np.array(P("0.5+0.2i, 0.5-0.2i, -0.1"))), atol=1e-3)


def test_all_reports_build_without_error():
    A = np.array([[-4, 3, 0], [-1, -4, -1], [0, 0, -2.]])
    B = np.array([[0], [0], [3.]])
    C = np.array([[1, 0, 0.]])
    reps = [rep_analise(A, B, C), rep_feedback(A, B, np.array(P("-2±1i, -10"))),
            rep_observer(A, C, np.array(P("-1±1i, -10"))), rep_discretize(A, B, 0.1),
            rep_canonical_tf([1, 3, 7, 5], [3, 4, -2]), rep_canonical_ss(A, B, C),
            rep_servo(A, B, C, np.array(P("-1±0.1i, -5, -12")))]
    for r in reps:
        assert r.error is None and r.steps
