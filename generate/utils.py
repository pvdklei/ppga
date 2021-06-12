import ganja
import math
import sympy as sp


def inner(a, b):
    gp = a * b
    out = a ^ b
    return gp - out


def apply(m, g):
    return m * g * ~m


def rotor(angle, line):
    return math.cos(angle / 2.0) + math.sin(angle / 2.0) * line.normalized()


def translator(dist, line):
    return 1.0 + dist / 2.0 * line


E0 = ganja.PGA3D(1.0, 1)           # ideal plane
E1 = ganja.PGA3D(1.0, 2)           # x=0 plane
E2 = ganja.PGA3D(1.0, 3)           # y=0 plane
E3 = ganja.PGA3D(1.0, 4)           # z=0 plane


def plane(a, b, c, d):
    return a * E1 + b * E2 + c * E3 + d * E0


E123 = E1 ^ E2 ^ E3
E032 = E0 ^ E3 ^ E2
E013 = E0 ^ E1 ^ E3
E021 = E0 ^ E2 ^ E1


def point(x, y, z):
    return E123 + -x * E032 + -y * E013 + -z * E021


def dir(x, y, z):
    return -x * E032 + -y * E013 + -z * E021


def inverse(self):
    rev = ~self
    fac = (self * rev)[0]
    return rev * (1/fac)


def normalized_inverse(v):
    return ~v


def sqrt(self):
    scalar = self[self._base.index("1")]
    left_num = 1 + self
    left_den = sp.sqrt(2 * (1 + scalar))
    right_num = self[self._base.index("e0123")]
    right_den = 2 * (1 + scalar)
    return left_num * (1 / left_den) * (1 - right_num * (1 / right_den))


def ssqrt(self):
    return (1 + self).normalized()


def linenorm(l):
    a = l[l._base.index("e01")]
    b = l[l._base.index("e02")]
    c = l[l._base.index("e03")]
    import numpy as np
    return np.sqrt(a**2+b**2+c**2)


def multivec(a: list) -> ganja.PGA3D:
    """Creates a multivec filled with sympy symbols or numbers"""

    multivec_a = ganja.PGA3D(0, 0)
    for (name, unite) in a:
        if name == 0 or name == "" or name == None or not name:
            continue
        varsym = 0.0
        if type(name) is str:
            varsym = sp.symbols(name)
        elif type(name) is int:
            varsym = name
        multivec_a += ganja.PGA3D(varsym, ganja.PGA3D._base.index(unite))
    return multivec_a


def trivec_sym(n): return multivec([(f"{n}[0]", "e123"), (f"{n}[1]", "e032"),
                                    (f"{n}[2]", "e013"), (f"{n}[3]", "e021")])


def point_sym(n): return multivec([(1, "e123"), (f"{n}[1]", "e032"),
                                   (f"{n}[2]", "e013"), (f"{n}[3]", "e021")])


def dir_sym(n): return multivec([(f"{n}[1]", "e032"),
                                 (f"{n}[2]", "e013"), (f"{n}[3]", "e021")])


def trivec_inf_sym(n): return multivec([(f"{n}[1]", "e032"),
                                        (f"{n}[2]", "e013"), (f"{n}[3]", "e021")])


def vec_sym(n): return multivec([(f"{n}[0]", "e0"), (f"{n}[1]", "e1"),
                                 (f"{n}[2]", "e2"), (f"{n}[3]", "e3")])


def e_bivec_sym(n): return multivec(
    [(f"{n}[0]", "e23"), (f"{n}[1]", "e31"), (f"{n}[2]", "e12")])


def v_bivec_sym(n): return multivec(
    [(f"{n}[0]", "e01"), (f"{n}[1]", "e02"), (f"{n}[2]", "e03")])


def motor_sym(s="s", ps="ps", e="e", v="v"):
    return scalar_sym(s) + pseudo_sym(ps) + v_bivec_sym(v) + e_bivec_sym(e)


def rotor_sym(s="s", e="e"):
    return scalar_sym(s) + e_bivec_sym(e)


def simple_motor_sym(s, e, v):
    return scalar_sym(s) + v_bivec_sym(v) + e_bivec_sym(e)


def line_sym(v="v", e="e"):
    return e_bivec_sym(e) + v_bivec_sym(v)


def translator_sym(s="s", v="v"):
    return scalar_sym(s) + v_bivec_sym(v)


def p1_sym(n="p1."):
    return multivec([(f"{n}x", "1"),
                     (f"{n}y", "e23"),
                     (f"{n}z", "e31"),
                     (f"{n}w", "e12")])


def p2_sym(n="p2."):
    return multivec([(f"{n}x", "e0123"),
                     (f"{n}y", "e01"),
                     (f"{n}z", "e02"),
                     (f"{n}w", "e03")])


def dir_sym_euc(n="d."): return -1 * multivec([(f"{n}x", "e032"),
                                               (f"{n}y", "e013"), (f"{n}z", "e021")])


def scalar_sym(n="s"):
    return multivec([(n, "1")])


def pseudo_sym(n="ps"):
    return multivec([(n, "e0123")])
