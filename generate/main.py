import utils
import motor
import numpy as np
import sympy


def exp_ganja():
    b = utils.line_sym(e="self.e_bivector")
    sdbb = utils.scalar_sym("sdbb")
    csdbb = utils.scalar_sym("csdbb")
    ssdbb = utils.scalar_sym("ssdbb")
    mbb = utils.pseudo_sym("mbb.0")
    left = csdbb
    middle = mbb * (1 / (2 * sdbb[0]) * ssdbb)
    right = (ssdbb - mbb * (1 / (2 * sdbb[0]))
             * csdbb) * b * utils.inverse(sdbb - mbb * (1 / 2 * sdbb[0]))

    exp = left + middle + right
    print(exp)


def log_4cs():
    # a = utils.scalar_sym("a")
    # bi = utils.pseudo_sym("bi")
    # l = utils.line_sym(e="le", v="lv")
    # print((a + bi) * l)

    w = utils.line_sym(e="we", v="wv")
    wdotwrev = utils.scalar_sym("wdotwrev")
    sqrt = utils.scalar_sym("sqrt_wdotwrev")
    wmeetwrev = utils.pseudo_sym("wmeetwrev.0")
    l = w * (1. / sqrt[0]) * (1 - 0.5 * wmeetwrev * (1 / wdotwrev[0]))
    a = utils.scalar_sym("a")
    bi = wmeetwrev * (1 / (2 * sqrt[0]))
    print((a+bi)*l)
    return


def main():

    t = utils.translator_sym("ts", "tv")
    r = utils.rotor_sym("rs", "re")
    m1 = utils.motor_sym("ms", "mps", "me", "mv")
    m2 = utils.motor_sym("s2", "ps2", "e2", "v2")
    p1 = utils.trivec_sym("p")
    p2 = utils.trivec_sym("p2")

    p0 = utils.point(0, 0, 0)
    e1 = utils.point(1, 0, 0)
    e2 = utils.point(0, 1, 0)
    p1 = np.array([2, 4, -3])
    p2 = np.array([9, 1, 0])

    n = np.linalg.norm(p2)
    p1 = utils.point(*(p1 / n))

    # p1 = utils.point(2, 4, -3)

    fromm = p0 & e1
    to = p0 & p1

    r = utils.sqrt(to * utils.inverse(fromm))

    print(r)


if __name__ == "__main__":
    main()
