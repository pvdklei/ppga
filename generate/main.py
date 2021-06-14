import utils
import motor
import numpy as np
import sympy


def main():
    p = utils.trivec_sym()
    print(p.normalized())
    return
    l1 = utils.point(3, -2, 7) & utils.point(7, 1, -1)
    l1 = l1.normalized()
    l2 = utils.point(7, 4, 2) & utils.point(8, -3, -6)
    l2 = l2.normalized()
    m = utils.move_to(l1, l2)
    print(l2)
    print(utils.apply(m, l1))
    return

    m = utils.motor_sym()
    l = utils.line_sym()
    print(utils.apply(m, l))


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


if __name__ == "__main__":
    main()
