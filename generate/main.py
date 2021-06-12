import utils
import motor
import numpy as np
import sympy


def main():

    # a = utils.scalar_sym("a")
    # bi = utils.pseudo_sym("bi")
    # l = utils.line_sym(e="le", v="lv")
    # print((a + bi) * l)

    # return

    # l1 = utils.line_sym(e="e1", v="v1")
    # l2 = utils.line_sym(e="e2", v="v2")
    # print(l1 * l2)

    # return

    w = utils.line_sym(e="we", v="wv")
    wdotwrev = utils.scalar_sym("wdotwrev")
    sqrt = utils.scalar_sym("sqrt_wdotwrev")
    wmeetwrev = utils.pseudo_sym("wmeetwrev.0")
    l = w * (1. / sqrt[0]) * (1 - 0.5 * wmeetwrev * (1 / wdotwrev[0]))
    a = utils.scalar_sym("a")
    bi = wmeetwrev * (1 / (2 * sqrt[0]))
    print((a+bi)*l)
    return

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

    return

    rot = utils.rotor(2.3, utils.point(2, 1, 8) & utils.point(8, 2, 3))
    trans = utils.translator(-5.6, utils.point(-4, 5, 4)
                             & utils.point(-4, -2, -1))
    mot = rot
    # for testing
    a = utils.point(2, -3.5, 1)
    a_ = utils.apply(mot, a)
    b = utils.point(-9, 4.3, 0)
    b_ = utils.apply(mot, b)
    c = utils.point(-7, 1.2, 0)
    c_ = utils.apply(mot, c)

    # p = (a & b).normalized()
    # p_ = (c_ & a_).normalized()
    # p = p * (1/utils.linenorm(p))
    # p_ = p_ * (1/utils.linenorm(p_))
    # p = utils.plane(3, 2, 4, 5)
    # p_ = utils.plane(-3, 9, 5, 5)
    # p = utils.point(3, 2, 4)
    # p_ = utils.point(-3, 9, 5)
    # m = utils.ssqrt(p_.normalized() * utils.inverse(p.normalized()))
    # res = m * p * ~m

    # print(p)
    # print(p.normalized())
    # print(p_)
    # print(res)

    # motor.from_point_correspondences(utils.trivec_inf_sym("a"),
    # utils.trivec_inf_sym("a_"),
    # utils.trivec_inf_sym("b"),
    # utils.trivec_inf_sym("b_"),
    # utils.trivec_inf_sym("c"),
    # utils.trivec_inf_sym("c_"))
    motor.from_point_correspondences(utils.point(1, 0, 0),
                                     utils.trivec_sym("a_"),
                                     utils.point(0, 1, 0),
                                     utils.trivec_sym("b_"),
                                     utils.point(0, 0, 1),
                                     utils.trivec_sym("c_"))
    # motor.from_point_correspondences(a, a_, b, b_, c, c_)


if __name__ == "__main__":
    main()
