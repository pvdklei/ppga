import utils
import sympy
import ganja


def main():
    outer_exp()


def bivec_decomp():
    b = utils.line_sym(e="be", v="bv")
    bdb = utils.scalar_sym("bdb")
    bmb = utils.pseudo_sym("bmb.0")
    eucl = b * (1 - 0.5 * bmb * (1/bdb[0]))
    van = 0.5 * b * bmb * (1 / bdb[0])
    print(eucl)
    print(van)


def outer_exp():
    a = utils.e_bivec_sym("a")
    b = utils.v_bivec_sym("b")
    s = a+b
    sds = utils.pseudo_sym("sds") # in fact 0.5 s ^ s
    sds = 0.5 * (s ^ s)
    exp = (1 + s + sds) * (1 / (1 - sds[0])) # euclidian part assumed normalized
                                   # so inverse == reverse
    print(exp)


def exp_4cs():
    eucl = utils.line_sym(e="eucl_e", v="eucl_v")
    van = utils.v_bivec_sym("van_v")
    chalfphi = utils.scalar_sym("cos_half_phi")
    shalfphi = utils.scalar_sym("sin_half_phi")
    exp = (1 + van) * (chalfphi - shalfphi * eucl)
    print(exp)


def exp_gunn():
    cu = utils.scalar_sym("cu")
    su = utils.scalar_sym("su")
    u = utils.scalar_sym("u")
    v = utils.pseudo_sym("v.0")
    b = utils.line_sym(e="be", v="bv")
    right = (v * cu + su) * utils.inverse(u+v) * b
    exp = cu-v*su + right
    print(exp)

    return
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
    w = utils.line_sym(e="we", v="wv")
    wdotwrev = utils.scalar_sym("wdotwrev")
    sqrt = utils.scalar_sym("sqrt_wdotwrev")
    wmeetwrev = utils.pseudo_sym("wmeetwrev.0")
    l = w * (1. / sqrt[0]) * (1 - 0.5 * wmeetwrev * (1 / wdotwrev[0]))
    a = utils.scalar_sym("a")
    bi = wmeetwrev * (1 / (2 * sqrt[0]))
    print((a+bi)*l)
    return


def log_gunn():
    u = utils.scalar_sym("u")
    v = utils.pseudo_sym("v")
    s = utils.scalar_sym("s2")
    p = utils.pseudo_sym("p2.0")
    b = utils.line_sym(e="be", v="bv")
    print((u+v)*b*utils.inverse(s+p))


if __name__ == "__main__":
    main()
