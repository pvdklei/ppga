import utils


def from_point_correspondences(a, a_, b, b_, c, c_):

    Va = utils.sqrt(a_ * utils.inverse(a))

    print("1")

    Ba = Va * b * ~Va
    Ba_sym = utils.simple_motor_sym("ba_s", "ba_e", "ba_v")
    to = (a_ & b_)
    from_ = (a_ & Ba_sym)
    Vb = utils.sqrt(to * utils.inverse(from_))

    print("2")

    Cba = Vb * Va * c * ~Va * ~Vb
    Cba_sym = utils.simple_motor_sym("cba_s", "cba_e", "cba_v")
    to = a_ & b_ & c_
    from_ = a_ & b_ & Cba_sym

    ratio = to * utils.inverse(from_)
    print("2.5")
    Vc = utils.sqrt(ratio)

    print("3")

    V = Vc * Vb * Va

    print("Va: ", Va)
    print("Vb: ", Vb)
    print("Vc: ", Vc)

    a__ = V * a * ~V
    b__ = V * b * ~V
    c__ = V * c * ~V
    # print(a_)
    # print(a__)
    # print(b_)
    # print(b__)
    # print(c_)
    # print(c__)

    print("\n\n V ======== \n\n")
    print(V)
