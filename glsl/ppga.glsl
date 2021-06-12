#ifndef PPGA_GLSL
#define PPGA_GLSL

// PGA4CS blocks-of-four memory layout
//
// p0 -> (e0, e1, e2, e3)
// p1 -> (1, e23, e31, e12)
// p2 -> (e0123, e01, e02, e03)
// p3 -> (e123, e032, e013, e021)

struct ppga_rotor
{
    vec4 p1;
};

struct ppga_motor
{
    vec4 p1;
    vec4 p2;
};

vec3 ppga_apply_rotor_to_direction(ppga_rotor r, vec3 d) {
    vec4 mask1 = vec4(1, 1, -1, -1);
    vec4 mask2 = vec4(1, 1, 1, -1);

    vec3 res = 2.0 * vec3(dot(d.yyzz * r.p1.wywx * r.p1.xzyz, mask2),
                          dot(d.xzzx * r.p1.ywxw * r.p1.zzyx, mask2),
                          dot(d.xxyy * r.p1.wxwx * r.p1.yzzy, mask2));
    vec4 rsq = r.p1 * r.p1;
    res += vec3(dot(mask1, d.x * rsq),
                dot(mask1, d.y * rsq.xzyw),
                dot(mask1, d.z * rsq.xwyz));
    return res;
}

vec3 ppga_apply_motor_to_direction(ppga_motor m, vec3 dir) {
    ppga_rotor r = ppga_rotor(m.p1);
    return ppga_apply_rotor_to_direction(r, dir);
}

vec3 ppga_apply_motor_to_origin(ppga_motor m) {
    vec3 res = m.p1.yzw * m.p2.x;
    res += m.p1.x * m.p2.yzw;
    res -= m.p1.zwy * m.p2.wyz;
    res += m.p1.wyz * m.p2.zwy;
    return 2.0 * res;
}

#endif