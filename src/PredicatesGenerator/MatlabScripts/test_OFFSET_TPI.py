import numpy as np
import time


def calculate_intersection_point_vector(a, b, c, o, p, q, r, s, t):
    # Calculate normals
    nor_opq = np.cross(p - o, q - o)
    nor_rst = np.cross(s - r, t - r)

    # Calculate m
    m0 = np.dot(o - a, nor_opq)
    m1 = np.dot(r - a, nor_rst)

    # Calculate L matrix
    L = np.array(
        [
            [np.dot(nor_opq, b - a), np.dot(nor_opq, c - a)],
            [np.dot(nor_rst, b - a), np.dot(nor_rst, c - a)],
        ]
    )

    # Calculate determinant of L matrix
    det_L = np.linalg.det(L)

    # Calculate lambda
    det_L_sub0 = np.array([[m0, np.dot(nor_opq, c - a)], [m1, np.dot(nor_rst, c - a)]])
    det_L_sub1 = np.array([[np.dot(nor_opq, b - a), m0], [np.dot(nor_rst, b - a), m1]])

    # Calculate lambda
    lambda_u = np.linalg.det(det_L_sub0)
    lambda_v = np.linalg.det(det_L_sub1)
    lambda_result = lambda_u * (b - a) + lambda_v * (c - a)

    return lambda_result, det_L


def calculate_intersection_point_number(a, b, c, o, p, q, r, s, t):
    xa, ya, za = a[0], a[1], a[2]
    xb, yb, zb = b[0], b[1], b[2]
    xc, yc, zc = c[0], c[1], c[2]
    xo, yo, zo = o[0], o[1], o[2]
    xp, yp, zp = p[0], p[1], p[2]
    xq, yq, zq = q[0], q[1], q[2]
    xr, yr, zr = r[0], r[1], r[2]
    xs, ys, zs = s[0], s[1], s[2]
    xt, yt, zt = t[0], t[1], t[2]

    xpo = xp - xo
    ypo = yp - yo
    zpo = zp - zo
    xqo = xq - xo
    yqo = yq - yo
    zqo = zq - zo

    xsr = xs - xr
    ysr = ys - yr
    zsr = zs - zr
    xtr = xt - xr
    ytr = yt - yr
    ztr = zt - zr

    ypo_zqo = ypo * zqo
    zpo_yqo = zpo * yqo
    zpo_xqo = zpo * xqo
    xpo_zqo = xpo * zqo
    xpo_yqo = xpo * yqo
    ypo_xqo = ypo * xqo

    x_nor_opq = ypo_zqo - zpo_yqo
    y_nor_opq = zpo_xqo - xpo_zqo
    z_nor_opq = xpo_yqo - ypo_xqo

    ysr_ztr = ysr * ztr
    zsr_ytr = zsr * ytr
    zsr_xtr = zsr * xtr
    xsr_ztr = xsr * ztr
    xsr_ytr = xsr * ytr
    ysr_xtr = ysr * xtr

    x_nor_rst = ysr_ztr - zsr_ytr
    y_nor_rst = zsr_xtr - xsr_ztr
    z_nor_rst = xsr_ytr - ysr_xtr

    xoa = xo - xa
    yoa = yo - ya
    zoa = zo - za

    xra = xr - xa
    yra = yr - ya
    zra = zr - za

    xba = xb - xa
    yba = yb - ya
    zba = zb - za

    xca = xc - xa
    yca = yc - ya
    zca = zc - za

    temp_0 = xoa * x_nor_opq
    temp_1 = yoa * y_nor_opq
    temp_2 = zoa * z_nor_opq
    temp_3 = temp_0 + temp_1
    oa_dot_nor_opq = temp_2 + temp_3

    temp_4 = xra * x_nor_rst
    temp_5 = yra * y_nor_rst
    temp_6 = zra * z_nor_rst
    temp_7 = temp_4 + temp_5
    ra_dot_nor_rst = temp_6 + temp_7

    temp_8 = xba * x_nor_opq
    temp_9 = yba * y_nor_opq
    temp_10 = zba * z_nor_opq
    temp_11 = temp_8 + temp_9
    ba_dot_nor_opq = temp_10 + temp_11

    temp_12 = xca * x_nor_opq
    temp_13 = yca * y_nor_opq
    temp_14 = zca * z_nor_opq
    temp_15 = temp_12 + temp_13
    ca_dot_nor_opq = temp_14 + temp_15

    temp_16 = xba * x_nor_rst
    temp_17 = yba * y_nor_rst
    temp_18 = zba * z_nor_rst
    temp_19 = temp_16 + temp_17
    ba_dot_nor_rst = temp_18 + temp_19

    temp_20 = xca * x_nor_rst
    temp_21 = yca * y_nor_rst
    temp_22 = zca * z_nor_rst
    temp_23 = temp_20 + temp_21
    ca_dot_nor_rst = temp_22 + temp_23

    temp_100 = ba_dot_nor_opq * ca_dot_nor_rst
    temp_101 = ca_dot_nor_opq * ba_dot_nor_rst
    lambda_d = temp_100 - temp_101

    temp_102 = oa_dot_nor_opq * ca_dot_nor_rst
    temp_103 = ca_dot_nor_opq * ra_dot_nor_rst
    det_sub0 = temp_102 - temp_103

    temp_105 = ba_dot_nor_opq * ra_dot_nor_rst
    temp_106 = oa_dot_nor_opq * ba_dot_nor_rst
    det_sub1 = temp_105 - temp_106

    xu = det_sub0 * xba
    yu = det_sub0 * yba
    zu = det_sub0 * zba
    xv = det_sub1 * xca
    yv = det_sub1 * yca
    zv = det_sub1 * zca

    lambda_x = xu + xv
    lambda_y = yu + yv
    lambda_z = zu + zv

    return np.array([lambda_x, lambda_y, lambda_z]), lambda_d


# Generate random sets of values for xa, ya, za, ..., yt
num_samples = 10000
samples = np.random.rand(
    num_samples, 27
)  # Each row corresponds to one set of parameters

# Call calculate_intersection_point for each set of parameters
time_vector = 0.0
time_number = 0.0

for i in range(num_samples):
    a = samples[i, :3]
    b = samples[i, 3:6]
    c = samples[i, 6:9]
    o = samples[i, 9:12]
    p = samples[i, 12:15]
    q = samples[i, 15:18]
    r = samples[i, 18:21]
    s = samples[i, 21:24]
    t = samples[i, 24:]
    time_start = time.time()
    v_lambda, v_lambda_d = calculate_intersection_point_vector(
        a, b, c, o, p, q, r, s, t
    )
    time_vector = time_vector + (time.time() - time_start)
    time_start = time.time()
    n_lambda, n_lambda_d = calculate_intersection_point_number(
        a, b, c, o, p, q, r, s, t
    )
    time_number = time_number + (time.time() - time_start)

    if (
        abs(v_lambda[0] - n_lambda[0]) > 1e-10
        or abs(v_lambda[1] - n_lambda[1]) > 1e-10
        or abs(v_lambda[2] - n_lambda[2]) > 1e-10
        or abs(v_lambda_d - n_lambda_d) > 1e-10
    ):
        print("detect difference.", v_lambda, v_lambda_d, n_lambda, n_lambda_d)

print(f"vector time: {time_vector}, number time: {time_number}")