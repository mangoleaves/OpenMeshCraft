import numpy as np
import time


def calculate_intersection_point_vector(a, b, o, p, q):
    nor_opq = np.cross(p - o, q - o)

    lambda_d = np.dot(b - a, nor_opq)

    t = -np.dot(a - o, nor_opq)

    lambda_result = t * (b - a)

    return lambda_result, lambda_d


def calculate_intersection_point_number(a, b, o, p, q):
    xa, ya, za = a[0], a[1], a[2]
    xb, yb, zb = b[0], b[1], b[2]
    xo, yo, zo = o[0], o[1], o[2]
    xp, yp, zp = p[0], p[1], p[2]
    xq, yq, zq = q[0], q[1], q[2]

    xba = xb - xa
    yba = yb - ya
    zba = zb - za

    xpo = xp - xo
    ypo = yp - yo
    zpo = zp - zo

    xqo = xq - xo
    yqo = yq - yo
    zqo = zq - zo

    xoa = xo - xa
    yoa = yo - ya
    zoa = zo - za

    ypo_zqo = ypo * zqo
    zpo_yqo = zpo * yqo
    zpo_xqo = zpo * xqo
    xpo_zqo = xpo * zqo
    xpo_yqo = xpo * yqo
    ypo_xqo = ypo * xqo
    xnor = ypo_zqo - zpo_yqo
    ynor = zpo_xqo - xpo_zqo
    znor = xpo_yqo - ypo_xqo

    x_ba_nor = xba * xnor
    y_ba_nor = yba * ynor
    z_ba_nor = zba * znor
    temp_0 = x_ba_nor + y_ba_nor
    lambda_n = temp_0 + z_ba_nor

    x_oa_nor = xoa * xnor
    y_oa_nor = yoa * ynor
    z_oa_nor = zoa * znor
    temp_1 = x_oa_nor + y_oa_nor
    lambda_d = temp_1 + z_oa_nor

    lambda_x = lambda_n * xba
    lambda_y = lambda_n * yba
    lambda_z = lambda_n * zba

    return np.array([lambda_x, lambda_y, lambda_z]), lambda_d


# Generate random sets of values for xa, ya, za, ..., yt
num_samples = 10000
samples = np.random.rand(
    num_samples, 15
)  # Each row corresponds to one set of parameters

# Call calculate_intersection_point for each set of parameters
time_vector = 0.0
time_number = 0.0

for i in range(num_samples):
    a = samples[i, :3]
    b = samples[i, 3:6]
    o = samples[i, 6:9]
    p = samples[i, 9:12]
    q = samples[i, 12:]
    time_start = time.time()
    v_lambda, v_lambda_d = calculate_intersection_point_vector(a, b, o, p, q)
    time_vector = time_vector + (time.time() - time_start)
    time_start = time.time()
    n_lambda, n_lambda_d = calculate_intersection_point_number(a, b, o, p, q)
    time_number = time_number + (time.time() - time_start)

    if (
        abs(v_lambda[0] - n_lambda[0]) > 1e-10
        or abs(v_lambda[1] - n_lambda[1]) > 1e-10
        or abs(v_lambda[2] - n_lambda[2]) > 1e-10
        or abs(v_lambda_d - n_lambda_d) > 1e-10
    ):
        print("detect difference.", v_lambda, v_lambda_d, n_lambda, n_lambda_d)

print(f"vector time: {time_vector}, number time: {time_number}")
