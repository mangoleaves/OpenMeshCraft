import numpy as np
import time


def calculate_intersection_point_vector(a, b, p, q):
    n = np.cross(a[0:2] - p[0:2], q[0:2] - p[0:2])

    lambda_d = -np.cross(b[0:2] - a[0:2], q[0:2] - p[0:2])

    lambda_result = n * (b - a)

    return lambda_result, lambda_d


def calculate_intersection_point_number(a, b, p, q):
    xa, ya, za = a[0], a[1], a[2]
    xb, yb, zb = b[0], b[1], b[2]
    xp, yp, zp = p[0], p[1], p[2]
    xq, yq, zq = q[0], q[1], q[2]

    xap = xa - xp
    yap = ya - yp

    yqp = yq - yp
    xqp = xq - xp

    xba = xb - xa
    yba = yb - ya
    zba = zb - za

    c1 = xap * yqp
    c2 = xqp * yap
    n = c1 - c2

    c3 = xba * yqp
    c4 = xqp * yba
    lambda_d = c4 - c3

    lambda_x = n * xba
    lambda_y = n * yba
    lambda_z = n * zba

    return np.array([lambda_x, lambda_y, lambda_z]), lambda_d


# Generate random sets of values for xa, ya, za, ..., yt
num_samples = 10000
samples = np.random.rand(
    num_samples, 12
)  # Each row corresponds to one set of parameters

# Call calculate_intersection_point for each set of parameters
time_vector = 0.0
time_number = 0.0

for i in range(num_samples):
    a = samples[i, :3]
    b = samples[i, 3:6]
    p = samples[i, 6:9]
    q = samples[i, 9:12]
    time_start = time.time()
    v_lambda, v_lambda_d = calculate_intersection_point_vector(a, b, p, q)
    time_vector = time_vector + (time.time() - time_start)
    time_start = time.time()
    n_lambda, n_lambda_d = calculate_intersection_point_number(a, b, p, q)
    time_number = time_number + (time.time() - time_start)

    if (
        abs(v_lambda[0] - n_lambda[0]) > 1e-10
        or abs(v_lambda[1] - n_lambda[1]) > 1e-10
        or abs(v_lambda[2] - n_lambda[2]) > 1e-10
        or abs(v_lambda_d - n_lambda_d) > 1e-10
    ):
        print("detect difference.", v_lambda, v_lambda_d, n_lambda, n_lambda_d)

print(f"vector time: {time_vector}, number time: {time_number}")
