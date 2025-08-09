import numba

@numba.njit
def rotate_pont(x, y, cx, cy, angle_rad):
    x_shifted = x - cx
    y_shifted = y - cy

    x_rot = x_shifted * np.cos(angle_rad) - y_shifted * np.sin(angle_rad)
    y_rot = x_shifted * np.sin(angle_rad) + y_shifted * np.cos(angle_rad)

    return x_rot + cx, y_rot + cy

@numba.njit
def inside_rotated_box(x, y, rec_points):
    p1, p2, p3, p4 = rec_points

    def sign(point, line_start, line_end):
        return (
            (line_end[0] - line_start[0]) * (point[1] - line_start[1]) -
            (line_end[1] - line_start[1]) * (point[0] - line_start[0])
        )

    b1 = sign((x, y), p1, p2) < 0.0
    b2 = sign((x, y), p2, p3) < 0.0
    b3 = sign((x, y), p3, p4) < 0.0
    b4 = sign((x, y), p4, p1) < 0.0

    return ( (b1 == b2) and (b2 == b3) and (b3 == b4)  )

@numba.njit
def compute_mean_intensity(image, box):
    img_h, img_w = image.shape[:2]

    x_min, x_max = min(box[:, 0]), max(box[:, 0])
    y_min, y_max = min(box[:, 1]), max(box[:, 1])

    x_min, x_max = max(x_min, 0), min(x_max, img_w - 1)
    y_min, y_max = max(y_min, 0), min(y_max, img_h - 1)

    pixel_sum = 0.0
    count = 0

    for x in range(x_min, x_max + 1):
        for y in range(y_min, y_max + 1):
            if inside_rotated_box(x, y, box):
                pixel_sum += image[y, x]
                count += 1

    if count == 0:
        return np.nan

    return pixel_sum / count