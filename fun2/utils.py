import numba
import cv2
import numpy as np

from typing import Tuple, Any
from scipy.spatial.distance import cdist


@numba.njit
def rotate_point(x: float, y: float, cx: float, cy: float, angle_rad: float) -> Tuple[float, float]:
    """Rotate a point (x, y) around center (cx, cy) by a given angle in radians.

    Args:
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        cx (float): x-coordinate of the center.
        cy (float): y-coordinate of the center.
        angle_rad (float): Rotation angle in radians.

    Returns:
        Tuple[float, float]: The rotated point's (x, y) coordinates.
    """
    x_shifted = x - cx
    y_shifted = y - cy

    x_rot = x_shifted * np.cos(angle_rad) - y_shifted * np.sin(angle_rad)
    y_rot = x_shifted * np.sin(angle_rad) + y_shifted * np.cos(angle_rad)

    return x_rot + cx, y_rot + cy

@numba.njit
def inside_rotated_box(x: float, y: float, rec_points: np.ndarray) -> bool:
    """Check if a point (x, y) is inside a rotated rectangle defined by 4 points.

    Args:
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        rec_points (np.ndarray): An array of shape (4,2) representing the rectangle's corners.

    Returns:
        bool: True if the point is inside the rectangle, False otherwise.
    """
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
def compute_mean_intensity(image: np.ndarray, box: np.ndarray) -> float:
    """Compute the mean intensity of pixels inside a rotated rectangle.

    Args:
        image (np.ndarray): 2D image array.
        box (np.ndarray): An array of shape (4,2) representing the rectangle's corners.

    Returns:
        float: The mean intensity within the box, or np.nan if no pixels are found.
    """
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

@numba.njit
def compute_sum_intensity(image: np.ndarray, box: np.ndarray) -> float:
    """Compute the total intensity (sum of pixel values) inside a rotated rectangle.

    Args:
        image (np.ndarray): 2D image array.
        box (np.ndarray): An array of shape (4,2) representing the rectangle's corners.

    Returns:
        float: The sum of pixel intensities within the box.
    """
    img_h, img_w = image.shape[:2]
    box[0,] = (box[0,] + box[1,]) / 2
    box[3,] = (box[2,] + box[3,]) / 2

    x_min, x_max = min(box[:, 0]), max(box[:, 0])
    y_min, y_max = min(box[:, 1]), max(box[:, 1])

    x_min, x_max = max(x_min, 0), min(x_max, img_w - 1)
    y_min, y_max = max(y_min, 0), min(y_max, img_h - 1)

    pixel_sum = 0.0

    for x in range(x_min, x_max + 1):
        for y in range(y_min, y_max + 1):
            if inside_rotated_box(x, y, box):
                val = image[y, x]
                if np.isnan(val):
                    continue
                pixel_sum += val

    return pixel_sum
    
def generate_box(center: tuple, width: float, height: float, angle: float) -> np.ndarray:
    """Generates a rotated rectangular box with given center, width, height, and angle."""
    rect = (
        (float(center[0]), float(center[1])),
        (float(width), float(height)), float(angle)
    )
    # only half of box
    box = cv2.boxPoints(rect).astype(np.intp)
    box[0] = (box[0] + box[1])//2
    box[3] = (box[2] + box[3])//2

    return box

def adjust_p_values(df):
    # extract p-values
    p_up, p_down = df['pval_upstream'].values, df['pval_downstream'].values
    try:
        p_up[np.isnan(p_up)] = 1
        p_down[np.isnan(p_down)] = 1

    except TypeError:
        p_up, p_down = np.asarray(p_up, dtype = np.float64), \
        np.asarray(p_down, dtype = np.float64)

    adj_p_up, adj_p_down = _adjust_p_values(p_up), _adjust_p_values(p_down)

    df['FDR_upstream'] = adj_p_up
    df['FDR_downstream'] = adj_p_down
    
    return df
    
def _adjust_p_values(p_vals):
    """adjust p-value using Benjamini-Hochberg (BH)"""

    n = len(p_vals)
    sorted_indices = np.argsort(p_vals)
    sorted_p_values = p_vals[sorted_indices]
    
    adjusted_p_values = np.zeros(n)
    for i in range(n):
        adjusted_p_values[i] = sorted_p_values[i] * n / (i + 1)
        
    adjusted_p_values = np.minimum.accumulate(adjusted_p_values[::-1])[::-1]
    adjusted_p_values = np.minimum(adjusted_p_values, 1)
    
    return adjusted_p_values[
        np.argsort(sorted_indices)
    ]
