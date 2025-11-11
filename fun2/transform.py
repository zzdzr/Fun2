"""
Module for transforming genuine length/width (kb) to pixel length/width.

This module provides a function to convert a given genuine length/width (in kilobases)
to its corresponding pixel length/width by accounting for a specified angle. The transformation
uses the formula:

    pixel_length = 2 * genuine_length(width) / cos(angle in radians)
"""

import numpy as np
from typing import Union
from fun2.planning.sampling_box import SamplingBoxAgent
# def transform(angle: float, resolution: int, genuine_length: Union[int, float]) -> float:
#     """Converts a genuine length/width (kb) to pixel length/width based on the given angle.

#     The conversion is performed using the formula:
    
#         pixel_length = 2 * genuine_length / cos(angle in radians)

#     Args:
#         angle (float): The angle in degrees.
#         genuine_length (int or float): The genuine length/width in kilobases.

#     Returns:
#         float: The corresponding pixel length/width.
#     """
#     print(f"genuine_length is {genuine_length}.")
#     # genuine_length = np. 

#     if angle == 90 or angle == 0:
#         pixel_length = 2 * genuine_length * 1000 / resolution
#     else:
#         rad = np.deg2rad(90-angle)
#         pixel_length = 2 * genuine_length * 1000 / (np.cos(rad) * resolution)

#     return pixel_length



# def pixel_to_genuine(angle: float, resolution: int, pixel_length: Union[int, float]) -> int:
#     """Converts a pixel length/width (kb) to genuine length/width based on the given angle.

#     The conversion is performed using the formula:
    
#         genuine_length = 0.5 * pixel_length * cos(angle in radians)

#     Args:
#         angle (float): The angle in degrees.
#         pixel_length (int or float): The corresponding pixel length/width.

#     Returns:
#         int: The genuine length/width in kilobases.
#     """
#     if angle == 90 or angle == 0:
#         genuine_length = 0.5 * pixel_length * resolution / 1000
#     else:
#         rad = np.deg2rad(90-angle)
#         genuine_length = 0.5 * pixel_length * np.cos(rad) * resolution / 1000

#     return genuine_length



# def pixel_to_genuine(angle: float, resolution: int, pixel_length: Union[int, float]) -> int:
#     """Converts a pixel length/width (kb) to genuine length/width based on the given angle.

#     The conversion is performed using the formula:
    
#         genuine_length = 0.5 * pixel_length * cos(angle in radians)

#     Args:
#         angle (float): The angle in degrees.
#         pixel_length (int or float): The corresponding pixel length/width.

#     Returns:
#         int: The genuine length/width in kilobases.
#     """
#     rad = np.deg2rad(angle)
    
#     genuine_length = 0.5 * pixel_length * (np.cos(rad) + np.sin(rad)) * resolution / 1000

#     return genuine_length


# def transform(angle: float, resolution: int, genuine_length: Union[int, float]) -> float:
#     """Converts a genuine length/width (kb) to pixel length/width based on the given angle.

#     The conversion is performed using the formula:
    
#         pixel_length = 2 * genuine_length / cos(angle in radians)

#     Args:
#         angle (float): The angle in degrees.
#         genuine_length (int or float): The genuine length/width in kilobases.

#     Returns:
#         float: The corresponding pixel length/width.
#     """
#     rad = np.deg2rad(angle)

#     pixel_length = 2 * genuine_length / (np.cos(rad) + np.sin(rad)) / resolution * 1000

#     return pixel_length


def mapping_coverage(agent: SamplingBoxAgent, bin_start_adjusted: int, resolution: int):
    
    box = agent.generate_box()
    h, _ = agent._image.shape[0], agent._image.shape[1]

    # four points of sampling box
    p1, p4 = box[0,], box[3,]
    p2, p3 = box[1,], box[2,]
    anchor = (p2 + p3) // 2

    elongation_down, elongation_up = min(h + bin_start_adjusted, anchor[0] + bin_start_adjusted), \
    max(0, anchor[1] + bin_start_adjusted)

    width_up, width_down = min(h + bin_start_adjusted - 1, max(0, (p1[0] + p1[1]) // 2 + bin_start_adjusted)), \
    min(h + bin_start_adjusted, (p4[0] + p4[1]) // 2 + bin_start_adjusted)

    return elongation_up * resolution, elongation_down * resolution, width_up * resolution, width_down * resolution

def transform(resolution: int, genuine_length: Union[int, float], type_: str = 'height') -> float:

    if type_ == 'height':
        pixel_length = genuine_length * 2 * np.sqrt(2) * 1000 / resolution # Should be genome length
    else:
        pixel_length = genuine_length * np.sqrt(2) * 1000 / resolution

    return pixel_length

def pixel_to_genuine(resolution: int, pixel_length: Union[int, float], type_: str = 'height') -> int:

    if type_ == 'height':
        genuine_length = pixel_length * resolution / 2 / 1000
    else:
        genuine_length = pixel_length * resolution / 1000

    return genuine_length