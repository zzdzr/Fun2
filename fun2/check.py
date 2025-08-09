import warnings
from typing import Tuple, Callable

# Decorator for checking and converting variable types
def check_box_type(func: Callable) -> Callable:
    def wrapper(*args, **kwargs):
        # Extract self and other arguments
        self = args[0]
        center = kwargs.get("center")
        height = kwargs.get("height")
        angle = kwargs.get("angle")
        layer_height = kwargs.get("layer_height")


        # Check and convert 'center'
        if not isinstance(center, tuple) or len(center) != 2:
            raise TypeError("center must be a Tuple[float, float] (for opencv support)")
        
        if not all(isinstance(coord, float) for coord in center):
            warnings.warn("center should be a Tuple[float, float]. Converting to float (for opencv support).")
            center = tuple(float(coord) for coord in center)
        
        # Check and convert 'height'
        if not isinstance(height, float):
            warnings.warn("height should be a float. Converting to float (for opencv support).")

        # Check and convert 'angle'
        if not isinstance(angle, float):
            warnings.warn("angle should be a float. Converting to float (for opencv support).")
            angle = float(angle)

        # Check layer height
        if not isinstance(layer_height, int):
            raise TypeError("layer height should be int.")
        
        return func(*args, **kwargs)
    
    return wrapper

def check_apply_transform(func: Callable) -> Callable:
    def wrapper(*args, **kwargs):

        self = args[0]
        transform_type = kwargs.get("transform_type")
        n_theta = kwargs.get("n_theta", 0)
        n_translocation = kwargs.get("n_translocation")
        n_extension = kwargs.get("n_extension")
        n_expand = kwargs.get("n_expand")

        if not isinstance(transform_type, str):
            raise TypeError("Invalid type, transform_type should be a string.")

        for param, name in [(n_theta, "n_theta"), (n_translocation, "n_translocation"),
                            (n_extension, "n_extension"), (n_expand, "n_expand")]:

            if not isinstance(param, (float, int)):
                warnings.warn(f"{name} should be a float or int. Converting to float.")
                param = float(param)
        
        return func(self, transform_type = transform_type, n_theta = n_theta, n_translocation = n_translocation,
            n_extension = n_extension, n_expand = n_expand)
    
    return wrapper