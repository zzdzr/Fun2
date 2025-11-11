import logging
import numpy as np
from scipy.optimize import minimize
from scipy.stats import gamma, wilcoxon
from typing import List, Dict, Tuple, Union, Any

from fun2.utils import rotate_point, compute_mean_intensity, compute_sum_intensity, inside_rotated_box # type: ignore


# Setup global logging configuration
logging.basicConfig(
    level = logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# -----------------------------------------------------------------------------
# SamplingBoxAgent: Extended class for sampling box operations.
# -----------------------------------------------------------------------------
class SamplingBoxAgent:
    """Agent for handling sampling box operations on images.

    This class extends a base SamplingBoxAgent and provides additional methods
    for computing intensities, performing geometric transformations, and calculating
    rewards based on signal intensities within a rotated sampling box.
    """
    def __init__(self, image: np.ndarray, center: Tuple[float, float], width: int,
                 height: int, angle: int = 0, edge_width: int = 8, layer_height: int = 10, random_seed = 1):
        """Initializes the SamplingBoxAgent.

        Args:
            image (np.ndarray): The input image (2D array).
            center (Tuple[float, float]): The center coordinates of the sampling box.
            width (int): The width of the sampling box.
            height (int): The height of the sampling box.
            angle (int, optional): The rotation angle in degrees. Defaults to 0.
            layer_height (int, optional): The height of each layer for reward computation. Defaults to 10.
        """
        self._image = image
        self.center = center
        self.width = width
        self.height = height
        self.angle = angle
        self.edge_width = edge_width
        self.layer_height = layer_height
        self.random_seed = random_seed

    @property
    def image(self) -> np.ndarray:
        """Returns the image associated with the sampling box."""
        return self._image
    
    def clone(self) -> 'SamplingBoxAgent':
        """Creates a clone of the current sampling box agent.

        Returns:
            SamplingBoxAgent: A new instance with the same parameters.
        """
        return SamplingBoxAgent(
            image=self._image,
            center=self.center,
            width=self.width,
            height=self.height,
            angle=self.angle,
            layer_height=self.layer_height,
            edge_width = self.edge_width
        )

    def _is_out_of_bounds(self) -> bool:

        box = self.generate_box()
        box[0, ] = (box[0,] + box[1,]) / 2
        box[3, ] = (box[2,] + box[3,]) / 2

        img_h, img_w = self.image.shape[:2]

        return not np.all(
            (0 <= box[:, 0]) & (box[:, 0] < img_h) & 
            (0 <= box[:, 1]) & (box[:, 1] < img_w)
        )


    def get_extension_length(self) -> float:
        """Calculates the extension length based on the sampling box geometry.

        Returns:
            float: The extension length.
        """
        box = self.generate_box()
        # Compute Euclidean distance between center and midpoint of two adjacent box edges.
        distance = np.linalg.norm(self.center - ((box[1] + box[2]) // 2), ord=2)
        return distance
    
    def get_state_vector(self) -> Dict[str, Any]:
        """Retrieves the state vector representing the agent's current parameters.

        Returns:
            Dict[str, Any]: A dictionary with keys 'center', 'width', 'height', 'angle'.
        """
        return {
            'center': self.center[0],
            'width': self.width,
            'height': self.height,
            'angle': self.angle
        }

    def set_state(self, key: str, value: Any) -> None:
        """Sets a state variable of the agent.

        Args:
            key (str): The name of the state variable.
            value (Any): The new value for the state variable.

        Raises:
            KeyError: If the specified state variable does not exist.
        """
        if key in self.get_state_vector():
            if key == 'center':
                self.__dict__[key] = (value, value)
            else:
                self.__dict__[key] = value
        else:
            raise KeyError(f"State variable {key} does not exist.")

    def generate_box(self, center: Tuple[float, float] = None, width: int = None,
                     height: int = None, angle: int = None) -> np.ndarray:
        """Generates the coordinates of the rotated sampling box.

        Args:
            center (Tuple[float, float], optional): Center of the box. Defaults to the agent's center.
            width (int, optional): Width of the box. Defaults to the agent's width.
            height (int, optional): Height of the box. Defaults to the agent's height.
            angle (int, optional): Rotation angle in degrees. Defaults to the agent's angle.

        Returns:
            np.ndarray: An array of shape (4, 2) representing the box corners.
        """
        import cv2
        center = center or self.center
        width = width or self.width
        height = height or self.height
        angle = angle or self.angle

        rect = ((float(center[0]), float(center[1])), (float(width), float(height)), float(angle))
        box = cv2.boxPoints(rect).astype(np.intp)
        if not np.all(np.isfinite(box)):
            raise ValueError("Invalid box coordinates: contains NaN or non-finite values.")
        return box

    def apply_transform(self, transform_type: str, n_expand: int = 0,
                        n_extension: int = 0, theta: float = 0,
                        n_translocation: int = 0) -> None:
        """Applies a transformation to the sampling box.

        Args:
            transform_type (str): Type of transformation to apply.
            n_expand (int, optional): Expansion value for width. Defaults to 0.
            n_extension (int, optional): Extension value for height. Defaults to 0.
            theta (float, optional): Angle change in degrees. Defaults to 0.
            n_translocation (int, optional): Translation value for the center. Defaults to 0.

        Raises:
            ValueError: If the transformation results in invalid dimensions.
        """
        if transform_type == 'comb':
            self.angle += theta
            self.center = (self.center[0] + n_translocation, self.center[1] + n_translocation)
            self.width += n_expand
            self.height += n_extension

        if not (np.isfinite(self.height) and np.isfinite(self.width)):
            raise ValueError("Invalid transformation: height or width is not finite.")

    def compute_mean_intensity(self, box: np.ndarray) -> float:
        """Computes the mean pixel intensity inside a given rotated box.

        Args:
            box (np.ndarray): The coordinates of the rotated box.

        Returns:
            float: The mean intensity within the box.
        """
        return compute_mean_intensity(self.image, box)

    def cal_rewards(self, eta: float = 0.5) -> float:
        """Calculates the reward based on signal intensities and gradients with decay weighting.

        This function splits the sampling box into multiple layers and calculates:
        - The signal strength in each layer.
        - The signal intensities in the upstream and downstream edge and inner regions.
        - The gradient computed as the minimum of (upstream_inner - upstream_edge) and 
            (downstream_inner - downstream_edge) for each layer.
        Then, it computes the differences between adjacent layers (signal differences) and applies
        an exponential decay weighting based on the layer index so that gradients at lower layers
        (closer to the root) are given higher importance. The final reward is the sum of the weighted 
        gradients minus the weighted signal differences.

        Args:
            eta (float, optional): Scaling factor for the penalty based on signal differences.
                Defaults to 0.5.
            lambda_decay (float, optional): Decay coefficient for gradient weighting. A higher value 
                implies a faster decay with increasing layer index. Defaults to 0.1.

        Returns:
            float: The computed reward value. If all signal strengths are NaN, returns np.nan.
        """
        # Split the sampling box into several regions for each layer.
        boxes, boxes_upstream_edge, boxes_downstream_edge, \
        boxes_upstream_inner, boxes_downstream_inner = self.split_box()

        # Compute the mean signal intensity in each region.
        signal_strengths = np.array(
            [self.compute_mean_intensity(box) for box in boxes]
        )
        signal_upstream_edge = np.array(
            [self.compute_mean_intensity(box) for box in boxes_upstream_edge]
        )
        signal_downstream_edge = np.array(
            [self.compute_mean_intensity(box) for box in boxes_downstream_edge]
        )

        signal_upstream_inner = np.array(
            [self.compute_mean_intensity(box) for box in boxes_upstream_inner]
        )

        signal_downstream_inner = np.array(
            [self.compute_mean_intensity(box) for box in boxes_downstream_inner]
        )

        # Calculate gradients: the difference between inner and edge signals.
        upstream_gradient = signal_upstream_inner - signal_upstream_edge
        downstream_gradient = signal_downstream_inner - signal_downstream_edge
        gradients = np.nanmean([upstream_gradient, downstream_gradient], axis=0)

        # Compute differences between adjacent layers in the main sampling box.
        signal_diffs = abs(np.diff(signal_strengths))

        # Sum the weighted gradients minus the weighted differences.
        delta = np.nansum(gradients[:-1] - eta * signal_diffs)

        if np.isnan(signal_strengths).all():
            logging.warning("All signal strengths are NaN; reward cannot be computed.")
            return np.nan
        
        # logging.warning(f"delta is {delta}, signal strength is {np.nansum(signal_strengths)}")
        return delta
  
    def split_box(self) -> Tuple[List[np.ndarray], List[np.ndarray],
                                                       List[np.ndarray], List[np.ndarray],
                                                       List[np.ndarray]]:
        """Generates multiple region boxes for reward computation based on the current sampling box.
        
        The following boxes are generated:
          - sampling_boxes: The main sampling boxes for each layer.
          - upstream_edge_boxes: The edge boxes on the upstream side.
          - downstream_edge_boxes: The edge boxes on the downstream side.
          - upstream_inner_boxes: The inner boxes near the upstream edge.
          - downstream_inner_boxes: The inner boxes near the downstream edge.
        
        Note:
          This method uses the object's properties (e.g., center, width, height, angle, layer_height).
        
        Returns:
            Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
                A tuple containing five lists of numpy arrays representing the coordinates of:
                (sampling_boxes, upstream_edge_boxes, downstream_edge_boxes,
                 upstream_inner_boxes, downstream_inner_boxes).
        
        Raises:
            ValueError: If height or layer_height is not finite or layer_height <= 0.
        """
        if not np.isfinite(self.height) or not np.isfinite(self.layer_height) or self.layer_height <= 0:
            raise ValueError("Invalid agent state: height or layer_height must be finite and > 0.")

        rad_angle = np.deg2rad(self.angle)  # Convert angle to radians
        num_layers = int(self.height // 2 // self.layer_height)  # Calculate number of layers

        sampling_boxes: List[np.ndarray] = []         # Main sampling boxes
        upstream_edge_boxes: List[np.ndarray] = []    # Upstream edge boxes
        downstream_edge_boxes: List[np.ndarray] = []  # Downstream edge boxes
        upstream_inner_boxes: List[np.ndarray] = []   # Upstream inner boxes
        downstream_inner_boxes: List[np.ndarray] = [] # Downstream inner boxes

        half_width = self.width / 2.0  # Half of the box width
        current_edge_width = min(self.edge_width, half_width)  # Adjust edge width dynamically


        half_edge_width = current_edge_width / 2.0  # Half of edge width

        for i in range(num_layers):
            offset = i * self.layer_height
            # Calculate the center of the current layer
            layer_center = (
                self.center[0] + offset * np.sin(rad_angle),
                self.center[1] - offset * np.cos(rad_angle)
            )

            # Calculate upstream and downstream centers by shifting by half the width
            upstream_center = rotate_point(
                layer_center[0] - half_width, layer_center[1],
                layer_center[0], layer_center[1], rad_angle
            )
            downstream_center = rotate_point(
                layer_center[0] + half_width, layer_center[1],
                layer_center[0], layer_center[1], rad_angle
            )

            # Calculate edge box centers by further shifting by half_edge_width
            upstream_edge_center = rotate_point(
                upstream_center[0] - half_edge_width, upstream_center[1],
                upstream_center[0], upstream_center[1], rad_angle
            )
            downstream_edge_center = rotate_point(
                downstream_center[0] + half_edge_width, downstream_center[1],
                downstream_center[0], downstream_center[1], rad_angle
            )
            # Calculate inner box centers by shifting inward
            upstream_inner_center = rotate_point(
                upstream_center[0] + half_edge_width, upstream_center[1],
                upstream_center[0], upstream_center[1], rad_angle
            )
            downstream_inner_center = rotate_point(
                downstream_center[0] - half_edge_width, downstream_center[1],
                downstream_center[0], downstream_center[1], rad_angle
            )

            # Generate boxes using the generate_box method (assumed defined elsewhere)
            sampling_boxes.append(self.generate_box(layer_center, self.width, self.layer_height))
            upstream_edge_boxes.append(self.generate_box(upstream_edge_center, current_edge_width, self.layer_height))
            downstream_edge_boxes.append(self.generate_box(downstream_edge_center, current_edge_width, self.layer_height))
            upstream_inner_boxes.append(self.generate_box(upstream_inner_center, current_edge_width, self.layer_height))
            downstream_inner_boxes.append(self.generate_box(downstream_inner_center, current_edge_width, self.layer_height))

        return sampling_boxes, upstream_edge_boxes, downstream_edge_boxes, upstream_inner_boxes, downstream_inner_boxes

    def _generate_padding_box(self):
        center = self.center
        width = self.width
        height = self.height
        angle = self.angle
        rad_angle = np.deg2rad(angle)

        upstream_center = rotate_point(
            center[0] - width, center[1],
            center[0], center[1], rad_angle
        )
        downstream_center = rotate_point(
            center[0] + width, center[1],
            center[0], center[1], rad_angle
        )
        
        upstream_backgroundbox = self.generate_box(
            center=upstream_center, 
            width=width,
            height=height, 
            angle=angle
        )
        
        downstream_backgroundbox = self.generate_box(
            center=downstream_center, 
            width=width,
            height=height, 
            angle=angle
        )
        
        return upstream_backgroundbox, downstream_backgroundbox

    def _rescale_array(self, arr, target_len=30):
        arr = np.asarray(arr)
        n = len(arr)
        if n == 0:

            return np.zeros(target_len)
        
        old_x = np.arange(n)
        new_x = np.linspace(0, n - 1, target_len)
        
        new_arr = np.interp(new_x, old_x, arr)
        return new_arr

    def _fit_zero_inflated_gamma(self, data):
        data = np.array(data)
        if len(data) == 0:
            return np.nan, np.nan, np.nan

        def negloglik(params):
            p, alpha, scale = params
            if not (0 <= p <= 1) or alpha <= 0 or scale <= 0:
                return np.inf

            zeros = data == 0
            nonzeros = ~zeros
            n0 = np.sum(zeros)
            n1 = np.sum(nonzeros)

            ll_zero = n0 * np.log(p + 1e-10) 
            ll_nonzero = np.sum(np.log((1 - p + 1e-10) * gamma.pdf(data[nonzeros], a=alpha, loc=0, scale=scale) + 1e-10))

            return -(ll_zero + ll_nonzero)

        zero_percent = np.mean(data == 0)
        bounds = [(0.001, 0.999), (1e-3, 1e3), (1e-3, 1e3)] 
        init_params = [zero_percent, 2.0, 1.0]

        result = minimize(negloglik, init_params, bounds=bounds, method='L-BFGS-B', options={'maxiter': 200})

        if result.success:
            return result.x
        else:
            return np.nan, np.nan, np.nan

    def _calculate_intensity_within_box(self, box):
        pixels = [] 
        #only half of a box is ok
        box[0,] = (box[0,] + box[1,]) / 2
        box[3,] = (box[2,] + box[3,]) / 2
        x_min, x_max = max(0, min(box[:, 0])), min(max(box[:, 0]), self.image.shape[0])
        y_min, y_max = max(0, min(box[:, 1])), min(max(box[:, 1]), self.image.shape[0])

        for x in range(x_min, x_max):
            for y in range(y_min, y_max):
                if inside_rotated_box(x, y, box):
                    pixels.append(self.image[y, x])
        try:
            p, alpha, scale = self._fit_zero_inflated_gamma(pixels)
        except:
            p, alpha, scale = np.nan, np.nan, np.nan
        intensity = (1 - p) * alpha * scale
        
        return intensity, p, alpha, scale, pixels

    def calculate_intensity(self) -> Tuple[float, float, float]:
        box = self.generate_box()
        sampling, _, _, _, _ = self._calculate_intensity_within_box(box)
        
        return sampling
    
    def compute_sum_interactions(self) -> float:
        """Computes the sum of pixel intensities inside a given rotated box.

        Args:
            box (np.ndarray): The coordinates of the rotated box.

        Returns:
            float: The sum of intensities within the box.
        """
        box = self.generate_box()

        return compute_sum_intensity(self.image, box)

    def calculate_quality(self) -> float:
        """Calculates the fountainess metric and total intensity from the sampling box.
        
        This function computes signal intensities for multiple regions (full box, edge and inner regions),
        calculates gradients, and then derives the fountainess value based on the median signal strength,
        gradient delta, and extension length.
        
        Returns:
            Tuple[float, float]: A tuple containing the fountainess value and the total intensity.
                Returns (np.nan, np.nan) if no valid boxes are available.
        """
        boxes, boxes_upstream_edge, boxes_downstream_edge, boxes_upstream_inner, boxes_downstream_inner = self.split_box()
        if not boxes:
            logging.warning("No boxes generated in split_box; returning NaN for fountainess and intensity.")
            return np.nan, np.nan
        
        signal_strengths = np.array([self.compute_mean_intensity(box) for box in boxes])
        signal_upstream_edge = np.array([self.compute_mean_intensity(box) for box in boxes_upstream_edge])
        signal_downstream_edge = np.array([self.compute_mean_intensity(box) for box in boxes_downstream_edge])
        signal_upstream_inner = np.array([self.compute_mean_intensity(box) for box in boxes_upstream_inner])
        signal_downstream_inner = np.array([self.compute_mean_intensity(box) for box in boxes_downstream_inner])
        

        # Calculate gradients for upstream and downstream regions
        upstream_gradient = signal_upstream_inner - signal_upstream_edge
        downstream_gradient = signal_downstream_inner - signal_downstream_edge
        
        # Use the minimum gradient as penalty
        min_gradient = np.nanmean([upstream_gradient, downstream_gradient], axis=0)
        signal_diffs = np.abs(np.diff(signal_strengths))  # differences between adjacent layers
        delta = np.nansum(min_gradient[:-1] - signal_diffs)
        
        if np.isnan(signal_strengths).all():
            logging.warning("All signal strengths are NaN; quality cannot be computed.")
            return np.nan, np.nan
        
        ext_length = self.get_extension_length()
        quality = delta / ext_length

        return quality


    def split_box_fdr(self):
        rad_angle = np.deg2rad(self.angle)
        num_layers = int(self.height // 2 // self.layer_height)
        layer_boxes, upstream_boxes, downstream_boxes = [], [], []

        for i in range(num_layers):

            offset = i * self.layer_height
            layer_center = (
                self.center[0] + offset * np.sin(rad_angle),
                self.center[1] - offset * np.cos(rad_angle)
            )

            upstream_center = rotate_point(
                layer_center[0] - self.width, layer_center[1],
                layer_center[0], layer_center[1], rad_angle
            )

            downstream_center = rotate_point(
                layer_center[0] + self.width, layer_center[1],
                layer_center[0], layer_center[1], rad_angle
            )
            if any(i < 0 for i in layer_center + upstream_center + downstream_center):
                break

            layer_boxes.append(self.generate_box(layer_center, self.width, self.layer_height))
            upstream_boxes.append(self.generate_box(upstream_center, self.width, self.layer_height))
            downstream_boxes.append(self.generate_box(downstream_center, self.width, self.layer_height))

        return layer_boxes, upstream_boxes, downstream_boxes

    def cal_foldchange(self) -> Tuple[float, float]:
        """Calculates the fold change of signal strengths in the upstream and downstream edge regions.
        
        Returns:
            Tuple[float, float]: The upstream and downstream fold change values.
        """
        boxes, boxes_upstream_edge, boxes_downstream_edge, boxes_upstream_inner, boxes_downstream_inner = self.split_box()
        signal_strengths = np.array([self.compute_mean_intensity(box) for box in boxes])
        signal_upstream_edge = np.array([self.compute_mean_intensity(box) for box in boxes_upstream_edge])
        signal_downstream_edge = np.array([self.compute_mean_intensity(box) for box in boxes_downstream_edge])
        
        # Compute fold change as the ratio between main sampling box and edge intensities.
        upstream_fold_change = np.nanmean(signal_strengths) / (np.nanmean(signal_upstream_edge) + 1e-10)
        downstream_fold_change = np.nanmean(signal_strengths) / (np.nanmean(signal_downstream_edge) + 1e-10)
        
        # logging.info("Calculated fold changes: upstream=%.4f, downstream=%.4f", upstream_fold_change, downstream_fold_change)
        return upstream_fold_change, downstream_fold_change

    def digest_extension(self) -> np.ndarray:
        """
        Retrieves the signal strengths for each layer within the sampling box.

        Returns:
            np.ndarray: An array containing the signal strengths for each layer.
        """
        boxes, _, _, _, _ = self.split_box()
        signal_strengths = np.array(
            [self.compute_mean_intensity(box) for box in boxes]
        )

        return signal_strengths

    def plot_box(self):
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list(
            'interaction', ['#FFFFFF', '#FFDFDF', '#FF7575', '#FF2626', "#F70000"]
        )

        box = self.generate_box()
        image_copy = np.copy(self.image)
        
        # adjust box
        box[0] = (box[0] + box[1])//2
        box[3] = (box[2] + box[3])//2

        patch = patches.Polygon(box, closed=True, edgecolor='blue', facecolor='none', linewidth=2)
        current_reward = self.cal_rewards()

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.add_patch(patch)
        ax.imshow(image_copy, cmap=cmap, vmax=8, vmin=0)

        plt.title(f"Rotated Sampling Box: center={self.center}, width={self.width}, height={self.height}, angle={self.angle}")
        ax.text(10, 10, f'SoN: {current_reward:.3f}', color='white', fontsize=15, bbox=dict(facecolor='black', alpha=0.5))

        plt.show()

    def perform_test(self):
        _, boxes_up_edge, boxes_dn_edge, boxes_up_inner, boxes_dn_inner = self.split_box()
  
        sig_up_edge = np.array([compute_mean_intensity(self.image, box) for box in boxes_up_edge])
        sig_dn_edge = np.array([compute_mean_intensity(self.image, box) for box in boxes_dn_edge])

        sig_up_inner = np.array([compute_mean_intensity(self.image, box) for box in boxes_up_inner])
        sig_dn_inner = np.array([compute_mean_intensity(self.image, box) for box in boxes_dn_inner])

        W_stat1, p_value1 = wilcoxon(sig_up_inner, sig_up_edge, alternative='greater')
        W_stat2, p_value2 = wilcoxon(sig_dn_inner, sig_dn_edge, alternative='greater')

        # d_up = paired_cohens_d(sig_up_inner, sig_up_edge)
        # d_dn = paired_cohens_d(sig_dn_inner, sig_dn_edge)

        n = len(sig_up_edge)
        r_rb1 = (2 * W_stat1) / (n * (n+1) / 2) - 1
        r_rb2 = (2 * W_stat2) / (n * (n+1) / 2) - 1

        return p_value1, p_value2, r_rb1, r_rb2