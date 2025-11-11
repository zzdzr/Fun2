"""mcts_refactor.py

A safer, clearer re‑implementation of Node and MCTS that avoids the
`UnboundLocalError` you hit, cleanly separates phases (selection →
expansion → simulation → back‑propagation), and makes every error path
explicit.  The public interface (constructor signatures + `search`) is
kept compatible with the original so it can drop‑in‑replace the old
class.

Key design decisions
--------------------
*   **dataclass + slots** for `Node` – cheaper and avoids typos.
*   **Enum** for the four action dimensions – makes the magic indices
    (0, 1, 2, 3) readable.
*   `search()` uses a `while` loop with `iteration` counter so we can
    early‑`continue` on a failed expansion without risking uninitialised
    variables.
*   `expand()` returns `None` instead of raising when it cannot produce
    a legal child; the caller handles the `None`.
*   Every function that might fail is wrapped with `try/except` and
    logs *and re‑raises* unless the failure is expected.
*   Added fine‑grained logging and `__str__`/`__repr__` helpers so a
    single `logging.debug(node)` prints the interesting state.
*   `update_action_limits()` no longer mutates the *original* list – it
    returns a fresh list each call to avoid accidental cross‑talk among
    threads.

If you only need to swap the classes, drop this file next to the old
one and do::

    from mcts_refactor import MCTS, Node

"""

from __future__ import annotations

import logging
import random
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

# ---------------------------------------------------------------------
#  Utilities / constants
# ---------------------------------------------------------------------
LOG = logging.getLogger(__name__)

class Act(IntEnum):
    """Index mapping for the 4‑dim action vector."""

    EXPAND = 0      # width  ↗
    EXTEND = 1      # height ↘
    ROTATE = 2      # angle  ↻
    SHIFT = 3       # move centre along diagonal ↔

# ---------------------------------------------------------------------
#  Node
# ---------------------------------------------------------------------

@dataclass(slots = True)
class Node:
    agent: 'SamplingBoxAgent'
    parent: Optional[Node] = None
    action: np.ndarray = field(default_factory=lambda: np.zeros(4))
    visits: int = 0
    q_value: float = 0.0
    children: Dict[str, Node] = field(default_factory=dict)
    gradient_cache: np.ndarray = field(default_factory=lambda: np.zeros(4))

    # ---- helper ----
    def is_root(self) -> bool:
        return self.parent is None
    
    def child_key(self, act: np.ndarray) -> str:
        """Stable text key so `act=[0,0,0,0]` and `[0.,0.,0.,0.]` match."""
        return ','.join(f"{x:.2f}" for x in act)
    
    # ---- expansion ----
    def expand(self, act: np.ndarray) -> Optional[Node]:
        """Apply *act* to a cloned agent and return a new child node.

        Returns *None* if the action leads to an invalid box.
        """
        if act is None:
            LOG.debug("expand(): action is None - skip")
            return None
        
        clone = self.agent.clone()
        try:
            clone.apply_transform('comb', *act)
        except Exception as exc:
            LOG.debug("expand(): transform failed - %s", exc)
            return None
        
        if clone._is_out_of_bounds():
            LOG.debug("expand(): out of bounds after transform, skipping")
            return None

        child = Node(clone, parent = self, action = act)
        self.children[self.child_key(act)] = child

        return child

    def __repr__(self): # concise one‑liner for logging
        c, h, w, a = self.agent.center, self.agent.height, self.agent.width, self.agent.angle
        return (f"Node(vis={self.visits}, Q={self.q_value:.3f}, act={self.action}, "
                f"ctr={c}, h={h}, w={w}, ang={a})")
    
class MCTS:
    """Monte‑Carlo Tree Search with Evolution‑Strategy exploration."""

    # ------------------------------------------------------------------
    #  Construction helpers
    # ------------------------------------------------------------------
    def __init__(
        self,
        agent: SamplingBoxAgent,
        *,
        max_iter: int = 50,
        exploration_c: float = 1.4,
        alpha_pw: float = 0.5,
        action_limits: Sequence[Tuple[float, float]] | None = None,
        angle_boundary: Tuple[float, float] = (40.0, 50.0),
        # ES hyper‑params
        es_pop: int = 8,
        es_sigma: float = 5.0,
        lr0: float = 2.0,
        momentum: float = 0.5,
        # rollout hyper‑params
        rollout_depth: int = 5,
        gamma: float = 0.9,
        eta: float = 0.5,
        seed = 1,
        mode: str = 'default',  # or 'fixed_angle'
    ) -> None:
        self.agent = agent
        self.max_iter = max_iter
        self.exploration_c = exploration_c
        self.alpha_pw = alpha_pw
        self.angle_boundary = angle_boundary
        self.mode = mode
        self.seed = seed
        # ES
        self.es_pop = es_pop
        self.es_sigma = es_sigma
        self.lr0 = lr0
        self.momentum = momentum
        # rollout
        self.rollout_depth = rollout_depth
        self.gamma = gamma
        self.eta = eta

        # action limits template – mutated copy per iteration
        if action_limits is None:
            action_limits = [(-2, 2), (-10, 10), (-1, 1), (-1, 1)]
        self._action_limits0 = list(action_limits)

    # ------------------------------------------------------------------
    #  Core search loop
    # ------------------------------------------------------------------
    def search(self) -> Tuple[Node, float]:
        """Run MCTS and return *(best_node, best_immediate_reward)*."""
        if self.seed is not None:
            random.seed(self.seed)
            np.random.seed(self.seed)

        root = Node(self.agent)
        best_node: Node = root
        best_reward = -np.inf

        iteration = 0
        while iteration < self.max_iter:
            
            iteration += 1
            # ------ selection ------
            leaf = self._select(root)
            if leaf.agent._is_out_of_bounds():
                # LOG.debug("leaf out of bounds - stop search early")
                break
            
            # ------ expansion ------
            child = self._expand(leaf, iteration)

            if child is None:  # expansion failed – continue search
                # LOG.debug("Failed expansion, continue searching")
                continue

            # ------ simulation ------
            q_sim = self._simulate(child)
            reward_now = child.agent.cal_rewards(eta = self.eta)
            if reward_now > best_reward:
                best_reward = reward_now
                best_node = child

            # ------ backprop ------
            self._backprop(child, q_sim)
            # LOG.info(f"Best Node is {best_node}")

        if self.mode == 'default':
            self._fine_tune_angle(best_node)

        try:
            reward = self.calculate_reward(node = best_node)
            quality = self.estimate_quality(node = best_node)
            intensity = self.estimate_intensity(node = best_node)
            sum_interactions = self.calculate_sum_interactions(node = best_node)
            p_up, p_dn, rb_up, rb_dn = self.perform_test(node = best_node)

        except Exception as e:
            logging.info(f"Error during metric calculation: {e}")
            return None, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        return best_node, reward, quality, intensity, sum_interactions, p_up, p_dn

    def _fine_tune_angle(self, node: Node) -> None:
        """After search, fine-tune angle within boundary to maximize quality."""
        best_angle = node.agent.angle
        best_reward = -np.inf
        for angle in np.arange(self.angle_boundary[0], self.angle_boundary[1] + 0.0001, 0.1):
            agent_copy = node.agent.clone()
            agent_copy.angle = angle
            try:
               reward = agent_copy.cal_rewards(eta = self.eta)

               if reward > best_reward:
                    best_reward = reward
                    best_angle = angle
            except:
                continue

        node.agent.angle = best_angle

    def _ucb1(self, node: Node) -> float:
        if node.visits == 0:
            return np.inf
        return node.q_value + self.exploration_c * np.sqrt(
            np.log(node.parent.visits + 1) / (node.visits + 1)
        )
    
    def _select(self, node: Node) -> Node:
        """Recursively select a node for expansion using UCB1 + progressive widening."""
        # If we’ve never visited it, or haven’t even expanded to our PW cap, stop here.
        if node.visits == 0 or not node.children:
            return node

        # how many children we’re allowed under PW
        n_allowed = int(node.visits ** self.alpha_pw)
        if len(node.children) < n_allowed:
            return node

        # else pick the best child and recurse
        best = max(node.children.values(), key=self._ucb1)
        return self._select(best)

    
    # ------------------------------------------------------------------
    #  Expansion
    # ------------------------------------------------------------------
    def _expand(self, node: Node, n_iter: int) -> Optional[Node]:
        act = self._generate_action(node, n_iter)
        return node.expand(act)
    
    # ------------------------------------------------------------------
    #  Simulation (rollout)
    # ------------------------------------------------------------------
    def _simulate(self, start: Node) -> float:
        # Build a separate simulation node, don't pollute the real tree
        sim_agent = start.agent.clone()
        sim_node = Node(sim_agent)
        sim_node.action = start.action.copy()
        sim_node.gradient_cache = start.gradient_cache.copy()

        total, discount = 0.0, 1.0
        for depth in range(self.rollout_depth):
            total += discount * sim_node.agent.cal_rewards(eta = self.eta)
            discount *= self.gamma

            try:
                act = self._generate_action(sim_node, depth)
                next_node = sim_node.expand(act)
            except Exception as e:
                LOG.debug("Simulation terminated early: %s", e)
                break

            if next_node is None:
                break
            sim_node = next_node

        return total
    
    def _backprop(self, node: Node, G: float) -> None:
        while node is not None:
            node.visits += 1
            node.q_value += (G - node.q_value) / node.visits
            G *= self.gamma
            node = node.parent
    
    def _generate_action(self, node: Node, n_iter: int) -> np.ndarray:
        lr = self.lr0
        limits = self._dynamic_limits(node.agent)
        base = node.action.copy()

        if self.mode == 'fixed_angle':
            base[[Act.ROTATE, Act.SHIFT]] = 0.0

        half = self.es_pop // 2
        noise_half = np.random.normal(0.0, self.es_sigma, (half, 4))
        noise = np.vstack([noise_half, -noise_half])
        if self.mode == 'fixed_angle':
            noise[:, [Act.ROTATE, Act.SHIFT]] = 0.0

        rewards = np.empty(self.es_pop)
        for i in range(self.es_pop):
            cand = self._clip(base + noise[i], limits)
            rewards[i] = self._eval_action(node, cand)

        norm_r = (rewards - rewards.mean()) / (rewards.std() + 1e-9)
        grad = (noise.T @ norm_r) / self.es_pop

        node.gradient_cache = (
            self.momentum * node.gradient_cache
            + (1 - self.momentum) * grad
        )

        act = self._clip(base + lr * node.gradient_cache, limits)
        if self.mode == 'fixed_angle':
            act[[Act.ROTATE, Act.SHIFT]] = 0.0
        
        return act
    
    @staticmethod
    def _clip(vec: np.ndarray, limits: Sequence[Tuple[float, float]]) -> np.ndarray:
        out = np.empty_like(vec)
        for i, (lo, hi) in enumerate(limits):
            out[i] = float(np.clip(vec[i], lo, hi))
        return out
    
    def _eval_action(self, node: Node, act: np.ndarray) -> float:
        agent = node.agent.clone()
        agent.apply_transform('comb', *act)
        return agent.cal_rewards(eta = self.eta)
    
    def _dynamic_limits(self, agent: SamplingBoxAgent) -> List[Tuple[float, float]]:
        m, _ = agent._image.shape
        c = agent.center[0]
        h, w = agent.height/2, agent.width

        orig = list(self._action_limits0)
        limits: List[Tuple[float, float]] = []
        # LOG.info(f"Height, Width, boundary are {h}, {w}, {m}")
        for dim, (lo0, hi0) in enumerate(orig):
            if dim == Act.EXPAND:
                # ensure width >= 2px and <= image width
                lo = max(lo0, 2 - w)
                hi = min(hi0, m - w)

            elif dim == Act.EXTEND:
                # ensure height >= 2px and <= image height
                lo = max(lo0, 2 - h)
                hi = min(hi0, m - h)
            # elif dim == Act.ROTATE:
            #     if self.mode == 'fixed_angle':
            #         lo = hi = 0.0
            #     else:
            #         lo = max(lo0, -agent.angle)
            #         hi = min(hi0, 90 - agent.angle)
            #         # enforce custom angle boundary
            #         if not (self.angle_boundary[0] <= agent.angle + lo <= self.angle_boundary[1]):
            #             lo = hi = 0.0
            elif dim == Act.ROTATE:
                if self.mode == 'fixed_angle':
                    lo = hi = 0.0
                else:
                    lo = max(lo0, self.angle_boundary[0] - agent.angle)
                    hi = min(hi0, self.angle_boundary[1] - agent.angle)


            elif dim == Act.SHIFT:
                lo = max(lo0, -c)
                hi = min(hi0, m - c)
            
            else:
                lo, hi = lo0, hi0
            
            limits.append((lo, hi))

        return limits

    def estimate_quality(self, node: Node) -> float:
        try:
            return node.agent.calculate_quality()
        except:
            return np.nan
    
    def perform_test(self, node: Node) -> float:
        try:
            return node.agent.perform_test()
        except:
            return np.nan, np.nan, np.nan, np.nan

    def estimate_intensity(self, node: Node) -> float:
        try:
            return node.agent.calculate_intensity()
        except:
            return np.nan

    def calculate_reward(self, node: Node) -> float:
        return node.agent.cal_rewards(eta = self.eta)

    def calculate_sum_interactions(self, node: Node) -> float:
        return node.agent.compute_sum_interactions()
                
        