"""
Position-specific quality score model for NEAT.

The model consists of an initial distribution giving the probability of observing each
quality score at the first position of a read, per-position marginals, and per-position
transition matrices.

If transition_distributions is None, the model independently samples from
per-position marginals.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

__all__ = ["MarkovQualityModel"]

@dataclass
class MarkovQualityModel:
    """
    Parameters
    ----------
    initial_distribution:
        A mapping from integer quality scores to their observed
        probabilities at position 0.
    position_distributions:
        A list of mappings from integer quality scores to their observed
        probabilities at each position. Element ``i`` in the list
        corresponds to position ``i`` in the read.
    max_quality:
        The maximum quality observed in the training data. Generated
        qualities will be clipped to not exceed the maximum quality.
    read_length:
        The length of reads used during training determines how many
        position-specific distributions exist.
    transition_distributions:
        Optional list trans[i][q_prev][q_next] = count/prob.
    """

    initial_distribution: Dict[int, float]
    position_distributions: List[Dict[int, float]]
    max_quality: int
    read_length: int
    transition_distributions: Optional[List[Dict[int, Dict[int, float]]]] = None

    def __post_init__(self) -> None:

        # Normalize the initial distribution
        total_init = float(sum(self.initial_distribution.values()))

        if total_init <= 0:
            raise ValueError(
                "Initial distribution must have positive probability mass."
            )

        self.initial_distribution = {
            int(k): float(v) / total_init for k, v in self.initial_distribution.items()
        }

        # Normalize per-position marginals
        norm_positions: List[Dict[int, float]] = []

        for dist in self.position_distributions:
            total = float(sum(dist.values()))

            if total <= 0:
                norm_positions.append({0: 1.0})
                continue

            norm_positions.append({int(k): float(v) / total for k, v in dist.items()})

        self.position_distributions = norm_positions

        # Ensure max_quality and read_length are integers
        self.max_quality = int(self.max_quality)
        self.read_length = int(self.read_length)

        # Precompute numpy arrays for fast sampling
        init_keys = list(self.initial_distribution.keys())
        init_vals = [self.initial_distribution[k] for k in init_keys]
        self._init_scores = np.asarray(init_keys, dtype=int)
        self._init_probs = np.asarray(init_vals, dtype=float)

        self._values_by_pos: List[np.ndarray] = []
        self._probs_by_pos: List[np.ndarray] = []

        for dist in self.position_distributions:
            keys = list(dist.keys())
            vals = [dist[k] for k in keys]
            self._values_by_pos.append(np.asarray(keys, dtype=int))
            self._probs_by_pos.append(np.asarray(vals, dtype=float))

        # Normalize and cache transitions
        self._has_transitions = False
        self._trans_values_by_pos: List[Dict[int, np.ndarray]] = []
        self._trans_probs_by_pos: List[Dict[int, np.ndarray]] = []

        if self.transition_distributions is not None:
            if len(self.transition_distributions) not in (0, max(0, self.read_length - 1)):
                raise ValueError(
                    "transition_distributions must have length read_length-1 "
                    f"(expected {max(0, self.read_length - 1)}, got {len(self.transition_distributions)})."
                )

            # Normalize each row q_prev -> distribution over q_next
            self._has_transitions = len(self.transition_distributions) > 0

            for pos_trans in self.transition_distributions:
                cache_vals: Dict[int, np.ndarray] = {}
                cache_probs: Dict[int, np.ndarray] = {}

                for q_prev, next_map in pos_trans.items():
                    total = float(sum(next_map.values()))
                    if total <= 0:
                        continue

                    keys = [int(k) for k in next_map.keys()]

                    keys = []
                    probs = []

                    for qn, c in next_map.items():
                        keys.append(int(qn))
                        probs.append(float(c) / total)

                    cache_vals[int(q_prev)] = np.asarray(keys, dtype=int)
                    cache_probs[int(q_prev)] = np.asarray(probs, dtype=float)

                self._trans_values_by_pos.append(cache_vals)
                self._trans_probs_by_pos.append(cache_probs)

    @property
    def quality_scores(self) -> List[int]:
        """List of supported quality scores."""

        return list(range(0, self.max_quality + 1))

    def _position_index_for_length(self, pos: int, length: int) -> int:
        """Map a position in a generated read onto a training position."""

        if length <= 1 or self.read_length <= 1:
            return 0

        idx = int(round((pos / max(1, length - 1)) * (self.read_length - 1)))

        if idx < 0:
            idx = 0

        elif idx > self.read_length - 1:
            idx = self.read_length - 1

        return idx

    def get_quality_scores(
        self,
        model_read_length: int,
        length: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Generate a synthetic quality score array.

        If transitions are provided, use the Markov-based model.

        If transitions are not provided, use an empirical version of the model.
        """

        _ = model_read_length

        if length <= 0:
            return np.zeros(0, dtype=int)

        qualities = np.zeros(length, dtype=int)

        # Sample initial quality from the starting distribution
        qualities[0] = int(rng.choice(self._init_scores, p=self._init_probs))

        for i in range(1, length):
            p_cur = self._position_index_for_length(i, length)

            q = None

            if self._has_transitions and self._trans_values_by_pos:

                p_prev = self._position_index_for_length(i - 1, length)

                # Transitions are defined for positions
                p_prev = min(p_prev, len(self._trans_values_by_pos) - 1)

                q_prev = int(qualities[i - 1])
                vals_map = self._trans_values_by_pos[p_prev]
                probs_map = self._trans_probs_by_pos[p_prev]

                if q_prev in vals_map:
                    vals = vals_map[q_prev]
                    probs = probs_map[q_prev]
                    q = int(rng.choice(vals, p=probs))

            # Use marginal method if no transition row
            if q is None:
                p_idx = min(p_cur, len(self._values_by_pos) - 1)
                vals = self._values_by_pos[p_idx]
                probs = self._probs_by_pos[p_idx]
                q = int(rng.choice(vals, p=probs))

            # Clip to valid range
            if q < 0:
                q = 0
            elif q > self.max_quality:
                q = self.max_quality

            qualities[i] = q

        return qualities
