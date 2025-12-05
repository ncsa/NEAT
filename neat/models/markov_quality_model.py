"""
Empirical position-specific quality score model for NEAT.

The model consists of an initial distribution giving the probability of observing each
quality score at the first position of a read and a sequence of per-position empirical
distributions giving the probability of observing each quality score at each subsequent
position.

At generation time, the model samples a quality score independently at
each position from the corresponding empirical distribution.

Although this model exists outside of ``neat.models.error_models``, it
parallels the ``TraditionalQualityModel`` class from NEAT’s core and can
be used interchangeably when constructing quality models.
"""

from dataclasses import dataclass
from typing import Dict, List

import numpy as np

__all__ = ["MarkovQualityModel"]

@dataclass
class MarkovQualityModel:
    """Empirical, position-specific quality score model.

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
    """

    initial_distribution: Dict[int, float]
    position_distributions: List[Dict[int, float]]
    max_quality: int
    read_length: int

    def __post_init__(self) -> None:

        # Normalize the initial distribution
        total_init = float(sum(self.initial_distribution.values()))

        if total_init <= 0:
            raise ValueError(
                "Initial distribution must have positive probability mass."
            )

        self.initial_distribution = {
            int(k): float(v) / total_init
            for k, v in self.initial_distribution.items()
        }

        # Normalize each position-specific distribution
        norm_positions: List[Dict[int, float]] = []

        for dist in self.position_distributions:
            total = float(sum(dist.values()))

            if total <= 0:
                # Fall back to a single mass at quality 0
                norm_positions.append({0: 1.0})
                continue

            norm_positions.append(
                {int(k): float(v) / total for k, v in dist.items()}
            )

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
        """Generate a synthetic quality score array using this empirical model.

        Parameters
        ----------
        model_read_length:
            The length of reads used during training.
        length:
            Desired length of the returned quality array.
        rng:
            A :class:`numpy.random.Generator` instance used for sampling.

        Returns
        -------
        numpy.ndarray
            An array of integers representing quality scores.
        """
        if length <= 0:
            return np.zeros(0, dtype=int)

        qualities = np.zeros(length, dtype=int)

        # Sample initial quality from the starting distribution
        qualities[0] = rng.choice(self._init_scores, p=self._init_probs)

        # Sample each subsequent position independently from its empirical per-position distribution
        for i in range(1, length):
            p_idx = self._position_index_for_length(i, length)
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
