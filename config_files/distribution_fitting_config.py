"""Parameters for fitting distributions."""

from dataclasses import dataclass


@dataclass
class DistributionFittingConfig:
    """Distribution fitting configuration."""

    power_law_fitting_minimum_value_data: float = 10.
    power_law_fitting_minimum_value_model: float = 10.