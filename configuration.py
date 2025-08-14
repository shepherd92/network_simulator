#!/usr/bin/env python
"""Configuration class."""

from datetime import datetime, UTC

from config_files.general_config import GeneralConfig
from config_files.data_set_config import DataSetConfig
from config_files.distribution_fitting_config import DistributionFittingConfig
from config_files.model_config import ModelConfig


class Configuration:
    """Class handling configuration parameters."""

    def __init__(self) -> None:
        """Construct a high level data node."""
        self.general = GeneralConfig()
        self.data_set_analysis = DataSetConfig()
        self.distribution_fitting = DistributionFittingConfig()
        self.model = ModelConfig()
        now = datetime.now(UTC).strftime('%Y%m%d_%H%M%S')
        self.general.output_dir = self.general.output_dir / now
