#!/usr/bin/env python
"""Configuration class."""

from datetime import datetime, UTC
import os

from config_files.general_config import GeneralConfig
from config_files.data_set_config import DataSetConfig
from config_files.model_config import ModelConfig


class Configuration:
    """Class handling configuration parameters."""

    def __init__(self) -> None:
        """Construct a high level data node."""
        self.general = GeneralConfig()
        self.data_set_analysis = DataSetConfig()
        self.model = ModelConfig()

        self.general.num_of_processes = self._determine_num_of_processes()

        now = datetime.now(UTC).strftime('%Y%m%d_%H%M%S')
        self.general.output_dir = self.general.output_dir / now

    def _determine_num_of_processes(self) -> int:

        cpu_count = os.cpu_count()
        if not cpu_count:
            num_of_processes = 1
        elif self.general.num_of_processes:
            num_of_processes = min(cpu_count, self.general.num_of_processes)
            num_of_processes = max(num_of_processes, 1)
        else:
            num_of_processes = max(cpu_count - 2, 1)

        return num_of_processes
