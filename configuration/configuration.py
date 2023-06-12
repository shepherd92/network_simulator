#!/usr/bin/env python
"""Configuration class."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
import os
from pathlib import Path
from typing import NamedTuple

from model.model import Model
from data_set.data_set import DataSet
from network.property import BaseNetworkProperty


class Directories(NamedTuple):
    """Directory configuration."""

    root: Path = Path()
    output: Path = Path()


class ModelFitting(NamedTuple):
    """Model fitting configuration."""


class ModelTesting(NamedTuple):
    """Model testing configuration for networks."""

    test_against_data_set: bool = False
    num_of_simulations: int = 0
    num_of_infinite_networks: int = 0


class ModelAnalysis(NamedTuple):
    """Model analysis configuration."""

    component_index_from_largest: int = 0
    num_of_infinite_networks: int = 0
    properties_to_calculate: list[BaseNetworkProperty.Type] = []


class ModelStatistics(NamedTuple):
    """Model statistics configuration."""

    num_of_simulations: int = 0


@dataclass
class Configuration:
    """Class handling configuration parameters."""

    class General(NamedTuple):
        """Configuration of the execution of the program."""

        log_level: str = 'WARNING'
        runtime_profiling: bool = False
        memory_profiling: bool = False
        tracing: bool = False
        num_of_processes: int = 1
        directories: Directories = Directories()

    class DataSetAnalysis(NamedTuple):
        """Data set configuration."""

        type_: DataSet.Type = DataSet.Type.INVALID
        plot: bool = False
        properties_to_calculate: list[BaseNetworkProperty.Type] = []

    class Model(NamedTuple):
        """Model configuration."""

        type_: Model.Type = Model.Type.INVALID
        fitting: ModelFitting = ModelFitting()
        network_testing: ModelTesting = ModelTesting()
        analysis: ModelAnalysis = ModelAnalysis()

    general = General()
    data_set_analysis = DataSetAnalysis()
    model = Model()

    def load(self, path: Path) -> None:
        """Construct a high level data node."""
        with open(path, encoding='utf-8') as config_file:
            params = json.load(config_file)

        now = datetime.utcnow().strftime('%Y%m%d_%H%M%S')
        num_of_processes = Configuration._determine_num_of_processes(int(params['general']['num_of_processes']))

        self.general = Configuration.General(
            log_level=str(params['general']['log_level']),
            runtime_profiling=bool(params['general']['runtime_profiling']),
            memory_profiling=bool(params['general']['memory_profiling']),
            tracing=bool(params['general']['tracing']),
            num_of_processes=num_of_processes,
            directories=Directories(
                root=Path(params['general']['directories']['root']),
                output=Path(params['general']['directories']['output']) / now,
            ),
        )

        self.data_set_analysis = Configuration.DataSetAnalysis(
            type_=DataSet.Type[params['data_set_analysis']['type']],
            plot=bool(params['data_set_analysis']['plot']),
            properties_to_calculate=[
                BaseNetworkProperty.Type[network_property]
                for network_property in params['data_set_analysis']['properties']
            ],
        )

        self.model = Configuration.Model(
            type_=Model.Type[str(params['model']['type'])],
            fitting=ModelFitting(
            ),
            network_testing=ModelTesting(
                test_against_data_set=bool(params['model']['testing']['test_against_data_set']),
                num_of_simulations=int(params['model']['testing']['num_of_simulations']),
                num_of_infinite_networks=int(params['model']['testing']['num_of_infinite_networks']),
            ),
            analysis=ModelAnalysis(
                component_index_from_largest=int(
                    params['model']['analysis']['component_index_from_largest']
                ),
                num_of_infinite_networks=int(params['model']['analysis']['num_of_infinite_networks']),
                properties_to_calculate=[
                    BaseNetworkProperty.Type[network_property]
                    for network_property in params['model']['analysis']['properties']
                ],
            ),
        )

    @staticmethod
    def _determine_num_of_processes(desired_num_of_processes: int | None) -> int:

        cpu_count = os.cpu_count()
        if not cpu_count:
            num_of_processes = 1
        elif desired_num_of_processes:
            num_of_processes = min(cpu_count, desired_num_of_processes)
            num_of_processes = max(num_of_processes, 1)
        else:
            num_of_processes = max(cpu_count - 2, 1)

        return num_of_processes