#!/usr/bin/env python3
"""This module is responsible for building a simplicial complex from a dataset."""

import argparse
import cProfile
import io
import logging
from logging import basicConfig, info
from pathlib import Path
from pstats import Stats, SortKey
from shutil import copytree, ignore_patterns
from subprocess import Popen, PIPE, call, check_output

import numpy as np

from configuration import Configuration
from config_files.properties.scalar_properties_to_fit import SCALAR_PROPERTY_PARAMS_TO_FIT
from config_files.properties.scalar_properties_to_test import SCALAR_PROPERTY_PARAMS_TO_TEST
from data_set.factory import load_data
from distribution.approximation import DistributionApproximation
from model.model import Model
from model.factory import create_model, load_default_parameters
from network.property import ScalarNetworkPropertyReport
from optimizer.model_optimizer import ModelOptimizer
from optimizer.factory import create_parameter_options
from reports.network_analysis import analyze_network
from reports.model_analysis import create_model_test_report
from tools.debugger import debugger_is_active


def main(configuration: Configuration) -> None:
    """Run program - main function."""
    num_of_processes = 1 \
        if debugger_is_active() or configuration.general.enable_profiling \
        else configuration.general.num_of_processes

    data_set_type = configuration.data_set.type_
    model_type = configuration.model.type_

    if data_set_required(configuration):
        data_set = load_data(data_set_type)

    if configuration.data_set.analysis.enable:

        (configuration.general.directories.output / 'data').mkdir(parents=True, exist_ok=True)
        analyze_network(
            data_set,
            configuration.data_set.analysis.properties_to_calculate,
            configuration.general.directories.output / 'data' /
            f'{data_set_type.name.lower()}_report.png',
        )

    model: Model = create_model(model_type)
    model.parameters = load_default_parameters(model_type)

    if configuration.model.fitting.enable:

        parameter_options = create_parameter_options(model_type, data_set)
        optimizer = ModelOptimizer(model, method='BFGS')

        target_values = [
            data_set.calc_scalar_property(one_scalar_property_params)
            for one_scalar_property_params in SCALAR_PROPERTY_PARAMS_TO_FIT
        ]

        info(f'Model {model_type} built with parameters {parameter_options}.')

        model.parameters = optimizer.fit(
            parameter_options=parameter_options,
            target_property_params=list(SCALAR_PROPERTY_PARAMS_TO_FIT),
            target_values=target_values,
            num_of_processes=num_of_processes,
            options={
                'maxiter': 10
            }
        )

    if configuration.model.network_testing.enable:
        info('Network testing started.')

        if configuration.model.network_testing.test_against_data_set:
            model.set_relevant_parameters_from_data_set(data_set)

        scalar_property_distributions = model.simulate(
            scalar_property_params_to_calculate=list(SCALAR_PROPERTY_PARAMS_TO_TEST),
            num_of_simulations=configuration.model.network_testing.num_of_simulations,
            num_of_processes=num_of_processes,
        )

        scalar_property_reports: list[ScalarNetworkPropertyReport] = []
        for property_params, distribution in zip(SCALAR_PROPERTY_PARAMS_TO_TEST, scalar_property_distributions):
            distribution_pair = DistributionApproximation(
                distribution,
                property_params.theoretical_approximation_type
            )
            distribution_pair.fit(property_params.fitting_parameters)

            if configuration.model.network_testing.test_against_data_set:
                data_set_value = data_set.calc_scalar_property(property_params)
            else:
                data_set_value = np.nan

            test_results = distribution_pair.run_test(data_set_value)

            scalar_property_report = ScalarNetworkPropertyReport(
                params=property_params,
                distributions=distribution_pair,
                test_results=test_results,
                data_point=data_set_value,
            )

            scalar_property_reports.append(scalar_property_report)

        model_test_save_dir = (configuration.general.directories.output / 'model_test')
        model_test_save_dir.mkdir(parents=True, exist_ok=True)

        for scalar_network_property_report in scalar_property_reports:
            scalar_property_save_dir = model_test_save_dir / scalar_network_property_report.params.name
            scalar_property_save_dir.mkdir(parents=True, exist_ok=True)

            pdfs = scalar_network_property_report.distributions.get_pdfs()
            pdfs.to_csv(scalar_property_save_dir / 'pdfs.csv', float_format='%.9f')
            value_sequence = scalar_network_property_report.distributions.empirical.value_sequence
            np.savetxt(scalar_property_save_dir / 'value_sequence.csv', value_sequence, delimiter=',')

            confidence_levels = [0.9, 0.95, 0.99]
            confidence_intervals = \
                scalar_network_property_report.distributions.get_confidence_intervals(confidence_levels)
            confidence_intervals.to_csv(scalar_property_save_dir / 'confidence_intervals.csv', float_format='%.4f')

        model_network_report_figure = create_model_test_report(scalar_property_reports)
        model_network_report_figure.savefig(model_test_save_dir / f'{model_type.name.lower()}_report.png')
        model_network_report_figure.clf()

        info('Network testing finished.')

    if configuration.model.analysis.enable:
        info('Model analysis started.')

        typical_network = model.generate_finite_network()
        analyze_network(
            typical_network,
            configuration.model.analysis.properties_to_calculate,
            configuration.general.directories.output /
            f'{model_type.name.lower()}_typical_network.png',
        )

        info('Model analysis finished.')


def data_set_required(configuration: Configuration) -> bool:
    """Check if loading of a data set is required."""
    return \
        configuration.data_set.analysis.enable or \
        configuration.model.fitting.enable or \
        (
            configuration.model.network_testing.enable and
            configuration.model.network_testing.test_against_data_set
        )


def create_parser() -> argparse.ArgumentParser:
    """Create parser for command line arguments."""
    parser = argparse.ArgumentParser(description='Load and prepare input data for training.')
    parser.add_argument('--config_dir', type=Path, default=Path('config_files'), help='Path to the config directory')
    return parser


def main_wrapper(params: argparse.Namespace) -> None:
    """Wrap the main function to separate profiling, logging, etc. from actual algorithm."""
    configuration = Configuration()
    configuration.load(params.config_dir / 'config.json')
    configuration.general.directories.output.mkdir(parents=True, exist_ok=True)
    copytree(
        params.config_dir,
        configuration.general.directories.output / 'config_files',
        ignore=ignore_patterns('__pycache__', '__init__.py')
    )

    basicConfig(
        filename=configuration.general.directories.output / 'log.txt',
        filemode='w',
        encoding='utf-8',
        level=getattr(logging, configuration.general.log_level)
    )
    logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
    logging.getLogger('PIL.PngImagePlugin').setLevel(logging.WARNING)
    logging.captureWarnings(True)

    if configuration.general.enable_profiling:
        profiler = cProfile.Profile()
        profiler.enable()

    main(configuration)

    if configuration.general.enable_profiling:
        profile_output_dir = configuration.general.directories.output
        profiler.disable()
        stream = io.StringIO()
        statistics = Stats(profiler, stream=stream).sort_stats(SortKey.CUMULATIVE)
        statistics.print_stats(str(configuration.general.directories.root), 100)
        statistics.dump_stats(profile_output_dir / 'profile_results.prof')
        profile_log_path = profile_output_dir / 'profile_results.txt'

        with open(profile_log_path, 'w', encoding='utf-8') as profile_log_file:
            profile_log_file.write(stream.getvalue())

        with Popen(
            ('gprof2dot', '-f', 'pstats', profile_output_dir / 'profile_results.prof'),
            stdout=PIPE
        ) as gprof_process:
            check_output(
                ('dot', '-Tpng', '-o', profile_output_dir / 'profile_results.png'),
                stdin=gprof_process.stdout
            )
            gprof_process.wait()

        call(('snakeviz', '--server', profile_output_dir / 'profile_results.prof'))


if __name__ == '__main__':
    argument_parser = create_parser()
    program_arguments = argument_parser.parse_args()
    main_wrapper(program_arguments)
