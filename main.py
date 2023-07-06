#!/usr/bin/env python3
"""This module is responsible for building a simplicial complex from a dataset."""

import argparse
import cProfile
from enum import Enum
import io
import logging
from logging import basicConfig, info
from pathlib import Path
from pstats import Stats, SortKey
from shutil import copytree, ignore_patterns
from subprocess import Popen, PIPE, call, check_output
from trace import Trace
import tracemalloc

import numpy as np
from tqdm import tqdm

from configuration.configuration import Configuration
from config_files.model_fitting import SCALAR_PROPERTY_PARAMS_TO_FIT
from config_files.properties_to_test import SCALAR_PROPERTY_PARAMS_TO_TEST
from data_set.factory import load_data
from distribution.approximation import DistributionApproximation
from model.model import Model
from model.factory import create_model, load_default_parameters
from network.property import ScalarNetworkPropertyReport
from optimizer.model_optimizer import ModelOptimizer
from optimizer.factory import create_parameter_options
from reports.network_analysis import analyze_finite_network, analyze_infinite_network_set
from reports.model_analysis import create_model_test_report
from tools.debugger import debugger_is_active


class Mode(Enum):
    """Mode in which the program runs."""

    ANALYSIS = 'analysis'  # to analyze a data set
    FITTING = 'fitting'  # to fit model parameters to data
    TESTING = 'testing'  # to analyze model properties
    EXAMPLE = 'example'  # to analyze an example network from the model


def main(mode: Mode, configuration: Configuration) -> None:
    """Run program - main function."""
    num_of_processes = 1 \
        if debugger_is_active() or configuration.general.runtime_profiling or configuration.general.memory_profiling \
        else configuration.general.num_of_processes

    data_set_type = configuration.data_set_analysis.type_
    model_type = configuration.model.type_

    if mode == Mode.ANALYSIS:
        data_set = load_data(data_set_type)
        (configuration.general.output_dir / 'data').mkdir(parents=True, exist_ok=True)
        analyze_finite_network(
            data_set,
            configuration.data_set_analysis.properties_to_calculate,
            configuration.data_set_analysis.plot,
            configuration.general.output_dir / 'data',
        )
    elif mode == Mode.FITTING:

        model: Model = create_model(model_type)
        model.parameters = load_default_parameters(model_type)

        data_set = load_data(data_set_type)
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
    elif mode == Mode.TESTING:
        info('Network testing started.')

        model: Model = create_model(model_type)
        model.parameters = load_default_parameters(model_type)

        if configuration.model.network_testing.test_against_data_set:
            data_set = load_data(data_set_type)
            model.set_relevant_parameters_from_data_set(data_set)

        scalar_property_distributions = model.simulate(
            scalar_property_params_to_calculate=list(SCALAR_PROPERTY_PARAMS_TO_TEST),
            num_of_simulations=configuration.model.network_testing.num_of_simulations,
            num_of_infinite_networks=configuration.model.network_testing.num_of_infinite_networks,
            num_of_processes=num_of_processes,
        )

        scalar_property_reports: list[ScalarNetworkPropertyReport] = []
        for property_params, distribution in tqdm(
            zip(SCALAR_PROPERTY_PARAMS_TO_TEST, scalar_property_distributions),
            total=len(SCALAR_PROPERTY_PARAMS_TO_TEST),
        ):
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

        model_test_save_dir = (configuration.general.output_dir / 'model_test')
        model_test_save_dir.mkdir(parents=True, exist_ok=True)

        for scalar_network_property_report in scalar_property_reports:
            scalar_property_save_dir = model_test_save_dir / scalar_network_property_report.params.name
            scalar_network_property_report.distributions.save(scalar_property_save_dir)

        model_network_report_figure = create_model_test_report(scalar_property_reports)
        model_network_report_figure.savefig(model_test_save_dir / f'{model_type.name.lower()}_report.png')
        model_network_report_figure.clf()

        info('Network testing finished.')
    elif mode == Mode.EXAMPLE:
        info('Model analysis started.')

        model: Model = create_model(model_type)
        model.parameters = load_default_parameters(model_type)

        model_analysis_save_dir = (configuration.general.output_dir / 'model_analysis')
        model_analysis_save_dir.mkdir(parents=True, exist_ok=True)

        typical_finite_network = model.generate_finite_network()
        analyze_finite_network(
            typical_finite_network,
            configuration.model.analysis.properties_to_calculate,
            model_analysis_save_dir,
        )
        typical_infinite_network_set = model.generate_infinite_network_set(
            configuration.model.analysis.num_of_infinite_networks,
            seed=0,
        )
        analyze_infinite_network_set(
            typical_infinite_network_set,
            configuration.model.analysis.properties_to_calculate,
            model_analysis_save_dir,
        )

        info('Model analysis finished.')
    else:
        raise NotImplementedError('Unknown mode.')


def data_set_required(mode: Mode, configuration: Configuration) -> bool:
    """Check if loading of a data set is required."""
    return \
        mode == Mode.ANALYSIS or \
        mode == Mode.FITTING or \
        (
            mode == Mode.TESTING and
            configuration.model.network_testing.test_against_data_set
        )


def create_parser() -> argparse.ArgumentParser:
    """Create parser for command line arguments."""
    parser = argparse.ArgumentParser(description='Load and prepare input data for training.')
    parser.add_argument('--mode', type=Mode, choices=list(Mode), help='Mode in which the program should be run.')
    parser.add_argument('--config_dir', type=Path, default=Path('config_files'), help='Path to the config directory.')
    return parser


def main_wrapper(params: argparse.Namespace) -> None:
    """Wrap the main function to separate profiling, logging, etc. from actual algorithm."""
    configuration = Configuration()
    configuration.general.output_dir.mkdir(parents=True, exist_ok=True)
    copytree(
        params.config_dir,
        configuration.general.output_dir / 'config_files',
        ignore=ignore_patterns('__pycache__', '__init__.py')
    )

    basicConfig(
        filename=configuration.general.output_dir / 'log.txt',
        filemode='w',
        encoding='utf-8',
        level=getattr(logging, configuration.general.log_level),
        format='%(asctime)s %(levelname)-8s %(filename)s.%(funcName)s.%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
    logging.getLogger('PIL.PngImagePlugin').setLevel(logging.WARNING)
    logging.captureWarnings(True)

    if configuration.general.runtime_profiling:
        runtime_profiler = cProfile.Profile()
        runtime_profiler.enable()

    if configuration.general.memory_profiling:
        tracemalloc.start()

    if configuration.general.tracing:
        tracer = Trace(
            trace=True,
            count=False,
        )
        tracer.runfunc(main, configuration)
        tracer_results = tracer.results()
        tracer_results.write_results(show_missing=True, coverdir=".")

    else:
        main(params.mode, configuration)

    if configuration.general.memory_profiling:
        memory_usage_snapshot = tracemalloc.take_snapshot()
        tracemalloc.stop()

        profile_output_dir = configuration.general.output_dir / 'profile_statistics'
        profile_output_dir.mkdir(parents=True, exist_ok=True)

        project_root_dir = Path(__file__).parent
        filter_ = tracemalloc.Filter(inclusive=True, filename_pattern=f'{str(project_root_dir)}/*')
        memory_statistics = memory_usage_snapshot.filter_traces([filter_]).statistics('lineno')
        memory_usage_snapshot.dump(profile_output_dir / 'memory_profile_results.prof')
        memory_profile_log_file = profile_output_dir / 'memory_profile_results.txt'
        with open(memory_profile_log_file, 'w', encoding='utf-8') as profile_log_file:
            profile_log_file.write('size,count,trace_back\n')
            for stat in memory_statistics:
                profile_log_file.write(f'{stat.size / 2**16},{stat.count},{stat.traceback}\n')
        tracemalloc.clear_traces()

    if configuration.general.runtime_profiling:
        runtime_profiler.disable()

        profile_output_dir = configuration.general.output_dir / 'profile_statistics'
        profile_output_dir.mkdir(parents=True, exist_ok=True)

        stream = io.StringIO()

        runtime_statistics = Stats(runtime_profiler, stream=stream).sort_stats(SortKey.CUMULATIVE)
        runtime_statistics.print_stats(str(configuration.general.root_dir), 100)
        runtime_statistics.dump_stats(profile_output_dir / 'runtime_profile_results.prof')
        runtime_profile_log_file = profile_output_dir / 'runtime_profile_results.txt'
        with open(runtime_profile_log_file, 'w', encoding='utf-8') as profile_log_file:
            profile_log_file.write(stream.getvalue())

        with Popen(
            ('gprof2dot', '-f', 'pstats', profile_output_dir / 'runtime_profile_results.prof'),
            stdout=PIPE
        ) as gprof_process:
            check_output(
                ('dot', '-Tpng', '-o', profile_output_dir / 'runtime_profile_results.png'),
                stdin=gprof_process.stdout
            )
            gprof_process.wait()

        call(('snakeviz', '--server', profile_output_dir / 'runtime_profile_results.prof'))


if __name__ == '__main__':
    argument_parser = create_parser()
    program_arguments = argument_parser.parse_args()
    main_wrapper(program_arguments)
