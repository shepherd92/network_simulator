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
import random
from shutil import copytree, ignore_patterns
from subprocess import Popen, PIPE, call, check_output
from trace import Trace
import tracemalloc

import numpy as np
from tqdm import tqdm

from configuration import Configuration
from config_files.properties_to_test import (
    get_finite_scalar_property_params,
    get_infinite_scalar_property_params,
)
from data_set.factory import load_data
from distribution.approximation import DistributionApproximation
from model.model import Model
from model.factory import create_model, load_default_parameters
from network.property import (
    DerivedNetworkProperty,
    ScalarNetworkPropertyReport,
)
from reports.data_set_network_analysis import analyze_data_set_network
from reports.model_example_network_analysis import (
    analyze_model_example_finite_network,
    analyze_model_example_infinite_network_set,
    create_infinite_network_plots,
)
from reports.model_testing import create_model_test_report


class Mode(Enum):
    """Mode in which the program runs."""

    ANALYSIS = 'analysis'  # to analyze a data set
    FITTING = 'fitting'  # to fit model parameters to data
    TESTING = 'testing'  # to analyze model properties
    EXAMPLE = 'example'  # to analyze an example network from the model


def main(mode: Mode, configuration: Configuration) -> None:
    """Run program - main function."""
    data_set_type = configuration.data_set_analysis.type_
    model_type = configuration.model.type_
    seed = configuration.general.seed if configuration.general.seed else random.randint(0, 2**31 - 1)

    if mode == Mode.ANALYSIS:
        data_set = load_data(data_set_type)
        (configuration.general.output_dir / 'data').mkdir(parents=True, exist_ok=True)
        analyze_data_set_network(
            data_set,
            configuration.data_set_analysis.properties_to_calculate,
            configuration.data_set_analysis.plot,
            configuration.general.output_dir / 'data',
        )
    elif mode == Mode.TESTING:
        info('Network testing started.')

        model: Model = create_model(model_type)
        model.parameters = load_default_parameters(model_type)

        if configuration.model.set_params_from_data_set:
            data_set = load_data(data_set_type)
            model.set_relevant_parameters_from_data_set(data_set)

        scalar_property_params = get_finite_scalar_property_params(model.parameters.gamma) \
            if configuration.model.network_testing.mode == Model.Mode.FINITE \
            else get_infinite_scalar_property_params(model.parameters.gamma)

        scalar_property_distributions = model.simulate(
            configuration.model.network_testing.mode,
            scalar_property_params_to_calculate=list(scalar_property_params),
            num_of_simulations=configuration.model.network_testing.num_of_simulations,
            initial_seed=seed,
            num_of_infinite_networks_per_simulation=configuration.model.network_testing.num_of_infinite_networks,
        )

        scalar_property_reports: list[ScalarNetworkPropertyReport] = []
        for property_params in tqdm(scalar_property_params):
            assert isinstance(property_params, DerivedNetworkProperty)
            distribution_pair = DistributionApproximation(
                scalar_property_distributions[property_params.name],
                property_params.theoretical_approximation_type
            )

            distribution_pair.fit(property_params.fitting_parameters)

            data_set_value = data_set.calc_scalar_property(property_params) \
                if configuration.model.set_params_from_data_set \
                else np.nan

            test_results = distribution_pair.run_test(data_set_value)

            scalar_property_report = ScalarNetworkPropertyReport(
                params=property_params,
                distributions=distribution_pair,
                test_results=test_results,
                data_point=data_set_value,
            )

            scalar_property_reports.append(scalar_property_report)

        model_test_save_dir = configuration.general.output_dir / 'model_test'
        model_test_save_dir.mkdir(parents=True, exist_ok=True)
        model.save_info(model_test_save_dir / 'model_info.csv')

        for scalar_network_property_report in scalar_property_reports:
            scalar_property_save_dir = model_test_save_dir / scalar_network_property_report.params.name
            scalar_network_property_report.distributions.save(scalar_property_save_dir)
            scalar_network_property_report.test_results.save(scalar_property_save_dir)

        model_network_report_figure = create_model_test_report(scalar_property_reports)
        model_network_report_figure.savefig(model_test_save_dir / f'{model_type.name.lower()}_report.png')
        model_network_report_figure.clf()

        info('Network testing finished.')
    elif mode == Mode.EXAMPLE:
        info('Model analysis started.')

        model: Model = create_model(model_type)
        model.parameters = load_default_parameters(model_type)

        if configuration.model.set_params_from_data_set:
            data_set = load_data(data_set_type)
            model.set_relevant_parameters_from_data_set(data_set)

        if configuration.model.analysis.finite.enable:
            model_analysis_save_dir = configuration.general.output_dir / 'model_analysis_finite'
            model_analysis_save_dir.mkdir(parents=True, exist_ok=True)
            typical_finite_network = model.generate_finite_network(seed)
            analyze_model_example_finite_network(
                typical_finite_network,
                configuration.model.analysis.finite.properties_to_calculate,
                configuration.model.analysis.finite.plot,
                model_analysis_save_dir,
            )

        if configuration.model.analysis.infinite.enable:
            model_analysis_save_dir = configuration.general.output_dir / 'model_analysis_infinite'
            model_analysis_save_dir.mkdir(parents=True, exist_ok=True)
            typical_infinite_network = model.generate_infinite_network(
                configuration.model.analysis.infinite.typical_mark,
                seed=seed
            )
            network_with_most_nodes = typical_infinite_network
            create_infinite_network_plots(network_with_most_nodes, model_analysis_save_dir)

        if configuration.model.analysis.infinite_set.enable:
            assert configuration.model.analysis.infinite_set.num_of_infinite_networks > 0
            model_analysis_save_dir = configuration.general.output_dir / 'model_analysis_infinite_set'
            model_analysis_save_dir.mkdir(parents=True, exist_ok=True)
            typical_infinite_network_set = model.generate_infinite_network_set(
                configuration.model.analysis.infinite_set.num_of_infinite_networks,
                seed=seed,
            )
            analyze_model_example_infinite_network_set(
                typical_infinite_network_set,
                configuration.model.analysis.properties_to_calculate_infinite,
                model_analysis_save_dir,
            )

        info('Model analysis finished.')
    else:
        raise NotImplementedError('Unknown mode.')

    print('\nDone.')


def data_set_required(mode: Mode, configuration: Configuration) -> bool:
    """Check if loading of a data set is required."""
    return \
        mode == Mode.ANALYSIS or \
        mode == Mode.FITTING or \
        (
            mode == Mode.TESTING and
            configuration.model.set_params_from_data_set
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
