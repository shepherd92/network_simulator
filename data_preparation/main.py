#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path
from shutil import copyfile, copytree

import numpy as np
import pandas as pd

from data_preparation.data_dirs.hypergraph import (
    betti_number_directories,
    data_analysis_directories,
    degree_distribution_directories,
    hypothesis_testing_directories,
    model_sample_for_data_directories,
    model_sample_directory,
    simplex_count_directories,
    output_path_root,
)
from data_preparation.latex_interface import (
    create_boxplots,
    create_dataset_parameter_estimates_table,
    create_dataset_properties_table,
    create_histograms_normal_plots,
    create_hypothesis_tests_plot,
    create_hypothesis_tests_table,
    create_persistence_diagrams,
    create_value_counts_log_plot,
)


def main() -> None:
    """Prepare all data for publication."""
    output_path_root.mkdir(parents=True, exist_ok=True)

    prepare_data_analysis_data(data_analysis_directories, hypothesis_testing_directories, output_path_root)
    prepare_model_sample_for_data_sets_data(model_sample_for_data_directories, output_path_root)
    # prepare_model_analysis_sample_plot(model_sample_directory_plot, output_path_root)
    prepare_model_sample_data(model_sample_directory, output_path_root)
    prepare_simulation_degree_distribution_data(degree_distribution_directories, output_path_root)
    prepare_simulation_betti_number_data(betti_number_directories, output_path_root)
    prepare_simulation_simplex_count_data(simplex_count_directories, output_path_root)
    prepare_hypothesis_testing_data(hypothesis_testing_directories, output_path_root)


def prepare_data_analysis_data(
    dataset_dirs: dict[str, Path],
    hypothesis_testing_dirs: dict[str, Path],
    output_dir: Path
) -> None:
    """Prepare all information related to the data analysis."""
    create_data_degree_distributions(dataset_dirs, output_dir)
    create_betti_numbers(dataset_dirs, output_dir)
    _copy_data_network_plots(dataset_dirs, output_dir)
    network_info = _merge_network_info(dataset_dirs, 'data', output_dir)
    create_dataset_properties_table(network_info, output_dir / 'data' / 'latex_tables' / 'network_info.tex')
    create_dataset_parameter_estimates_table(
        dataset_dirs,
        hypothesis_testing_dirs,
        output_dir / 'data' / 'latex_tables' / 'parameter_estimates.tex'
    )
    copy_persistence_data(dataset_dirs, output_dir / 'data' / 'persistence_diagrams')
    create_persistence_diagrams(dataset_dirs, output_dir / 'data' / 'persistence_diagrams')


def copy_persistence_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Copy persistence data."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for dataset_name, directory in directories.items():
        persistence_file_name = directory / 'data' / 'persistence.csv'
        if not persistence_file_name.is_file():
            continue

        persistence = pd.read_csv(persistence_file_name)
        persistence_finite_death = persistence.replace([np.inf], np.nan)
        persistence_finite_death.dropna(subset=['birth', 'death'], inplace=True)
        persistence_finite_death = persistence_finite_death.astype('int32')
        inf_value = int(1.05 * max(
            persistence['birth'].max(),
            persistence_finite_death['death'].max(),
        )) + 1
        persistence.replace([np.inf], inf_value, inplace=True)
        persistence = persistence.astype('int32')
        persistence.drop_duplicates(inplace=True, ignore_index=True)
        out_file_name = output_dir / '_'.join([dataset_name, 'persistence.csv'])
        persistence.to_csv(out_file_name, index=False)


def prepare_model_sample_for_data_sets_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    _copy_model_samples_for_data_sets_network_plots(directories, output_dir)


def prepare_model_sample_data(directory: Path, output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    (output_dir / 'model_sample').mkdir(parents=True, exist_ok=True)

    captions = {
        'total_degree_distribution': '$0$-coface degree distribution',
        'interaction_degree_distribution': r'$\PP$-vertex degree distribution',
        'interaction_dimension_distribution': r"$\PP'$-vertex distribution",
        'ho_degree_distribution_1': '$1$-coface degree distribution',
        # 'ho_degree_distribution_2': 'Triangle--tetrahedron degree distribution',
    }

    for property_name in [
        'total_degree_distribution',
        'interaction_degree_distribution',
        'interaction_dimension_distribution',
        'ho_degree_distribution_1',
        # 'ho_degree_distribution_2',
    ]:
        directory_name = directory / 'model_analysis_finite' / property_name
        if directory_name.is_dir():
            copytree(directory_name, output_dir / 'model_sample' / property_name, dirs_exist_ok=True)
            create_value_counts_log_plot(
                in_dir=directory_name,
                out_file_name=output_dir / 'model_sample' / property_name / 'distribution_latex_figure.tex',
                csv_file_name=Path('data/model_sample') / property_name / 'value_counts.csv',
                x_column=0,
                y_column=1,
                caption=captions[property_name]
            )

    for image_name in [
        'network.png',
        'network_fixed_vertex_positions.png',
    ]:
        file_name = directory / 'model_analysis_finite' / image_name
        if file_name.is_file():
            copyfile(file_name, output_dir / 'model_sample' / image_name)


def prepare_simulation_degree_distribution_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the higher-order degree distributions."""
    _merge_degree_exponents(directories, output_dir)
    _create_boxplots_degree_distributions(directories, output_dir)


def prepare_simulation_simplex_count_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the Betti numbers."""
    for property_name in [
        'num_of_edges_normal_mle',
        'num_of_edges_stable',
        'num_of_triangles_normal_mle',
        'num_of_triangles_stable',
    ]:
        _merge_property_tables(
            directories,
            'model_test',
            property_name,
            output_dir / 'model_test',
            tables_to_merge=[
                'distribution_infos', 'linear_histograms', 'theoretical_pdfs', 'qq_plots',
            ],
        )
        create_histograms_normal_plots(directories, 'model_test', property_name, output_dir / 'model_test')


def prepare_simulation_betti_number_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the Betti numbers."""
    for property_name in [
        'betti_number_0_normal',
        'betti_number_0_stable',
        'betti_number_1_normal',
        'betti_number_1_stable',
        # 'betti_number_2_normal',
        # 'betti_number_2_stable',
    ]:
        _merge_property_tables(
            directories,
            'model_test',
            property_name,
            output_dir / 'model_test',
            tables_to_merge=[
                'distribution_infos', 'linear_histograms', 'theoretical_pdfs', 'qq_plots',
            ],
        )
        create_histograms_normal_plots(directories, 'model_test', property_name, output_dir / 'model_test')


def prepare_hypothesis_testing_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare hypothesis testing data."""
    (output_dir / 'hypothesis_test').mkdir(parents=True, exist_ok=True)
    _merge_model_info(directories, 'model_test', output_dir / 'hypothesis_test')

    figure_captions = {
        'computer_science': 'cs',
        'engineering': 'eess',
        'mathematics': 'math',
        'statistics': 'stat',
    }

    table_captions = {
        'num_of_edges_stable': 'Results of the hypothesis test for the number of edges',
        'num_of_triangles_stable': 'Results of the hypothesis test for the number of triangles',
        'betti_number_1_stable': 'Results of the hypothesis test for the first Betti numbers',
    }

    table_labels = {
        'num_of_edges_stable': 'tab:num_of_edges_hypothesis_tests',
        'num_of_triangles_stable': 'tab:num_of_triangles_hypothesis_tests',
        'betti_number_1_stable': 'tab:betti_number_1_hypothesis_tests',
    }

    for property_name in [
        # 'vertex_degree_exponent',
        # 'edge_degree_exponent',
        # 'triangle_degree_exponent',
        # 'num_of_edges_normal_mle',
        'num_of_edges_stable',
        # 'num_of_triangles_normal_mle',
        'num_of_triangles_stable',
        # 'betti_number_0_normal',
        # 'betti_number_0_stable',
        # 'betti_number_1_normal',
        'betti_number_1_stable',
        # 'betti_number_2_normal',
        # 'betti_number_2_stable',
    ]:
        _merge_property_tables(
            directories,
            'model_test',
            property_name,
            output_dir / 'hypothesis_test',
            tables_to_merge=[
                'distribution_infos', 'test_results',
                'linear_histograms', 'theoretical_pdfs',
            ],
        )

        for dataset_name, directory in directories.items():
            create_hypothesis_tests_plot(
                in_dir=directory / 'model_test' / property_name,
                out_file_name=output_dir / 'hypothesis_test' / property_name /
                f'hypothesis_test_{dataset_name}_latex_figure.tex',
                csv_file_dir=Path('data/hypothesis_test') / property_name,
                dataset_name=dataset_name,
                caption=figure_captions[dataset_name],
            )

        create_hypothesis_tests_table(
            directories,
            property_name,
            table_captions[property_name],
            table_labels[property_name],
            output_dir / 'hypothesis_test' / 'latex_tables' / f'{property_name}_hypothesis_test_table.tex'
        )


def _merge_degree_exponents(directories: dict[str, Path], output_dir: Path):
    for property_name in [
        'vertex_degree_exponent',
        'edge_degree_exponent',
        # 'triangle_degree_exponent',
        'interaction_vertex_degree_exponent',
        'vertex_interaction_degree_exponent',
        'edge_interaction_degree_exponent',
        # 'triangle_interaction_degree_exponent',
    ]:
        _merge_property_tables(
            directories,
            'model_test',
            property_name,
            output_dir / 'model_test',
            tables_to_merge=[
                'distribution_infos', 'quantiles',
            ],
        )


def _create_boxplots_degree_distributions(directories: dict[str, Path], output_dir: Path):

    GAMMA = 0.7
    GAMMA_PRIME = 0.2

    for property_name in [
        'vertex_degree_exponent',
        'edge_degree_exponent',
        # 'triangle_degree_exponent',
        'interaction_vertex_degree_exponent',
        'vertex_interaction_degree_exponent',
        'edge_interaction_degree_exponent',
        # 'triangle_interaction_degree_exponent',
    ]:
        caption = ''
        theoretical_value = 0.
        if property_name == 'vertex_degree_exponent':
            caption = r'$0$-coface degree exponent distribution'
            theoretical_value = np.nan
        elif property_name == 'edge_degree_exponent':
            caption = r'$1$-coface degree exponent distribution'
            theoretical_value = np.nan
        elif property_name == 'triangle_degree_exponent':
            caption = r'$2$-coface degree exponent distribution'
            theoretical_value = np.nan
        elif property_name == 'interaction_vertex_degree_exponent':
            caption = r"$\PP'$-vertex degree exponent distribution"
            theoretical_value = 1. + 1. / GAMMA_PRIME
        elif property_name == 'vertex_interaction_degree_exponent':
            caption = r'$\PP$-vertex degree exponent distribution'
            theoretical_value = 1. + 1. / GAMMA
        elif property_name == 'edge_interaction_degree_exponent':
            caption = r'$1$-coface degree exponent distribution'
            theoretical_value = 2. / GAMMA
        elif property_name == 'triangle_interaction_degree_exponent':
            caption = r'$2$-coface degree exponent distribution'
            theoretical_value = 3. / GAMMA - 1.

        minimum_degrees = set([key[1] for key in directories.keys()])
        for minimum_degree in minimum_degrees:
            directories_filtered = {
                key: directory for key, directory in directories.items() if key[1] == minimum_degree
            }
            create_boxplots(
                directories_filtered,
                'model_test',
                property_name,
                output_dir / 'model_test',
                caption,
                theoretical_value,
                suffix=f'{minimum_degree}'
            )


def create_data_degree_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create higher order degree distributions."""
    for property_name in [
        'total_degree_distribution',
        'interaction_degree_distribution',
        'ho_degree_distribution_1',
        # 'ho_degree_distribution_2',
        # 'ho_degree_distribution_3',
        # 'simplex_dimension_distribution',
        'interaction_dimension_distribution',
        # 'facet_dimension_distribution',
    ]:
        _merge_property_tables(
            directories,
            'data',
            property_name,
            output_dir / 'data',
            tables_to_merge=[
                'distribution_infos', 'value_counts',
            ],
        )
        for dataset_name, directory in directories.items():
            create_value_counts_log_plot(
                in_dir=directory / 'data' / property_name,
                out_file_name=output_dir / 'data' / property_name / f'{dataset_name}_distribution_latex_figure.tex',
                csv_file_name=Path('data/data') / property_name / 'value_counts.csv',
                x_column='value',
                y_column=dataset_name,
                caption='',
            )


def create_betti_numbers(directories: dict[str, Path], output_dir: Path) -> None:
    """Create Betti numbers."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        in_file_name = directory / 'data' / 'betti_numbers.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name, index_col=0, dtype=int)
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_dir / 'data' / 'betti_numbers.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.index.name = 'dimension'
        merged_data_frame.fillna(0, inplace=True)
        merged_data_frame.astype(int).to_csv(out_file_name)


def _merge_network_info(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> pd.DataFrame:
    """Create and merge info."""
    (output_path / 'data').mkdir(parents=True, exist_ok=True)

    data_frames: list[pd.DataFrame] = []
    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        current_info = pd.read_csv(directory / subdirectory_name / 'network_info.csv')
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.to_csv(output_path / 'data' / 'network_info.csv')

    return merged_data_frame


def _copy_data_network_plots(directories: dict[str, Path], output_path: Path) -> None:
    """Copy network plots to the output path."""
    (output_path / 'data').mkdir(parents=True, exist_ok=True)
    for dataset_name, directory in directories.items():

        file_name = directory / 'data' / 'network.png'
        if file_name.is_file():
            copyfile(file_name, output_path / 'data' / f'{dataset_name}_network.png')


def _copy_model_samples_for_data_sets_network_plots(
    directories: dict[str, Path],
    output_path: Path
) -> None:
    (output_path / 'data').mkdir(parents=True, exist_ok=True)
    for dataset_name, directory in directories.items():
        file_name_original = directory / 'model_analysis_finite' / 'network.png'
        if file_name_original.is_file():
            copyfile(file_name_original, output_path / 'data' / f'{dataset_name}_model_sample_network.png')

        file_name_fixed = directory / 'model_analysis_finite' / 'network_fixed_vertex_positions.png'
        if file_name_fixed.is_file():
            copyfile(file_name_fixed, output_path / 'data' / f'{dataset_name}_network_fixed_vertex_positions.png')


def _merge_property_tables(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
    tables_to_merge: list[str],
) -> None:
    """Merge property tables."""
    (output_dir / property_name).mkdir(parents=True, exist_ok=True)
    if 'distribution_infos' in tables_to_merge:
        _merge_distribution_info(directories, input_subdir, property_name, output_dir)
    if 'linear_histograms' in tables_to_merge:
        _merge_histograms(directories, input_subdir, property_name, output_dir)
    if 'value_sequences' in tables_to_merge:
        _merge_value_sequences(directories, input_subdir, property_name, output_dir)
    if 'value_counts' in tables_to_merge:
        _merge_value_counts(directories, input_subdir, property_name, output_dir)
    if 'qq_plots' in tables_to_merge:
        _merge_qq_plots(directories, input_subdir, property_name, output_dir)
    if 'confidence_intervals' in tables_to_merge:
        _merge_confidence_intervals(directories, input_subdir, property_name, output_dir)
    if 'quantiles' in tables_to_merge:
        _merge_quantiles(directories, input_subdir, property_name, output_dir)
    if 'theoretical_pdfs' in tables_to_merge:
        _merge_theoretical_pdfs(directories, input_subdir, property_name, output_dir)
    if 'test_results' in tables_to_merge:
        _merge_test_results(directories, input_subdir, property_name, output_dir)


def _merge_model_info(
    directories: dict[str, Path],
    input_subdir: str,
    output_path: Path,
) -> None:
    """Create and merge model info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        in_file_name = directory / input_subdir / 'model_info.csv'
        if not in_file_name.is_file():
            continue

        current_info = pd.read_csv(in_file_name)
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        out_file_name = output_path / 'model_info.csv'
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.index.name = 'name'
        merged_data_frame.to_csv(out_file_name)


def _merge_distribution_info(
    directories: dict[str | tuple[int, ...], Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Create and merge network info."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'distribution_info.csv'
        if not in_file_name.is_file():
            continue

        current_info = pd.read_csv(in_file_name)
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        out_file_name = output_path / property_name / 'distribution_info.csv'
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.index.name = 'name'
        merged_data_frame.to_csv(out_file_name)


def _merge_histograms(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Merge linearly binned histograms."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'histogram_linear.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name)
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_bin_left_limit', f'{dataset_name}_value']
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'linear_histograms.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(out_file_name, index=False)


def _merge_value_sequences(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'value_sequence.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name, header=None,)
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'value_sequences.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(out_file_name, index=False)


def _merge_value_counts(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'value_counts.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name, index_col=0, header=None, dtype=int)
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'value_counts.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.index.name = 'value'
        merged_data_frame.fillna(0, inplace=True)
        merged_data_frame.sort_index(inplace=True)
        merged_data_frame.astype(int).to_csv(out_file_name)


def _merge_qq_plots(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'qq_plot_points.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name,)
        if current_distribution.empty:
            continue

        current_distribution.columns = ['_'.join([dataset_name, column_name])
                                        for column_name in current_distribution.columns]
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'qq_plot.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(out_file_name, index=False)


def _merge_confidence_intervals(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'confidence_intervals.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name, index_col=0)
        if current_distribution.empty:
            continue

        current_distribution.columns = [
            f'{dataset_name}_empirical_lower',
            f'{dataset_name}_empirical_upper',
            f'{dataset_name}_theoretical_lower',
            f'{dataset_name}_theoretical_upper',
        ]
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'confidence_intervals.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.index.name = 'name'
        merged_data_frame.to_csv(out_file_name)


def _merge_quantiles(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'quantiles.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name, index_col=0)
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_empirical', f'{dataset_name}_theoretical']
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'quantiles.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.index.name = 'name'
        merged_data_frame.to_csv(out_file_name)


def _merge_theoretical_pdfs(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Merge theoretical probability density functions."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'pdfs.csv'
        if not in_file_name.is_file():
            continue

        current_distribution = pd.read_csv(in_file_name, usecols=[0, 2])
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_value', f'{dataset_name}_pdf']
        data_frames.append(current_distribution)

    if data_frames:
        out_file_name = output_path / property_name / 'theoretical_pdf.csv'
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(out_file_name, index=False)


def _merge_test_results(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Merge test results table."""
    data_frames: list[pd.DataFrame] = []

    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(map(str, key))
        in_file_name = directory / input_subdir / property_name / 'test_results.csv'
        if not in_file_name.is_file():
            continue

        current_info = pd.read_csv(in_file_name)
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        out_file_name = output_path / property_name / 'test_results.csv'
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.index.name = 'name'
        merged_data_frame.to_csv(out_file_name)


if __name__ == '__main__':
    main()
