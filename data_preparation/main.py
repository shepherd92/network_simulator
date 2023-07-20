#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path
from shutil import copyfile

import pandas as pd

input_base_path = Path('../output')
output_path = Path('../output/_prepared_data')

model_sample_directory = input_base_path / '20230713_204802'

data_analysis_directories = {
    'computer_science': input_base_path / '20230714_131929',
    'biology':          input_base_path / '20230714_124427',
    'economics':        input_base_path / '20230714_124535',
    'engineering':      input_base_path / '20230714_131441',
    'finance':          input_base_path / '20230714_124618',
    'mathematics':      input_base_path / '20230714_125222',
    'statistics':       input_base_path / '20230714_131150',
}

model_sample_for_data_directories = {
    'computer_science': input_base_path / '20230714_173606',
    'biology':          input_base_path / '20230714_185135',
    'economics':        input_base_path / '20230714_190142',
    'engineering':      input_base_path / '20230714_190208',
    'finance':          input_base_path / '20230714_190228',
    'mathematics':      input_base_path / '20230714_190308',
    'statistics':       input_base_path / '20230714_190334',
}

degree_distribution_directories = {
    '10':       input_base_path / '20230614_165042',
    '100':      input_base_path / '20230614_165118',
    '1000':     input_base_path / '20230614_165434',
    '10000':    input_base_path / '20230614_170118',
    '100000':   input_base_path / '20230614_185624',
    'infinite': input_base_path / '20230615_074249',
}

simplex_count_directories = {
    # the key is the value of gamma * 100
    '25': input_base_path / '20230615_103808',
    '40': input_base_path / '20230703_101414',
    '45': input_base_path / '20230703_130440',
    '50': input_base_path / '20230616_061230',
    '60': input_base_path / '20230616_193506',
    '75': input_base_path / '20230615_204859',
}

betti_number_directories = {
    # the key is the value of gamma * 100
    # '10': input_base_path / '20230608_191114',
    # '20': input_base_path / '20230608_191120',
    '25': input_base_path / '20230609_212258',
    # '30': input_base_path / '20230608_191124',
    # '40': input_base_path / '20230608_191128',
    '50': input_base_path / '20230609_212245',
    # '60': input_base_path / '20230608_205831',
    # '70': input_base_path / '20230608_205839',
    '75': input_base_path / '20230609_142002',
    # '80': input_base_path / '20230608_205845',
    # '90': input_base_path / '20230608_205850',
}

hypothesis_testing_directories = {
    'computer_science': input_base_path / '20230712_075521',
    'biology':          input_base_path / '20230710_164059',
    'economics':        input_base_path / '20230710_162805',
    'engineering':      input_base_path / '20230711_110334',
    'finance':          input_base_path / '20230710_162239',
    'mathematics':      input_base_path / '20230710_200618',
    'statistics':       input_base_path / '20230710_143554',
}


def main() -> None:
    """Prepare all data for publication."""
    output_path.mkdir(parents=True, exist_ok=True)

    prepare_data_analysis_data(data_analysis_directories, output_path)
    prepare_model_sample_for_data_sets_data(model_sample_for_data_directories, output_path)
    prepare_model_analysis_data(model_sample_directory, output_path)
    prepare_simulation_degree_distribution_data(degree_distribution_directories, output_path)
    prepare_simulation_betti_number_data(betti_number_directories, output_path)
    prepare_simulation_simplex_count_data(simplex_count_directories, output_path)
    prepare_hypothesis_testing_data(hypothesis_testing_directories, output_path)


def prepare_data_analysis_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    create_degree_distributions(directories, output_dir)
    create_dimension_distributions(directories, output_dir)
    create_betti_numbers(directories, output_dir)
    _copy_data_network_plots(directories, output_dir)
    _merge_network_info(directories, 'data', output_dir)


def prepare_model_sample_for_data_sets_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    _copy_model_samples_for_data_sets_network_plots(directories, output_dir)


def prepare_model_analysis_data(directory: Path, output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    _copy_model_sample_network_plots(directory, output_dir)


def prepare_simulation_degree_distribution_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the higher-order degree distributions."""
    _merge_degree_exponents(directories, output_dir)


def prepare_simulation_simplex_count_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the Betti numbers."""
    for property_name in [
        'num_of_edges_normal_mle',
        'num_of_edges_stable',
        'num_of_triangles_normal_mle',
        'num_of_triangles_stable',
    ]:
        _merge_property_tables(directories, 'model_test', property_name, output_dir / 'model_test')


def prepare_simulation_betti_number_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the Betti numbers."""
    for property_name in [
        'betti_number_0_normal',
        'betti_number_0_stable',
        'betti_number_1_normal',
        'betti_number_1_stable',
        'betti_number_2_normal',
        'betti_number_2_stable',
    ]:
        _merge_property_tables(directories, 'model_test', property_name, output_dir / 'model_test')


def prepare_hypothesis_testing_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare hypothesis testing data."""
    (output_dir / 'hypothesis_test').mkdir(parents=True, exist_ok=True)
    _merge_model_info(directories, 'model_test', output_dir / 'hypothesis_test')

    for property_name in [
        'vertex_degree_exponent',
        'edge_degree_exponent',
        'triangle_degree_exponent',
        'num_of_edges_normal_mle',
        'num_of_edges_stable',
        'num_of_triangles_normal_mle',
        'num_of_triangles_stable',
        'betti_number_0_normal',
        'betti_number_0_stable',
        'betti_number_1_normal',
        'betti_number_1_stable',
        'betti_number_2_normal',
        'betti_number_2_stable',
    ]:
        _merge_property_tables(directories, 'model_test', property_name, output_dir / 'hypothesis_test')


def _merge_degree_exponents(directories: dict[str, Path], output_dir: Path):
    for property_name in [
        'vertex_degree_exponent',
        'edge_degree_exponent',
        'triangle_degree_exponent',
    ]:
        _merge_property_tables(directories, 'model_test', property_name, output_dir / 'model_test')


def create_degree_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create higher order degree distributions."""
    for property_name in [
        'total_degree_distribution',
        'ho_degree_distribution_1',
        'ho_degree_distribution_2',
        'ho_degree_distribution_3',
    ]:
        _merge_property_tables(directories, 'data', property_name, output_dir / 'data')


def create_dimension_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create dimension distributions."""
    for property_name in [
        'simplex_dimension_distribution',
        'interaction_dimension_distribution',
        'facet_dimension_distribution',
    ]:
        _merge_property_tables(directories, 'data', property_name, output_dir / 'data')


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


def _merge_network_info(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
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


def _copy_data_network_plots(directories: dict[str, Path], output_path: Path) -> None:
    """Copy network plots to the output path."""
    (output_path / 'data').mkdir(parents=True, exist_ok=True)
    for dataset_name, directory in directories.items():
        file_name = directory / 'data' / 'network.png'
        if not file_name.is_file():
            continue
        copyfile(file_name, output_path / 'data' / f'{dataset_name}_network.png')


def _copy_model_sample_network_plots(directory: Path, output_path: Path) -> None:
    """Copy network plots to the output path."""
    (output_path / 'model_sample').mkdir(parents=True, exist_ok=True)
    file_name = directory / 'model_analysis_finite' / 'network.png'
    if file_name.is_file():
        copyfile(file_name, output_path / 'model_sample' / 'finite_network.png')


def _copy_model_samples_for_data_sets_network_plots(
    directories: dict[str, Path],
    output_path: Path
) -> None:
    (output_path / 'data').mkdir(parents=True, exist_ok=True)
    for dataset_name, directory in directories.items():
        file_name = directory / 'model_analysis_finite' / 'network.png'
        if not file_name.is_file():
            continue
        copyfile(file_name, output_path / 'data' / f'{dataset_name}_model_sample_network.png')


def _merge_property_tables(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
) -> None:

    (output_dir / property_name).mkdir(parents=True, exist_ok=True)

    _merge_distribution_info(directories, input_subdir, property_name, output_dir)
    _merge_histograms(directories, input_subdir, property_name, output_dir)
    # _merge_value_sequences(directories, input_subdir, property_name, output_dir)
    _merge_value_counts(directories, input_subdir, property_name, output_dir)
    _merge_qq_plots(directories, input_subdir, property_name, output_dir)
    _merge_confidence_intervals(directories, input_subdir, property_name, output_dir)
    _merge_quantiles(directories, input_subdir, property_name, output_dir)
    _merge_theoretical_pdfs(directories, input_subdir, property_name, output_dir)
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
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Create and merge network info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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
        merged_data_frame.astype(int).to_csv(out_file_name)


def _merge_qq_plots(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_path: Path,
) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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

    for dataset_name, directory in directories.items():
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
