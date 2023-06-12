#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path
from shutil import copyfile

import pandas as pd


def main() -> None:
    """Prepare all dat for publication."""
    input_base_path = Path('../output')
    output_path = Path('../output/_prepared_data')
    output_path.mkdir(parents=True, exist_ok=True)

    data_analysis_directories = {
        # 'finance':     input_base_path / '20230606_074218',
        # 'biology':     input_base_path / '20230606_074555',
        # 'statistics':  input_base_path / '20230606_075103',
        # 'mathematics': input_base_path / '20230606_075318',
        # 'economics':   input_base_path / '20230606_080916',
        # 'engineering': input_base_path / '20230606_081030',
        'finance':     input_base_path / '20230606_184218',
        'biology':     input_base_path / '20230606_184236',
        'statistics':  input_base_path / '20230606_184248',
        'mathematics': input_base_path / '20230606_185539',
        'economics':   input_base_path / '20230606_184303',
        'engineering': input_base_path / '20230606_185553',
    }

    degree_distribution_directories = {
        '10':       input_base_path / '20230610_092933',
        '100':      input_base_path / '20230610_092951',
        '1000':     input_base_path / '20230610_093131',
        '10000':    input_base_path / '20230610_093202',
        '100000':   input_base_path / '20230610_095550',
        'infinite': input_base_path / '20230610_130000',
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

    prepare_data_analysis_data(data_analysis_directories, output_path)
    prepare_simulation_degree_distribution_data(degree_distribution_directories, output_path)
    prepare_simulation_betti_number_data(betti_number_directories, output_path)


def prepare_data_analysis_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    create_ordinary_degree_distributions(directories, output_dir)
    create_higher_order_degree_distributions(directories, output_dir)
    create_dimension_distributions(directories, output_dir)
    copy_network_plots(directories, output_dir)
    _merge_network_info(directories, 'data', output_dir / 'network_info.csv')


def prepare_simulation_degree_distribution_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the higher-order degree distributions."""
    _merge_ho_degree_exponents(directories, output_dir)


def prepare_simulation_betti_number_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the Betti numbers."""
    for dimension in range(3):
        _merge_histograms(
            directories,
            f'model_test/betti_number_{dimension}_normal',
            output_dir / f'betti_number_{dimension}_histograms.csv',
        )
        _merge_value_sequence(
            directories,
            f'model_test/betti_number_{dimension}_normal',
            output_dir / f'betti_number_{dimension}_value_sequences.csv',
        )
        _merge_distribution_info(
            directories,
            f'model_test/betti_number_{dimension}_normal',
            output_dir / f'betti_number_{dimension}_normal_distribution_info.csv',
        )
        _merge_distribution_info(
            directories,
            f'model_test/betti_number_{dimension}_stable',
            output_dir / f'betti_number_{dimension}_stable_distribution_info.csv',
        )
        _merge_qq_plots(
            directories,
            f'model_test/betti_number_{dimension}_stable',
            output_dir / f'betti_number_{dimension}_stable_qq_plot.csv',
        )
        _merge_qq_plots(
            directories,
            f'model_test/betti_number_{dimension}_normal',
            output_dir / f'betti_number_{dimension}_normal_qq_plot.csv',
        )


def _merge_ho_degree_exponents(directories: dict[str, Path], output_dir: Path):
    for dimension in ['vertex', 'edge', 'triangle']:
        merged_exponent_data_frame_rows: list[dict[str, float]] = []
        for network_size, directory in directories.items():
            subdirectory_name = f'{dimension}_degree_exponent'
            subdirectory_path = directory / 'model_test' / subdirectory_name

            if not subdirectory_path.is_dir():
                continue

            confidence_intervals = pd.read_csv(subdirectory_path / 'confidence_intervals.csv', index_col=0)
            quantiles = pd.read_csv(subdirectory_path / 'quantiles.csv', index_col=0)
            distribution_info = pd.read_csv(subdirectory_path / 'distribution_info.csv')

            merged_data_frame_row = {
                'network_size':         network_size,
                'num_of_values':        distribution_info['empirical_num_of_values'].iloc[0],
                'mean':                 distribution_info['empirical_mean'].iloc[0],
                'stddev':               distribution_info['empirical_std_dev'].iloc[0],
                'confidence_99_lower':  confidence_intervals['empirical_lower'][0.99],
                'confidence_99_upper':  confidence_intervals['empirical_upper'][0.99],
                'confidence_95_lower':  confidence_intervals['empirical_lower'][0.95],
                'confidence_95_upper':  confidence_intervals['empirical_upper'][0.95],
                'confidence_90_lower':  confidence_intervals['empirical_lower'][0.90],
                'confidence_90_upper':  confidence_intervals['empirical_upper'][0.90],
                'quantile_0':           quantiles['empirical'][0],
                'quantile_25':          quantiles['empirical'][25],
                'quantile_50':          quantiles['empirical'][50],
                'quantile_75':          quantiles['empirical'][75],
                'quantile_100':         quantiles['empirical'][100],
            }

            merged_exponent_data_frame_rows.append(merged_data_frame_row)

        merged_data_frame = pd.DataFrame(merged_exponent_data_frame_rows)
        merged_data_frame.to_csv(output_dir / f'{dimension}_degree_exponents.csv', index=False)


def create_ordinary_degree_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create ordinary degree distributions."""
    _merge_value_counts(
        directories,
        'data/total_degree_distribution',
        output_dir / 'total_degree_distribution_value_counts.csv',
    )
    _merge_distribution_info(
        directories,
        'data/total_degree_distribution',
        output_dir / 'total_degree_distribution_info.csv',
    )


def create_higher_order_degree_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create higher order degree distributions."""
    for order in range(1, 4):

        _merge_value_counts(
            directories,
            f'data/ho_degree_distribution_{order}',
            output_dir / f'ho_degree_distribution_{order}_value_counts.csv',
        )
        _merge_distribution_info(
            directories,
            f'data/ho_degree_distribution_{order}',
            output_dir / f'ho_degree_distribution_{order}_info.csv',
        )


def create_dimension_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create dimension distributions."""
    for distribution_type in ['simplex', 'interaction', 'facet']:
        _merge_value_counts(
            directories,
            f'data/{distribution_type}_dimension_distribution',
            output_dir / f'{distribution_type}_dimension_distribution_value_counts.csv',
        )
        _merge_distribution_info(
            directories,
            f'data/{distribution_type}_dimension_distribution',
            output_dir / f'{distribution_type}_dimension_distribution_info.csv',
        )


def copy_network_plots(directories: dict[str, Path], output_path: Path) -> None:
    """Copy network plots to the output path."""
    for dataset_name, directory in directories.items():
        file_name = directory / 'data' / 'network.png'
        if not file_name.is_file():
            continue

        copyfile(directory / 'data' / 'network.png', output_path / f'{dataset_name}_network.png')


def _merge_value_counts(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory_name / 'value_counts.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(
            file_name,
            index_col=0,
            header=None,
            dtype=int
        )
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.index.name = 'value'
        merged_data_frame.fillna(0, inplace=True)
        merged_data_frame.astype(int).to_csv(output_path)


def _merge_histograms(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Merge linearly binned histograms."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory_name / 'histogram_linear.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name)
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_bin_left_limit', f'{dataset_name}_value']
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path, index=False)


def _merge_value_sequence(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory_name / 'value_sequence.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, header=None,)
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path, index=False)


def _merge_qq_plots(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory_name / 'qq_plot_points.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name,)
        if current_distribution.empty:
            continue

        current_distribution.columns = ['_'.join([dataset_name, column_name])
                                        for column_name in current_distribution.columns]
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path, index=False)


def _merge_empirical_pdfs(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory_name / 'pdfs.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, usecols=[0, 1])
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_value', f'{dataset_name}_pdf']
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path, index=False)


def _merge_distribution_info(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Create and merge info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / subdirectory_name).is_dir():
            continue

        current_info = pd.read_csv(directory / subdirectory_name / 'distribution_info.csv')
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.to_csv(output_path)


def _merge_network_info(directories: dict[str, Path], subdirectory_name: str, output_path: Path) -> None:
    """Create and merge info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        current_info = pd.read_csv(directory / subdirectory_name / 'network_info.csv')
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.to_csv(output_path)


if __name__ == '__main__':
    main()
