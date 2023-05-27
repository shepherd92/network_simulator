#!/usr/bin/env python

from pathlib import Path
from shutil import copyfile

import pandas as pd


def main() -> None:
    """Prepare all dat for publication."""
    input_base_path = Path('../output')
    output_path = Path('../output/_prepared_data')
    output_path.mkdir(parents=True, exist_ok=True)

    directories = {
        'finance':     input_base_path / '20230525_104118',
        'biology':     input_base_path / '20230525_103656',
        'statistics':  input_base_path / '20230525_101139',
        'mathematics': input_base_path / '20230525_101354',
        'economics':   input_base_path / '20230525_104241',
        'engineering': input_base_path / '20230525_104359',
    }

    prepare_data_analysis_data(directories, output_path)


def prepare_data_analysis_data(directories: dict[str, Path], output_path: Path) -> None:
    """Prepare all information related to the data analysis."""
    create_ordinary_degree_distributions(directories, output_path)
    create_higher_order_degree_distributions(directories, output_path)
    create_dimension_distributions(directories, output_path)
    copy_network_plots(directories, output_path)
    create_network_info_table(directories, output_path)


def create_ordinary_degree_distributions(directories: dict[str, Path], output_path: Path) -> None:
    """Create ordinary degree distributions."""
    _merge_distribution_integer_histograms(
        directories,
        'total_degree_distribution.csv',
        output_path
    )
    _merge_data_frames(
        directories,
        'total_degree_distribution_info.csv',
        output_path
    )


def create_higher_order_degree_distributions(directories: dict[str, Path], output_path: Path) -> None:
    """Create higher order degree distributions."""
    for order in range(1, 4):
        _merge_distribution_integer_histograms(
            directories,
            f'ho_degree_distribution_{order}.csv',
            output_path
        )
        _merge_data_frames(
            directories,
            f'ho_degree_distribution_{order}_info.csv',
            output_path
        )


def create_dimension_distributions(directories: dict[str, Path], output_path: Path) -> None:
    """Create dimension distributions."""
    for distribution_type in ['simplex', 'interaction', 'facet']:
        _merge_distribution_integer_histograms(
            directories,
            f'{distribution_type}_dimension_distribution.csv',
            output_path
        )
        _merge_data_frames(
            directories,
            f'{distribution_type}_dimension_distribution_info.csv',
            output_path
        )


def create_network_info_table(directories: dict[str, Path], output_path: Path) -> None:
    """Create and merge info."""
    _merge_data_frames(directories, 'network_info.csv', output_path)


def copy_network_plots(directories: dict[str, Path], output_path: Path) -> None:
    """Copy network plots to the output path."""
    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        copyfile(directory / 'data' / 'network.png', output_path / f'{dataset_name}_network.png')


def _merge_distribution_integer_histograms(directories: dict[str, Path], file_name: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        current_ho_degree_distribution = pd.read_csv(
            directory / 'data' / file_name,
        )
        if current_ho_degree_distribution.empty:
            continue

        current_ho_degree_distribution.index = range(len(current_ho_degree_distribution))
        current_ho_degree_distribution.index.name = 'value'
        current_ho_degree_distribution.drop('bin_left_limit', axis=1, inplace=True)
        current_ho_degree_distribution.columns = [f'{dataset_name}']
        data_frames.append(current_ho_degree_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path / file_name)


def _merge_data_frames(directories: dict[str, Path], file_name: str, output_path: Path) -> None:
    """Create and merge info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        current_info = pd.read_csv(directory / 'data' / file_name)
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.to_csv(output_path / file_name)


if __name__ == '__main__':
    main()
