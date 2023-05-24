#!/usr/bin/env python

from json import load
from pathlib import Path
from shutil import copyfile

import pandas as pd


def main() -> None:
    """Prepare all dat for publication."""
    input_base_path = Path('../output')
    output_path = Path('../output/_prepared_data')
    output_path.mkdir(parents=True, exist_ok=True)

    directories = {
        'finance':     input_base_path / '20230524_150453',
        'biology':     input_base_path / '20230524_150207',
        'statistics':  input_base_path / '20230524_150716',
        'mathematics': input_base_path / '',
    }

    prepare_data_analysis_data(directories, output_path)


def prepare_data_analysis_data(directories: dict[str, Path], output_path: Path) -> None:
    """Prepare all information related to the data analysis."""

    create_higher_order_degree_distributions(directories, output_path)
    copy_network_plots(directories, output_path)
    create_info_table(directories, output_path)


def create_higher_order_degree_distributions(directories: dict[str, Path], output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        for order in range(1, 4):
            current_ho_degree_distribution = pd.read_csv(
                directory / 'data' / f'ho_degree_distribution_{order}.csv',
            )
            if current_ho_degree_distribution.empty:
                continue

            current_ho_degree_distribution.index = range(len(current_ho_degree_distribution))
            current_ho_degree_distribution.index.name = 'degree'
            current_ho_degree_distribution.drop('bin_left_limit', axis=1, inplace=True)
            current_ho_degree_distribution.columns = [f'{dataset_name}_{order}']
            data_frames.append(current_ho_degree_distribution)

    merged_data_frame = pd.concat(data_frames, axis=1)
    merged_data_frame.to_csv(output_path / 'ho_degree_distributions.csv')


def create_info_table(directories: dict[str, Path], output_path: Path) -> None:
    """Create and merge info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        with open(directory / 'data' / 'info.json') as file:
            current_info = pd.DataFrame(load(file))
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    merged_data_frame = pd.concat(data_frames, axis=0)
    merged_data_frame.to_csv(output_path / 'info.csv')


def copy_network_plots(directories: dict[str, Path], output_path: Path) -> None:
    """Copy network plots to the output path."""
    for dataset_name, directory in directories.items():
        if not (directory / 'data').is_dir():
            continue

        copyfile(directory / 'data' / 'network.png', output_path / f'{dataset_name}_network.png')


if __name__ == '__main__':
    main()
