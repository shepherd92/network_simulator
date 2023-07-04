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
        'finance':     input_base_path / '20230606_184218',  # average degree: 3.2621
        'biology':     input_base_path / '20230606_184236',  # average degree: 6.2481
        'statistics':  input_base_path / '20230606_184248',  # average degree: 5.0948
        'mathematics': input_base_path / '20230606_185539',  # average degree: 4.5596
        'economics':   input_base_path / '20230606_184303',  # average degree: 3.0593
        'engineering': input_base_path / '20230606_185553',  # average degree: 7.0446
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

    prepare_data_analysis_data(data_analysis_directories, output_path)
    prepare_simulation_degree_distribution_data(degree_distribution_directories, output_path)
    prepare_simulation_betti_number_data(betti_number_directories, output_path)
    prepare_simulation_simplex_count_data(simplex_count_directories, output_path)


def prepare_data_analysis_data(directories: dict[str, Path], output_dir: Path) -> None:
    """Prepare all information related to the data analysis."""
    create_degree_distributions(directories, output_dir)
    create_dimension_distributions(directories, output_dir)
    copy_network_plots(directories, output_dir)
    _merge_network_info(directories, 'data', output_dir / 'network_info.csv')


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
        _merge_property_tables(directories, f'model_test/{property_name}', output_dir)


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
        _merge_property_tables(directories, f'model_test/{property_name}', output_dir)


def _merge_degree_exponents(directories: dict[str, Path], output_dir: Path):
    for property_name in [
        'vertex_degree_exponent',
        'edge_degree_exponent',
        'triangle_degree_exponent',
    ]:
        _merge_property_tables(directories, f'model_test/{property_name}', output_dir)


def create_degree_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create higher order degree distributions."""
    for property_name in [
        'total_degree_distribution',
        'ho_degree_distribution_1',
        'ho_degree_distribution_2',
        'ho_degree_distribution_3',
    ]:
        _merge_property_tables(directories, f'data/{property_name}', output_dir)


def create_dimension_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create dimension distributions."""
    for property_name in [
        'simplex_dimension_distribution',
        'interaction_dimension_distribution',
        'facet_dimension_distribution',
    ]:
        _merge_property_tables(directories, f'data/{property_name}', output_dir)


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


def copy_network_plots(directories: dict[str, Path], output_path: Path) -> None:
    """Copy network plots to the output path."""
    for dataset_name, directory in directories.items():
        file_name = directory / 'data' / 'network.png'
        if not file_name.is_file():
            continue
        copyfile(directory / 'data' / 'network.png', output_path / f'{dataset_name}_network.png')


def _merge_property_tables(directories: dict[str, Path], subdirectory: str, output_dir: Path) -> None:

    (output_dir / subdirectory).mkdir(parents=True, exist_ok=True)

    _merge_distribution_info(directories, subdirectory, output_dir)
    _merge_histograms(directories, subdirectory, output_dir)
    _merge_value_sequence(directories, subdirectory, output_dir)
    _merge_value_counts(directories, subdirectory, output_dir)
    _merge_qq_plots(directories, subdirectory, output_dir)
    _merge_confidence_intervals(directories, subdirectory, output_dir)
    _merge_quantiles(directories, subdirectory, output_dir)
    _merge_theoretical_pdfs(directories, subdirectory, output_dir)


def _merge_value_counts(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'value_counts.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, index_col=0, header=None, dtype=int)
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.index.name = 'value'
        merged_data_frame.fillna(0, inplace=True)
        merged_data_frame.astype(int).to_csv(output_path / subdirectory / 'value_counts.csv')


def _merge_histograms(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Merge linearly binned histograms."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'histogram_linear.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name)
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_bin_left_limit', f'{dataset_name}_value']
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path / subdirectory / 'linear_histograms.csv', index=False)


def _merge_value_sequence(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'value_sequence.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, header=None,)
        if current_distribution.empty:
            continue

        current_distribution.columns = [dataset_name]
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path / subdirectory / 'value_sequences.csv', index=False)


def _merge_qq_plots(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'qq_plot_points.csv'
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
        merged_data_frame.to_csv(output_path / subdirectory / 'qq_plot.csv', index=False)


def _merge_confidence_intervals(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'confidence_intervals.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, index_col=0)
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
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path / subdirectory / 'confidence_intervals.csv')


def _merge_quantiles(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'quantiles.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, index_col=0)
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_empirical', f'{dataset_name}_theoretical']
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path / subdirectory / 'quantiles.csv')


def _merge_theoretical_pdfs(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Crete higher order degree distributions."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        file_name = directory / subdirectory / 'pdfs.csv'
        if not file_name.is_file():
            continue

        current_distribution = pd.read_csv(file_name, usecols=[0, 2])
        if current_distribution.empty:
            continue

        current_distribution.columns = [f'{dataset_name}_value', f'{dataset_name}_pdf']
        data_frames.append(current_distribution)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=1)
        merged_data_frame.to_csv(output_path / subdirectory / 'theoretical_pdf.csv', index=False)


def _merge_distribution_info(directories: dict[str, Path], subdirectory: str, output_path: Path) -> None:
    """Create and merge info."""
    data_frames: list[pd.DataFrame] = []

    for dataset_name, directory in directories.items():
        if not (directory / subdirectory).is_dir():
            continue

        current_info = pd.read_csv(directory / subdirectory / 'distribution_info.csv')
        current_info.index = [dataset_name]
        data_frames.append(current_info)

    if data_frames:
        merged_data_frame = pd.concat(data_frames, axis=0)
        merged_data_frame.to_csv(output_path / subdirectory / 'distribution_info.csv')


if __name__ == '__main__':
    main()
