#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path
from shutil import copyfile, copytree

import pandas as pd

input_base_path = Path('../output')
output_path_root = Path('../output/_prepared_data')

model_sample_directory = input_base_path / '20230724_112745'
model_sample_directory_plot = input_base_path / '20230713_204802'

parameter_a_directories = {
    # the key is the value of "a" * 100
    '50':  input_base_path / '20230805_172644',
    '60':  input_base_path / '20230805_172653',
    '70':  input_base_path / '20230805_172704',
    '80':  input_base_path / '20230805_172714',
    '90':  input_base_path / '20230805_172723',
    '100': input_base_path / '20230805_172733',
    '150': input_base_path / '20230805_172912',
    '200': input_base_path / '20230805_172921',
}

data_analysis_directories = {
    # 'computer_science': input_base_path / '20231208_080033',
    # 'biology':          input_base_path / '20231208_083044',
    # 'economics':        input_base_path / '20231208_083202',
    # 'engineering':      input_base_path / '20231208_083229',
    # 'finance':          input_base_path / '20231208_083259',
    # 'mathematics':      input_base_path / '20231208_075316',
    # 'statistics':       input_base_path / '20231208_074547',
    'computer_science': input_base_path / '20240312_193303',
    'biology':          input_base_path / '20240311_161350',
    'economics':        input_base_path / '20240311_161414',
    'engineering':      input_base_path / '20240312_191910',
    'finance':          input_base_path / '20240311_160315',
    'mathematics':      input_base_path / '20240312_190315',
    'statistics':       input_base_path / '20240311_161506',
}

# model_sample_for_data_directories = {
#     'computer_science': input_base_path / '20230714_173606',
#     'biology':          input_base_path / '20230714_185135',
#     'economics':        input_base_path / '20230714_190142',
#     'engineering':      input_base_path / '20230714_190208',
#     'finance':          input_base_path / '20230714_190228',
#     'mathematics':      input_base_path / '20230714_190308',
#     'statistics':       input_base_path / '20230714_190334',
# }

model_sample_for_data_directories = {
    'computer_science': input_base_path / '20240330_120256',
    'biology':          input_base_path / '20240330_115643',
    'economics':        input_base_path / '20240330_115023',
    'engineering':      input_base_path / '20240330_113433',
    'finance':          input_base_path / '20240330_113104',
    'mathematics':      input_base_path / '20240330_110429',  # '20240316_210724',
    'statistics':       input_base_path / '20240330_094048',
}

degree_distribution_directories = {
    # '10':       input_base_path / '20230614_165042',
    # '100':      input_base_path / '20230614_165118',
    # '1000':     input_base_path / '20230614_165434',
    # '10000':    input_base_path / '20230614_170118',
    # '100000':   input_base_path / '20230614_185624',
    # 'infinite': input_base_path / '20230615_074249',
    # '10':       input_base_path / '20240405_132454',
    '100':      input_base_path / '20240405_132519',
    '1000':     input_base_path / '20240405_132526',
    '10000':    input_base_path / '20240405_132536',
    '100000':   input_base_path / '20240405_132544',
    # '1000000':  input_base_path / '20240323_185549',
}

# simplex_count_directories = {
#     # the key is the value of gamma * 100
#     '25': input_base_path / '20230801_200720',
#     '40': input_base_path / '20230703_101414',
#     '45': input_base_path / '20230703_130440',
#     '50': input_base_path / '20230801_100128',
#     '60': input_base_path / '20230731_124944',
#     '75': input_base_path / '20230615_204859',
# }

simplex_count_directories = {
    # the key is <gamma_prime*100>_<gamma*100>
    '05_25': input_base_path / '20240404_120532',
    '05_50': input_base_path / '20240404_131443',
    '05_75': input_base_path / '20240404_120519',
    '10_25': input_base_path / '20240404_120510',
    '10_50': input_base_path / '20240404_120505',
    '10_75': input_base_path / '20240404_120459',
    '15_25': input_base_path / '20240404_112232',
    '15_50': input_base_path / '20240404_112223',
    '15_75': input_base_path / '20240404_112216',
    '20_25': input_base_path / '20240404_112203',
    '20_50': input_base_path / '20240404_112159',
    '20_75': input_base_path / '20240404_112153',
    '25_25': input_base_path / '20240404_102317',
    '25_50': input_base_path / '20240404_102313',
    '25_75': input_base_path / '20240404_102305',
    '30_25': input_base_path / '20240404_102255',
    '30_50': input_base_path / '20240404_102249',
    '30_75': input_base_path / '20240404_092900',
}

# betti_number_directories = {
#     # the key is the value of gamma * 100
#     # '10': input_base_path / '20230608_191114',
#     # '20': input_base_path / '20230608_191120',
#     '25': input_base_path / '20230609_212258',  # 20230804_115431
#     # '30': input_base_path / '20230608_191124',
#     # '40': input_base_path / '20230608_191128',
#     '50': input_base_path / '20230609_212245',  # 20230803_105044
#     # '60': input_base_path / '20230608_205831',  # 20230803_100554
#     '67': input_base_path / '20230809_173951',
#     # '70': input_base_path / '20230608_205839',
#     # '75': input_base_path / '20230609_142002',  # 20230804_115418
#     # '80': input_base_path / '20230608_205845',
#     # '90': input_base_path / '20230608_205850',
# }

betti_number_directories = {
    # the key is <gamma_prime*100>_<gamma*100>
    '05_25': input_base_path / '20240404_120532',
    '05_50': input_base_path / '20240404_131443',
    '05_75': input_base_path / '20240404_120519',
    '10_25': input_base_path / '20240404_120510',
    '10_50': input_base_path / '20240404_120505',
    '10_75': input_base_path / '20240404_120459',
    '15_25': input_base_path / '20240404_112232',
    '15_50': input_base_path / '20240404_112223',
    '15_75': input_base_path / '20240404_112216',
    '20_25': input_base_path / '20240404_112203',
    '20_50': input_base_path / '20240404_112159',
    '20_75': input_base_path / '20240404_112153',
    '25_25': input_base_path / '20240404_102317',
    '25_50': input_base_path / '20240404_102313',
    '25_75': input_base_path / '20240404_102305',
    '30_25': input_base_path / '20240404_102255',
    '30_50': input_base_path / '20240404_102249',
    '30_75': input_base_path / '20240404_092900',
}

# hypothesis_testing_directories = {
#     'computer_science': input_base_path / '20230807_181159',  # 20230712_075521
#     'biology':          input_base_path / '20230710_164059',  # 20230710_164059
#     'economics':        input_base_path / '20230710_162805',  # 20230710_162805
#     'engineering':      input_base_path / '20230806_210326',  # 20230711_110334
#     'finance':          input_base_path / '20230710_162239',  # 20230710_162239
#     'mathematics':      input_base_path / '20230806_200515',  # 20230710_200618
#     'statistics':       input_base_path / '20230806_205545',  # 20230710_143554
# }


hypothesis_testing_directories = {
    'computer_science': input_base_path / '20240322_200625',
    'engineering':      input_base_path / '20240322_092931',
    'mathematics':      input_base_path / '20240321_221229',
    'statistics':       input_base_path / '20240322_085600',
    # 'finance':          input_base_path / '20240322_090205',
}


def main() -> None:
    """Prepare all data for publication."""
    output_path_root.mkdir(parents=True, exist_ok=True)

    prepare_data_analysis_data(data_analysis_directories, output_path_root)
    prepare_model_sample_for_data_sets_data(model_sample_for_data_directories, output_path_root)
    # prepare_model_analysis_sample_plot(model_sample_directory_plot, output_path_root)
    # prepare_model_analysis_data(model_sample_directory, output_path_root)
    prepare_simulation_degree_distribution_data(degree_distribution_directories, output_path_root)
    prepare_simulation_betti_number_data(betti_number_directories, output_path_root)
    prepare_simulation_simplex_count_data(simplex_count_directories, output_path_root)
    prepare_hypothesis_testing_data(hypothesis_testing_directories, output_path_root)


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
    (output_dir / 'model_sample').mkdir(parents=True, exist_ok=True)

    directory_name = directory / 'model_analysis_finite' / 'total_degree_distribution'
    if directory_name.is_dir():
        copytree(directory_name, output_dir / 'model_sample' / 'total_degree_distribution', dirs_exist_ok=True)

    directory_name = directory / 'model_analysis_finite' / 'ho_degree_distribution_1'
    if directory_name.is_dir():
        copytree(directory_name, output_dir / 'model_sample' / 'ho_degree_distribution_1', dirs_exist_ok=True)

    directory_name = directory / 'model_analysis_finite' / 'ho_degree_distribution_2'
    if directory_name.is_dir():
        copytree(directory_name, output_dir / 'model_sample' / 'ho_degree_distribution_2', dirs_exist_ok=True)


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
        _merge_property_tables(directories, 'model_test', property_name, output_dir / 'model_test')
        _create_histograms_normal_plots(directories, 'model_test', property_name, output_dir / 'model_test')


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
        _create_histograms_normal_plots(directories, 'model_test', property_name, output_dir / 'model_test')


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
        'interaction_vertex_degree_exponent',
        'vertex_interaction_degree_exponent',
    ]:
        _merge_property_tables(directories, 'model_test', property_name, output_dir / 'model_test')


def _create_boxplots_degree_distributions(directories: dict[str, Path], output_dir: Path):

    GAMMA = 0.7
    GAMMA_PRIME = 0.2

    for property_name in [
        'vertex_degree_exponent',
        'edge_degree_exponent',
        # 'triangle_degree_exponent',
        'interaction_vertex_degree_exponent',
        'vertex_interaction_degree_exponent',
    ]:
        caption = ''
        theoretical_value = 0.
        if property_name == 'vertex_degree_exponent':
            caption = 'Vertex--edge degree exponent distribution'
            theoretical_value = 1. + 1. / GAMMA
        elif property_name == 'edge_degree_exponent':
            caption = 'Edge--triangle degree exponent distribution'
            theoretical_value = 2. / GAMMA
        elif property_name == 'triangle_degree_exponent':
            caption = 'Triangle--tetrahedron degree exponent distribution'
            theoretical_value = 3. / GAMMA - 1.
        elif property_name == 'interaction_vertex_degree_exponent':
            caption = 'Interaction--vertex degree exponent distribution'
            theoretical_value = 1. + 1. / GAMMA_PRIME
        elif property_name == 'vertex_interaction_degree_exponent':
            caption = 'Vertex--interaction degree exponent distribution'
            theoretical_value = 1. + 1. / GAMMA

        _create_boxplots(directories, 'model_test', property_name, output_dir / 'model_test', caption, theoretical_value)


def _create_boxplots(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
    caption: str,
    theoretical_value: float
) -> None:
    """Create boxplots."""
    # read quantiles
    quantiles = pd.concat([
        pd.read_csv(directory / input_subdir / property_name / 'quantiles.csv', index_col=0)['empirical']
        for directory in directories.values()
    ], axis=1)
    quantiles.columns = directories.keys()
    quantiles.index.name = 'quantile'

    out_file_name = output_dir / property_name / 'boxplot_latex_figure.tex'
    out_file_name.parent.mkdir(parents=True, exist_ok=True)

    y_range = quantiles.to_numpy().max() - quantiles.to_numpy().min()
    y_min = min(quantiles.to_numpy().min(), theoretical_value) - 0.05 * y_range
    y_max = min(max(quantiles.to_numpy().max(), theoretical_value) + 0.05 * y_range, 2 * theoretical_value)

    with open(out_file_name, 'w') as out_file:
        print(r'\begin{subfigure}[t]{0.22\textwidth} \resizebox{\textwidth}{!}{', file=out_file)
        print(r'    \begin{tikzpicture}', file=out_file)
        print(r'        \begin{axis}[', file=out_file)
        print(f'                xmin=0, xmax={len(directories) + 1}, ymin={y_min:.4f}, ymax={y_max:.4f},', file=out_file)
        print(r'                xlabel=network size,', file=out_file)
        print(r'                axis x line=middle, axis y line=middle,', file=out_file)
        print(r'                xtick={', end='', file=out_file)
        for i in range(len(directories) - 1):
            print(f'{i + 1}, ', end='', file=out_file)
        print(f'{len(directories)}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                ytick={', end='', file=out_file)
        for i in range(int(y_min) + 1, int(y_max)):
            print(f'{i}, ', end='', file=out_file)
        print(f'{int(y_max)}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                xticklabels={', end='', file=out_file)
        for dataset_name in list(directories.keys())[:-1]:
            print(f'{dataset_name}, ', end='', file=out_file)
        print(f'{list(directories.keys())[-1]}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'        ]', file=out_file)
        for i, dataset_name in enumerate(directories.keys()):
            print(r'            \addplot+[boxplot, boxplot/draw direction=y, black, solid, boxplot prepared={', end='', file=out_file)
            print(f'draw position={i + 1}, ', end='', file=out_file)
            print(f'lower whisker={quantiles[dataset_name][0]}, ', end='', file=out_file)
            print(f'lower quartile={quantiles[dataset_name][25]}, ', end='', file=out_file)
            print(f'median={quantiles[dataset_name][50]}, ', end='', file=out_file)
            print(f'upper quartile={quantiles[dataset_name][75]}, ', end='', file=out_file)
            print(f'upper whisker={quantiles[dataset_name][100]}, ', end='', file=out_file)
            print(f'sample size={100}', end='', file=out_file)
            print(r'}] coordinates {};', file=out_file)
        print(r'            \addplot', end='', file=out_file)
        print(f'[red, thick, domain=0:{len(directories) + 1}] ', end='', file=out_file)
        print(r'{', end='', file=out_file)
        print(f'{theoretical_value:.4f}', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'      \end{axis}', file=out_file)
        print(r'  \end{tikzpicture}} \caption{', end='', file=out_file)
        print(caption, end='', file=out_file)
        print(r'}', file=out_file)
        print(r'\end{subfigure}', file=out_file)


def create_degree_distributions(directories: dict[str, Path], output_dir: Path) -> None:
    """Create higher order degree distributions."""
    for property_name in [
        'total_degree_distribution',
        'interaction_degree_distribution',
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


def prepare_model_analysis_sample_plot(directory: Path, output_path: Path) -> None:
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


def _create_value_counts_log_plots(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
) -> None:
    for dataset_name, directory in directories.items():
        histogram_file_name = directory / input_subdir / property_name / 'histogram_linear.csv'
        if not histogram_file_name.is_file():
            continue

        current_distribution = pd.read_csv(histogram_file_name)
        if current_distribution['bin_left_limit'].isna().iloc[0]:
            continue

        info = pd.read_csv(directory / input_subdir / property_name / 'distribution_info.csv')
        mean = info['empirical_mean'].iloc[0]
        std = info['empirical_std_dev'].iloc[0]

        histogram_x_values = current_distribution['bin_left_limit']
        histogram_y_values = current_distribution['value']
        histogram_bin_size = histogram_x_values[1] - histogram_x_values[0]

        histogram_x_range = [
            histogram_x_values.min(),
            histogram_x_values.max() + histogram_bin_size,
        ]

        histogram_x_limits = [
            histogram_x_range[0] - (histogram_x_range[1] - histogram_x_range[0]) * 0.05,
            histogram_x_range[1] + (histogram_x_range[1] - histogram_x_range[0]) * 0.05,
        ]
        histogram_y_limits = [0, histogram_y_values.max() * 1.1,]

        out_file_name = output_dir / property_name / 'distribution_latex_figures' / f'figure_{dataset_name}.tex'
        out_file_name.parent.mkdir(parents=True, exist_ok=True)
    '''
  \begin{subfigure}[b]{0.24\textwidth} \resizebox{\textwidth}{!}{
    \begin{tikzpicture}
      \begin{axis} [
          xmin=1, xmax=220, ymin=10^(-0.2), ymax=10^(4.5),
          clip mode=individual,
          smooth, xmode=log, ymode=log,
          xlabel=vertex degree,
          ylabel style={
            yshift=-0.6cm
          },
          legend style={at={(0,0)}, anchor=south west, at={(axis description cs:0.025,0.025)}}, legend cell align={left},
        ]
        \addplot[only marks, mark options={scale=0.4}] table [x=value, y=statistics, col sep=comma] {data_total_degree_distribution_value_counts.csv};
        \addplot[red, ultra thick, no marks, domain={10:3000}] {10^(5.85)*x^(-2.9635)};
      \end{axis}
    \end{tikzpicture}}
  \end{subfigure}
  '''
    pass


def _create_histograms_normal_plots(
    directories: dict[str, Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
) -> None:
    for dataset_name, directory in directories.items():
        histogram_file_name = directory / input_subdir / property_name / 'histogram_linear.csv'
        if not histogram_file_name.is_file():
            continue

        current_distribution = pd.read_csv(histogram_file_name)
        if current_distribution['bin_left_limit'].isna().iloc[0]:
            continue

        info = pd.read_csv(directory / input_subdir / property_name / 'distribution_info.csv')
        mean = info['empirical_mean'].iloc[0]
        std = info['empirical_std_dev'].iloc[0]

        histogram_x_values = current_distribution['bin_left_limit']
        histogram_y_values = current_distribution['value']
        histogram_bin_size = histogram_x_values[1] - histogram_x_values[0]

        histogram_x_range = [
            histogram_x_values.min(),
            histogram_x_values.max() + histogram_bin_size,
        ]

        histogram_x_limits = [
            histogram_x_range[0] - (histogram_x_range[1] - histogram_x_range[0]) * 0.05,
            histogram_x_range[1] + (histogram_x_range[1] - histogram_x_range[0]) * 0.05,
        ]
        histogram_y_limits = [0, histogram_y_values.max() * 1.1,]

        out_file_name = output_dir / property_name / 'distribution_latex_figures' / f'figure_{dataset_name}.tex'
        out_file_name.parent.mkdir(parents=True, exist_ok=True)
        with open(out_file_name, 'w') as out_file:
            print(r'\begin{subfigure}[b]{0.28\textwidth}', file=out_file)
            print(r'    \resizebox{0.48\textwidth}{!}{', file=out_file)
            print(r'        \begin{tikzpicture}', file=out_file)
            print(r'            \begin{axis} [', file=out_file)
            print(r'                    grid=none, ticks=none, smooth,', file=out_file)
            print(f'                    xmin={histogram_x_limits[0]:.4f}, xmax={histogram_x_limits[1]:.4f}, ymin={histogram_y_limits[0]:.10f}, ymax={histogram_y_limits[1]:.10f},', file=out_file)
            print(r'                ]', file=out_file)
            print(r'                \addplot[ybar interval, mark=no, fill=gray!50!white, ', file=out_file, end='')
            print(f'restrict x to domain={histogram_x_limits[0]:.4f}:{histogram_x_limits[1]:.4f}]', file=out_file)
            print(f'                table[x={dataset_name}_bin_left_limit, y={dataset_name}_value, col sep=comma]', file=out_file)
            print(r'                    {', file=out_file, end='')
            print(f'data/{input_subdir}/{property_name}/linear_histograms.csv', file=out_file, end='')
            print(r'};', file=out_file)
            print(rf'               \addplot[ultra thick, red, samples=100, domain={histogram_x_limits[0]:.4f}:{histogram_x_limits[1]:.4f}]', file=out_file, end='')
            print(r' {', file=out_file, end='')
            print(f'normalpdf(x, {mean:.4f}, {std:.4f})', file=out_file, end='')
            print(r'};', file=out_file)
            print(r'            \end{axis}', file=out_file)
            print(r'        \end{tikzpicture}', file=out_file)
            print(r'    }', file=out_file)
            print(r'    \hfill', file=out_file)
            print(r'    \resizebox{0.48\textwidth}{!}{', file=out_file)
            print(r'        \begin{tikzpicture}', file=out_file)
            print(r'            \begin{axis} [', file=out_file)
            print(r'                    grid=none, ticks=none, smooth,', file=out_file)
            print(r'                    xmin=-3, xmax=3, ymin=-3, ymax=3,', file=out_file)
            print(r'                    clip mode=individual,', file=out_file)
            print(r'                    scaled y ticks=false, % remove axis multiplier', file=out_file)
            print(r'                    yticklabel style={text width=2em,align=right},', file=out_file)
            print(r'                    y tick label style={/pgf/number format/fixed} % suppress scientific notation', file=out_file)
            print(r'                ]', file=out_file)
            print(r'                \addplot[black, only marks, mark=*, mark size=1.5]', file=out_file)
            print(f'                table [x={dataset_name}_theoretical, y={dataset_name}_empirical, col sep=comma]', file=out_file)
            print(r'                    {', file=out_file, end='')
            print(f'data/{input_subdir}/{property_name}/qq_plot.csv', file=out_file, end='')
            print(r'};', file=out_file)
            print(r'                \addplot[ultra thick, red, samples=100, domain=-3:3]{x};', file=out_file)
            print(r'            \end{axis}', file=out_file)
            print(r'        \end{tikzpicture}', file=out_file)
            print(r'    }', file=out_file)
            print(r'\end{subfigure}', file=out_file)


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
