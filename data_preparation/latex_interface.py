#!/usr/bin/env python3
"""Create latex plots."""

from pathlib import Path

import pandas as pd
import numpy as np


def create_boxplots(
    directories: dict[tuple[int, int], Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
    caption: str,
    theoretical_value: float,
    suffix: str = '',
) -> None:
    """Create boxplots."""
    # read quantiles
    quantiles = pd.concat([
        pd.read_csv(directory / input_subdir / property_name / 'quantiles.csv', index_col=0)['empirical']
        for directory in directories.values()
    ], axis=1)

    quantiles.columns = [
        f'{gamma_prime}_{gamma}'
        for gamma_prime, gamma in directories.keys()
    ]
    quantiles.index.name = 'quantile'

    out_file_name = output_dir / property_name / f'boxplot_{suffix}_latex_figure.tex'
    out_file_name.parent.mkdir(parents=True, exist_ok=True)

    y_range = quantiles.to_numpy().max() - quantiles.to_numpy().min()
    y_min = min(quantiles.to_numpy().min(), theoretical_value) - 0.05 * y_range
    y_max = min(max(quantiles.to_numpy().max(), theoretical_value) + 0.05 * y_range, 2 * theoretical_value)

    with open(out_file_name, 'w') as out_file:
        print(r'\begin{subfigure}[t]{0.2\textwidth} \resizebox{\textwidth}{!}{', file=out_file)
        print(r'    \begin{tikzpicture}', file=out_file)
        print(r'        \begin{axis}[', file=out_file)
        print(f'                xmin=0, xmax={len(directories) + 1}, ymin={y_min}, ymax={y_max},', file=out_file)
        # print(r'                xlabel=network size,', file=out_file)
        # print(r'                axis x line=middle, axis y line=middle,', file=out_file)
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
        for network_size, _ in list(directories.keys())[:-1]:
            print(f'{network_size}, ', end='', file=out_file)
        print(f'{list(directories.keys())[-1][0]}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'        ]', file=out_file)
        for i, key in enumerate(directories.keys()):
            column_name = f'{key[0]}_{key[1]}'
            print(r'            \addplot+[boxplot, boxplot/draw direction=y, black, solid, boxplot prepared={', end='', file=out_file)
            print(f'draw position={i + 1}, ', end='', file=out_file)
            print(f'lower whisker={quantiles[column_name][0]}, ', end='', file=out_file)
            print(f'lower quartile={quantiles[column_name][25]}, ', end='', file=out_file)
            print(f'median={quantiles[column_name][50]}, ', end='', file=out_file)
            print(f'upper quartile={quantiles[column_name][75]}, ', end='', file=out_file)
            print(f'upper whisker={quantiles[column_name][100]}, ', end='', file=out_file)
            print(f'sample size={100}', end='', file=out_file)
            print(r'}] coordinates {};', file=out_file)
        print(r'            \addplot', end='', file=out_file)
        print(f'[red, thick, domain=0:{len(directories) + 1}] ', end='', file=out_file)
        print(r'{', end='', file=out_file)
        print(f'{theoretical_value:.4f}', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'      \end{axis}', file=out_file)
        print(r'  \end{tikzpicture}}', file=out_file)
        # print(r'\caption{', end='', file=out_file)
        # print(caption, end='', file=out_file)
        # print(r'}', file=out_file)
        print(r'\end{subfigure}', file=out_file)


def create_histograms_normal_plots(
    directories: dict[str | tuple[int, ...], Path],
    input_subdir: str,
    property_name: str,
    output_dir: Path,
) -> None:
    for key, directory in directories.items():
        dataset_name = key if isinstance(key, str) else '_'.join(str(i) for i in key)
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


def create_value_counts_log_plot(
    in_dir: Path,
    out_file_name: Path,
    csv_file_name: Path,
    x_column: str | int,
    y_column: str | int,
    caption: str,
) -> None:
    """Create value counts log plots."""
    value_counts_file_name = in_dir / 'value_counts.csv'
    if not value_counts_file_name.is_file():
        return

    info = pd.read_csv(
        in_dir / 'distribution_info.csv',
        converters={'theoretical_domain': pd.eval}
    )

    value_counts = pd.read_csv(value_counts_file_name, dtype=int, header=None)
    value_counts.columns = ['value', 'count']
    x_min = max(value_counts['value'].min(), 1.)**1.05 / value_counts['value'].max()**0.05
    x_max = value_counts['value'].max()**1.05 / max(value_counts['value'].min(), 1.)**0.05
    y_min = max(value_counts['count'].min(), 1.)**1.05 / value_counts['count'].max()**0.05
    y_max = value_counts['count'].max()**1.05 / max(value_counts['count'].min(), 1.)**0.05

    pdf = pd.read_csv(in_dir / 'pdfs.csv', dtype=float)
    pdf.columns = ['x', 'empirical', 'theoretical']
    exponent = info['theoretical_exponent'][0]

    theoretical_x_min = info['theoretical_domain'][0][0] \
        if 'theoretical_domain' in info \
        else info['theoretical_domain_min'][0]

    theoretical_integral = np.trapz(pdf['theoretical'].to_numpy(), pdf['x'].to_numpy())
    theoretical_integral = 1. / (1. - exponent) * (x_max**(1. - exponent) - theoretical_x_min**(1. - exponent))
    empirical_integral = value_counts[value_counts['value'] >= theoretical_x_min]['count'].sum()
    coefficient = empirical_integral / theoretical_integral

    out_file_name.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file_name, 'w') as out_file:
        print(r'\begin{subfigure}[b]{0.2\textwidth}', file=out_file)
        print(r'    \resizebox{\textwidth}{!}{', file=out_file)
        print(r'        \begin{tikzpicture}', file=out_file)
        print(r'            \begin{axis} [', file=out_file)
        print(f'                    xmin={x_min:.4f}, xmax={x_max:.4f}, ymin={y_min:.4f}, ymax={y_max},', file=out_file)
        print(r'                    clip mode=individual, smooth, xmode=log, ymode=log,', file=out_file)
        print(r'                    legend style={', file=out_file)
        print(r'                            at={(axis description cs:0.025,0.025)},', file=out_file)
        print(r'                            anchor=south west,', file=out_file)
        print(r'                            legend cell align={left}', file=out_file)
        print(r'                        }', file=out_file)
        print(r'                ]', file=out_file)
        print(r'                \addplot[only marks, mark options={scale=0.4}]', file=out_file)
        x_column_name = f'x index={x_column}' if isinstance(x_column, int) else f'x={x_column}'
        y_column_name = f'y index={y_column}' if isinstance(y_column, int) else f'y={y_column}'
        print(f'                table [{x_column_name}, {y_column_name}, col sep=comma]', file=out_file)
        print(r'                    {', end='', file=out_file)
        print(f'{csv_file_name}', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'                \addplot[red, ultra thick, no marks, domain={', end='', file=out_file)
        print(f'{theoretical_x_min:.4f}:{x_max:.4f}', end='', file=out_file)
        print(r'}] {', end='', file=out_file)
        print(f'{coefficient:.4f}*x^(-{exponent:.4f})', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'            \end{axis}', file=out_file)
        print(r'        \end{tikzpicture}', file=out_file)
        print(r'    }', file=out_file)
        print(r'    \caption{', end='', file=out_file)
        print(caption, end='', file=out_file)
        print(r'}', file=out_file)
        print(r'\end{subfigure}', file=out_file)


def create_dataset_properties_table(
    info_table: pd.DataFrame,
    out_file_name: Path,
) -> None:
    """Create info table."""
    out_file_name.parent.mkdir(parents=True, exist_ok=True)
    dataset_name_translation = {
        'computer_science': 'cs',
        'engineering': 'eess',
        'mathematics': 'math',
        'statistics': 'stat',
    }
    with open(out_file_name, 'w') as out_file:
        print(r'\begin{table} [h] \centering \caption{Main properties of the datasets} \label{tab:dataset_properties}', file=out_file)
        print(r'    \begin{tabular}{|l|r|r|r|r|} \hline', file=out_file)
        print(r'        dataset &  authors & documents & components & size of largest component \\', file=out_file)
        print(r'        \hline', file=out_file)
        for name_in_table, name_to_print in dataset_name_translation.items():
            info_row = info_table.loc[name_in_table]
            print(f'        {name_to_print: <7} & ', end='', file=out_file)
            print(f'${info_row['num_of_vertices']:,}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info_row['num_of_interactions']:,}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info_row['num_of_components']:,}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info_row['num_of_vertices_in_component_0']:,}$ '.replace(',', r'\,'), end='', file=out_file)
            print(r'\\', file=out_file)
        print(r'        \hline', file=out_file)
        print(r'    \end{tabular}', file=out_file)
        print(r'\end{table}', file=out_file)


def create_dataset_parameter_estimates_table(
    dataset_directories: dict[str, Path],
    hypothesis_testing_directories: dict[str, Path],
    out_file_name: Path,
) -> None:
    """Create parameter estimates table."""
    out_file_name.parent.mkdir(parents=True, exist_ok=True)
    dataset_name_translation = {
        'computer_science': 'cs',
        'engineering': 'eess',
        'mathematics': 'math',
        'statistics': 'stat',
    }
    with open(out_file_name, 'w') as out_file:
        print(r"\begin{table} [h] \centering \caption{Fitted power-law exponents and the inferred model parameters} \label{tab:dataset_parameter_estimates}", file=out_file)
        print(r'    \begin{tabular}{|l|r|r|r|r|r|r|r|} \hline', file=out_file)
        print(r"        \multirow{2}{*}{dataset} & \multicolumn{2}{|c|}{vertex--interaction degree} & \multicolumn{2}{|c|}{interaction--vertex degree} & \multirow{2}{*}{$\b$} & \multirow{2}{*}{$\la$} & \multirow{2}{*}{$\la'$} \\", file=out_file)
        print(r"            & exponent & $\g$ & exponent & $\g'$ & & & \\ \hline", file=out_file)
        for name_in_table, name_to_print in dataset_name_translation.items():
            vertex_interaction_exponent = \
                pd.read_csv(dataset_directories[name_in_table] / 'data' / 'interaction_degree_distribution' / 'distribution_info.csv')['theoretical_exponent'][0]
            interaction_vertex_exponent = \
                pd.read_csv(dataset_directories[name_in_table] / 'data' / 'interaction_dimension_distribution' / 'distribution_info.csv')['theoretical_exponent'][0]

            model_info = pd.read_csv(hypothesis_testing_directories[name_in_table] / 'model_test' / 'model_info.csv')
            lambda_ = model_info['network_size'][0]
            lambda_prime = model_info['interaction_intensity'][0]
            beta = model_info['parameter_beta'][0]

            print(f'        {name_to_print: <7} & ', end='', file=out_file)
            print(f'${vertex_interaction_exponent:.2f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${1. / (vertex_interaction_exponent - 1.):.2f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${interaction_vertex_exponent:.2f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${1. / (interaction_vertex_exponent - 1.):.2f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${beta:.2E}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${lambda_:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${lambda_prime:,.0f}$ '.replace(',', r'\,'), end='', file=out_file)
            print(r'\\', file=out_file)
        print(r'        \hline', file=out_file)
        print(r'    \end{tabular}', file=out_file)
        print(r'\end{table}', file=out_file)


def create_hypothesis_tests_table(
    directories: dict[str, Path],
    property_name: str,
    caption: str,
    out_file_name: Path,
) -> None:
    """Create hypothesis tests table."""
    out_file_name.parent.mkdir(parents=True, exist_ok=True)
    dataset_name_translation = {
        'computer_science': 'cs',
        'engineering': 'eess',
        'mathematics': 'math',
        'statistics': 'stat',
    }
    with open(out_file_name, 'w') as out_file:
        print(r'\begin{table} [h] \centering \caption{', file=out_file, end='')
        print(caption, file=out_file, end='')
        print(r'} \label{tab:edge_counts_hypothesis_tests}', file=out_file)
        print(r'    \begin{tabular}{|l|r|r|r|r|r|r|} \hline', file=out_file)
        print(r'        dataset & dataset value & $\hat\alpha$ & $\hat\b$ & location & scale & $p$-value \\ \hline', file=out_file)
        for name_in_table, name_to_print in dataset_name_translation.items():
            info = pd.read_csv(directories[name_in_table] / 'model_test' / property_name / 'distribution_info.csv')
            test = pd.read_csv(directories[name_in_table] / 'model_test' / property_name / 'test_results.csv')
            print(f'        {name_to_print: <7} & ', end='', file=out_file)
            print(f'${test['point_value'][0]:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info['theoretical_alpha'][0]:.4f}$ & ', end='', file=out_file)
            print(f'${info['theoretical_beta'][0]:.1f}$ & ', end='', file=out_file)
            print(f'${info['theoretical_location'][0]:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info['theoretical_scale'][0]:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${test['point_p_value'][0]:.4f}$', end='', file=out_file)
            print(r'\\', file=out_file)
        print(r'        \hline', file=out_file)
        print(r'    \end{tabular}', file=out_file)
        print(r'\end{table}', file=out_file)


def create_triangle_counts_hypothesis_tests_table(
    directories: dict[str, Path],
    out_file_name: Path,
) -> None:
    out_file_name.parent.mkdir(parents=True, exist_ok=True)
    dataset_name_translation = {
        'computer_science': 'cs',
        'engineering': 'eess',
        'mathematics': 'math',
        'statistics': 'stat',
    }
    with open(out_file_name, 'w') as out_file:
        print(r'\begin{table} [h] \centering \caption{Hypothesis tests for the edge counts} \label{tab:edge_counts_hypothesis_tests}', file=out_file)
        print(r'    \begin{tabular}{|l|r|r|r|r|r|r|} \hline', file=out_file)
        print(r'        dataset & edge count & $\hat\alpha$ & $\hat\b$ & location & scale & $p$-value \\ \hline', file=out_file)
        for name_in_table, name_to_print in dataset_name_translation.items():
            info = pd.read_csv(directories[name_in_table] / 'model_test' / 'num_of_triangles_stable' / 'distribution_info.csv')
            test = pd.read_csv(directories[name_in_table] / 'model_test' / 'num_of_triangles_stable' / 'test_results.csv')
            print(f'        {name_to_print: <7} & ', end='', file=out_file)
            print(f'${test['point_value'][0]:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info['theoretical_alpha'][0]:.4f}$ & ', end='', file=out_file)
            print(f'${info['theoretical_beta'][0]:.1f}$ & ', end='', file=out_file)
            print(f'${info['theoretical_location'][0]:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${info['theoretical_scale'][0]:,.0f}$ & '.replace(',', r'\,'), end='', file=out_file)
            print(f'${test['point_p_value'][0]:.4f}$', end='', file=out_file)
            print(r'\\', file=out_file)
        print(r'        \hline', file=out_file)
        print(r'    \end{tabular}', file=out_file)
        print(r'\end{table}', file=out_file)


# \begin{table}[!htb]
#   \begin{minipage}{.5\linewidth} \centering
#     \centering \caption{Number of simplices of different dimensions in the datasets} \label{simplex_counts_table}
#   \begin{tabular}{|l|r|r|r|} \hline
#     dataset & vertices &       edges &   triangles \\ \hline
#     cs      & 433\,244 & 2\,073\,235 & 4\,055\,220 \\
#     eess    &  77\,686 &    276\,947 &    562\,382 \\
#     math    & 198\,601 &    455\,130 &    321\,406 \\
#     stat    &  44\,380 &    114\,003 &    135\,800 \\ \hline
#   \end{tabular}
#   \end{minipage}
#   \qquad \qquad
#   \begin{minipage}{.5\linewidth} \centering
#     \caption{Estimated parameters of the stable distributions for triangle counts} \label{triangle_count_stable_distribution_parameters}
#     \begin{tabular}{|l|r|r|r|r|r|r|} \hline
#       dataset & $\hat\alpha$ & $\hat\b$ &     location &    scale \\ \hline
#       cs      &         1.39 &      1.0 & 18\,785\,263 & 504\,582 \\
#       eess    &         1.98 &      1.0 &  1\,911\,396 &  38\,527 \\
#       math    &         1.79 &      1.0 &  2\,027\,542 &  28\,774 \\
#       stat    &         1.96 &      1.0 &     566\,665 &  15\,352 \\ \hline
#     \end{tabular}
#   \end{minipage}
# \end{table}

# \begin{table}[!htb]
#   \begin{minipage}[t]{.5\linewidth}
#     \centering \caption{Betti numbers of the datasets} \label{betti_numbers_table}
#       \begin{tabular}{|l|r|r|} \hline
#         dataset & Betti-0 &  Betti-1 \\ \hline
#         cs      & 22\,576 & 168\,770 \\
#         eess    &  5\,533 &   7\,419 \\
#         math    & 26\,197 &  78\,009 \\
#         stat    &  4\,049 &   7\,275 \\ \hline
#       \end{tabular}
#     \end{minipage}
# \qquad \qquad
#   \begin{minipage}[t]{.5\linewidth}
#     \centering \caption{Parameter estimates of the stable distributions for Betti-1} \label{betti_number_stable_distribution_parameters}
#     \begin{tabular}{|l|r|r|r|r|} \hline
#       dataset & $\hat\a$ & $\hat\b$ & location & scale \\ \hline
#       cs      &     1.39 &     -1.0 &       37 & 12.83 \\
#       eess    &     1.98 &     -1.0 &      105 &  9.66 \\
#       math    &     1.79 &     -1.0 &      490 & 19.16 \\
#       stat    &     1.96 &     -1.0 &      126 &  8.99 \\ \hline
#     \end{tabular}
#   \end{minipage}
# \end{table}


def create_hypothesis_tests_plot(
    in_dir: Path,
    out_file_name: Path,
    csv_file_dir: Path,
    dataset_name: str | int,
    caption: str,
) -> None:
    """Create value counts log plots."""
    value_counts_file_name = in_dir / 'value_counts.csv'
    if not value_counts_file_name.is_file():
        return

    info = pd.read_csv(
        in_dir / 'distribution_info.csv',
        converters={'theoretical_domain': pd.eval}
    )

    # read data
    value_counts = pd.read_csv(value_counts_file_name, dtype=int, header=None)
    value_counts.columns = ['value', 'count']
    histogram = pd.read_csv(in_dir / 'histogram_linear.csv')
    pdfs = pd.read_csv(in_dir / 'pdfs.csv', dtype=float)
    pdfs.columns = ['x', 'empirical', 'theoretical']

    test_info = pd.read_csv(in_dir / 'test_results.csv')
    point_value = test_info['point_value'].iloc[0]
    point_p_value = test_info['point_p_value'].iloc[0]

    # determine x domain
    x_min = info['empirical_domain_min'].iloc[0]
    x_max = info['empirical_domain_max'].iloc[0]

    if point_value < x_min:
        point_value_position = 'left_from_distribution'
    elif point_value > x_max:
        point_value_position = 'right_from_distribution'
    else:
        point_value_position = 'inside_distribution'

    x_range = x_max - x_min
    x_min -= 0.05 * x_range
    x_max += 0.05 * x_range

    # determine_y_domain
    y_min = 0
    y_max = max(pdfs['theoretical'].max(), histogram['value'].max()) * 1.1

    def get_xticks() -> list[int]:
        digits_to_round_xticks = int(np.log10(x_max - x_min)) - 1
        num_of_x_ticks = 2
        xticks = [
            int(round(x_min + (i + 1) * (x_max - x_min) / (num_of_x_ticks + 1), -digits_to_round_xticks))
            for i in range(num_of_x_ticks)
        ]
        return xticks

    def print_histogram_group_plot(side: str) -> None:
        xticks = get_xticks()
        print(r'                \nextgroupplot[', file=out_file)
        if side == 'right':
            print(r'                    axis x discontinuity=parallel,', file=out_file)
            print(f'                    xmin={x_min - 0.1 * x_range:.4f}, xmax={x_max:.4f},', file=out_file)
        elif side == 'left':
            print(f'                    xmin={x_min:.4f}, xmax={x_max:.4f},', file=out_file)
        print(r'                    xtick={', end='', file=out_file)
        for xtick in xticks[:-1]:
            print(f'{xtick}, ', end='', file=out_file)
        print(f'{xticks[-1]}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                    xticklabels={', end='', file=out_file)
        for xtick in xticks[:-1]:
            print(f'${xtick:,}$'.replace(',', r'\,'), end='', file=out_file)
            print(', ', end='', file=out_file)
        print(f'${xticks[-1]:,}$'.replace(',', r'\,'), end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                    ytick=\empty,', file=out_file)
        print(r'                    yticklabel=\empty,', file=out_file)
        if side == 'right':
            print(r'                    axis y line=right,', file=out_file)
        elif side == 'left':
            print(r'                    axis y line=left,', file=out_file)
        print(r'                    scaled x ticks=false,', file=out_file)
        print(r'                    scaled y ticks=false,', file=out_file)
        print(r'                    xtick scale label code/.code={},', file=out_file)
        print(r'                    xlabel=\phantom a,', file=out_file)
        print(r'                    ylabel=\empty,', file=out_file)
        print(r'                    width=7.0cm', file=out_file)
        print(r'                ]', file=out_file)
        print(r'                \addplot[ybar interval, mark=no, fill=gray!50!white,]', file=out_file)
        print(f'                table[x={dataset_name}_bin_left_limit, y={dataset_name}_value, col sep=comma]', file=out_file)
        print(r'                    {', end='', file=out_file)
        print(f'{csv_file_dir}/linear_histograms.csv', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'                \addplot[ultra thick, red, no marks, restrict x to domain={', end='', file=out_file)
        print(f'{x_min}:{x_max}', end='', file=out_file)
        print(r'}]', file=out_file)
        print(f'                table [x={dataset_name}_value, y={dataset_name}_pdf, col sep=comma]', file=out_file)
        print(r'                    {', end='', file=out_file)
        print(f'{csv_file_dir}/theoretical_pdf.csv', end='', file=out_file)
        print(r'};', file=out_file)

    def print_point_value_groupplot(side: str) -> None:
        print(r'                \nextgroupplot[', file=out_file)
        if side == 'right':
            print(r'                    axis x discontinuity=parallel,', file=out_file)
        print(r'                    scaled x ticks=false,', file=out_file)
        # print(r'                    scaled y ticks=false,', file=out_file)
        print(f'                    xmin={point_value - 10:.4f}, xmax={point_value + 10:.4f},', file=out_file)
        print(r'                    xtick={', end='', file=out_file)
        print(f'{point_value}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                    xticklabels={', end='', file=out_file)
        print(f'${point_value:,}$'.replace(',', r'\,'), end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                    ytick=\empty,', file=out_file)
        print(r'                    yticklabel=\empty,', file=out_file)
        if side == 'right':
            print(r'                    axis y line=right,', file=out_file)
        elif side == 'left':
            print(r'                    axis y line=left,', file=out_file)
        print(r'                    axis line style={-},', file=out_file)
        print(r'                    width=3.5cm', file=out_file)
        print(r'                ]', file=out_file)
        print(r'                \addplot[ultra thick, green!70!black, mark=none, dashed]', file=out_file)
        print(r'                coordinates {', end='', file=out_file)
        print(f'({point_value}, 0) ({point_value}, {y_max})', end='', file=out_file)
        print(r'};', file=out_file)

    def print_joint_group_plot() -> None:
        xticks = get_xticks()
        print(r'                \nextgroupplot[', file=out_file)
        print(f'                    xmin={x_min:.4f}, xmax={x_max:.4f},', file=out_file)
        print(r'                    xtick={', end='', file=out_file)
        for xtick in xticks[:-1]:
            print(f'{xtick}, ', end='', file=out_file)
        print(f'{xticks[-1]}', end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                    xticklabels={', end='', file=out_file)
        for xtick in xticks[:-1]:
            print(f'${xtick:,}$'.replace(',', r'\,'), end='', file=out_file)
            print(', ', end='', file=out_file)
        print(f'${xticks[-1]:,}$'.replace(',', r'\,'), end='', file=out_file)
        print(r'},', file=out_file)
        print(r'                    ylabel=\empty,', file=out_file)
        print(r'                    ytick=\empty,', file=out_file)
        print(r'                    yticklabel=\empty,', file=out_file)
        print(r'                    axis line style={-},', file=out_file)
        print(r'                    width=8.375cm,', file=out_file)
        print(r'                    xtick scale label code/.code={},', file=out_file)
        print(r'                    scaled x ticks=false,', file=out_file)
        print(r'                    scaled y ticks=false', file=out_file)
        print(r'                ]', file=out_file)
        print(r'                \addplot[ybar interval, mark=no, fill=gray!50!white,]', file=out_file)
        print(f'                table[x={dataset_name}_bin_left_limit, y={dataset_name}_value, col sep=comma]', file=out_file)
        print(r'                    {', end='', file=out_file)
        print(f'{csv_file_dir}/linear_histograms.csv', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'                \addplot[ultra thick, red, no marks]', file=out_file)
        print(f'                table [x={dataset_name}_value, y={dataset_name}_pdf, col sep=comma]', file=out_file)
        print(r'                    {', end='', file=out_file)
        print(f'{csv_file_dir}/theoretical_pdf.csv', end='', file=out_file)
        print(r'};', file=out_file)
        print(r'                \addplot[ultra thick, green!70!black, mark=none, dashed]', file=out_file)
        print(r'                coordinates {', end='', file=out_file)
        print(f'({point_value}, 0) ({point_value}, {y_max})', end='', file=out_file)
        print(r'};', file=out_file)

    # create plot
    out_file_name.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file_name, 'w') as out_file:
        print(r'\begin{subfigure}[b]{0.22\textwidth}', file=out_file)
        print(r'    \resizebox{\textwidth}{!}{', file=out_file)
        print(r'        \begin{tikzpicture}', file=out_file)
        print(r'            \begin{groupplot} [', file=out_file)
        print(f'                    ymin={y_min}, ymax={y_max},', file=out_file)
        if point_value_position in ['right_from_distribution', 'left_from_distribution']:
            print(r'                    height=7.5cm,', file=out_file)
        elif point_value_position == 'inside_distribution':
            print(r'                    height=7.0cm,', file=out_file)
        print(r'                    group style={', file=out_file)
        print(r'                            xticklabels at=edge bottom,', file=out_file)
        print(r'                            horizontal sep=0pt,', file=out_file)
        if point_value_position in ['right_from_distribution', 'left_from_distribution']:
            print(r'                            group size=2 by 1', file=out_file)
        elif point_value_position == 'inside_distribution':
            print(r'                            group size=1 by 1', file=out_file)
        print(r'                        }', file=out_file)
        print(r'            ]', file=out_file)

        if point_value_position in ['right_from_distribution', 'left_from_distribution']:
            if point_value_position == 'left_from_distribution':
                print_point_value_groupplot('left')
                print_histogram_group_plot('right')
            elif point_value_position == 'right_from_distribution':
                print_histogram_group_plot('left')
                print_point_value_groupplot('right')
        else:
            print_joint_group_plot()

        print(r'            \end{groupplot}', file=out_file)
        print(r'            \node[anchor=north east, draw] at ', end='', file=out_file)
        if point_value_position == 'right_from_distribution':
            print('(4.5, 5.5)', file=out_file)
        if point_value_position == 'inside_distribution':
            print('(4.5, 5.0)', file=out_file)
        if point_value_position == 'left_from_distribution':
            print('(4.5, 5.5)', file=out_file)
        print(r'            {\small{\textbf{p-value:} $\mathbf{', end='', file=out_file)
        print(f'{point_p_value:.4f}', end='', file=out_file)
        print(r'}$}};', file=out_file)
        print(r'        \end{tikzpicture}', file=out_file)
        print(r'    }', file=out_file)
        print(r'    \caption{', end='', file=out_file)
        print(f'{caption}', end='', file=out_file)
        print(r'}', file=out_file)
        print(r'\end{subfigure}', file=out_file)
