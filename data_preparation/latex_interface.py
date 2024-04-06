#!/usr/bin/env python3
"""Create latex plots."""

from pathlib import Path

import pandas as pd
import numpy as np


def create_boxplots(
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


def create_histograms_normal_plots(
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
        print(f'                    xmin={x_min:.4f}, xmax={x_max:.4f}, ymin={y_min:.4f}, ymax={y_max:.4f},', file=out_file)
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
