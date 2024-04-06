#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path

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
    'computer_science': input_base_path / '20231208_080033',
    'biology':          input_base_path / '20231208_083044',
    'economics':        input_base_path / '20231208_083202',
    'engineering':      input_base_path / '20231208_083229',
    'finance':          input_base_path / '20231208_083259',
    'mathematics':      input_base_path / '20231208_075316',
    'statistics':       input_base_path / '20231208_074547',
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
    '25': input_base_path / '20230801_200720',
    '40': input_base_path / '20230703_101414',
    '45': input_base_path / '20230703_130440',
    '50': input_base_path / '20230801_100128',
    '60': input_base_path / '20230731_124944',
    '75': input_base_path / '20230615_204859',
}

betti_number_directories = {
    # the key is the value of gamma * 100
    # '10': input_base_path / '20230608_191114',
    # '20': input_base_path / '20230608_191120',
    '25': input_base_path / '20230609_212258',  # 20230804_115431
    # '30': input_base_path / '20230608_191124',
    # '40': input_base_path / '20230608_191128',
    '50': input_base_path / '20230609_212245',  # 20230803_105044
    # '60': input_base_path / '20230608_205831',  # 20230803_100554
    '67': input_base_path / '20230809_173951',
    # '70': input_base_path / '20230608_205839',
    # '75': input_base_path / '20230609_142002',  # 20230804_115418
    # '80': input_base_path / '20230608_205845',
    # '90': input_base_path / '20230608_205850',
}

hypothesis_testing_directories = {
    'computer_science': input_base_path / '20230807_181159',  # 20230712_075521
    'biology':          input_base_path / '20230710_164059',  # 20230710_164059
    'economics':        input_base_path / '20230710_162805',  # 20230710_162805
    'engineering':      input_base_path / '20230806_210326',  # 20230711_110334
    'finance':          input_base_path / '20230710_162239',  # 20230710_162239
    'mathematics':      input_base_path / '20230806_200515',  # 20230710_200618
    'statistics':       input_base_path / '20230806_205545',  # 20230710_143554
}
