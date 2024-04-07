#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path

input_base_path = Path('../output')
output_path_root = Path('../output/_prepared_data')

model_sample_directory = input_base_path / '20240406_094950'

data_analysis_directories = {
    'computer_science': input_base_path / '20240407_073645',
    'biology':          input_base_path / '20240407_053735',
    'economics':        input_base_path / '20240407_053830',
    'engineering':      input_base_path / '20240407_054143',
    'finance':          input_base_path / '20240407_053843',
    'mathematics':      input_base_path / '20240407_054707',
    'statistics':       input_base_path / '20240407_054022',
}

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
    # '10':       input_base_path / '20240405_132454',
    '100':      input_base_path / '20240405_132519',
    '1000':     input_base_path / '20240405_132526',
    '10000':    input_base_path / '20240405_132536',
    '100000':   input_base_path / '20240405_132544',
    # '1000000':  input_base_path / '20240323_185549',
}

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

hypothesis_testing_directories = {
    'computer_science': input_base_path / '20240322_200625',
    'engineering':      input_base_path / '20240322_092931',
    'mathematics':      input_base_path / '20240321_221229',
    'statistics':       input_base_path / '20240322_085600',
    # 'finance':          input_base_path / '20240322_090205',
}
