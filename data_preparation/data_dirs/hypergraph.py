#!/usr/bin/env python
"""Merge results from different runs for easier comparison and analysis."""

from pathlib import Path

input_base_path = Path('../output')
output_path_root = Path('../output/_prepared_data')

model_sample_directory = input_base_path / '20240406_094950'

data_analysis_directories = {
    'computer_science': input_base_path / '20240407_073645',
    # 'biology':          input_base_path / '20240407_053735',
    # 'economics':        input_base_path / '20240407_053830',
    'engineering':      input_base_path / '20240407_054143',
    # 'finance':          input_base_path / '20240407_053843',
    'mathematics':      input_base_path / '20240407_054707',
    'statistics':       input_base_path / '20240407_054022',
}

model_sample_for_data_directories = {
    'computer_science': input_base_path / '20240330_120256',
    # 'biology':          input_base_path / '20240330_115643',
    # 'economics':        input_base_path / '20240330_115023',
    'engineering':      input_base_path / '20240330_113433',
    # 'finance':          input_base_path / '20240330_113104',
    'mathematics':      input_base_path / '20240330_110429',  # '20240316_210724',
    'statistics':       input_base_path / '20240330_094048',
}

degree_distribution_directories = {
    # key is (network size, minimum degree)
    # (10, 10):       input_base_path / '20240405_132454',
    (100, 5):      input_base_path / '20240409_143839',
    (1000, 5):     input_base_path / '20240409_143905',
    (10000, 5):    input_base_path / '20240409_143916',
    (100000, 5):   input_base_path / '20240409_144046',
    (100, 8):      input_base_path / '20240409_160840',
    (1000, 8):     input_base_path / '20240409_160848',
    (10000, 8):    input_base_path / '20240409_160859',
    (100000, 8):   input_base_path / '20240409_160916',
    (100, 10):      input_base_path / '20240405_132519',
    (1000, 10):     input_base_path / '20240405_132526',
    (10000, 10):    input_base_path / '20240405_132536',
    (100000, 10):   input_base_path / '20240405_132544',
    # (1000000, 10):  input_base_path / '20240323_185549',
    (100, 12):      input_base_path / '20240409_161132',
    (1000, 12):     input_base_path / '20240409_161141',
    (10000, 12):    input_base_path / '20240409_161156',
    (100000, 12):   input_base_path / '20240409_161210',
    (100, 15):      input_base_path / '20240409_161644',
    (1000, 15):     input_base_path / '20240409_161658',
    (10000, 15):    input_base_path / '20240409_161711',
    (100000, 15):   input_base_path / '20240409_161718',
}

simplex_count_directories = {
    # the key is <gamma_prime*100>_<gamma*100>
    (5, 25): input_base_path / '20240404_120532',
    (5, 50): input_base_path / '20240404_131443',
    (5, 75): input_base_path / '20240404_120519',
    (10, 25): input_base_path / '20240404_120510',
    (10, 50): input_base_path / '20240404_120505',
    (10, 75): input_base_path / '20240404_120459',
    (15, 25): input_base_path / '20240404_112232',
    (15, 50): input_base_path / '20240404_112223',
    (15, 75): input_base_path / '20240404_112216',
    (20, 25): input_base_path / '20240404_112203',
    (20, 50): input_base_path / '20240404_112159',
    (20, 75): input_base_path / '20240404_112153',
    (25, 25): input_base_path / '20240404_102317',
    (25, 50): input_base_path / '20240404_102313',
    (25, 75): input_base_path / '20240404_102305',
    (30, 25): input_base_path / '20240404_102255',
    (30, 50): input_base_path / '20240404_102249',
    (30, 75): input_base_path / '20240404_092900',
}

betti_number_directories = {
    # the key is (<gamma_prime*100>, <gamma*100>)
    (5, 25): input_base_path / '20240404_120532',
    (5, 50): input_base_path / '20240404_131443',
    (5, 75): input_base_path / '20240404_120519',
    (10, 25): input_base_path / '20240404_120510',
    (10, 50): input_base_path / '20240404_120505',
    (10, 75): input_base_path / '20240404_120459',
    (15, 25): input_base_path / '20240404_112232',
    (15, 50): input_base_path / '20240404_112223',
    (15, 75): input_base_path / '20240404_112216',
    (20, 25): input_base_path / '20240404_112203',
    (20, 50): input_base_path / '20240404_112159',
    (20, 75): input_base_path / '20240404_112153',
    (25, 25): input_base_path / '20240404_102317',
    (25, 50): input_base_path / '20240404_102313',
    (25, 75): input_base_path / '20240404_102305',
    (30, 25): input_base_path / '20240404_102255',
    (30, 50): input_base_path / '20240404_102249',
    (30, 75): input_base_path / '20240404_092900',
}

hypothesis_testing_directories = {
    'computer_science': input_base_path / '20240322_200625',
    'engineering':      input_base_path / '20240322_092931',
    'mathematics':      input_base_path / '20240321_221229',
    'statistics':       input_base_path / '20240408_185157',
}
