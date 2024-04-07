from os import listdir
from os.path import isfile, join
import re
from pathlib import Path
from shutil import rmtree

python_files_with_data_dirs = Path('data_preparation/data_dirs')
data_dirs = Path('../output')

all_file_content = ''
for filename in listdir(python_files_with_data_dirs):
    full_path = join(str(python_files_with_data_dirs), filename)
    if isfile(full_path) and full_path.endswith(".py"):
        with open(full_path) as file:
            all_file_content += (file.read())

directories_to_keep = set(re.findall(r'202[0-9]{5}_[0-9]{6}', all_file_content))
directories_to_keep.add('_prepared_data')

directories_to_remove = set(listdir(data_dirs)) - directories_to_keep
print(f'Directories to remove: {directories_to_remove}')
for filename in sorted(directories_to_remove):
    rmtree(join(str(data_dirs), filename))
