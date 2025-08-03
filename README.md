# Installation

python -mvenv ~/venv/network_simulator

. ~/venv/network_simulator/bin/activate

pip install -r requirements.txt

sudo apt update

sudo apt install libgudhi-dev python3-gudhi pybind11-dev

# Usage

. ~/venv/network_simulator/bin/activate

## Dataset Analysis

python3 main.py --mode analysis --config_dir config_files

## Model Testing

python3 main.py --mode testing --config_dir config_files

## Model Example

python3 main.py --mode example --config_dir config_files
