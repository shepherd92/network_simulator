#!/usr/bin/env python
"""Configuration class."""

from dataclasses import dataclass

from pathlib import Path


@dataclass
class GeneralConfig:
    """Configuration of the execution of the program."""

    log_level = 'DEBUG'
    runtime_profiling = False
    memory_profiling = False
    tracing = False
    seed = 1733852529  # 0 for random seed
    root_dir = Path('/home/au725389/save/research/higher_order_adrcm/code/')
    output_dir = Path('../output')
