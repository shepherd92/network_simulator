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
    num_of_processes: int = 4
    root_dir = Path('/home/au725389/save/research/tda/code/')
    output_dir = Path('../output')
