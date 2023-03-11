#!/usr/bin/env python3
"""Parameter options for model fitting / optimization."""

from logging import warning
from typing import NamedTuple

import numpy as np

from model.model import Model


class ParameterOption(NamedTuple):
    """Represent a single model parameter option for fitting."""

    parameter_name: str = ''
    initial_guess: float | int = 0.
    lower_bound: float | int = 0.
    upper_bound: float | int = 0.
    is_integer: bool = False
    optimize: bool = True
    step_for_gradient: float = 1. if is_integer else 1e-3


class ModelParameterOptions:
    """Represent a base class for model parameter option sets."""

    def __init__(self) -> None:
        """Create default model parameter options."""
        self._current_index = 0
        self.options: tuple[ParameterOption, ...] = ()

    def get_parameter_names_to_optimize(self) -> list[str]:
        """Extract parameter names to be optimized."""
        return [
            parameter_option.parameter_name
            for parameter_option in self.options
            if parameter_option.optimize
        ]

    def create_initial_guess(self) -> tuple[float, ...]:
        """Convert parameter options to initial guess tuple."""
        initial_guess = tuple(
            parameter_option.initial_guess
            for parameter_option in self.options
            if parameter_option.optimize
        )
        return initial_guess

    def create_bounds_tuple(self) -> tuple[tuple[float, float], ...]:
        """Convert parameter options to parameter bounds tuple."""
        bounds = (
            (parameter_option.lower_bound, parameter_option.upper_bound)
            for parameter_option in self.options
            if parameter_option.optimize
        )
        return tuple(bounds)

    def create_model_parameters_from_initial_guess(self, model: Model) -> Model.Parameters:
        """Convert parameter options to Parameters object of initial guesses."""
        model_parameters = model.Parameters()

        for parameter_option in self.options:
            setattr(
                model_parameters,
                parameter_option.parameter_name,
                parameter_option.initial_guess
            )

        return model_parameters

    def create_model_parameters_from_guess(self, model: Model, guess: np.ndarray,) -> Model.Parameters:
        """Convert an array of parmeters to Parameters object."""
        model_parameters = model.Parameters()

        guess_current_index = 0

        for parameter_option in self.options:
            if parameter_option.optimize:
                value = \
                    type(getattr(model_parameters, parameter_option.parameter_name))(
                        guess[guess_current_index]
                    )
                if value < parameter_option.lower_bound:
                    warning(
                        f'Guessed parameter {parameter_option.parameter_name} ' +
                        f'has value {value}, which smaller than lower bound: ' +
                        f'{parameter_option.lower_bound}.'
                    )
                    # pull the parameter inside the bounds and add a small random number
                    value = parameter_option.lower_bound + \
                        np.random.random() * 1e-2 * np.abs(parameter_option.lower_bound)
                if value > parameter_option.upper_bound:
                    warning(
                        f'Guessed parameter {parameter_option.parameter_name} ' +
                        f'has value {value}, which is greater than upper bound: ' +
                        f'{parameter_option.upper_bound}.'
                    )
                    # pull the parameter inside the bounds and subtract a small random number
                    value = parameter_option.upper_bound - \
                        np.random.random() * 1e-2 * np.abs(parameter_option.upper_bound)
                setattr(
                    model_parameters,
                    parameter_option.parameter_name,
                    value
                )
                guess_current_index += 1
            else:
                setattr(
                    model_parameters,
                    parameter_option.parameter_name,
                    parameter_option.initial_guess
                )

        return model_parameters
