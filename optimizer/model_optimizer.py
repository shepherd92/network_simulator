#!/usr/bin/env python3
"""This module is responsible for fitting a module to a dataset."""

from logging import info, warning
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import minimize, OptimizeResult

from model.model import Model
from network.property import DerivedNetworkProperty
from optimizer.parameter_option import ModelParameterOptions


class ModelOptimizer:
    """Optimize a model with the given method.

    The main purpose is to give a target network, and fit a model to some summaries of this target
    network by adjusting some given parameters within the given bounds.
    """

    def __init__(self, model: Model, method: str) -> None:
        """Create a model optimizer."""
        self._model = model
        self._method = method

    def fit(
        self,
        parameter_options: ModelParameterOptions,
        target_property_params: list[DerivedNetworkProperty],
        target_values: list[float | int],
        options: dict[str, Any]
    ) -> Model.Parameters:
        """Fit model parameters to data.

        The initial guess for the parameters is the same as the values set for the model.

        Args:
            parameter_options (ModelParameterOptions): a tuple of ParameterOptions specifying the
                optimization options; should contain the same number of members as the number of
                parameters to be optimized
            target_property_params: (list[ScalarNetworkProperty]): targets with calculators
            target_values (list[float | int, ...]): target values to fit the model to
            options: (dict[str, Any]): options for the optimization algorithm
        Returns:
            optimized model parameters
        """
        info(
            'Model fitting started: ' +
            f'model type {type(self.model)}, ' +
            f'parameter_options: {parameter_options.options}'
        )

        initial_guess = parameter_options.create_initial_guess()
        if not any((option.optimize for option in parameter_options.options)):
            warning('No parameters set to be optimized, exiting optimization.')
            model_parameters = parameter_options.create_model_parameters_from_guess(
                self.model,
                initial_guess,
            )
            return model_parameters

        bounds = parameter_options.create_bounds_tuple()
        optimization_options = self._create_optimization_options(parameter_options)
        options |= optimization_options

        result: OptimizeResult = minimize(
            ModelOptimizer._objective_function,
            x0=initial_guess,
            args=(
                self.model,
                target_property_params,
                target_values,
                parameter_options,
            ),
            method=self.method,
            bounds=bounds,
            options=options,
        )

        if result.success:
            info('Network fitting succeeded.')
        else:
            warning(f'Network fitting failed with message {result.message}.')

        optimized_model_parameters = parameter_options.create_model_parameters_from_guess(self.model, result.x)

        info(f'Model fitting finished: fitted parameters are {optimized_model_parameters}.')

        return optimized_model_parameters

    @staticmethod
    def _objective_function(
        guessed_parameters: npt.NDArray[np.float_],
        model: Model,
        target_property_params: list[DerivedNetworkProperty],
        target_values: list[float | int],
        parameter_options: ModelParameterOptions,
    ) -> float:
        """Target function to be used for optimization."""
        model.parameters = parameter_options.create_model_parameters_from_guess(
            model,
            guessed_parameters
        )

        num_of_simulations, network_properties_for_estimation = \
            ModelOptimizer._estimate_required_number_of_simulations(
                model,
                target_property_params,
            )
        num_of_additional_simulations = \
            max(num_of_simulations - len(network_properties_for_estimation), 0)

        if num_of_additional_simulations != 0:
            additional_network_properties = model.simulate(
                target_property_params,
                num_of_additional_simulations,
            )
            network_properties = \
                pd.concat([network_properties_for_estimation, additional_network_properties])
        else:
            network_properties = network_properties_for_estimation

        summary_means = network_properties.mean()

        # sum all the relative differences of the targets
        loss = 1. / len(target_property_params) * np.sqrt(np.sum([
            (summary - target)**2 / target
            for summary, target in zip(summary_means, target_values)
        ]))

        info(f'Guessing parameters {model.parameters}; loss: {loss:.6f}')

        return loss

    @staticmethod
    def _estimate_required_number_of_simulations(
        model: Model,
        target_property_params: list[DerivedNetworkProperty],
    ) -> tuple[int, pd.DataFrame]:
        """Estmate how many simulations need to be executed.

        The goal is to reduce the relative standard deviation of the calculated summaries below a
        threshold.
        """
        num_of_simulations_to_estimate_mean_and_std = 100
        max_relative_std = 0.01

        network_properties_to_estimate_mean_and_std = model.simulate(
            target_property_params,
            num_of_simulations_to_estimate_mean_and_std,
        )
        means_estimate = network_properties_to_estimate_mean_and_std.mean().to_numpy()
        vars_estimate = network_properties_to_estimate_mean_and_std.var().to_numpy()

        num_of_simulations = int(np.ceil(np.nanmax(
            1 + vars_estimate / (max_relative_std * means_estimate)**2,
            initial=100
        )))

        return num_of_simulations, network_properties_to_estimate_mean_and_std

    @staticmethod
    def _create_optimization_options(parameter_options: ModelParameterOptions) -> dict[str, Any]:
        optimization_options = {'eps': np.array([
            parameter_option.step_for_gradient
            for parameter_option in parameter_options.options
            if parameter_option.optimize
        ])}
        return optimization_options

    @property
    def model(self) -> Model:
        """Get the network model to be optimized."""
        return self._model

    @property
    def model_parameters_type(self) -> type:
        """Return the model parameters type."""
        return self.model.Parameters

    @property
    def method(self) -> str:
        """Return the optimization method."""
        return self._method
