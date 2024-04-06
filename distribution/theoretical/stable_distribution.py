#!/usr/bin/env python3
"""Levy alpha-stable theoretical distribution."""

from __future__ import annotations

from dataclasses import dataclass, astuple
from enum import Enum, auto
from logging import warning

import numpy as np
import numpy.typing as npt
from scipy.optimize import OptimizeResult, minimize
from scipy.stats import levy_stable
from scipy.stats._warnings_errors import FitError

from distribution.distribution import Distribution
from distribution.empirical_distribution import EmpiricalDistribution
from distribution.theoretical.theoretical_distribution import TheoreticalDistribution


class StableDistribution(TheoreticalDistribution):
    """Levy-stable theoretical distribution."""

    @dataclass
    class DistributionParameters(TheoreticalDistribution.DistributionParameters):
        """Parameters of the stable distribution."""

        alpha: float = np.nan
        beta: float = np.nan
        location: float = np.nan
        scale: float = np.nan

    @dataclass
    class DomainCalculation(TheoreticalDistribution.FittingParameters.DomainCalculation):
        """Parameters of domain calculation."""

    @dataclass
    class ParameterFitting(TheoreticalDistribution.FittingParameters.ParameterFitting):
        """Parameters of fitting."""

        class Method(Enum):
            """Method used for fitting the power law distribution."""

            QUICK_FIT = auto()
            MLE_SCIPY = auto()
            OPTIMIZATION = auto()

        method: Method
        fixed_parameters: StableDistribution.DistributionParameters

    def __init__(self) -> None:
        """Create a default stable distribution."""
        super().__init__()
        self._parameters = StableDistribution.DistributionParameters()

    def calc_quantiles(self, quantiles_to_calculate: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        assert ((quantiles_to_calculate >= 0.) & (quantiles_to_calculate <= 1.)).all(), \
            f'Quntiles to calculate must be in [0, 1], but they are {quantiles_to_calculate}'
        quantiles = levy_stable.ppf(quantiles_to_calculate, *astuple(self._parameters))
        return quantiles

    def info(self) -> dict[str, int | float]:
        """Return a dict representation based on the distribution properties."""
        result = super().info()
        result.update({
            'distribution_type': 'stable',
            'alpha': self._parameters.alpha,
            'beta': self._parameters.beta,
            'location': self._parameters.location,
            'scale': self._parameters.scale,
        })
        return result

    def _fit_domain(
        self,
        empirical_distribution: EmpiricalDistribution,
        domain_calculation_parameters: DomainCalculation,
    ) -> None:
        self._domain = Distribution.Domain(-np.inf, np.inf)

    def _fit_parameters(
        self,
        empirical_distribution: EmpiricalDistribution,
        parameter_fitting_parameters: ParameterFitting
    ) -> None:
        Method = StableDistribution.ParameterFitting.Method
        if parameter_fitting_parameters.method == Method.MLE_SCIPY:
            self._parameters = self._estimate_parameters_mle_scipy(
                empirical_distribution,
                parameter_fitting_parameters.fixed_parameters,
            )
        elif parameter_fitting_parameters.method == Method.QUICK_FIT:
            warning(f'{Method.QUICK_FIT} cannot handle fixed parameters.')
            self._parameters = self._estimate_parameters_quick_fit(
                empirical_distribution,
            )
        elif parameter_fitting_parameters.method == Method.OPTIMIZATION:
            warning(f'{Method.OPTIMIZATION} cannot handle fixed parameters.')
            self._parameters = self._estimate_parameters_optimization(
                empirical_distribution,
            )
        else:
            assert False, f'Unknown fitting method: {parameter_fitting_parameters.method}.'

    def _estimate_parameters_mle_scipy(
        self,
        empirical_distribution: EmpiricalDistribution,
        fixed_parameters: DistributionParameters,
    ) -> DistributionParameters:
        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)

        kwargs: dict[str, np.float_] = {}
        if not np.isnan(fixed_parameters.alpha):
            kwargs['f0'] = fixed_parameters.alpha
        if not np.isnan(fixed_parameters.beta):
            kwargs['f1'] = fixed_parameters.beta
        if not np.isnan(fixed_parameters.location):
            kwargs['floc'] = fixed_parameters.location
        if not np.isnan(fixed_parameters.scale):
            kwargs['fscale'] = fixed_parameters.scale

        try:
            params = levy_stable.fit(value_sequence, **kwargs)
        except FitError:
            # sometimes the fitting converges to parameters outside the domain
            params = (np.nan, np.nan, np.nan, np.nan)

        return StableDistribution.DistributionParameters(
            params[0],
            params[1],
            params[2],
            params[3]
        )

    def _estimate_parameters_quick_fit(
        self,
        empirical_distribution: EmpiricalDistribution
    ) -> DistributionParameters:
        def pconv(
            alpha: float,
            beta: float,
            location: float,
            sigma: float
        ) -> StableDistribution.DistributionParameters:
            result = StableDistribution.DistributionParameters(
                alpha=alpha,
                beta=beta,
                location=location - sigma * beta * np.tan(np.pi * alpha / 2.),
                scale=sigma
            )
            return result

        value_sequence = empirical_distribution.get_value_sequence_in_domain(self.domain)
        estimated_parameters = pconv(*levy_stable._fitstart(value_sequence))  # pylint: disable=protected-access
        return estimated_parameters

    def _estimate_parameters_optimization(
        self,
        empirical_distribution: EmpiricalDistribution
    ) -> DistributionParameters:

        initial_guess = (
            self.parameters.alpha,
            self.parameters.beta,
            self.parameters.location,
            self.parameters.scale
        )

        bounds: tuple[tuple[float, float], ...] = (
            (1., 2.),  # alpha bounds
            (0.999, 1.001),  # beta bounds
            (0., np.inf),  # location bounds
            (0.001, np.inf),  # scale bounds
        )

        eps = np.array([
            5e-3,  # alpha eps
            5e-3,  # beta eps
            1e-2 * self.parameters.location,  # location eps
            1e-2 * self.parameters.scale,  # scale eps
        ])

        # set validity to True so that the Kolmogorov-Smirnov can be calculated
        self._valid = True

        result: OptimizeResult = minimize(
            StableDistribution._objective_function_for_fitting,
            x0=initial_guess,
            args=(
                self,
                empirical_distribution,
            ),
            method='L-BFGS-B',
            bounds=bounds,
            options={
                'eps': eps
            }
        )
        self._valid = False

        fitted_prameters = StableDistribution.DistributionParameters(*result.x)
        return fitted_prameters

    @staticmethod
    def _objective_function_for_fitting(
        guess: npt.NDArray[np.float_],
        theoretical_distribution: StableDistribution,
        empirical_distribution: EmpiricalDistribution
    ) -> float:
        parameters = StableDistribution.DistributionParameters(*tuple(guess.tolist()))
        print(f'Guessing: {parameters}')
        theoretical_distribution.parameters = parameters
        kolmogorov_smirnov_statistic = theoretical_distribution.kolmogorov_smirnov(
            empirical_distribution,
            empirical_distribution.natural_x_values
        )
        print(f'KS-statistic: {kolmogorov_smirnov_statistic}')

        if np.isnan(kolmogorov_smirnov_statistic):
            kolmogorov_smirnov_statistic = np.inf

        return kolmogorov_smirnov_statistic

    def _pdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the PDF of the distribution evaluted at the given x_values."""
        pdf_values = levy_stable.pdf(x_values, *astuple(self._parameters))
        return pdf_values

    def _cdf_in_domain(self, x_values: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Return the CDF of the distribution evaluted at the given x_values."""
        cdf_values = levy_stable.cdf(x_values, *astuple(self._parameters))
        return cdf_values

    @property
    def parameters(self) -> StableDistribution.DistributionParameters:
        """Return the parameters of the distribution."""
        return self._parameters
