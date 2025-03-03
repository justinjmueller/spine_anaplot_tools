import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
import pandas as pd
import matplotlib.pyplot as plt

from artists import SpineArtist
from style import Style
from variable import Variable

class SpineSpectra(SpineArtist):
    """
    A base class designed to encapsulate spectra for multiple variables
    in an ensemble of samples. Though intended for use with a
    collection of variables, this class is meant to be an abstraction
    representing a single plot instance or plot content on an axis. 

    Attributes
    ----------
    _variables : list
        The list of Variable objects for the spectra.
    _categories : dict
        A dictionary of the categories for the spectra. This serves as
        a map between the category label in the input TTree and the
        category label for the aggregated data (and therefore what is
        shown in a single legend entry).
    _plotdata : dict
        A dictionary of the data for the spectra. This is a map between
        the category label for the spectra and the histogram data for
        that category.
    """
    def __init__(self, variables, categories, colors) -> None:
        """
        Initializes the SpineSpectra object with the given kwargs.

        Parameters
        ----------
        variables : list
            The list of Variable objects for the spectra.
        categories : dict
            A dictionary of the categories for the spectra. This serves
            as a map between the category label in the input TTree and
            the category label for the aggregated data (and therefore
            what is shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the spectra.
            This serves as a map between the category label for the
            spectra (value in the `_categories` dictionary) and the color
            to use for the histogram. The color can be any valid matplotlib
            color string or a cycle indicator (e.g. 'C0', 'C1', etc.).

        Returns
        -------
        None.
        """
        self._variables = variables
        self._categories = categories
        self._colors = colors
        self._plotdata = None
        self._binedges = None
        self._onebincount = None

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra object. The most basic form
        of this method is to save the exposure information from the
        ordinate sample.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra1D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)

    def mark_pot(self, ax) -> None:
        """
        Add the POT information to the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to add the POT information to.

        Returns
        -------
        None.
        """
        yrange = ax.get_ylim()
        usey = yrange[1] + 0.02*(yrange[1] - yrange[0])
        xrange = ax.get_xlim()
        usex = xrange[1] - 0.02*(xrange[1] - xrange[0])
        mag = int(np.floor(np.log10(self._exposure)))
        usepot = self._exposure/10**mag
        s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
        ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

    def mark_preliminary(self, ax, label) -> None:
        """
        Add a preliminary label to the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to add the preliminary label to.
        label : str
            The label to add to the plot to indicate that the plot is
            preliminary.

        Returns
        -------
        None.
        """
        yrange = ax.get_ylim()
        usey = yrange[1] + 0.025*(yrange[1] - yrange[0])
        xrange = ax.get_xlim()
        usex = xrange[0] + 0.025*(xrange[1] - xrange[0])
        ax.text(x=usex, y=usey, s=label, fontsize=14, color='#d67a11')

    def fit_with_function(self, ax, bin_centers, data, bin_edges, fit_type, range=(-1,1)) -> None:
        """
        Fit the data with a given function and plot the fit on the axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to plot the fit on.
        bin_centers : list
            The list of bin centers for the data.
        data : list
            The list of data values for the bins.
        bin_edges : list
            The list of bin edges for the data.
        fit_type : str
            The type of fit to perform on the data. The default is
            None, which will not perform any fit. The options are:
                'crystal_ball' - Perform a Crystal Ball fit on the data.
                'gaussian'     - Perform a Gaussian fit on the data.
        range : tuple, optional
            The range to fit the data over. The default is (-1,1).

        Returns
        -------
        None.
        """
        if fit_type == 'crystal_ball':
            # Crystal Ball fit
            # First, estimate the parameters "reasonably" well to avoid
            # the fit getting stuck in a weird place.
            cb_mean = np.average(bin_centers, weights=data)
            cb_sigma = np.sqrt(np.average((bin_centers - cb_mean)**2, weights=data))
            cb_norm = data * np.diff(bin_edges)
            cb_alpha = 1.5
            cb_n = 5
            initial_guess = [cb_alpha, cb_n, cb_mean, cb_sigma, np.sum(cb_norm)]
            popt, pcov = curve_fit(self.crystal_ball, bin_centers[0], data[0], p0=initial_guess)
            
            # Label with estimated parameters and +/- 1 sigma
            cb_label = f'Crystal Ball Fit\n'
            cb_label += f'$\\mu$={popt[2]:.2f}$\\pm${np.sqrt(pcov[2,2]):.2f}\n'
            cb_label += f'$\\sigma$={popt[3]:.2f}$\\pm${np.sqrt(pcov[3,3]):.2f}\n'
            cb_label += f'$\\alpha$={popt[0]:.2f}$\\pm${np.sqrt(pcov[0,0]):.2f}\n'
            cb_label += f'n={popt[1]:.2f}$\\pm${np.sqrt(pcov[1,1]):.2f}'
            
            # Plot the fit
            x = np.linspace(*range, 1000)
            ax.plot(x, self.crystal_ball(x, *popt), 'r-', label=cb_label)
        
        elif fit_type == 'gaussian':
            # Gaussian fit
            # First, estimate the parameters "reasonably" well to avoid
            # the fit getting stuck in a weird place. This fit can use
            # the "standard" Gaussian fit function within scipy.
            g_mean = np.average(bin_centers, weights=data)
            g_sigma = np.sqrt(np.average((bin_centers - g_mean)**2, weights=data))
            g_norm = data * np.diff(bin_edges)
            initial_guess = [g_mean, g_sigma, np.sum(g_norm)]
            popt, pcov = curve_fit(self.gaussian, bin_centers[0], data[0], p0=initial_guess)
            
            # Label with estimated parameters and +/- 1 sigma
            g_label = f'Gaussian Fit\n'
            g_label += f'$\\mu$={popt[0]:.2f}$\\pm${np.sqrt(pcov[0,0]):.2f}\n'
            g_label += f'$\\sigma$={popt[1]:.2f}$\\pm${np.sqrt(pcov[1,1]):.2f}'

            # Plot the fit
            x = np.linspace(*range, 1000)
            ax.plot(x, self.gaussian(x, *popt), 'r-', label=g_label)

    @staticmethod
    def crystal_ball(x, alpha, n, mean, sigma, N):
        """
        Crystal Ball probability density function.

        Parameters:
        -----------
        x : array-like
            Input values where the function is evaluated.
        alpha : float
            Parameter defining the point where the Gaussian core transitions to the power-law tail.
            (alpha > 0 typically)
        n : float
            Parameter that defines the power of the tail.
        mean : float
            Mean (peak location) of the Gaussian core.
        sigma : float
            Standard deviation (width) of the Gaussian core.
        N : float
            Normalization factor (area under the curve).

        Returns:
        --------
        array-like
            Function values evaluated at x.
        """
        x = np.array(x)
        z = (x - mean) / sigma
        
        # Define constants
        abs_alpha = np.abs(alpha)

        # Guard against invalid values
        if abs_alpha == 0 or n <= 0 or sigma <= 0:
            return np.full_like(x, np.nan)
        
        # Normalization factor
        A = (n / abs_alpha)**n * np.exp(-abs_alpha**2 / 2)
        B = n / abs_alpha - abs_alpha

        with np.errstate(divide='ignore', invalid='ignore'):
            # Gaussian core
            gaussian_core = np.exp(-z**2 / 2)

            # Power-law tail
            tail_argument = B - z
            tail = np.where(tail_argument > 0, A * tail_argument ** (-n), 0)

        # Piecewise function
        condition = z > -abs_alpha
        result = np.where(condition, gaussian_core, tail)

        # Avoid NaNs in the output explicitly
        result = np.nan_to_num(result, nan=0.0)

        # Normalize to N
        return N * result

    @staticmethod
    def gaussian(x, mean, sigma, N):
        """
        Gaussian probability density function.

        Parameters:
        -----------
        x : array_like
            Independent variable.
        mean : float
            Mean (peak location) of the Gaussian.
        sigma : float
            Standard deviation (width) of the Gaussian.
        N : float
            Normalization factor (area under the curve).

        Returns:
        --------
        array_like
            Gaussian function evaluated at x.
        """
        prefactor = N / (sigma * np.sqrt(2 * np.pi))
        exponent = -0.5 * ((x - mean) / sigma)**2
        return np.exp(exponent) * N / (sigma * np.sqrt(2 * np.pi))