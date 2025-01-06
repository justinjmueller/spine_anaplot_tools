import numpy as np
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