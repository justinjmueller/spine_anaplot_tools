import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectra import SpineSpectra
from style import Style
from variable import Variable

class SpineSpectra2D(SpineSpectra):
    """
    A class designed to encapsulate a pair of variables' spectrum for
    an ensemble of samples. The class method add_sample() can be used
    to add a sample to the SpineSpectra2D.

    Attributes
    ----------
    _variables : list
        The list of Variable objects for the spectrum.
    _categories : dict
        A dictionary of the categories for the spectrum. This serves as
        a map between the category label in the input TTree and the
        category label for the spectrum (and therefore what is shown
        in a single legend entry).
    _colors : dict
        A dictionary of the colors for the categories in the spectrum.
        This serves as a map between the category label for the
        spectrum (value in the `_categories` dictionary) and the color
        to use for the histogram. The color can be any valid matplotlib
        color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
    _plotdata : dict
        A dictionary of the data for the spectrum. This is a map between
        the category label for the spectrum and the histogram data for
        that category.
    """
    def __init__(self, variables, categories, colors, category_types) -> None:
        """
        Initializes the SpineSpectra2D object.

        Parameters
        ----------
        variables : list
            The list of Variable objects for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This serves
            as a map between the category label in the input TTree and
            the category label for the spectrum (and therefore what is
            shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the spectrum.
            This serves as a map between the category label for the
            spectrum (value in the `_categories` dictionary) and the color
            to use for the histogram. The color can be any valid matplotlib
            color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the spectrum.
            This serves as a map between the category label for the spectrum
            (value in the `_categories` dictionary) and the type of plot to
            use for the histogram. The type should be either 'histogram' or
            'scatter' to correspond to a stacked histogram or scatter plot,
            respectively.

        Returns
        -------
        None.
        """
        super().__init__(variables, categories, colors)
        self._category_types = category_types
        self._plotdata_diagonal = None
        self._binedges_diagonal = None

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra2D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra2D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
        if self._plotdata_diagonal is None:
            self._plotdata_diagonal = {}
            self._binedges_diagonal = {}

        data, weights = sample.get_data([self._variables[0]._key, self._variables[1]._key])        
        for category, values in data.items():
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros((self._variables[0]._nbins, self._variables[1]._nbins))
            h = np.histogram2d(values[0], values[1], bins=(self._variables[0]._nbins, self._variables[1]._nbins), range=(self._variables[0]._range, self._variables[1]._range), weights=weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

            if self._categories[category] not in self._plotdata_diagonal:
                self._plotdata_diagonal[self._categories[category]] = np.zeros(self._variables[0]._nbins)
            diag = np.divide(values[1] - values[0], values[0])
            h = np.histogram(diag, bins=self._variables[0]._nbins, range=(-1,1), weights=weights[category])
            self._plotdata_diagonal[self._categories[category]] += h[0]
            self._binedges_diagonal[self._categories[category]] = h[1]

    def draw(self, ax, style, show_option='2d') -> None:
        """
        Plots the data for the SpineSpectra2D object.

        Parameters
        ----------
        style : Style
            The Style object to use for the plot.
        path : str
            The path to the output directory
        name : str
            The name of the output image file.
        
        Returns
        -------
        None.
        """
        if show_option == '2d' and self._plotdata is not None:
            values = np.sum([v for v in self._plotdata.values()], axis=0)
            binedges = self._binedges[list(self._plotdata.keys())[0]]
            ax.imshow(values.T, extent=(binedges[0], binedges[-1], binedges[0], binedges[-1]), aspect='auto', origin='lower')
            ax.set_xlabel(self._variables[0]._xlabel)
            ax.set_ylabel(self._variables[1]._xlabel)

        if show_option == 'projection' and self._plotdata_diagonal is not None:
            labels, data = zip(*self._plotdata_diagonal.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges_diagonal[l][:-1] + np.diff(self._binedges_diagonal[l]) / 2 for l in labels]

            ax.hist(bincenters, weights=data, bins=self._variables[0]._nbins, range=(-1,1), histtype='barstacked', label=labels, color=colors, stacked=True)
            ax.set_xlabel('(Y-X)/X')
            ax.set_ylabel('Entries')
        
        if style.get_mark_pot():
            self.mark_pot(ax)
        if style.get_mark_preliminary() is not None:
            self.mark_preliminary(ax, style.get_mark_preliminary())
        if style.get_title() is not None:
            ax.set_title(style.get_title())