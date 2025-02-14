import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectra import SpineSpectra
from style import Style
from variable import Variable

class SpineSpectra1D(SpineSpectra):
    """
    A class designed to encapsulate a single variable's spectrum for an
    ensemble of samples. This is a specialization of the SpineSpectra
    class that is intended to be used for 1D spectra. At its core, this
    is a simple histogram plot of a single variable.

    Attributes
    ----------
    _variable : Variable
        The Variable object for the spectrum.
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
    def __init__(self, variable, categories, colors, category_types) -> None:
        """
        Initializes the SpineSpectra1D.

        Parameters
        ----------
        variable : Variable
            The Variable object for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This
            serves as a map between the category label in the input
            TTree and the category label for the spectrum (and
            therefore what is shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the color to use for the histogram. The color can be
            any valid matplotlib color string or a cycle indicator
            (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the type of plot to use for the histogram. The type
            should be either 'histogram' or 'scatter' to correspond to
            a stacked histogram or scatter plot, respectively.

        Returns
        -------
        None.
        """
        super().__init__([variable,], categories, colors)
        self._variable = self._variables[0]
        self._category_types = category_types

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra1D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

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

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
            self._onebincount = {}
        data, weights = sample.get_data([self._variable._key,])
        for category, values in data.items():
            values = values[0]
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros(self._variable._nbins)
                self._onebincount[self._categories[category]] = 0
            h = np.histogram(values, bins=self._variable._nbins, range=self._variable._range, weights=weights[category])
            self._onebincount[self._categories[category]] += np.sum(weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

    def draw(self, ax, style) -> None:
        """
        Plots the data for the SpineSpectra1D object.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        style : Style
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.

        Returns
        -------
        None.
        """
        ax.set_xlabel(self._variable._xlabel)
        ax.set_ylabel('Candidates')
        ax.set_xlim(*self._variable._range)

        if self._plotdata is not None:
            labels, data = zip(*self._plotdata.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in labels]

            histogram_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'histogram']
            scatter_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'scatter']

            denominator = np.sum([self._onebincount[labels[i]] for i in histogram_mask])
            counts = [x for x in self._onebincount.values()]
            if style.get_show_component_number() and style.get_show_component_percentage():
                hlabel = lambda x : f'{np.sum(x):.1f}, {np.sum(x)/denominator:.2%}'
                slabel = lambda x : f'{np.sum(x):.1f}'
                labels = [f'{label} ({hlabel(d) if li in histogram_mask else slabel(d)})' for li, (label, d) in enumerate(zip(labels, counts))]
            elif style.get_show_component_number():
                labels = [f'{label} ({np.sum(d):.1f})' for label, d in zip(labels, counts)]
            elif style.get_show_component_percentage():
                labels = [f'{label} ({np.sum(d)/denominator:.2%})' if li in histogram_mask else label for li, (label, d) in enumerate(zip(labels, counts))]

            if style.get_invert_stack_order():
                reduce = lambda x : [x[i] for i in histogram_mask[::-1]]
            else:
                reduce = lambda x : [x[i] for i in histogram_mask]
            
            ax.hist(reduce(bincenters), weights=reduce(data), bins=self._variable._nbins, range=self._variable._range, label=reduce(labels), color=reduce(colors), **style.plot_kwargs)

            reduce = lambda x : [x[i] for i in scatter_mask]
            for i, label in enumerate(reduce(labels)):
                ax.errorbar(bincenters[scatter_mask[i]], data[scatter_mask[i]], yerr=np.sqrt(data[scatter_mask[i]]), fmt='o', label=label, color=colors[scatter_mask[i]])
        
        if style.get_invert_stack_order():
            h, l = ax.get_legend_handles_labels()
            ax.legend(h[::-1], l[::-1])
        else:
            ax.legend()
        if style.get_mark_pot():
            self.mark_pot(ax)
        if style.get_mark_preliminary() is not None:
            self.mark_preliminary(ax, style.get_mark_preliminary())
        if style.get_title() is not None:
            ax.set_title(style.get_title())