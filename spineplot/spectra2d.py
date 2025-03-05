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
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
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
    def __init__(self, variables, categories, colors, category_types, title=None) -> None:
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
        title : str, optional
            The title of the artist. This will be placed at the top of
            the axis assigned to the artist. The default is None.

        Returns
        -------
        None.
        """
        super().__init__(variables, categories, colors, title)
        self._category_types = category_types
        self._plotdata_diagonal = None
        self._binedges_diagonal = None
        self._plotdata_absdiag = None
        self._binedges_absdiag = None

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
        if self._plotdata_absdiag is None:
            self._plotdata_absdiag = {}
            self._binedges_absdiag = {}
        

        # Check if a mask is present for the variables. If so, we need
        # to combine the masks for the two variables.
        if self._variables[0].mask is not None and self._variables[1].mask is not None:
            joint_mask = f'{self._variables[0].mask} and {self._variables[1].mask}'
        elif self._variables[0].mask is not None:
            joint_mask = self._variables[0].mask
        elif self._variables[1].mask is not None:
            joint_mask = self._variables[1].mask
        else:
            joint_mask = None

        data, weights = sample.get_data([self._variables[0]._key, self._variables[1]._key], joint_mask)        
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

            if self._categories[category] not in self._plotdata_absdiag:
                self._plotdata_absdiag[self._categories[category]] = np.zeros(self._variables[0]._nbins)
            absdiag = values[1] - values[0]
            h = np.histogram(absdiag, bins=self._variables[0]._nbins, range=(-100,100), weights=weights[category])
            self._plotdata_absdiag[self._categories[category]] += h[0]
            self._binedges_absdiag[self._categories[category]] = h[1]

    def draw(self, ax, style, show_option='2d', draw_identity=True,
             override_xlabel=None, invert_stack_order=False, fit_type=None,
             logx=False, logy=False) -> None:
        """
        Plots the data for the SpineSpectra2D object.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object to draw the plot on.
        style : Style
            The Style object to use for the plot.
        show_option : str
            The option to use for the plot. This can be one of a few
            options (default is '2d'):
                '2d'         - Draw a 2D histogram of the data.
                'projection' - Draw a projection of the data about the
                               diagonal
                'absdiff'    - Draw a absolute difference of the data
        draw_identity : bool
            A flag to indicate if the identity line should be drawn on
            the plot. The default is True.
        override_xlabel : str
            An optional override for the x-axis label.
        invert_stack_order : bool
            A flag to indicate if the stack order in the legend should
            be inverted. The default is False.
        fit_type : str
            The type of fit to perform on the data. The default is
            None, which will not perform any fit. The options are:
                'crystal_ball' - Perform a Crystal Ball fit on the data.
                'gaussian'     - Perform a Gaussian fit on the data.
        logx : bool
            A flag to indicate if the x-axis should be logarithmic.
            The default is False.
        logy : bool
            A flag to indicate if the y-axis should be logarithmic.
            The default is False.
        
        Returns
        -------
        None.
        """
        ax.set_title(self._title)
        
        if show_option == '2d' and self._plotdata is not None:
            values = np.sum([v for v in self._plotdata.values()], axis=0)
            binedges = self._binedges[list(self._plotdata.keys())[0]]
            ax.imshow(values.T, extent=(binedges[0], binedges[-1], binedges[0], binedges[-1]), aspect='auto', origin='lower')
            ax.set_xlabel(self._variables[0]._xlabel if override_xlabel is None else override_xlabel)
            ax.set_ylabel(self._variables[1]._xlabel)
            
            # Draw the identity line. This must span the full range
            # of the plot, so we need to find the minimum and maximum
            # of the range for the plot.
            if draw_identity:
                min_range = min([binedges[0], binedges[-1]])
                max_range = max([binedges[0], binedges[-1]])
                ax.plot([min_range, max_range], [min_range, max_range], 'k--')

        if show_option == 'projection' and self._plotdata_diagonal is not None:
            labels, data = zip(*self._plotdata_diagonal.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges_diagonal[l][:-1] + np.diff(self._binedges_diagonal[l]) / 2 for l in labels]

            ax.hist(bincenters, weights=data, bins=self._variables[0]._nbins, range=(-1,1), histtype='barstacked', label=labels, color=colors, stacked=True)
            ax.set_xlabel('(Y-X)/X' if override_xlabel is None else override_xlabel)
            ax.set_ylabel('Entries')

            if fit_type is not None:
                super().fit_with_function(ax, bincenters[0], np.sum(data, axis=0), self._binedges_diagonal[labels[0]], fit_type)

            if invert_stack_order:
                h, l = ax.get_legend_handles_labels()
                ax.legend(h[::-1], l[::-1])
            else:
                ax.legend()

        if show_option == 'absdiff' and self._plotdata_absdiag is not None:
            labels, data = zip(*self._plotdata_absdiag.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges_absdiag[l][:-1] + np.diff(self._binedges_absdiag[l]) / 2 for l in labels]

            ax.hist(bincenters, weights=data, bins=self._variables[0]._nbins, range=(-100,100), histtype='barstacked', label=labels, color=colors, stacked=True)
            ax.set_xlabel('(Y-X)' if override_xlabel is None else override_xlabel)
            ax.set_ylabel('Entries')

            if invert_stack_order:
                h, l = ax.get_legend_handles_labels()
                ax.legend(h[::-1], l[::-1])
            else:
                ax.legend()
        
        if style.mark_pot:
            self.mark_pot(ax)
        if style.mark_preliminary is not None:
            self.mark_preliminary(ax, style.mark_preliminary)

        # Set the axis to be logarithmic if requested.
        if logx:
            # Modify the x-axis limits to ensure that the lower limit
            # is greater than zero. The lower edge needs to be at least
            # 3 orders of magnitude less than the maximum value in the
            # plot.
            if self._variable._range[0] == 0:    
                xhigh_exporder = np.floor(np.log10(self._variable._range[1]))
                xlow = xhigh_exporder - 3
                ax.set_xlim(10**xlow, self._variable._range[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
