import numpy as np
from scipy.stats import binom
import pandas as pd
import matplotlib.pyplot as plt

from artists import SpineArtist
from style import Style
from variable import Variable

class SpineEfficiency(SpineArtist):
    """
    A class designed to encapsulate the calculation of the efficiency
    of the selection as a function of the specified variable. The method
    employed is a Bayesian approach, where the true efficiency is assumed
    to be a binomial random variable with a uniform prior. The posterior
    distribution is then calculated using the binomial likelihood and the
    prior. 

    Attributes
    ----------
    variable : Variable
        The variable to calculate the efficiency with respect to.
    samples : list
        A list of samples to use in the calculation of the efficiency.
    """
    def __init__(self, variable, categories, cuts):
        """
        Parameters
        ----------
        variable : Variable
            The variable to calculate the efficiency with respect to.
        categories : dict
            A dictionary mapping the category key to the category name.
        cuts : dict
            A dictionary mapping the cut key to the cut label.
        """
        self._variable = variable
        self._samples = list()
        self._categories = categories
        self._cuts = cuts
        self._posteriors = dict()

    def draw(self, ax, show_option, show_seqeff=True, show_unseqeff=True, groups=['All',], style=None):
        """
        Draw the artist on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        show_option : str, optional
            The option to use when showing the artist. The default is
            None. This is intended to be used in cases where the artist
            can be shown in different ways (e.g. 2D vs projection of 2D
            down to 1D).
        show_seqeff : bool, optional
            A flag to indicate if the sequential (cumulative)
            efficiency should be shown. The default is True.
        show_unseqeff : bool, optional
            A flag to indicate if the unsequential (single-cut)
            efficiency should be shown. The default is True.
        groups : list, optional
            A list of groups to show in the plot. This configures the
            the artist to show only the specified groups. The default
            is ['All',].
        style : Style, optional
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.

        Returns
        -------
        None.
        """
        # Reduce the efficiency calculation to the final result
        final_posteriors = self.reduce()

        # Lambda function to extract the cv, msigma, and psigma
        # values from the final_posteriors dictionary for a given
        # key.
        get_values = lambda key : [final_posteriors[i][key] for i in [1,2,3]]

        # Lambda formatter to round the values to two decimal
        # places, display as percentage, and add super and 
        # subscripts for the error values.
        formatter = lambda x,y,z: f'${100*x:.2f}^{{\ +{100*y:.2f}}}_{{\ -{100*z:.2f}}}$'

        if show_option == 'table':

            # Clear up the axis because we are going to draw a table
            # on it (no need for any other plot elements).
            ax.axis('tight')
            ax.axis('off')

            # Create the table data.
            results = pd.DataFrame({r'   ': list(), r'Cut': list(), r'Efficiency [%]': list(), r'Cumulative [%]': list()})
            group_endpoint = dict()

            # Loop over the groups requested in the plot.
            for group in groups:
                if group not in final_posteriors[0].keys():
                    continue

                # Extract the cv, msigma, and psigma values for the
                # group.
                cv, msigma, psigma = [final_posteriors[i][group] for i in [1,2,3]]

                seq = lambda x : [v for k, v in x.items() if 'seq_' in k and 'unseq' not in k]
                unseq = lambda x : [v for k, v in x.items() if 'unseq_' in k]

                entry = {   r'   ': [group,] + [r'' for _ in range(1, len(self._cuts))],
                            r'Cut': self._cuts.values(),
                            'Differential\nEfficiency [%]': [formatter(x,y,z) for x,y,z in zip(unseq(cv), unseq(msigma), unseq(psigma))],
                            'Cumulative\nEfficiency [%]': [formatter(x,y,z) for x,y,z in zip(seq(cv), seq(msigma), seq(psigma))] }
                results = pd.concat([results, pd.DataFrame(entry)])
                group_endpoint[group] = len(results)

            # Trim the table to only show the columns that are
            # requested by the user.
            cols = [r'   ', r'Cut',] if len(groups) > 1 else [r'Cut',]
            if show_unseqeff:
                cols.append('Differential\nEfficiency [%]')
            if show_seqeff:
                cols.append('Cumulative\nEfficiency [%]')
            results = results[cols]

            # Rename "Cumulative" to "Efficiency" if it is the only
            # efficiency to show.
            if show_seqeff and not show_unseqeff:
                results.rename(columns={r'Cumulative Efficiency [%]': r'Efficiency [%]'}, inplace=True)

            table_data = [results.columns.to_list()] + results.values.tolist()
            table = ax.table(cellText=table_data, colLabels=None, loc='center', cellLoc='center', edges='T')
            table.scale(1, 2.75)
            for i in range(2, len(table_data)):
                if i == len(table_data) - 1:
                    for j in range(len(table_data[i])):
                        table[i, j].visible_edges = 'B'
                        #table[i, j].set_height(0.1)
                else:
                    for j in range(len(table_data[i])):
                        table[i, j].visible_edges = 'open'
                    if i in group_endpoint.values():
                        table[i, 0].visible_edges = 'B'


    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the efficiency calculation. The calculation of 
        the efficiency uses a Bayesian approach, where the true
        efficiency is assumed to be a binomial random variable with a
        uniform prior. The posterior distribution is then calculated
        using the binomial likelihood and the prior. The nature of this
        calculation means that it can be performed on a sample-by-sample
        (and category-by-category) basis, which is why the sample is
        processed individually.

        Parameters
        ----------
        sample : Sample
            The sample to add to the efficiency calculation.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)
        self._samples.append(sample)

        # Calculate the efficiency for the sample
        self.calculate(sample)

    def bin(self):
        """
        Bin the variable for each sample.

        Returns
        -------
        None.
        """
        range_mask = (self.df[bin_var] > bin_range[0]) & (self.df[bin_var] < bin_range[1])
        bins = np.percentile(self.df[bin_var][range_mask], np.linspace(0, 100, nbins))
        bin_centers = (bins[1:] + bins[:-1])/2
        bin_error = np.diff(bins)/2
        bin_assignment = np.digitize(self.df[bin_var], bins)

        success = [np.sum(self.df[cutname][bin_assignment == i]) for i in range(1,nbins)]
        total = [np.sum(bin_assignment == i) for i in range(1,nbins)]
        return success, total, bin_centers, bin_error

    def calculate(self, sample, significance=0.6827):
        """
        Calculate the efficiency of the selection on the sample as a
        function of the variable of interest. The method employed is a
        Bayesian approach, where the true efficiency is assumed to be a
        binomial random variable with a uniform prior. The posterior
        distribution is then calculated using the binomial likelihood
        and the prior. The nature of this calculation means that it can
        be performed on a sample-by-sample (and category-by-category)
        basis, which is why each category within the sample is processed
        individually.

        Parameters
        ----------
        sample : Sample
            The sample to calculate the efficiency for.
        significance : float, optional
            The significance level to use when calculating the
            efficiency. The default is 0.6827.

        Returns
        -------
        None.
        """
        efficiencies = np.linspace(0.0, 1, int(1e6))
        data, _ = sample.get_data([self._variable._key, *self._cuts.keys()])
        #print(data)
        for category, values in data.items():
            if category not in self._categories:
                continue
            if category not in self._posteriors:
                self._posteriors[category] = {f'seq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()}
                self._posteriors[category].update({f'unseq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()})
            total = len(values[0])
            for ci, (cut, cutname) in enumerate(self._cuts.items()):
                # Sequential cuts
                success = np.sum(np.all(values[1:ci+2], axis=0))
                self._posteriors[category][f'seq_{cut}'] *= binom.pmf(success, total, efficiencies)
                self._posteriors[category][f'seq_{cut}'] /= np.sum(self._posteriors[category][f'seq_{cut}'])

                # Non-sequential cuts
                success = np.sum(values[ci+1].to_numpy(dtype=np.bool))
                self._posteriors[category][f'unseq_{cut}'] *= binom.pmf(success, total, efficiencies)
                self._posteriors[category][f'unseq_{cut}'] /= np.sum(self._posteriors[category][f'unseq_{cut}'])

    @staticmethod
    def multiply_posteriors(pos0, pos1, group0, group1):
        """
        Update the first posterior with the product of the two
        posterior distributions. This is a simple process, but there
        are some subtleties to consider. First and foremost, the
        resulting posterior is not, in general, normalized. It is also
        important to note that the two posteriors may not overlap
        within the precision of the calculation. This happens when two
        groups / categories prefer very different values of the true
        efficiency. In this case, it does not make much sense to
        combine the two posteriors. Instead, the function raises an
        exception to indicate that the two posteriors are not
        compatible.

        Parameters
        ----------
        pos0 : np.ndarray
            The first posterior distribution.
        pos1 : np.ndarray
            The second posterior distribution.
        group0 : str
            The name of the first group, which is used to identify the
            posterior distribution if an exception is raised.
        group1 : str
            The name of the second group, which is used to identify the
            posterior distribution if an exception is raised.
        
        Returns
        -------
        pos : np.ndarray
            The product of the two posterior distributions.

        Raises
        ------
        ValueError
            If the two posteriors are not compatible.
        """
        pos = pos0*pos1
        if np.sum(pos) == 0:
            raise ValueError(f"Posteriors for groups '{group0}' and '{group1}' are non-overlapping within machine-precision.")
        pos /= np.sum(pos)
        return pos

    def reduce(self, significance=0.6827):
        """
        Reduce the efficiency calculation to the final result.

        Parameters
        ----------
        significance : float, optional
            The significance level to use when calculating the
            efficiency. The default is 0.6827.

        Returns
        -------
        final_posteriors : dict
            A dictionary containing the final efficiency posteriors.
        cv : dict
            A dictionary containing the central value of the efficiency.
        msigma : dict
            A dictionary containing the negative error on the efficiency.
        psigma : dict
            A dictionary containing the positive error on the efficiency.
        """
        final_posteriors = {}

        # Reduce by grouped categories
        for category in self._posteriors.keys():
            if self._categories[category] not in final_posteriors.keys():
                final_posteriors[self._categories[category]] = {f'seq_{c}' : np.ones(self._posteriors[category][f'seq_{c}'].shape) for c in self._cuts.keys()}
                final_posteriors[self._categories[category]].update({f'unseq_{c}' : np.ones(self._posteriors[category][f'unseq_{c}'].shape) for c in self._cuts.keys()})
            for cut in self._cuts.keys():
                final_posteriors[self._categories[category]][f'seq_{cut}'] = SpineEfficiency.multiply_posteriors(final_posteriors[self._categories[category]][f'seq_{cut}'], self._posteriors[category][f'seq_{cut}'], self._categories[category], category)
                final_posteriors[self._categories[category]][f'unseq_{cut}'] = SpineEfficiency.multiply_posteriors(final_posteriors[self._categories[category]][f'unseq_{cut}'], self._posteriors[category][f'unseq_{cut}'], self._categories[category], category)

        # Reduce across groups
        final_posteriors['All'] = {f'seq_{c}' : np.ones(list(final_posteriors.values())[0][f'seq_{c}'].shape) for c in self._cuts.keys()}
        final_posteriors['All'].update({f'unseq_{c}' : np.ones(list(final_posteriors.values())[0][f'unseq_{c}'].shape) for c in self._cuts.keys()})
        for category in final_posteriors.keys():
            for cut in self._cuts.keys():
                try:
                    final_posteriors['All'][f'seq_{cut}'] = SpineEfficiency.multiply_posteriors(final_posteriors['All'][f'seq_{cut}'], final_posteriors[category][f'seq_{cut}'], 'All', category)
                    final_posteriors['All'][f'unseq_{cut}'] = SpineEfficiency.multiply_posteriors(final_posteriors['All'][f'unseq_{cut}'], final_posteriors[category][f'unseq_{cut}'], 'All', category)
                except ValueError as e:
                    print(e)
                
        cv = {}
        msigma = {}
        psigma = {}
        efficiencies = np.linspace(0.0, 1, int(1e6))
        for category in final_posteriors.keys():
            cv[category] = {c : efficiencies[int(np.argmax(p))] for c, p in final_posteriors[category].items()}
            sig = [0.5-significance/2.0, 0.5+significance/2.0]
            cumulative = {c : np.cumsum(p) for c, p in final_posteriors[category].items()}
            msigma[category] = {c : cv[category][c]-efficiencies[int(np.argmax(cumulative[c] > sig[0]))] for c in final_posteriors[category].keys()}
            psigma[category] = {c : efficiencies[int(np.argmax(cumulative[c] > sig[1]))]-cv[category][c] for c in final_posteriors[category].keys()}
        
        return final_posteriors, cv, msigma, psigma