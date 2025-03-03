import numpy as np

class Variable:
    """
    A class designed to encapsulate the configuration of a single
    variable for the analysis.

    Attributes
    ----------
    _name : str
        The name of the variable.
    _key : str
        The key/name of the branch in the input TTree containing the
        variable data.
    _range : tuple
        The range of the variable.
    _nbins : int
        The number of bins for the variable.
    _xlabel : str
        The x-axis label for the variable.
    _bin_edges : numpy.ndarray
        The bin edges for the variable.
    _bin_centers : numpy.ndarray
        The bin centers for the variable.
    _bin_widths : numpy.ndarray
        The bin widths for the variable.
    """
    def __init__(self, name, key, range, nbins, xlabel) -> None:
        """
        Initializes the Variable object with the given kwargs.

        Parameters
        ----------
        name : str
            The name of the variable.
        key : str
            The key/name of the branch in the input TTree containing the
            variable data.
        range : tuple
            The range of the variable.
        nbins : int
            The number of bins for the variable.
        xlabel : str
            The x-axis label for the variable.

        Returns
        -------
        None
        """
        self._name = name
        self._key = key
        self._range = range
        self._nbins = nbins
        self._xlabel = xlabel
        self._bin_edges = np.linspace(self._range[0], self._range[1], self._nbins+1)
        self._bin_centers = 0.5*(self._bin_edges[1:] + self._bin_edges[:-1])
        self._bin_widths = self._bin_edges[1:] - self._bin_edges[:-1]