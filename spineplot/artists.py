from abc import ABC, abstractmethod
import matplotlib.pyplot as plt

class SpineArtist(ABC):
    """
    A class to represent an artist that will be added to a figure. This
    class serves as a base class for any artist that will be plotted
    with the package. All other artist classes will inherit from this
    class; this class will not be used directly and only contains
    abstract methods that must be implemented by any derived classes.

    """
    def __init__(self):
        pass

    @abstractmethod
    def draw(self, ax, style=None):
        """
        Draw the artist on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        style : Style, optional
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.

        Returns
        -------
        None.
        """
        pass

    @abstractmethod
    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the artist. Each artist is assumed to render a
        collection of samples, so this method is intended to enforce
        this assumption.

        Parameters
        ----------
        sample : Sample
            The sample to add to the artist.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        if is_ordinate:
            self._exposure_type = sample._exposure_type
            if self._exposure_type == 'pot':
                self._exposure = sample._exposure_pot
            else:
                self._exposure = sample._exposure_livetime