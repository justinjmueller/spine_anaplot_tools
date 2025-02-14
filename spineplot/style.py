import matplotlib.pyplot as plt

class Style:
    """
    A class designed to encapsulate the configuration of a plotting
    style for a given plot.

    Attributes
    ----------
    _style : str
        The name of the style sheet to use for the plot.
    _default_figsize : tuple
        The default size of the figure to create.
    _title : str
        The title to place at the top of the plot.
    _mark_pot : bool
        A flag toggling the display of the total POT (exposure) at the
        top of the plot above the axis and below the title.
    _mark_preliminary : str
        A string to be used a label to indicate that the plot is
        preliminary. If None, no label is added.
    _show_component_number : bool
        A flag toggling the display of the total number of selected
        interactions per component in the legend.
    _show_component_percentage : bool
        A flag toggling the display of the percentage of selected
        interactions per component in the legend.
    _invert_stack_order : bool
        A flag toggling the inversion of the stack order for the
        components in the histogram.
    _plot_kwargs : dict
        A dictionary containing the keyword arguments to be passed
        to the plotting function.
    """
    def __init__(self, style_sheet, default_figsize, title, mark_pot, mark_preliminary, show_component_number, show_component_percentage, invert_stack_order, plot_kwargs) -> None:
        """
        Initializes the Style object with the given kwargs.

        Parameters
        ----------
        style_sheet : str
            The name of the style sheet to use for the plot.
        default_figsize : tuple
            The default size of the figure to create.
        title : str
            The title to place at the top of the plot.
        mark_pot : bool
            A flag toggling the display of the total POT (exposure) at the
            top of the plot above the axis and below the title.
        mark_preliminary : str
            A string to be used a label to indicate that the plot is
            preliminary. If None, no label is added.
        show_component_number : bool
            A flag toggling the display of the total number of selected
            interactions per component in the legend.
        show_component_percentage : bool
            A flag toggling the display of the percentage of selected
            interactions per component in the legend.
        invert_stack_order : bool
            A flag toggling the inversion of the stack order for the
            components in the histogram.
        plot_kwargs : dict
            A dictionary containing the keyword arguments to be passed
            to the plotting function.

        Returns
        -------
        None
        """
        self._style = style_sheet
        self._default_figsize = default_figsize
        self._title = None if title == 'none' else title
        self._mark_pot = mark_pot
        self._mark_preliminary = None if mark_preliminary == 'none' else mark_preliminary
        self._show_component_number = show_component_number
        self._show_component_percentage = show_component_percentage
        self._invert_stack_order = invert_stack_order
        self._plot_kwargs = plot_kwargs

    def __enter__(self):
        """
        Sets the style for the plot.

        Returns
        -------
        self : Style
            The created Style object.
        """
        plt.style.use(self._style)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Resets the style for the plot.

        Parameters
        ----------
        exc_type : type
            The exception type.
        exc_val : Exception
            The exception value.
        exc_tb : traceback
            The traceback.

        Returns
        -------
        None
        """
        plt.style.use('default')

    def get_style(self) -> str:
        """
        Returns the value of the style attribute.

        Returns
        -------
        str
            The value of the style attribute.
        """
        return self._style
    
    @property
    def default_figsize(self):
        """
        Returns the value of the default_figsize attribute.

        Returns
        -------
        tuple
            The value of the default_figsize attribute.
        """
        return self._default_figsize

    def get_title(self) -> str:
        """
        Returns the value of the title attribute.

        Returns
        -------
        str
            The value of the title attribute.
        """
        return self._title
    
    def get_mark_pot(self) -> bool:
        """
        Returns the value of the mark_pot attribute.

        Returns
        -------
        bool
            The value of the mark_pot attribute.
        """
        return self._mark_pot
    
    def get_mark_preliminary(self) -> str:
        """
        Returns the value of the mark_preliminary attribute.

        Returns
        -------
        str
            The value of the mark_preliminary attribute.
        """
        return self._mark_preliminary

    def get_show_component_number(self) -> bool:
        """
        Returns the value of the show_component_number attribute.

        Returns
        -------
        bool
            The value of the show_component_number attribute.
        """
        return self._show_component_number
    
    def get_show_component_percentage(self) -> bool:
        """
        Returns the value of the show_component_percentage attribute.

        Returns
        -------
        bool
            The value of the show_component_percentage attribute.
        """
        return self._show_component_percentage
    
    def get_invert_stack_order(self) -> bool:
        """
        Returns the value of the invert_stack_order attribute.

        Returns
        -------
        bool
            The value of the invert_stack_order attribute.
        """
        return self._invert_stack_order

    @property
    def plot_kwargs(self) -> dict:
        """
        Returns the value of the plot_kwargs attribute.

        Returns
        -------
        dict
            The value of the plot_kwargs attribute.
        """
        return self._plot_kwargs