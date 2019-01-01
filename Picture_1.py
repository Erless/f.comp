#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt


"""
method that creates a graphic space were the data, given as the sample, is represented in a figure.
Sample:
{
    'items': [
    {
        'x': [1, 2, 3, 4],
        'y': [4, 7, 3, 9],
        'legend': 'caption',
        'linestyle': ':'
    },
    {
        'x': [1, 2, 3, 4],
        'y': [5, 5, 5, 5],
        'legend': caption 2',
        'linestyle': '-'
    }
    ],
    'title': 'Title of the graph',
    'xlabel': 'Tiempo',
    'ylabel': 'Posicion',
}
"""
class Picture:

    def __init__(self, params=None):
        """
        Here we create an empty space, with a grid.
        """
        plt.figure()
        plt.grid()

        """
        Now the class gets all the other attributes other than items, which are the title and the label for the axis.
        If they are not specified, they are set with predefined names (title, xlabel, ylabel respectively.
        Lastly it adds the parameters in 'items', which are the data to plot and some additional information about it.)
        """
        if params.get('title'):
            plt.title(params['title'])
        
        if params.get('xlabel'):
            plt.xlabel(params['xlabel'])
            
        if params.get('ylabel'):
            plt.ylabel(params['ylabel'])
            
        for item in params.get('items', []):
            self.add_item(item)
                
    
    def add_item(self, item):
        """
        This method adds the points to plot to the figure and sets the caption and the linestyle to default values
        (empty for the caption and '-' for the linestyle).
        """
        plt.plot(
                item['x'],
                item['y'],
                label='' if not item.get('legend') else item['legend'],
                linestyle='-' if not item.get('linestyle') else item['linestyle']
            )
    
    def show_plot(self):
        """
        Plots the data added to the figure and makes the legend visible.
        """
        plt.legend()
        plt.show()