# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 11:01:55 2019

@author: UNC
"""

def plot():

    import math
    import numpy as np

    # set up the meshgrid
    # floors and ceilings are taken because the range functions need
    # whole numbers
    # the ranges are taken from the data to make sure the grid covers
    # the values of the data
    x_min = math.floor(X[:,0].min())
    x_max = math.ceil(X[:,0].max())
    y_min = math.floor(X[:,1].min())
    y_max = math.ceil(X[:,1].max())
    x_,y_ = np.meshgrid(range(x_min,x_max), range(y_min,y_max))


    import plotly.plotly as py
    import plotly.graph_objs as go

    # for random color generation
    np.random.seed(10)
    color = np.random.randn(len(y))

    # the hypothesis function to generate the plane
    hypothesis = linreg.intercept_ + x_*linreg.coef_[0]  + y_*linreg.coef_[1]

    trace1 = go.Scatter3d(
        x=X[:,0],
        y=X[:,1],
        z=y,
        mode='markers',
        name='Dummy Data',
        marker=dict(
            size=20,
            color=color,               
            colorscale='Rainbow',   
            opacity=0.9
        )
    )

    trace2 = go.Surface(
            x = x_,
            y = y_,
            z = hypothesis,
            name='Plane of Best Fit',
            opacity = 0.7,
            colorscale='Greys',
            showscale= False
        )

    data = [trace1, trace2]
    layout = go.Layout(
            title='3D Plane of Best Fit Through Generated Dummy Data',
            margin=dict(
                l=0,
                r=0,
                b=10,
                t=100  # the title is obscured if the top margin is not adjusted
        )
    )
    fig = go.Figure(data=data, layout=layout)
    return py.iplot(fig, filename='3d-scatter-with-plane')
