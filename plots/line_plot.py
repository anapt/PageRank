import plotly.plotly as py
import plotly.graph_objs as go
import plotly

x = []

trace = go.Scatter(
    y = x
)

data = [trace]
bandxaxis = go.XAxis(
    title="Iteration",
    showgrid=True,
    showline=True,
    ticks="",
    showticklabels=True,
    mirror=True,
    linewidth=2
)

layout = go.Layout(title = 'Error until convergence // 6012 nodes',
                xaxis = bandxaxis,
                yaxis = dict(title = 'Error')
              )

fig = go.Figure(data=data, layout = layout)
#plotly.offline.plot(fig, filename='Imperative-2')
#plotly.offline.plot(fig, filename='scatter-plot-png', image='svg')
plotly.offline.plot(fig, filename='scatter-plot-png', image='png')
