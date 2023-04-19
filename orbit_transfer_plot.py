import numpy as np
from astropy import units as u
from poliastro.bodies import *
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
import plotly.graph_objs as go
from plotly.subplots import make_subplots

class OrbitTransferPlot:
    def plot(self):
        fig = self._create_figure(zoom=False)
        return fig

    def plot_zoomed(self):
        fig = self._create_figure(zoom=True)
        return fig

    def _create_figure(self, zoom=False,attractor=Earth):
        R = np.linspace(2, 75, num=200)
        Rstar = [15.58, 40, 60, 100, 200, np.inf]

        hohmann_data = np.zeros_like(R)
        bielliptic_data = np.zeros((len(R), len(Rstar)))

        ss_i = Orbit.circular(attractor, 1.8 * u.km)
        r_i = ss_i.a
        v_i = np.sqrt(ss_i.v @ ss_i.v)
        for ii, r in enumerate(R):
            r_f = r * r_i
            man = Maneuver.hohmann(ss_i, r_f)
            hohmann_data[ii] = (man.get_total_cost() / v_i).decompose().value
            for jj, rstar in enumerate(Rstar):
                r_b = rstar * r_i
                man = Maneuver.bielliptic(ss_i, r_b, r_f)
                bielliptic_data[ii, jj] = (
                    (man.get_total_cost() / v_i).decompose().value
                )

        fig = make_subplots()

        fig.add_trace(go.Scatter(x=R, y=hohmann_data, name="Hohmann"))

        color_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
        legend_labels = [
            "Bielliptic (R* = 15.58)",
            "Bielliptic (R* = 40)",
            "Bielliptic (R* = 60)",
            "Bielliptic (R* = 100)",
            "Bielliptic (R* = 200)",
            "Bielliptic (R* = inf)",
        ]

        for jj in range(len(Rstar)):
            fig.add_trace(
                go.Scatter(
                    x=R,
                    y=bielliptic_data[:, jj],
                    name=legend_labels[jj],
                    line=dict(color=color_list[jj]),
                )
            )

        if zoom:
            x_range = [11.0, 16.0]
            y_range = [0.52, 0.545]
        else:
            x_range = [2, 75]
            y_range = [0.35, 0.6]

        fig.update_layout(
            title="Hohmann vs bielliptic transfers",
            xaxis_title="R",
            yaxis_title="Relative change in velocity",
            xaxis=dict(range=x_range),
            yaxis=dict(range=y_range),
            legend_title="Transfer Type",
        )

        return fig
