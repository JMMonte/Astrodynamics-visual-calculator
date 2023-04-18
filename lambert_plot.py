import numpy as np
import plotly.graph_objects as go
from poliastro.core import iod

class LambertPlot:
    '''Plot the Lambert problem in Python with Plotly.
    The Lambert problem is a special case of the inverse optimal control problem. This means that the problem can be solved by minimizing the cost function:
    .. math::
        T = \int_{0}^{1} \sqrt{1 + \dot{r}^2} \, \mathrm{d}t
    where :math:`r` is the position vector and :math:`\dot{r}` is the velocity vector. The solution to this problem is given by the Lambert's problem.
    The Lambert's problem is a special case of the inverse optimal control problem. This means that the problem can be solved by minimizing the cost function:
    .. math::
        T = \int_{0}^{1} \sqrt{1 + \dot{r}^2} \, \mathrm{d}t
    where :math:`r` is the position vector and :math:`\dot{r}` is the velocity vector. The solution to this problem is given by the Lambert's problem.
    
    In this plot, the Lambert's problem is solved for different values of the Lambert parameter :math:`\lambda` and the number of revolutions :math:`M`.
    The Lambert parameter is defined as:
    .. math::
        \lambda = \frac{r_f - r_i}{r_f + r_i}
    where :math:`r_i` and :math:`r_f` are the initial and final position vectors, respectively.
    The number of revolutions :math:`M` is the number of times the spacecraft crosses the plane of the orbit.

    References
    ----------
    p. 3, https://www.researchgate.net/publication/228988510_The_Lambert_Problem
    poliastro/src/poliastro/core/iod.py

    '''
    def __init__(self):
        self.M_list = 0, 1, 2, 3 # M = 0 is the special case
        self.ll_list = 1, 0.9, 0.7, 0, -0.7, -0.9, -1 # ll is the Lambert parameter and this is the range of values to plot
        self.x = np.linspace(-1, 2, num=1000)

    def create_plot(self):
        fig = go.Figure()

        line_styles = ["solid", "dash"]
        colors = ["white"] * len(self.ll_list)

        for M, linestyle in zip(self.M_list, line_styles):
            for ll, color in zip(self.ll_list, colors):
                T_x0 = np.zeros_like(self.x)
                for ii in range(len(self.x)):
                    y = iod._compute_y(self.x[ii], ll)
                    T_x0[ii] = iod._tof_equation_y(self.x[ii], y, 0.0, ll, M)
                if M == 0 and ll == 1:
                    T_x0[self.x > 0] = np.nan
                elif M > 0:
                    T_x0[self.x > 1] = np.nan
                fig.add_trace(go.Scatter(x=self.x, y=T_x0, mode="lines", line=dict(color=color, dash=linestyle), name=f"M={M}, λ={ll}"))

        fig.add_shape(type="line", x0=1, x1=1, y0=0, y1=10, yref="y", xref="x", line=dict(color="white", width=1))

        fig.update_xaxes(title_text="x (Normalized Position Vector)", tickvals=(-1, 0, 1, 2), ticktext=["-1", "0", "1", "2"])
        fig.update_yaxes(title_text="Time of Flight (Normalized)", tickvals=(0, np.pi, 2 * np.pi, 3 * np.pi), ticktext=["0", "1\pi", "2 \pi", "3 \pi"])

        fig.update_layout(title="Revisiting Lambert's Problem in Python with Plotly", yaxis_range=[0, 10])

        return fig

