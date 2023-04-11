import streamlit as st
import poliastro
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.core.spheroid_location import cartesian_to_ellipsoidal
import plotly.graph_objects as go
import numpy as np
from astropy import units as u


class OrbitalManeuversCalculator:
    def __init__(self):
        self.R_EARTH = Earth.R.to(u.km).value
        self.MU_EARTH = Earth.k.to(u.km**3 / u.s**2).value

    @staticmethod
    def to_unit(value, output_unit):
        return value << output_unit

    @staticmethod
    def sample_orbit(orbit, num_points=100):
        time_delta = orbit.period / num_points
        positions = []

        for i in range(num_points):
            propagated_orbit = orbit.propagate(time_delta * i)
            pos = propagated_orbit.r
            positions.append(pos)

        return positions

    def get_orbit_parameters(self, title):
        st.sidebar.subheader(title)
        a = st.sidebar.number_input(f"{title} Semi-major axis (km)", min_value=self.R_EARTH, value=7000.0, max_value=self.R_EARTH * 50, step=1.0)
        e = st.sidebar.number_input(f"{title} Eccentricity", min_value=0.0, max_value=1.0, value=0.0, step=0.01)
        inc = st.sidebar.number_input(f"{title} Inclination (deg)", min_value=0.0, max_value=180.0, value=28.5, step=0.1)
        raan = st.sidebar.number_input(f"{title} Right ascension of the ascending node (deg)", min_value=0.0, max_value=360.0, value=0.0, step=0.1)
        argp = st.sidebar.number_input(f"{title} Argument of perigee (deg)", min_value=0.0, max_value=360.0, value=0.0, step=0.1)
        nu = st.sidebar.number_input(f"{title} True anomaly (deg)", min_value=0.0, max_value=360.0, value=0.0, step=0.1)

        return Orbit.from_classical(Earth, a << u.km, e * u.one, inc << u.deg, raan << u.deg, argp << u.deg, nu << u.deg)

    def calculate_maneuver(self, orb1, orb2):
        man = Maneuver.hohmann(orb1, orb2.a)
        orb_int, orb_fin = orb1.apply_maneuver(man, intermediate=True)

        return man, orb_int, orb_fin

    def display_maneuver_parameters(self, man):
        st.header("Maneuver parameters")
        st.write(f"Transfer time: {man.get_total_time().to(u.h)}")
        st.write(f"Delta-v at perigee: {man[0][1].to(u.km / u.s)}")
        st.write(f"Delta-v at apogee: {man[1][1].to(u.km / u.s)}")
        st.write(f"Total delta-v: {man.get_total_cost().to(u.km / u.s)}")

    def plot_orbit(self, orbit, label, color, fig_2d, fig_3d):
        positions = self.sample_orbit(orbit, num_points=100)
        x = [pos[0].to(u.km).value for pos in positions]
        y = [pos[1].to(u.km).value for pos in positions]
        z = [pos[2].to(u.km).value for pos in positions]

        ellipsoidal_coords = [cartesian_to_ellipsoidal(Earth.R.value, Earth.R_polar.value, pos_x, pos_y, pos_z) for pos_x, pos_y, pos_z in zip(x, y, z)]
        lat = [coords[0] for coords in ellipsoidal_coords]
        lon = [coords[1] for coords in ellipsoidal_coords]
        h = [coords[2] for coords in ellipsoidal_coords]

        fig_2d.add_trace(go.Scatter(x=lon,
                                    y=lat,
                                    mode="lines",
                                    name=label,
                                    line=dict(color=color)))

        fig_3d.add_trace(go.Scatter3d(x=x,
                                    y=y,
                                    z=z,
                                    mode="lines",
                                    name=label,
                                    line=dict(color=color)))

    def run(self):
        st.sidebar.title("Orbital Maneuvers Calculator")
        st.sidebar.markdown("Enter the orbital parameters of two orbits, and the system will calculate the resulting Hohmann maneuver, displaying its parameters, the delta-v required, and the ground tracks on Earth in 2D and 3D.")

        orb1 = self.get_orbit_parameters("Initial orbit")
        orb2 = self.get_orbit_parameters("Final orbit")

        man, orb_int, orb_fin = self.calculate_maneuver(orb1, orb2)

        self.display_maneuver_parameters(man)

        st.header("Orbits and ground tracks")
        fig_2d = go.Figure()
        fig_3d = go.Figure()

        fig_3d.add_trace(go.Surface(x=self.R_EARTH * np.cos(np.linspace(0, 2 * np.pi, 100)) * np.cos(np.linspace(-np.pi / 2, np.pi / 2, 100))[:, None],
                                    y=self.R_EARTH * np.sin(np.linspace(0, 2 * np.pi, 100)) * np.cos(np.linspace(-np.pi / 2, np.pi / 2, 100))[:, None],
                                    z=self.R_EARTH * np.sin(np.linspace(-np.pi / 2, np.pi / 2, 100))[:, None],
                                    colorscale="Blues",
                                    showscale=False))

        self.plot_orbit(orb1, "Initial orbit", "green", fig_2d, fig_3d)
        self.plot_orbit(orb_int, "Intermediate orbit", "orange", fig_2d, fig_3d)
        self.plot_orbit(orb_fin, "Final orbit", "red", fig_2d, fig_3d)

        fig_2d.update_layout(title="Ground tracks on Earth",
                             xaxis_title="Longitude (deg)",
                             yaxis_title="Latitude (deg)",
                             width=800,
                             height=600)

        fig_3d.update_layout(title="Orbits in space",
                             scene=dict(xaxis_title="x (km)",
                                        yaxis_title="y (km)",
                                        zaxis_title="z (km)"),
                             width=800,
                             height=600)

        st.plotly_chart(fig_2d)
        st.plotly_chart(fig_3d)


if __name__ == "__main__":
    calculator = OrbitalManeuversCalculator()
    calculator.run()