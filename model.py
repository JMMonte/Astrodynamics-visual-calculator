import datetime
import numpy as np
import plotly.graph_objs as go
import streamlit as st
from astropy.time import Time
from astropy import units as u
from plotly.subplots import make_subplots
from poliastro.bodies import *
from poliastro.maneuver import Maneuver
from poliastro.earth.plotting import GroundtrackPlotter
from poliastro.earth import EarthSatellite
from typing import List


def configure_page():
    """Configures the page layout and title."""
    st.set_page_config(
        page_title="Orbit maneuver calculator",
        page_icon="🚀",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def define_main_attractors():
    """Returns a dictionary of the main attractors in the solar system."""
    return {
        "Earth": Earth,
        "Moon": Moon,
        "Mars": Mars,
        "Jupiter": Jupiter,
        "Saturn": Saturn,
        "Uranus": Uranus,
        "Neptune": Neptune,
        "Pluto": Pluto,
        "Mercury": Mercury,
        "Venus": Venus,
        "Sun": Sun,
    }

def planet_texture(planet):
    """Returns the texture of a planet."""
    planet_textures = {
        "Mercury": "https://www.solarsystemscope.com/textures/download/8k_mercury.jpg",
        "Venus": "https://www.solarsystemscope.com/textures/download/2k_venus_atmosphere.jpg",
        "Earth": "https://www.solarsystemscope.com/textures/download/8k_earth_daymap.jpg",
        "Moon": "https://www.solarsystemscope.com/textures/download/8k_moon.jpg",
        "Mars": "https://www.solarsystemscope.com/textures/download/8k_mars.jpg",
        "Jupiter": "https://www.solarsystemscope.com/textures/download/8k_jupiter.jpg",
        "Saturn": "https://www.solarsystemscope.com/textures/download/8k_saturn.jpg",
        "Uranus": "https://www.solarsystemscope.com/textures/download/2k_uranus.jpg",
        "Neptune": "https://www.solarsystemscope.com/textures/download/2k_neptune.jpg",
        "Pluto": "https://www.solarsystemscope.com/textures/download/4k_ceres_fictional.jpg",
        "Sun": "https://www.solarsystemscope.com/textures/download/2k_sun.jpg"
    }
    # Match planet name with texture
    return planet_textures[planet]

def to_unit(value, unit):
    """Converts a value to a unit."""
    if unit == "seconds":
        return value * u.s
    elif unit == "minutes":
        return value * u.min
    elif unit == "hours":
        return value * u.h
    else:  # days
        return value * u.day


def plotly_orbit_plotter(orbit_list, attractor, intersections=None, labels=None):
    """Plots a list of orbits in 3D using plotly.
    Parameters:
    orbit_list: list of poliastro.twobody.orbit.Orbit
        List of orbits to plot
    attractor: poliastro.bodies.Body
        Main attractor of the orbits
    intersections: list of poliastro.twobody.orbit.Orbit
        List of orbits to plot as intersections
    labels: list of str
        List of labels for the orbits
    Returns:
    fig: plotly.graph_objects.Figure
    """
    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "scatter3d"}]])

    if labels is None:
        labels = ["Orbit"] * len(orbit_list)

    for orbit, label in zip(orbit_list, labels):
        r = orbit.sample().xyz.T
        x, y, z = r[:, 0].to(u.km).value, r[:, 1].to(u.km).value, r[:, 2].to(u.km).value
        # get color of orbit from orbit name
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                name=label,
            )
        )
    # Add attractor
    u_rad = u.km
    thetas = u.Quantity(np.linspace(0, 2 * np.pi, 100), u.rad)
    phis = u.Quantity(np.linspace(0, np.pi, 100), u.rad)
    radius_equatorial = attractor.R.to(u_rad).value
    radius_polar = attractor.R_polar.to(u_rad).value
    x_center = radius_equatorial * np.outer(np.cos(thetas), np.sin(phis))
    y_center = radius_equatorial * np.outer(np.sin(thetas), np.sin(phis))
    z_center = radius_polar * np.outer(np.ones_like(thetas), np.cos(phis))

    # add a streamlit dropdown to select the color scale

    fig.add_trace(
        go.Surface(
            x=x_center, y=y_center, z=z_center, colorscale="Viridis", showscale=False
        )
    )

    # Add intersections
    if intersections is not None:
        for intersection in intersections:
            x, y, z = intersection["point"].to(u.km).value
            scatter = go.Scatter3d(
                x=[x],
                y=[y],
                z=[z],
                mode="markers",
                marker=dict(size=8, color="purple", opacity=1),
            )
            fig.add_trace(scatter)

    fig.update_layout(scene=dict(aspectmode="data"))
    # set the height of the plot to 600px
    fig.update_layout(height=800)

    return fig

def get_orbit_parameters(mission_epoch, attractor, orbit_name, altitude=500.0, ecc=0.0, inclination=135.0, raan=0.0, argp=0.0, nu=0.0):
    """Returns the orbit parameters for a given orbit name."""
    # define any orbit based on streamlit inputs for all orbit parameters
    orbit_altitude = st.number_input(
        "Orbit altitude [km]",
        min_value=0.0,
        value=altitude,
        step=50.0,
        key=f"{orbit_name}_altitude",
        help="Altitude of the orbit above the Planet's surface (counted from equatorial radius).",
    )
    orbit_altitude = orbit_altitude + attractor.R.to(u.km).value
    orbit_ecc = st.number_input(
        "Orbit eccentricity",
        min_value=0.0,
        max_value=1.0,
        value=ecc,
        step=0.01,
        key=f"{orbit_name}_ecc",
        help="Eccentricity of the orbit is the ratio of the distance between the two foci of the ellipse and the distance between the center of the ellipse and one of the foci. For more information, see https://en.wikipedia.org/wiki/Eccentricity_(mathematics)",
    )
    orbit_inclination = st.number_input(
        "Orbit inclination [deg]",
        min_value=0.0,
        max_value=180.0,
        value=inclination,
        step=1.0,
        key=f"{orbit_name}_inclination",
        help="Inclination of the orbit with respect to the equatorial plane of the Planet. 0.0 for equatorial orbits, 90.0 for polar orbits.",
    )
    orbit_raan = st.number_input(
        "Orbit RAAN [deg]",
        min_value=0.0,
        max_value=360.0,
        value=raan,
        step=1.0,
        key=f"{orbit_name}_raan",
        help="The Right Ascension of the Ascending Node of the trajectory orbit is the angle between the ascending node and the vernal equinox. For more information, see https://en.wikipedia.org/wiki/Right_ascension_of_the_ascending_node.",
    )
    orbit_argp = st.number_input(
        "Orbit argument of perigee [deg]",
        min_value=0.0,
        max_value=360.0,
        value=argp,
        step=1.0,
        key=f"{orbit_name}_argp",
        help="The argument of perigee is the angle between the ascending node and the perigee. For more information, see https://en.wikipedia.org/wiki/Argument_of_perigee.",
    )
    orbit_nu = st.number_input(
        "Orbit true anomaly [deg]",
        min_value=0.0,
        max_value=360.0,
        value=nu,
        step=1.0,
        key=f"{orbit_name}_nu",
        help="The true anomaly is the angle between the perigee and the spacecraft. For more information, see https://en.wikipedia.org/wiki/True_anomaly.",
    )
    orbit_name = (
        attractor,
        orbit_altitude * u.km,
        orbit_ecc * u.one,
        orbit_inclination * u.deg,
        orbit_raan * u.deg,
        orbit_argp * u.deg,
        orbit_nu * u.deg,
        mission_epoch,
    )

    # return orbit
    return orbit_name


def get_mission_epoch():
    """Returns the mission epoch."""
    now = datetime.datetime.now()
    date = st.sidebar.date_input("Select the mission epoch", value=now)
    time = st.sidebar.time_input("Select the mission epoch time", value=now.time())
    date = datetime.datetime(
        date.year, date.month, date.day, time.hour, time.minute, time.second
    )
    return Time(date, scale="utc")

def define_maneuver_types():
    return {
        "Hohmann transfer": Maneuver.hohmann,
        "Bielliptic transfer": Maneuver.bielliptic,
        "Lambert transfer": Maneuver.lambert,
    }

def show_maneuver_data(maneuver_type):
    if maneuver_type == "Hohmann transfer":
        with st.sidebar.expander("Hohmann transfer equation"):
            st.latex(
                r"\Delta v = \sqrt{\frac{2\mu}{r_1} - \frac{\mu}{a}} - \sqrt{\frac{\mu}{r_1}}"
            )
        
            # explain the hohmann equation in markdown and each of it's variables
            st.markdown(
                r"""
            The Hohmann transfer equation is used to calculate the delta-v required to transfer from one orbit to another. The equation is given by:
            - $\Delta v$ is the change in velocity required to transfer from one orbit to another
            - $\mu$ is the gravitational parameter of the main attractor
            - $r_1$ is the radius of the initial orbit
            - $a$ is the semi-major axis of the target orbit
            """
            )
    elif maneuver_type == "Bielliptic transfer":
        with st.sidebar.expander("Bielliptic transfer equation"):
            st.latex(
                r"\Delta v = \sqrt{\frac{2\mu}{r_1} - \frac{\mu}{a}} - \sqrt{\frac{\mu}{r_1}} + \sqrt{\frac{2\mu}{r_1} - \frac{\mu}{r_b}} - \sqrt{\frac{\mu}{r_b}}"
            )
            # explain the bielliptic equation in markdown and each of it's variables
            st.markdown(
                r"""
            A bielliptic transfer is a two-burn transfer orbit that uses two burns to transfer from one orbit to another.
            The Bielliptic transfer equation is used to calculate the delta-v required to transfer from one orbit to another. The equation is given by:
            - $\Delta v$ is the change in velocity required to transfer from one orbit to another
            - $\mu$ is the gravitational parameter of the main attractor
            - $r_1$ is the radius of the initial orbit
            - $a$ is the semi-major axis of the target orbit
            - $r_b$ is the radius of the intermediate orbit
            """
            )
    elif maneuver_type == "Lambert transfer":
        # write the lambert maneuver equation in latex format
        with st.sidebar.expander("Lambert transfer equation"):
            st.latex(
                r"\Delta v = \sqrt{\frac{2\mu}{r_1} - \frac{\mu}{a}} - \sqrt{\frac{\mu}{r_1}} + \sqrt{\frac{2\mu}{r_1} - \frac{\mu}{r_b}} - \sqrt{\frac{\mu}{r_b}}"
            )
            # explain the lambert equation in markdown and each of it's variables
            st.markdown(
                r"""
            A Lambert transfer is a type of transfer orbit that uses a single burn to transfer from one orbit to another in a minimum amount of time.
            The Lambert transfer equation is used to calculate the delta-v required to transfer from one orbit to another. The equation is given by:
            - $\Delta v$ is the change in velocity required to transfer from one orbit to another
            - $\mu$ is the gravitational parameter of the main attractor
            - $r_1$ is the radius of the initial orbit
            - $a$ is the semi-major axis of the target orbit
            - $r_b$ is the radius of the intermediate orbit
            """
            )

def plot_groundtracks(orbits, attractor, texture_file: str):
    groundtrack_data = []

    for orbit_name, orbit in orbits.items():
        groundtrack = GroundtrackPlotter(orbit)
        lons, lats, _ = groundtrack.to_latlon()

        groundtrack_data.append(
            go.Scattergeo(
                lon=lons,
                lat=lats,
                mode="lines",
                name=orbit_name,
                line=dict(width=2)
            )
        )

    fig = go.Figure(
        data=groundtrack_data,
        layout=go.Layout(
            title=dict(text=f"Groundtracks for {attractor.name}"),
            showlegend=True,
            geo=dict(
                projection_type="equirectangular",
                showland=True,
                landcolor="rgb(243, 243, 243)",
                countrycolor="rgb(204, 204, 204)",
                bgcolor="rgba(0, 0, 0, 0)",
                showocean=True,
                oceancolor="rgb(0, 0, 255)",
                showcountries=True,
                showframe=False,
                showcoastlines=True,
                coastlinecolor="rgb(255, 255, 255)",
                lonaxis=dict(range=[-180, 180]),
                lataxis=dict(range=[-90, 90]),
                showlakes=True,
                lakecolor="rgb(0, 0, 255)",
                showrivers=True,
                rivercolor="rgb(0, 0, 255)",
                projection_rotation_lon=-90,
                center=dict(lat=0, lon=0),
            ),
        ),
    )

    fig.update_geos(
        visible=True,
        showcountries=True,
        showcoastlines=True,
        showlakes=True,
        showocean=True,
        oceancolor="rgb(0, 0, 255)",
        resolution=50,
        bgcolor="rgba(0, 0, 0, 0)",
        projection_rotation_lon=-90,
        center=dict(lat=0, lon=0),
    )

    fig.update_geos(projection_type="equirectangular")

    # Add the planet texture to the background
    fig.update_geos(
        showland=True,
        landcolor="rgba(0, 0, 0, 0)",
        showocean=True,
        oceancolor="rgba(0, 0, 0, 0)",
        showlakes=True,
        lakecolor="rgba(0, 0, 0, 0)",
        showrivers=True,
        rivercolor="rgba(0, 0, 0, 0)",
        bgcolor="rgba(0, 0, 0, 0)",
        fitbounds="locations",
        landimage=texture_file
    )

    return fig

def plot_groundplots(orbits, t_span, projection="equirectangular"):
    # Generate an instance of the plotter, add title and show latlon grid
    gp = GroundtrackPlotter()
    gp.update_layout(title="Groundtrack of Orbits")

    colors = ["red", "blue", "green", "purple", "orange", "yellow"]

    # Plot each orbit in the dictionary
    for label, orbit in orbits.items():
        spacecraft = EarthSatellite(orbit, None)
        # get a random color for the orbit
        random = np.random.RandomState(1234)
        color = random.choice(colors)
        print(color)

        gp.plot(
            spacecraft,
            t_span,
            label=label,
            color=color,
            marker={
                "size": 10,
                "symbol": "triangle-right",
                "line": {"width": 1, "color": color},
            },
        )
    
    gp.update_layout(
        showlegend=True,
        legend_title_text="Orbit",
        legend=dict(x=0, y=1),
        height=700,
        geo=dict(
            showland=True,
            showocean=True,
            showlakes=True,
            showcoastlines=True,
            showcountries=True,
            showsubunits=True,
            projection_type=projection,
            showframe=True,
            resolution=50,
            bgcolor="rgba(0, 0, 0,0)",
        ),
    )
    # Show the plot
    return gp.fig