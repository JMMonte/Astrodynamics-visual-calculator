import numpy as np
import streamlit as st
from astropy import units as u
import pandas as pd
from poliastro.bodies import *
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.util import time_range
from model import *

configure_page()
planets = define_main_attractors()
maneuver = define_maneuver_types()
mission_epoch = get_mission_epoch()

st.title("Orbit maneuver calculator")
st.write(f"Mission epoch {mission_epoch}")
with st.sidebar:
    st.subheader("Maneuver planner")
    attractor = st.selectbox("Select the main attractor", list(planets.keys()))
    attractor = planets[attractor]
    maneuver_type = st.selectbox("Select the maneuver type", list(maneuver.keys()))
    with st.expander("Initial orbit details"):
        initial_orbit = get_orbit_parameters(
            mission_epoch, attractor, "initial", altitude=342.0
        )
        initial_orbit = Orbit.from_classical(
            initial_orbit[0], # attractor
            initial_orbit[1], # semi-major axis
            initial_orbit[2], # eccentricity
            initial_orbit[3], # inclination
            initial_orbit[4], # raan
            initial_orbit[5], # argp
            initial_orbit[6], # nu
            initial_orbit[7], # epoch
        )

    # Set the final orbit
    with st.expander("Target orbit details"):
        if maneuver_type == "Hohmann transfer":
            target_orbit_altitude = st.number_input(
                "Target orbit apogee altitude (km)", min_value=0.0, value=541.0
            )
            target_orbit_radius = (target_orbit_altitude * u.km) + attractor.R.to(u.km)
            maneuver = Maneuver.hohmann(initial_orbit, target_orbit_radius)
            transfer_orbit = initial_orbit.apply_maneuver(maneuver, intermediate=True)[0]
            # get target orbit from the target orbit radius and the remaining initial orbit parameters. must be circular.
            target_orbit = Orbit.from_classical(
                attractor,
                target_orbit_radius,
                0.0 * u.dimensionless_unscaled,
                initial_orbit.inc,
                initial_orbit.raan,
                initial_orbit.argp,
                initial_orbit.nu,
                epoch=mission_epoch + maneuver.get_total_time()
            )
            orbits = {"Initial Orbit" : initial_orbit,"Transfer Orbit" : transfer_orbit,"Target Orbit" : target_orbit}
        elif maneuver_type == "Lambert transfer":
            target_orbit_epoch_offset = st.sidebar.number_input("Target orbit epoch offset (seconds)", min_value=0.0, value=3600.0)
            target_orbit = get_orbit_parameters(
                mission_epoch + target_orbit_epoch_offset * u.s, # Add the time offset here
                attractor,
                "target",
                altitude=1221.0,
                raan=0.0,
                argp=29.0,
                nu=4.0,
                inclination=120.0
            )
            target_orbit = Orbit.from_classical(
                target_orbit[0], # attractor
                target_orbit[1], # semi-major axis
                target_orbit[2] * u.dimensionless_unscaled, # eccentricity
                target_orbit[3], # inclination
                target_orbit[4], # raan
                target_orbit[5], # argp
                target_orbit[6], # nu
                target_orbit[7], # epoch
            )
            
            # Add a small perturbation to the initial position vector
            perturbation_vector = np.random.normal(0, 1e-6, 3) * u.km
            initial_orbit = initial_orbit.propagate(1 * u.s)
            initial_orbit.r = initial_orbit.r + perturbation_vector

            maneuver = Maneuver.lambert(initial_orbit, target_orbit)
            transfer_orbit = initial_orbit.apply_maneuver(maneuver, intermediate=True)[0]
            orbits = {"Initial Orbit" : initial_orbit,"Transfer Orbit" : transfer_orbit,"Target Orbit" : target_orbit}


        elif maneuver_type == "Bielliptic transfer":
            # define Altitude of the intermediate orbit
            intermediate_orbit_altitude = st.number_input(
                "Intermediate orbit altitude (km)", min_value=0.0, value=8000.0
            )
            # Define Final orbital radius
            target_orbit_altitude = st.number_input(
                "Target orbit altitude (km)", min_value=0.0, value=2541.0
            )
            target_orbit_radius = (target_orbit_altitude * u.km) + attractor.R.to(u.km)
            # bielliptic maneuver preformes 2 impulses to transfer from one orbit to another. the first impulse is to transfer to an intermediate orbit and the second impulse is to transfer to the target orbit.
            # let's first calculate the intermediate orbit
            # get the intermediate orbit radius
            intermediate_orbit_radius = (
                intermediate_orbit_altitude * u.km
            ) + attractor.R.to(u.km)
            maneuver = Maneuver.bielliptic(
                initial_orbit, intermediate_orbit_radius, target_orbit_radius
            )
            # get the 1st transfer orbit from the first impulse returned from the bielliptci maneuver
            intermediate_orbit_1 = initial_orbit.apply_maneuver(
                maneuver, intermediate=True
            )[0]
            # get the 2nd transfer orbit from the second impulse returned from the bielliptic maneuver
            intermediate_orbit_2 = initial_orbit.apply_maneuver(
                maneuver, intermediate=True
            )[1]
            # get the target orbit from the target orbit radius and the remaining initial orbit parameters. must be circular.
            target_orbit = Orbit.from_classical(
                attractor,
                target_orbit_radius,
                0.0 * u.dimensionless_unscaled,
                initial_orbit.inc,
                initial_orbit.raan,
                initial_orbit.argp,
                initial_orbit.nu,
                epoch=mission_epoch + maneuver.get_total_time(),
            )
            orbits = {"Initial Orbit" : initial_orbit,"Intermediat Orbit 1" : intermediate_orbit_1,"Intermediate Orbit 2" : intermediate_orbit_2,"Target Orbit" : target_orbit}

# Convert orbits into pandas dataframe
orbits_df = pd.DataFrame(
    {
        "Semi-major axis (km)": [orbit.a.to(u.km).value for orbit in orbits.values()],
        "Eccentricity": [orbit.ecc.value for orbit in orbits.values()],
        "Inclination (deg)": [orbit.inc.to(u.deg).value for orbit in orbits.values()],
        "RAAN (deg)": [orbit.raan.to(u.deg).value for orbit in orbits.values()],
        "Argument of periapsis (deg)": [
            orbit.argp.to(u.deg).value for orbit in orbits.values()
        ],
        "True anomaly (deg)": [orbit.nu.to(u.deg).value for orbit in orbits.values()],
        "Epoch": [orbit.epoch.iso for orbit in orbits.values()],
    },
    index=orbits.keys(),
)
# Display the orbits dataframe
st.dataframe(orbits_df, use_container_width=True)
show_maneuver_data(maneuver_type)

if maneuver_type == "Hohmann transfer" or maneuver_type == "Lambert transfer":
    st.subheader("Maneuver:")
    col1, col2, col3, col4 = st.columns(4)
    # display the delta-v of the first impulse
    impulse1 = maneuver.impulses[0][1]
    col1.metric(value=f"{impulse1[0]:.2f}", label="Delta-V1")
    # display the delta-v of the second impulse
    impulse2 = maneuver.impulses[1][1]
    col2.metric(value=f"{impulse2[0]:.2f}", label="Delta-V2")
    col3.metric(
        value=f"{np.linalg.norm(maneuver.impulses[0][1]) + np.linalg.norm(maneuver.impulses[1][1]):.2f}",
        label="Total Delta-V",
    )
    col4.metric(
        value=f"{maneuver.get_total_time().to(u.s).value:.2f} s", label="Total time"
    )

elif maneuver_type == "Bielliptic transfer":
    st.subheader("Maneuver:")
    col1, col2 = st.columns(2)
    col1.metric(value=f"{maneuver.impulses[0][1].value[0]:.4f} m/s", label="Delta-V1")
    col2.metric(value=f"{maneuver.impulses[1][1].value[0]:.4f} m/s", label="Delta-V2")
    col1.metric(value=f"{maneuver.impulses[2][1].value[0]:.4f} m/s", label="Delta-V3")
    col2.metric(value=f"{maneuver.get_total_cost():.4f} m/s", label="Total Delta-V")

# plot the orbits using plotly
st.subheader("Orbits projection:")
fig1 = plotly_orbit_plotter(
    orbits.values(),
    attractor,
    labels=orbits.keys(),
)

st.plotly_chart(fig1, use_container_width=True)

# plot the orbits using the groundtrack plotter native to poliastro
if attractor == Earth:
    # Define the time span for the orbits
    t_span = time_range(initial_orbit.epoch - 1.5 * u.h, periods=150, end=target_orbit.epoch + 1.5 * u.h)

    # Call the function with the orbits dictionary and the time span
    st.plotly_chart(plot_groundplots(orbits, t_span), use_container_width=True)
    st.plotly_chart(plot_groundplots(orbits, t_span, projection="orthographic"), use_container_width=True)
