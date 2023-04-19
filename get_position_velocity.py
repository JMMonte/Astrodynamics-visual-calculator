from astropy import units as u
from poliastro.bodies import *
from poliastro.core.elements import coe2rv
from poliastro.twobody import Orbit

class GetPositionVectors:
    def __init__(self, orbit):
        self.a = orbit.a
        self.ecc = orbit.ecc
        self.inc = orbit.inc
        self.raan = orbit.raan
        self.argp = orbit.argp
        self.body = orbit.attractor
        self.k = orbit.attractor.k.to(u.km**3 / u.s**2).value  # Gravitational parameter in km^3/s^2
        self.p = (self.a * (1 - self.ecc**2)).to(u.km).value  # Semi-latus rectum in km

    def get_position_velocity(self, nu):
        '''
        Return the position and velocity vector for any orbit.
        '''
        r_ijk, v_ijk = coe2rv(self.k, self.p, self.ecc.value, self.inc.to(u.rad).value,
                              self.raan.to(u.rad).value, self.argp.to(u.rad).value, nu.to(u.rad).value)

        position = r_ijk * u.km
        velocity = v_ijk * u.m / u.s

        return position, velocity
    
    def get_periapsis_apoapsis_positions(self):
        '''
        Return the position and velocity vector for periapsis and apoapsis of any given orbit.
        '''
        nu_periapsis = 0 * u.deg
        position_periapsis, _ = self.get_position_velocity(nu_periapsis)

        nu_apoapsis = 180 * u.deg
        position_apoapsis, _ = self.get_position_velocity(nu_apoapsis)

        return position_periapsis, position_apoapsis
