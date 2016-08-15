import math
import numpy as np
from abc import ABCMeta, abstractmethod

import scipy.constants as sconst
from scipy.stats import maxwell

from profilehooks import profile


class ParticleAbstract(object):
    """
    Abstract class which represents particle.
    Every implementation of potential should inherit this class.
    """
    __metaclass__ = ABCMeta

    def __init__(self, box, mass):
        self.box = box
        self.resultant_force = np.array([0, 0, 0], dtype=float)
        self.old_resultant_force = np.array([0, 0, 0], dtype=float)

        self.mass = mass

        # Variable contains information about distances to other particles
        self.distances = None

        self.set_interactions_and_forces()

    def __str__(self):
        return 'r={position}, v={velocity}'.format(**{
            'position': self.r,
            'velocity': self.v,
        })

    @abstractmethod
    def set_interactions_and_forces(self):
        """
        Sets interactions and forces describing potential.
        It must be implemented.
        """
        pass

    def set_position(self, r):
        """
        Method sets position of particle.
        param r: Array3D
        returns: self
        """
        self.r = r

        return self

    def set_velocity(self, v):
        """
        Method sets velocity of particle.
        param v: Array3D
        returns: self
        """
        self.v = v

        return self

    def adjust_position_to_periodic_boundary_conditions(self):
        """
        Adjusts position of particles.
        Take into account periodic boundary conditions.
        """
        self.r = np.where(self.r >= self.box.dim, self.r - self.box.dim, self.r)
        self.r = np.where(self.r < 0, self.box.dim + self.r, self.r)

    def simulate_position_change(self):
        """
        Simulates the movement of the particle.
        Calls adjust_position_to_periodic_boundary_conditions to make sure
        that particle is not outside the box.
        """
        self.r += self.v * self.box.dt + 0.5 * self.resultant_force/self.mass * self.box.dt2

        if self.r[0] < 0:
            self.box.sum_velocity_change -= 2*self.v[0]

        self.adjust_position_to_periodic_boundary_conditions()

    def simulate_velocity_change(self):
        """
        Simulates changing velocity in response to the resultant force.
        """
        self.v += 0.5 * (self.old_resultant_force + self.resultant_force)/self.mass * self.box.dt

    def adjust_forces(self):
        """
        Adjusts forces according to calculations.
        It also make sure that the old resultant force is stored.
        """
        self.old_resultant_force = self.resultant_force
        self.calculate_force_from_patricles()

    def get_distances_from_particles(self, particles):
        """
        It calculates distances between particular partcile and rest particles:
        Takes into account method of images.
        param particles: Array of ParticleAbstract objects
        returns: Array of floats with distances
        """
        if self.distances is not None:
            return self.distances
        if not particles:
            return []
        delta = np.abs(np.array([particle.r for particle in particles], dtype=float) - self.r)
        delta = np.where(delta > 0.5 * self.box.dim, self.box.dim - delta, delta)
        self.distances = np.sqrt((delta ** 2).sum(axis=-1))
        return self.distances

    def calculate_force_from_patricles(self):
        """
        Calculates resultant force acting on the particle.
        """
        particles = [particle for particle in self.box.particles if particle is not self]
        distances = self.get_distances_from_particles(particles)

        self.resultant_force = 0
        for i, particle in enumerate(particles):
            distance = distances[i]
            if distance < self.box.r_c or distance > self.box.r_m:
                continue
            for interaction in self.forces:
                delta = particle.r - self.r
                self.resultant_force += np.where(
                    np.abs(delta) > 0.5 * self.box.dim,
                    delta - self.box.dim * (np.sign(delta)),
                    delta
                ) * interaction(distance) / distance

    def clear_distances(self):
        """
        Clears distances variable
        """
        self.distances = None

    def get_kinetic_energy(self):
        """
        Calculates kinetic energy of particle.
        """
        return 0.5 * self.mass * np.linalg.norm(self.v)**2

    @abstractmethod
    def get_potential_energy(self, particles):
        """
        Calculates potential energy of particle.
        It must be implemented.
        params particles: Array of ParticleAbstract objects
        """
        pass


class ParticleLennard(ParticleAbstract):

    def __init__(self, box, mass, eps, sigma):
        self.eps = eps
        self.sigma = sigma

        super().__init__(box, mass)

    def set_interactions_and_forces(self):
        LENNARD_JONES_INTERACTION = lambda r: 4*self.eps * ((self.sigma/r)**12 - (self.sigma/r)**6)
        LENNARD_JONES_FORCE = lambda r: -48*self.eps/self.sigma * ((self.sigma/r)**13 - 0.5*(self.sigma/r)**7)
        self.interactions = [LENNARD_JONES_INTERACTION]
        self.forces = [LENNARD_JONES_FORCE]

    def get_potential_energy(self, particles):
        particles = [particle for particle in particles if particle is not self]
        distances = self.get_distances_from_particles(particles)

        potential_energy = 0
        for distance in distances:
            for interaction in self.interactions:
                potential_energy += interaction(distance)

        return potential_energy


class Box:
    N_STEPS_CUT_OFF = 10

    def __init__(self, dim, temperature, num_particles, particles_class, dt, duration, file_name, r_c=0):
        self.dim = dim
        self.temperature = temperature
        self.num_particles = num_particles
        self.particles_class = particles_class

        self.file_name = file_name

        self.t = 0.0
        self.dt = dt
        self.dt2 = dt**2
        self.duration = duration

        self.sum_kinetic_energy = 0.0
        self.sum_velocity_change = 0.0

        self.r_c = r_c
        self.r_m = r_c + Box.N_STEPS_CUT_OFF * r_c

        self.particles = []
        self.create_particles()

        self.set_initial_conditions()

        for p in self.particles:
            print(p.r, p.v)

    def set_initial_conditions(self):
        """
        Sets initial conditions inside the box.
        It sets positions and velocities of each particles.
        """
        i = 0
        for x in np.arange(0, self.dim[0], 2.45*self.r_m):
            for y in np.arange(0, self.dim[1], 2.45*self.r_m):
                for z in np.arange(0, self.dim[2], 2.45*self.r_m):
                    if i >= len(self.particles):
                        break
                    self.particles[i].set_position(np.array([x, y, z], dtype=float))
                    i += 1

        # positions = np.random.uniform(0, self.dim, size=(self.num_particles, 3))
        for i, particle in enumerate(self.particles):
            scale = math.sqrt(sconst.k * self.temperature / particle.mass)
            # particle.set_position(positions[i])
            particle.set_velocity(maxwell.rvs(size=3, scale=scale)/math.sqrt(3)*np.random.choice([-1, 1], 3))

    def get_summary_kinetic_energy(self):
        """
        It retrieves information about summary kinetic energy inside the box.
        returns: float
        """
        kinetic_energy = 0
        for particle in self.particles:
            kinetic_energy += particle.get_kinetic_energy()

        return kinetic_energy

    def get_summary_potential_energy(self):
        """
        It retrieves information about summary potential energy inside the box.
        returns: float
        """
        potential_energy = 0
        for particle in self.particles:
            potential_energy += particle.get_potential_energy(self.particles)

        return potential_energy/2

    def get_temperature(self):
        """
        Calculates temperature of gas.
        """
        steps = self.t / self.dt
        return 2/3 * self.sum_kinetic_energy / (self.num_particles * sconst.k * steps)

    def get_pressure(self):
        """
        Calculates pressure of gas.
        """
        return self.particles[0].mass * self.sum_velocity_change / (self.t * self.dim[1] * self.dim[2])


    def create_particles(self):
        """
        Creates particles.
        Adds all particles to array self.particles
        """
        self.particles = []
        args = []
        if self.particles_class == ParticleLennard:
            eps = 125.7 * sconst.k
            sigma = 0.3345 * 10**-9
            mass = 39.95 * 10**-3 / sconst.N_A
            args = [mass, eps, sigma]

        for i in range(self.num_particles):
            self.particles.append(self.particles_class(self, *args))

    @profile
    def simulate(self):
        """
        Simulates movements and interactions between all particles.
        Saves results into chosen file.
        """
        for particle in box.particles:
            particle.calculate_force_from_patricles()
            particle.clear_distances()

        with open(self.file_name, 'w+') as f:
            while self.t < self.duration:
                self.t += self.dt

                for particle in box.particles:
                    particle.simulate_position_change()

                for particle in box.particles:
                    particle.adjust_forces()
                    particle.simulate_velocity_change()

                # particle = box.particles[0]
                kinetic_energy = box.get_summary_kinetic_energy()
                potential_energy = box.get_summary_potential_energy()
                summary_energy = kinetic_energy + potential_energy
                box.sum_kinetic_energy += kinetic_energy
                f.write('%s %s %s %s %s %s\n' % (
                # f.write('%s %s %s %s %s %s %s\n' % (
                    # self.t, particle.r[0], particle.r[1], particle.r[2], kinetic_energy, potential_energy, summary_energy,
                    self.t, kinetic_energy, potential_energy, summary_energy, box.get_temperature(), box.get_pressure(),
                ))

                for particle in box.particles:
                    particle.clear_distances()

        print('Temperature: ', box.get_temperature())
        print('Pressure: ', box.get_pressure())


dim = 1e-8
r_c = 0.25*(0.3345 * 10**-9)
box = Box(np.array([dim, dim, dim]), 500, 100, ParticleLennard, 1e-14, 1e-10, 'simulation-param.dat', r_c)
# box.particles[0].set_position(np.array([r_c, 0, 0], dtype=float)).set_velocity(np.array([1000, 0, 0], dtype=float))
# box.particles[1].set_position(np.array([6*r_c, 0, 0], dtype=float)).set_velocity(np.array([-1000, 0, 0], dtype=float))

# print(box.particles[0])
# print(box.particles[1])

box.simulate()
