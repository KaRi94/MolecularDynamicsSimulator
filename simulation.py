import math
import numpy as np
from abc import ABCMeta, abstractmethod

from profilehooks import profile


class ParticleAbstract(object):

    N_A = 6.022140857 * 10**23
    BOLTZMANN_CONST = 1.38 * 10**-23

    __metaclass__ = ABCMeta

    def __init__(self, box):
        self.box = box
        self.set_initial_conditions()
        self.resultant_force = np.array([0, 0, 0], dtype=float)

        # Variable contains information about distances to other particles
        self.distances = []

        self.set_interactions_and_forces()

    @abstractmethod
    def set_interactions_and_forces(self):
        pass

    def set_position(self, r):
        self.r = r

        return self

    def set_velocity(self, v):
        self.v = v

        return self

    def set_initial_conditions(self):
        # rozklad jednorodny dla polozenia + dla predkosci
        pass

    def adjust_position_to_periodic_boundary_conditions(self):
        self.r = np.where(self.r >= self.box.dim, self.r - self.box.dim, self.r)
        self.r = np.where(self.r < 0, self.box.dim + self.r, self.r)

    def simulate_position_change(self):
        self.r += self.v * self.box.dt

        self.adjust_position_to_periodic_boundary_conditions()

    def simulate_velocity_change(self):
        self.v += self.resultant_force/self.mass * self.box.dt

    def get_distances_from_particles(self, particles):
        if self.distances:
            return self.distances
        delta = np.abs(np.array([particle.r for particle in particles], dtype=float) - self.r)
        delta = np.where(delta > 0.5 * self.box.dim, self.box.dim - delta, delta)
        self.distances = np.sqrt((delta ** 2).sum(axis=-1))
        return self.distances

    def calculate_force_from_patricles(self):
        particles = [particle for particle in self.box.particles if particle is not self]
        distances = self.get_distances_from_particles(particles)

        self.resultant_force = 0
        for i, particle in enumerate(particles):
            distance = distances[i]
            # if distance < self.box.r_c or distances > self.box.r_m:
            #     print(distance)
            #     continue
            for interaction in self.forces:
                delta = particle.r - self.r
                self.resultant_force += np.where(
                    np.abs(delta) > 0.5 * self.box.dim,
                    delta - self.box.dim * (np.sign(delta)),
                    delta
                ) * interaction(distance) / distance

        # Remove distances array
        self.distances = []

    def get_kinetic_energy(self):
        return 0.5 * self.mass * np.linalg.norm(self.v)**2

    @abstractmethod
    def get_potential_energy(self, particles):
        pass


class ParticleLennard(ParticleAbstract):

    def __init__(self, box, mass, eps, sigma):
        self.mass = mass
        self.eps = eps
        self.sigma = sigma

        super().__init__(box)

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

    N_STEPS_CUT_OFF = 20

    def __init__(self, dim, num_particles, particles_class, dt, duration, r_c=0):
        self.dim = dim
        self.num_particles = num_particles
        self.particles_class = particles_class

        self.dt = dt
        self.duration = duration

        self.r_c = r_c
        self.r_m = r_c + Box.N_STEPS_CUT_OFF * r_c

        self.particles = []
        self.create_particles()

    def get_summary_kinetic_energy(self):
        kinetic_energy = 0
        for particle in self.particles:
            kinetic_energy += particle.get_kinetic_energy()

        return kinetic_energy

    def get_summary_potential_energy(self):
        potential_energy = 0
        for particle in self.particles:
            potential_energy += particle.get_potential_energy(self.particles)

        return potential_energy/2

    def create_particles(self):
        self.particles = []
        args = []
        if self.particles_class == ParticleLennard:
            eps = 125.7 * ParticleAbstract.BOLTZMANN_CONST
            sigma = 0.3345 * 10**-9
            mass = 39.95 * 10**-3 / ParticleAbstract.N_A
            args = [mass, eps, sigma]

        for i in range(self.num_particles):
            self.particles.append(self.particles_class(self, *args))

    @profile
    def simulate(self):
        t = 0.0

        with open('particle.dat', 'w+') as f:
            while t < self.duration:
                t += self.dt
                particle = box.particles[0]
                kinetic_energy = box.get_summary_kinetic_energy()
                potential_energy = box.get_summary_potential_energy()
                summary_energy = kinetic_energy + potential_energy
                f.write('%s %s %s %s %s %s %s\n' % (
                    t, particle.r[0], particle.r[1], particle.r[2], kinetic_energy, potential_energy, summary_energy
                ))
                for particle in box.particles:
                    particle.calculate_force_from_patricles()

                for particle in box.particles:
                    particle.simulate_velocity_change()
                    particle.simulate_position_change()


dim = 1e-8
box = Box(np.array([dim, dim, dim]), 2, ParticleLennard, 1e-15, 5e-12, 2.5*(0.3345 * 10**-9))
box.particles[0].set_position(np.array([0, dim/2, dim/2 + 4*(0.3345 * 10**-9)], dtype=float)).set_velocity(np.array([0, 0, -1000], dtype=float))
box.particles[1].set_position(np.array([0, dim/2, dim/2 - 4*(0.3345 * 10**-9)], dtype=float)).set_velocity(np.array([0, 0, 1000], dtype=float))

box.simulate()
