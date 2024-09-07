#include <iostream>
#include <vector>
#include <array>
#include <random>

class MDSimulation {
  public:
    MDSimulation(double box_size_in, int nparticles_in);
    void run(int nsteps, double dt);
  private:
    double cutoff;                 // Cutoff radius for Lennard-Jones potential
    double cutoff2;                // Square of the cutoff radius
    double lj_potential_at_cutoff; // Value of the Lennard-Jones potential at the cutoff
    double box_size;               // Length of each side of the periodic simulation cell, which is cubic.
    int nparticles;                // Number of particles in the simulation
    std::vector<std::array<double, 3>> positions;  // Position of the particles
    std::vector<std::array<double, 3>> velocities; // Velocities of the particles
    std::vector<std::array<double, 3>> forces;     // Forces on the particles

    double minimum_image(double dist);
    double lj_potential(double r2);
    double lj_potential_with_cutoff(double r2);
    double lj_force(double r2);
    double lj_force_with_cutoff(double r2);
};

/*! \brief Initialize a molecular dynamics simulation
 *
 * \param [in]  box_size_in
 *                   Length of each side of the periodic simulation cell, which is cubic.
 * \param [in]  nparticles_in
 *                   Number of particles in the simulation.
 */
MDSimulation::MDSimulation(double box_size_in, int nparticles_in) {
    box_size = box_size_in;
    nparticles = nparticles_in;
    cutoff = 2.5;
    cutoff2 = cutoff * cutoff;

    // Initialize the particles on a rough grid
    int particles_per_side = std::ceil( std::pow(nparticles, 1.0/3.0) );
    double particle_spacing = box_size / (particles_per_side + 1);
    for (int iparticle = 0; iparticle < nparticles; iparticle++) {
        int ix = iparticle % particles_per_side;
        int iy = (iparticle / particles_per_side) % particles_per_side;
        int iz = iparticle / (particles_per_side * particles_per_side);
        positions.push_back({
            particle_spacing * ix + ( 0.5 * particle_spacing ),
            particle_spacing * iy + ( 0.5 * particle_spacing ),
            particle_spacing * iz + ( 0.5 * particle_spacing )
            });
    }

    // Initialize the velocities randomly
    for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
        /* The random number generators here use the particle index as the seed.
           This isn't something you would normally do, but it is quite helpful
           in this case for the purpose of ensuring that the velocities are 
           reproducible with respect to parallelization. */
        std::mt19937 gen(iparticle);
        std::uniform_real_distribution<double> random_vel(-0.5, 0.5);
        velocities.push_back({random_vel(gen), random_vel(gen), random_vel(gen)});
    }

    // Initialize the forces
    for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
        forces.push_back({0.0, 0.0, 0.0});
    }

    // Determine the Lennard-Jones potential at the cutoff
    lj_potential_at_cutoff = lj_potential(cutoff2);
}

/*! \brief Evaluate the Lennard-Jones potential associated with a specific particle separation.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double MDSimulation::lj_potential(double r2) {
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    return 4.0 * (inv_r6 * inv_r6 - inv_r6);
}

/*! \brief Evaluate the Lennard-Jones potential associated with a specific particle separation, with a cutoff.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double MDSimulation::lj_potential_with_cutoff(double r2) {
    if (r2 < cutoff2) {
        return lj_potential(r2) - lj_potential_at_cutoff;
    }
    else {
        return 0.0;
    }
}

/*! \brief Evaluate the Lennard-Jones force for a specific particle separation.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double MDSimulation::lj_force(double r2) {
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    return 24.0 * inv_r2 * (2.0 * inv_r6 * inv_r6 - inv_r6);
}

/*! \brief Evaluate the Lennard-Jones force for a specific particle separation, with a cutoff.
 *
 * \param [in]  r2
 *                   Square of the distance between two particles.
 */
double MDSimulation::lj_force_with_cutoff(double r2) {
    if (r2 < cutoff2) {
        return lj_force(r2);
    }
    else {
        return 0.0;
    }
}

/*! \brief Account for periodic boundary conditions by returning the smallest magnitude distance between two particles.
 *
 * \param [in]  dist
 *                   Distance between the particles
 */
double MDSimulation::minimum_image(double dist) {
    if (dist > 0.5 * box_size) dist -= box_size;
    if (dist < -0.5 * box_size) dist += box_size;
    return dist;
}

/*! \brief Run a molecular dynamics simulation.
 *
 * \param [in]  nsteps
 *                   Number of time integration steps to perform.
 * \param [in]  dt
 *                   Size of the timestep (reduced Lennard-Jones units).
 */
void MDSimulation::run(int nsteps, double dt) {

    // Main simulation loop
    for (int istep = 0; istep < nsteps; ++istep) {

        // Update the particle velocities and positions
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {

            // Update the positions
            positions[iparticle][0] += velocities[iparticle][0] * dt;
            positions[iparticle][1] += velocities[iparticle][1] * dt;
            positions[iparticle][2] += velocities[iparticle][2] * dt;

            // Apply periodic boundary conditions; ensure that particles outside the box wrap to the other side
            for (int idimension = 0; idimension < 3; ++idimension) {
                if (positions[iparticle][idimension] < 0.0) positions[iparticle][idimension] += box_size;
                if (positions[iparticle][idimension] >= box_size) positions[iparticle][idimension] -= box_size;
            }

        }

        // Zero the energy and forces
        double potential_energy = 0.0;
        double kinetic_energy = 0.0;

        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            forces[iparticle] = {0.0, 0.0, 0.0};
        }

        /* Evaluate the energy and forces
           Note: This code actually evaluates the interaction between each pair of atoms twice.
           In many molecular dynamics codes, it would be typical for the inner loop over jparticle
           to be evaluated over the range [iparticle + 1, nparticles] instead of [0, nparticles].
           That would ensure that each pair of iparticle and jparticle is only evaluated once,
           and would nearly cut the amount of work required to evaluate forces by a factor of two.
           Unfortunately, it would also make it a bit trickier to parallelize, so you can be glad
           that this loop doesn't do that.  You should still be sure to handle the condition
           inside the inner loop that verifies that the two particles are different.
        */
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            for (int jparticle = 0; jparticle < nparticles; ++jparticle) {
                if ( iparticle != jparticle ) { // Only compute the interactions between different particles
                    double dx = minimum_image(positions[iparticle][0] - positions[jparticle][0]);
                    double dy = minimum_image(positions[iparticle][1] - positions[jparticle][1]);
                    double dz = minimum_image(positions[iparticle][2] - positions[jparticle][2]);
                    double r2 = (dx * dx) + (dy * dy) + (dz * dz);

                    double f = lj_force_with_cutoff(r2);
                    forces[iparticle][0] += f * dx;
                    forces[iparticle][1] += f * dy;
                    forces[iparticle][2] += f * dz;
                
                    potential_energy += 0.5 * lj_potential_with_cutoff(r2);
                }
            }
        }

        // Compute the kinetic energy
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            kinetic_energy += 0.5 * velocities[iparticle][0] * velocities[iparticle][0];
            kinetic_energy += 0.5 * velocities[iparticle][1] * velocities[iparticle][1];
            kinetic_energy += 0.5 * velocities[iparticle][2] * velocities[iparticle][2];
        }

        // Update the particle velocities
        for (int iparticle = 0; iparticle < nparticles; ++iparticle) {
            velocities[iparticle][0] += forces[iparticle][0] * dt;
            velocities[iparticle][1] += forces[iparticle][1] * dt;
            velocities[iparticle][2] += forces[iparticle][2] * dt;
        }

        // Print output
        std::cout << "Iteration " << istep << std::endl;
        std::cout << "    Potential Energy: " << potential_energy << std::endl;
        std::cout << "    Kinetic Energy:   " << kinetic_energy << std::endl;
        std::cout << "    Total Energy:     " << potential_energy + kinetic_energy << std::endl << std::endl;
    }

    std::cout << "Simulation completed." << std::endl;
}

int main(int argc, char** argv) {
    MDSimulation mysimulation(20.0, 1000);
    mysimulation.run(100, 0.005);
    return 0;
}