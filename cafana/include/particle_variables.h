/**
 * @file particle_variables.h
 * @brief Header file for definitions of variables which act on single
 * particles.
 * @details This file contains definitions of variables which act on single
 * particles. Each variable is implemented as a function which takes a particle
 * object as an argument and returns a double. These variables are intended to
 * be used to define more complex variables which act on interactions.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_VARIABLES_H
#define PARTICLE_VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

/**
 * @namespace pvars
 * @brief Namespace for organizing generic variables which act on single
 * particles.
 * @details This namespace is intended to be used for organizing variables which
 * act on single particles. Each variable is implemented as a function which
 * takes a particle object as an argument and returns a double. The function
 * should be templated on the type of particle object if the variable is
 * intended to be used on both true and reconstructed particles.
 * @note The namespace is intended to be used in conjunction with the
 * vars namespace, which is used for organizing variables which act on
 * interactions.
 */
namespace pvars
{
    /**
     * @brief Variable for the best estimate of the particle energy.
     * @details At the most basic decision level, this is based on the
     * shower/track designation. Showers can only be reconstructed
     * calorimetrically, while tracks can be reconstructed calorimetrically,
     * by range (if contained), or by multiple scattering (if exiting).
     * @tparam T the type of particle.
     * @param p the particle to apply the variable on.
     * @return the best estimate of the particle energy.
     */
    template<class T>
        double energy(const T & p)
        {
            double energy = 0;
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            {
                energy += p.energy_deposit;
            }
            else
            {
                // Check if the particle is a shower.
                if(p.pid < 2) energy += p.calo_ke;
                else
                {
                    if(p.is_contained) energy += p.csda_ke;
                    else energy += p.mcs_ke;
                }
            }
            return energy;
        }

    /**
     * @brief Variable for true particle starting kinetic energy.
     * @details The starting kinetic energy is defined as the total energy
     * minus the rest mass energy of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the starting kinetic energy of the particle.
     */
    template<class T>
        double ke(const T & p)
        {
            double energy(0);
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            {
                energy = p.energy_init;
                switch(p.pid)
                {
                    case 1:
                        energy -= ELECTRON_MASS;
                        break;
                    case 2:
                        energy -= MUON_MASS;
                        break;
                    case 3:
                        energy -= PION_MASS;
                        break;
                    case 4:
                        energy -= PROTON_MASS;
                        break;
                    default:
                        break;
                }
            }
            else
            {
                energy = pvars::energy(p);
            }
            return energy;
        }
    
    /**
     * @brief Variable for the transverse momentum of a particle.
     * @details The transverse momentum is defined as the square root of the
     * sum of the squares of the x and y components of the momentum. This
     * variable is useful for identifying particles which are produced in a
     * transverse direction to the beam.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the transverse momentum of the particle.
     */
    template<class T>
        double transverse_momentum(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::sqrt(std::pow(p.truth_momentum[0], 2) + std::pow(p.truth_momentum[1], 2));
            else
                return std::sqrt(std::pow(p.momentum[0], 2) + std::pow(p.momentum[1], 2));
        }

    /**
     * @brief Variable for the polar angle (w.r.t the z-axis) of the particle.
     * @details The polar angle is defined as the arccosine of the z-component
     * of the momentum vector. This variable is useful for identifying particles
     * which are produced transversely to the beam.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the polar angle of the particle.
     */
    template<class T>
        double polar_angle(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(p.truth_start_dir[2]);
            else
                return std::acos(p.start_dir[2]);
        }

    /**
     * @brief Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @details The azimuthal angle is defined as the arccosine of the x-component
     * of the momentum vector divided by the square root of the sum of the squares
     * of the x and y components of the momentum vector.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the azimuthal angle of the particle.
     */
    template<class T>
        double azimuthal_angle(const T & p)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::acos(p.truth_start_dir[0] / std::sqrt(std::pow(p.truth_start_dir[0], 2) + std::pow(p.truth_start_dir[1], 2)));
            else
                return std::acos(p.start_dir[0] / std::sqrt(std::pow(p.start_dir[0], 2) + std::pow(p.start_dir[1], 2)));
        }

    /**
     * @brief Variable for the x-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-coordinate of the particle end point.
     */
    template<class T>
        double end_x(const T & p)
        {
            return p.end_point[0];
        }

    /**
     * @brief Variable for the y-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-coordinate of the particle end point.
     */
    template<class T>
        double end_y(const T & p)
        {
            return p.end_point[1];
        }
    
    /**
     * @brief Variable for the z-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-coordinate of the particle end point.
     */
    template<class T>
        double end_z(const T & p)
        {
            return p.end_point[2];
        }
    
    /**
     * @brief Variable for the muon softmax score of the particle.
     * @details The muon softmax score represents the confidence that the
     * network has in the particle being a muon. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a muon.
     * @param p the particle to apply the variable on.
     * @return the muon softmax score of the particle.
     */
    double muon_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[2];
    }

    /**
     * @brief Variable for the pion softmax score of the particle.
     * @details The pion softmax score represents the confidence that the
     * network has in the particle being a pion. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a pion.
     * @param p the particle to apply the variable on.
     * @return the pion softmax score of the particle.
     */
    double pion_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[3];
    }

    /**
     * @brief Variable for the proton softmax score of the particle.
     * @details The proton softmax score represents the confidence that the
     * network has in the particle being a proton. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a proton.
     * @param p the particle to apply the variable on.
     * @return the proton softmax score of the particle.
     */
    double proton_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[4];
    }

    /**
     * @brief Variable for the "MIP" softmax score of the particle.
     * @details The "MIP" softmax score is calculated as the sum of the softmax
     * scores for the muon and pion. The score represents the confidence that
     * the network has in the particle being a minimum ionizing particle.
     * @param p the particle to apply the variable on.
     * @return the "MIP" softmax score of the particle.
     */
    double mip_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[2] + p.pid_scores[3];
    }

    /**
     * @brief Variable for the "hadron" softmax score of the particle.
     * @details The "hadron" softmax score is calculated as the sum of the softmax
     * scores for the pion and proton. The score represents the confidence that
     * the network has in the particle being a hadron.
     * @param p the particle to apply the variable on.
     * @return the "hadron" softmax score of the particle.
     */
    double hadron_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[3] + p.pid_scores[4];
    }

    /**
     * @brief Variable for the primary softmax score of the particle.
     * @details The primary softmax score represents the confidence that the
     * network has in the particle being a primary particle. The score is between
     * 0 and 1, with 1 being the most confident that the particle is a primary
     * particle.
     * @param p the particle to apply the variable on.
     * @return the primary softmax score of the particle.
     */
    double primary_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.primary_scores[1];
    }

    /**
     * @brief Variable for the secondary softmax score of the particle.
     * @details The secondary softmax score represents the confidence that the
     * network has in the particle being a secondary particle. The score is between
     * 0 and 1, with 1 being the most confident that the particle is a secondary
     * particle.
     * @param p the particle to apply the variable on.
     * @return the secondary softmax score of the particle.
     */
    double secondary_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.primary_scores[0];
    }
}
#endif // PARTICLE_VARIABLES_H