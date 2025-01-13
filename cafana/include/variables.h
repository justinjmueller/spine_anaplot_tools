/**
 * @file variables.h
 * @brief Header file for definitions of analysis variables.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions. Each variable is implemented
 * as a function which takes an interaction object as an argument and returns a
 * double. These are the building blocks for producing high-level plots of the
 * selected interactions.
 * @author mueller@fnal.gov
*/
#ifndef VARIABLES_H
#define VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/Proxy/EpilogFwd.h"

#include "include/particle_variables.h"
#include "include/particle_cuts.h"
#include "include/cuts.h"
#include "include/utilities.h"

/**
 * @namespace vars
 * @brief Namespace for organizing generic variables which act on interactions.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions. Each variable is implemented as a function which takes
 * an interaction object as an argument and returns a double. The function
 * should be templated on the type of interaction object if the variable is
 * intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the
 * pvars namespace, which is used for organizing variables which act on single
 * particles.
 */
namespace vars
{
    /**
     * @brief Variable for the neutrino ID of the interaction.
     * @details This variable is intended to provide a unique identifier for
     * each parent neutrino within the event record. This number is assigned
     * starting at 0 for the first neutrino in the event and is incremented
     * for each subsequent neutrino. Non-neutrino interactions are assigned
     * a value of -1.
     * @param obj the interaction to apply the variable on.
     * @return the neutrino ID.
     */
    double neutrino_id(const caf::SRInteractionTruthDLPProxy & obj) { return obj.nu_id; }

    /**
     * @brief Variable for a basic enumeration of interaction categories by
     * the interaction mode.
     * @details This variable is based on the interaction mode and is intended
     * to provide a simple classification of interactions. 
     * @param obj the interaction to apply the variable on.
     * @return the interaction category.
     */
    double neutrino_interaction_mode(const caf::SRInteractionTruthDLPProxy & obj)
    {
        double cat(-1);
        if(cuts::neutrino(obj))
            cat = obj.nu_interaction_mode;
        return cat;
    }

    /**
     * @brief Variable for the true neutrino energy.
     * @units GeV
     * @tparam T the type of object.
     * @param obj the interaction to apply the variable on.
     * @return the true neutrino energy.
     */
    template<class T>
        double true_neutrino_energy(const T & obj) { return obj.nu_energy_init; }

    /**
     * @brief Variable for the true neutrino baseline. Currently, if the
     * distance is not available, the value is set to 585.0 meters.
     * @units meters
     * @tparam T the type of object.
     * @param obj interaction to apply the variable on.
     * @return the true neutrino baseline.
     */
    template<class T>
        double true_neutrino_baseline(const T & obj) { return obj.nu_distance_travel > 0 ? double(obj.nu_distance_travel) : 585.0; }

    /**
     * @brief Variable for the true neutrino PDG code.
     * @units none
     * @tparam T the type of object.
     * @param obj interaction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    template<class T>
        double true_neutrino_pdg(const T & obj) { return obj.nu_pdg_code; }

    /**
     * @brief Variable for the true neutrino current value.
     * @units none
     * @tparam T the type of object.
     * @param obj interaction to apply the variable on.
     * @return the true neutrino current value.
     */
    template<class T>
        double true_neutrino_cc(const T & obj) { return obj.nu_current_type; }

    /**
     * @brief Variable for the containment status of the interaction.
     * @details The containment status is determined upstream in the SPINE
     * reconstruction and is based on the set of all points in the interaction,
     * which must be contained within the volume of the TPC that created them.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the containment status of the interaction.
     */
    template<class T>
        double containment(const T & obj) { return obj.is_contained; }

    /**
     * @brief Variable for the fiducial volume status of the interaction.
     * @details The fiducial volume status is determined upstream in the SPINE
     * reconstruction and is a requirement that the interaction vertex is within
     * the fiducial volume of the TPC.
     */
    template<class T>
        double fiducial(const T & obj) { return obj.is_fiducial; }

    /**
     * @brief Variable for total visible energy of interaction.
     * @units GeV
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
        double visible_energy(const T & obj)
        {
            double energy(0);
            for(const auto & p : obj.particles)
            {
                if(p.is_primary)
                {
                    energy += pvars::energy(p);
                    if(p.pid == 2) energy += MUON_MASS;
                    else if(p.pid == 3) energy += PION_MASS;
                }
            }
            return energy/1000.0;
        }

    /**
     * @brief Variable for the flash time of the interaction.
     * @details The flash time is the time of the flash observed in the PMTs
     * and associated with the charge deposition in the interaction using the
     * OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash time of the interaction.
     */
    template<class T>
        double flash_time(const T & obj) { return obj.flash_time; }

    /**
     * @brief Variable for the flash total photoelectron count of the
     * interaction.
     * @details The flash total photoelectron count is the total number of
     * photoelectrons observed in the PMTs and associated with the charge
     * deposition in the interaction using the OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash total photoelectron count of the interaction.
     */
    template<class T>
        double flash_total_pe(const T & obj) { return obj.flash_total_pe; }

    /**
     * @brief Variable for the flash hypothesis total photoelectron count of
     * the interaction.
     * @details The flash hypothesis total photoelectron count is the total
     * number of photoelectrons predicted by OpT0Finder for the interaction in
     * the flash associated with the charge deposition in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash hypothesis total photoelectron count of the interaction.
     */
    template<class T>
        double flash_hypothesis(const T & obj) { return obj.flash_hypothesis; }

    /**
     * @brief Variable for the x-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the x-coordinate of the interaction vertex.
     */
    template<class T>
        double vertex_x(const T & obj)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return obj.truth_vertex[0];
            else
                return obj.vertex[0];
        }

    /**
     * @brief Variable for the y-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the y-coordinate of the interaction vertex.
     */
    template<class T>
        double vertex_y(const T & obj)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return obj.truth_vertex[1];
            else
                return obj.vertex[1];
        }

    /**
     * @brief Variable for the z-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the z-coordinate of the interaction vertex.
     */
    template<class T>
        double vertex_z(const T & obj)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return obj.truth_vertex[2];
            else
                return obj.vertex[2];
        }

    /**
     * @brief Variable for the x-coordinate of the leading muon end point.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The end point is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the x-coordinate of the leading muon end point.
     */
    template<class T>
        double leading_muon_end_x(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            return m.end_point[0];
        }

    /**
     * @brief Variable for the y-coordinate of the leading muon end point.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The end point is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the y-coordinate of the leading muon end point.
     */
    template<class T>
        double leading_muon_end_y(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            return m.end_point[1];
        }

    /**
     * @brief Variable for the z-coordinate of the leading muon end point.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The end point is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the z-coordinate of the leading muon end point.
     */
    template<class T>
        double leading_muon_end_z(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            return m.end_point[2];
        }

    /**
     * @brief Variable for the x-coordinate of the leading proton end point.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. The end point is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the x-coordinate of the leading proton end point.
     */
    template<class T>
        double leading_proton_end_x(const T & obj)
        {
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            return p.end_point[0];
        }

    /**
     * @brief Variable for the y-coordinate of the leading proton end point.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. The end point is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the y-coordinate of the leading proton end point.
     */
    template<class T>
        double leading_proton_end_y(const T & obj)
        {
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            return p.end_point[1];
        }
    
    /**
     * @brief Variable for the z-coordinate of the leading proton end point.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. The end point is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the z-coordinate of the leading proton end point.
     */
    template<class T>
        double leading_proton_end_z(const T & obj)
        {
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            return p.end_point[2];
        }

    /**
     * Variable for the muon softmax score for the leading muon of the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the muon softmax score of the leading muon.
    */
    
    template<class T>
        double electron_softmax(const T & interaction)
        {
            sauto & m(obj.particles[utilities::leading_particle_index(obj, 1)]);
            return m.pid_scores[1];
        }

    /**
     * @brief Variable for the muon softmax score of the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The softmax score can be thought of as a "confidence"
     * value for the particle identification.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the muon softmax score of the leading muon.
     */
    template<class T>
        double leading_muon_softmax(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            return m.pid_scores[2];
        }

    /**
     * @brief Variable for the proton softmax score of the leading proton.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. The softmax score can be thought of as a "confidence"
     * value for the particle identification.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the proton softmax score of the leading proton.
     */
    template<class T>
        double leading_proton_softmax(const T & obj)
        {
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            return p.pid_scores[4];
        }

    /**
     * @brief Variable for the "MIP" softmax score of the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The softmax score can be thought of as a "confidence"
     * value for the particle identification. In this case, the "MIP" score
     * is calculated as the sum of the softmax scores for the muon and pion.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the "MIP" softmax score of the leading muon.
     */
    template<class T>
        double leading_muon_mip_softmax(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            return m.pid_scores[2] + m.pid_scores[3];
        }

    /**
     * @brief Variable for finding the leading muon kinetic energy.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
     */
    template<class T>
        double leading_muon_ke(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 2));
            double energy(pvars::energy(obj.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = pvars::ke_init(obj.particles[i]);
            return energy;
        }

    /**
     * @brief Variable for finding the leading muon kinetic energy.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
     */
    template<class T>
        double leading_electron_ke(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 2));
            double energy(pvars::energy(obj.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = pvars::ke_init(obj.particles[i]);
            return energy;
        }
    

    /**
     * @brief Variable for finding the leading proton kinetic energy.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the kinetic energy of the leading proton.
     */
    template<class T>
        double leading_proton_ke(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 4));
            double energy(pvars::energy(obj.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = pvars::ke_init(obj.particles[i]);
            return energy;
        }
    
    /**
     * @brief Variable for the transverse momentum of the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The transverse momentum is defined as the square root
     * of the sum of the squares of the x and y components of the momentum.
     * This variable is useful for identifying muons which are produced in a
     * transverse direction to the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the leading muon.
     */
    template<class T>
        double leading_muon_pt(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 2));
            return pvars::transverse_momentum(obj.particles[i]);
        }

    /**
     * @brief Variable for the transverse momentum of the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. The transverse momentum is defined as the square root
     * of the sum of the squares of the x and y components of the momentum.
     * This variable is useful for identifying muons which are produced in a
     * transverse direction to the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the leading muon.
     */
    template<class T>
        double leading_electron_pt(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 1));
            return pvars::transverse_momentum(obj.particles[i]);
        }

    /**
     * @brief Variable for the transverse momentum of the leading proton.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. The transverse momentum is defined as the square root
     * of the sum of the squares of the x and y components of the momentum.
     * This variable is useful for identifying protons which are produced in a
     * transverse direction to the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the leading proton.
     */
    template<class T>
        double leading_proton_pt(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 4));
            return pvars::transverse_momentum(obj.particles[i]);
        }

    /**
     * @brief Variable for the muon polar angle.
     * @details The polar angle is defined as the arccosine of the z-component
     * of the momentum vector. This variable is useful for identifying muons
     * which are produced transversely to the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the polar angle of the leading muon.
     */
    template<class T>
        double muon_polar_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 2));
            return pvars::polar_angle(obj.particles[i]);
        }

    /**
     * @brief Variable for the muon azimuthal angle.
     * @details The azimuthal angle is defined as the arccosine of the x-component
     * of the momentum vector divided by the square root of the sum of the squares
     * of the x and y components of the momentum vector.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the azimuthal angle of the leading muon.
     */
    template<class T>
        double muon_azimuthal_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 2));
            return pvars::azimuthal_angle(obj.particles[i]);
        }

    /**
     * @brief Variable for the muon polar angle.
     * @details The polar angle is defined as the arccosine of the z-component
     * of the momentum vector. This variable is useful for identifying muons
     * which are produced transversely to the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the polar angle of the leading muon.
     */
    template<class T>
        double electron_polar_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 1));
            return pvars::polar_angle(obj.particles[i]);
        }

    /**
     * @brief Variable for the muon azimuthal angle.
     * @details The azimuthal angle is defined as the arccosine of the x-component
     * of the momentum vector divided by the square root of the sum of the squares
     * of the x and y components of the momentum vector.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the azimuthal angle of the leading muon.
     */
    template<class T>
        double electron_azimuthal_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 1));
            return pvars::azimuthal_angle(obj.particles[i]);
        }

    /**
     * @brief Variable for the transverse momentum of the interaction.
     * @details The transverse momentum is defined as the square root of the
     * sum of the squares of the x and y components of the momentum. This
     * variable is useful for identifying interactions that have a significant
     * transverse component (e.g. missing momentum).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
     */
    template<class T>
        double interaction_pt(const T & obj)
        {
            double px(0), py(0);
            for(const auto & p : obj.particles)
                if(p.is_primary)
                {
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        px += p.truth_momentum[0];
                        py += p.truth_momentum[1];
                    }
                    else
                    {
                        px += p.momentum[0];
                        py += p.momentum[1];
                    }
                }
            return std::sqrt(std::pow(px, 2) + std::pow(py, 2));
        }

    /**
     * @brief Variable for phi_T of the interaction.
     * @details phi_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the leading muon and the total hadronic
     * system. This variable is sensitive to the presence of F.S.I.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the phi_T of the interaction.
     */
    template<class T>
        double phiT(const T & obj)
        {
            double lpx(0), lpy(0), hpx(0), hpy(0);
            for(const auto & p : obj.particles)
                if(pcuts::final_state_signal(p))
                {
                    if(p.pid > 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            hpx += p.truth_momentum[0];
                            hpy += p.truth_momentum[1];
                        }
                        else
                        {
                            hpx += p.momentum[0];
                            hpy += p.momentum[1];
                        }
                    }
                    else if(p.pid == 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            lpx += p.truth_momentum[0];
                            lpy += p.truth_momentum[1];
                        }
                        else
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                    }
                }
            return std::acos((-hpx * lpx - hpy * lpy) / (std::sqrt(std::pow(hpx, 2) + std::pow(hpy, 2)) * std::sqrt(std::pow(lpx, 2) + std::pow(lpy, 2))));
        }

    /**
     * @brief Variable for alpha_T of the interaction.
     * @details alpha_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the total hadronic system and the muon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the alpha_T of the interaction.
     */
    template<class T>
        double alphaT(const T & obj)
        {
            double lpx(0), lpy(0), px(0), py(0);
            for(const auto & p : obj.particles)
                if(pcuts::final_state_signal(p))
                {
                    if(p.pid <= 2)
                    {
                        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                        {
                            lpx += p.truth_momentum[0];
                            lpy += p.truth_momentum[1];
                        }
                        else
                        {
                            lpx += p.momentum[0];
                            lpy += p.momentum[1];
                        }
                    }
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        px += p.truth_momentum[0];
                        py += p.truth_momentum[1];
                    }
                    else
                    {
                        px += p.momentum[0];
                        py += p.momentum[1];
                    }
                }
            return std::acos((-px * lpx - py * lpy) / (std::sqrt(std::pow(px, 2) + std::pow(py, 2)) * std::sqrt(std::pow(lpx, 2) + std::pow(lpy, 2))));
        }
}
#endif // VARIABLES_H
