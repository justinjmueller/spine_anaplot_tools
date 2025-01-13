/**
 * @file vars_muon2024.h
 * @brief Header file for definitions of analysis variables specific to the
 * muon2024 analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions specific to the muon2024
 * analysis. Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author mueller@fnal.gov
 */
#ifndef VARS_NUE2024_H
#define VARS_NUE2024_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/utilities.h"
#include "include/cuts.h"
#include "include/nue2024/cuts_nue2024.h"

/**
 * @namespace vars::muon2024
 * @brief Namespace for organizing variables specific to the muon2024 analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the muon2024 analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::nue2024
{
    /**
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 0: 1e1p (contained and fiducial)
     * 1: 1e1p (not contained or not fiducial)
     * 2: 1eNp (N > 1, contained and fiducial)
     * 3: 1eNp (N > 1, not contained or fiducial)
     * 4: 1eX (not 1eNp, contained and fiducial)
     * 5: 1eX (not 1eNp, not contained or fiducial)
     * 6: Other nu
     * 7: cosmic
     * @param obj The interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    double category(const caf::SRInteractionTruthDLPProxy & obj)
    {
        double cat(7);
        if(cuts::nue2024::signal_1e1p(obj)) cat = 0;
        else if(cuts::nue2024::nonsignal_1e1p(obj)) cat = 1;
        else if(cuts::nue2024::signal_1eNp(obj)) cat = 2;
        else if(cuts::nue2024::nonsignal_1eNp(obj)) cat = 3;
        else if(cuts::nue2024::signal_1eX(obj)) cat = 4;
        else if(cuts::nue2024::nonsignal_1eX(obj)) cat = 5;
        else if(cuts::neutrino(obj)) cat = 6;
        return cat;
    }

    /**
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 0: 1e1p (contained and fiducial)
     * 1: 1e1p (not contained or not fiducial)
     * 2: 1eNp (N > 1, contained and fiducial)
     * 3: 1eNp (N > 1, not contained or fiducial)
     * 4: 1eX (not 1eNp, contained and fiducial)
     * 5: 1eX (not 1eNp, not contained or fiducial)
     * 6: Other nu
     * 7: cosmic
     * @param obj The interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    double category_topology(const caf::SRInteractionTruthDLPProxy & obj)
    {   
        double cat(7);
        if(interaction.nu_id >= 0){
            std::vector<uint32_t> counts(utilities::count_primaries(obj));
            if(counts[1] == 1 && counts[2] == 0)
                {
                    if(counts[0] == 0 && counts[3] == 0 && counts[4] == 1 && cuts::containment_cut(interaction)) cat = 0; //&& interaction.is_fiducial
                            else if(counts[0] == 0 && counts[3] == 0 && counts[4] == 1) cat = 8;
                            else if(counts[0] == 0 && counts[3] == 0 && counts[4] == 0) cat = 1;
                            else if(counts[0] == 0 && counts[3] == 0 && counts[4] > 1 && cuts::containment_cut(interaction)) cat = 2; //&& interaction.is_fiducial
                            else if(counts[0] == 0 && counts[3] == 0 && counts[4] > 1) cat = 8;
                            else if(counts[0] == 0 && counts[3] == 1 && counts[4] == 1) cat = 3;
                            else if(interaction.current_type == 0) cat = 4;
                        }
                        else if(interaction.current_type == 0 && counts[2] == 1) cat = 7;
                        else if(interaction.current_type == 0 && interaction.pdg_code == 12) cat = 4;
                        else if(interaction.current_type == 0 && interaction.pdg_code == 14) cat = 7;
                        else if(interaction.current_type == 1) cat = 5;

        }
        
        if(cuts::muon2024::signal_1e1p(obj)) cat = 0;
        else if(cuts::muon2024::nonsignal_1e1p(obj)) cat = 1;
        else if(cuts::muon2024::signal_1eNp(obj)) cat = 2;
        else if(cuts::muon2024::nonsignal_1eNp(obj)) cat = 3;
        else if(cuts::muon2024::signal_1eX(obj)) cat = 4;
        else if(cuts::muon2024::nonsignal_1eX(obj)) cat = 5;
        else if(cuts::neutrino(obj)) cat = 6;
        return cat;
    }
    /**
     * @brief Variable for the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
     */
    template<class T>
        double opening_angle(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                return std::acos(m.truth_start_dir[0] * p.truth_start_dir[0] + m.truth_start_dir[1] * p.truth_start_dir[1] + m.truth_start_dir[2] * p.truth_start_dir[2]);
            else
                return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }

    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double NuMI_angle(const T & p)
        {
            double x;
            double y;
            double z;
            double r;
            x = (31512.0380) - p.start_point[0];
            y = (3364.4912) - p.start_point[1];
            z = (73363.2532) - p.start_point[2];
            r = std::sqrt(std::pow(x, 2)+std::pow(y, 2)+std::pow(z, 2));
            x = x/r;
            y = y/r;
            z = z/r;
            return std::acos(x *particle.start_dir[0] + y *particle.start_dir[1]+z *particle.start_dir[2]);
            
        }
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    /*
    template<class T>
        double azimuthal_angle(const T & p)
        {
            if(pa.start_dir[1] >0)
                return std::acos(p.start_dir[0] / std::sqrt(std::pow(p.start_dir[0], 2) + std::pow(p.start_dir[1], 2)));
            else
                return -std::acos(p.start_dir[0] / std::sqrt(std::pow(p.start_dir[0], 2) + std::pow(p.start_dir[1], 2)));
            
        }
    template<class T>
        double polar_angle(const T & p)
        {
            return std::acos(p.start_dir[2]);
        }
    */
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double NuMI_polar_angle(const T & p)
        {
            double r;
            std::vector<double> dir_vector(3,0);
               
            dir_vector[0] = (31512.0380) + interaction.vertex[0];
            dir_vector[1] = (3364.4912) + interaction.vertex[1];
            dir_vector[2] = (73363.2532) + interaction.vertex[2];
            r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
            dir_vector[0] = dir_vector[0]/r;
            dir_vector[1] = dir_vector[1]/r;
            dir_vector[2] = dir_vector[2]/r;                                      
            return std::acos(dir_vector[2]);
        }
    
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double NuMI_azimuthal_angle(const T & p)
        {
            double r;
            std::vector<double> dir_vector(3,0);
               
            dir_vector[0] = (31512.0380) + interaction.vertex[0];
            dir_vector[1] = (3364.4912) + interaction.vertex[1];
            dir_vector[2] = (73363.2532) + interaction.vertex[2];
            r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
            dir_vector[0] = dir_vector[0]/r;
            dir_vector[1] = dir_vector[1]/r;
            dir_vector[2] = dir_vector[2]/r;                                      
            if(dir_vector[1] >0)
                return std::acos(dir_vector[0] / std::sqrt(std::pow(dir_vector[0], 2) + std::pow(dir_vector[1], 2)));
            else
                return -std::acos(dir_vector[0] / std::sqrt(std::pow(dir_vector[0], 2) + std::pow(dir_vector[1], 2)));
        }
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double NuMI_transverse_momentum(const T & obj, int pid)
        {
            //TVector3 beamdir(0, 0, 1); // BNB
        std::vector<double> dir_vector(3,0);                    
        dir_vector[0] = (31512.0380) + interaction.vertex[0];
        dir_vector[1] = (3364.4912) + interaction.vertex[1];
        dir_vector[2] = (73363.2532) + interaction.vertex[2];
        double r = std::sqrt(std::pow(dir_vector[0], 2)+std::pow(dir_vector[1], 2)+std::pow(dir_vector[2], 2));
        dir_vector[0] = dir_vector[0]/r;
        dir_vector[1] = dir_vector[1]/r;
        dir_vector[2] = dir_vector[2]/r;                                                                                                                                                                                              
        TVector3 beamdir(dir_vector[0], dir_vector[1], dir_vector[2]); // NuMI                                                                                                                                                                                                 

        // Output                                                                                                                                                                                                                                                     
        double pT0(0), pT1(0), pT2(0);

        // Loop over particles                                                                                                                                                                                                                                        
        for(auto & part : interaction.particles)
          {

            if(!part.is_primary or part.pid != 4 ) continue;

            // pT = p - pL                                                                                                                                                                                                                                            
            //    = p-(p dot beamdir) * beamdir                                                                                                                                                                                                                         
            TVector3 p;
            TVector3 pL;
            TVector3 pT;

            p.SetX(part.momentum[0]);
            p.SetY(part.momentum[1]);
            p.SetZ(part.momentum[2]);

            pL = p.Dot(beamdir) * beamdir;
            pT = p - pL;
            pT0 += pT[0];
            pT1 += pT[1];
            pT2 += pT[2];

                                                                                                                                                                                                                                          
        }
        TVector3 ppT(pT0,pT1,pT2);
        return ppT;
        }
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double delta_pt(const T & interaction)
        {
            TVector3 plT(NuMI_transverse_momentum(interaction,1));
            TVector3 ppT(NuMI_transverse_momentum(interaction,4));
            TVector3 delta_p;
            delta_p.SetX(plT[0]+ppT[0]);
            delta_p.SetY(plT[1]+ppT[1]);
            delta_p.SetZ(plT[2]+ppT[2]);
            return delta_p.Mag();
        }
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double delta_alphaT(const T & interaction)
        {
            TVector3 plT(NuMI_transverse_momentum(interaction,1));
            TVector3 ppT(NuMI_transverse_momentum(interaction,4));
            TVector3 delta_p;
            delta_p.SetX(plT[0]+ppT[0]);
            delta_p.SetY(plT[1]+ppT[1]);
            delta_p.SetZ(plT[2]+ppT[2]);
            double delta_a = std::acos(delta_p.Dot(-plT)/(plT.Mag() * delta_p.Mag()));
            return delta_a;
        }
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double delta_phiT(const T & interaction)
        {
            TVector3 plT(NuMI_transverse_momentum(interaction,1));
            TVector3 ppT(NuMI_transverse_momentum(interaction,4));
            double delta_phi = std::acos(-plT.Dot(ppT)/(plT.Mag() * ppT.Mag()));
            return delta_phi;
        }
    /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double leading_proton_softmax(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return interaction.particles[i].pid_scores[4];
        }
        /**
     * @brief Varianble for the angle with respect to NuMI beam line.
     * @details The NuMI beam angle is ~23° from the BNB beam line and the particle
     * angle is defined at the vector from the particle startpoint and 
     * (31512.0380,3364.4912,73363.2532).
     * @tparam T the type of interaction (true or reco).
     * @param p the particle to apply the variable on.
     * @return the particle angle with respect to NuMI beam.
     */
    template<class T>
        double leading_electron_softmax(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            return interaction.particles[i].pid_scores[1];
        }
    template<class T>
        double electron_transverse_momentum_mag(const T & interaction)
        {
            TVector3 plT(NuMI_transverse_momentum(interaction, 1));
            

            return plT.Mag();
        }
    template<class T>
        double proton_transverse_momentum_mag(const T & interaction)
        {
            
            TVector3 ppT(NuMI_transverse_momentum(interaction, 4));
            
            return ppT.Mag();
        }
    template<class T>
        double leading_electron_NuMI_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            double angle(NuMI_angle(interaction.particles[i]));
            return angle;
        }
    template<class T>
        double leading_proton_NuMI_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            double angle(NuMI_angle(interaction.particles[i]));
            return angle;
        }
    template<class T>
        double leading_electron_NuMI_polar_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            double angle(NuMI_polar_angle(interaction.particles[i]));
            return angle;
        }
    template<class T>
        double leading_proton_NuMI_polar_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            double angle(NuMI_polar_angle(interaction.particles[i]));
            return angle;
        }
    template<class T>
        double leading_electron_NuMI_azimuthal_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 1));
            double angle(NuMI_azimuthal_angle(interaction.particles[i]));
            return angle;
        }
    template<class T>
        double leading_proton_NuMI_azimuthal_angle(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            double angle(NuMI_azimuthal_angle(interaction.particles[i]));
            return angle;
        }
}
#endif // VARS_NUE2024_H