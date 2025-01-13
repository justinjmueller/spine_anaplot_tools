/**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select particles. Each cut is implemented as a function which takes an
 * particle object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections/interaction cuts.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_CUTS_H
#define PARTICLE_CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "particle_variables.h"

/**
 * @namespace pcuts
 * @brief Namespace for organizing generic cuts which act on particles.
 * @details This namespace is intended to be used for organizing cuts which act
 * on particles. Each cut is implemented as a function which takes a
 * particle object as an argument and returns a boolean. The function should
 * be templated on the type of particle object if the cut is intended to be
 * used on both true and reconstructed particles.
 */
namespace pcuts
{
    /**
     * @brief Check if the particle meets final state signal requirements.
     * @details must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(p.is_primary)
            {
                double energy(pvars::energy(p));
                if((p.pid == 2 && energy > 143.425) || (p.pid != 2 && p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 50))
                    passes = true;
            }
            return passes;
        }
    /**
     * @brief Check if the particle meets final state signal requirements.
     * @details must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
        bool final_state_signal_nue(const T & p)
        {
            bool passes(false);
            if(p.is_primary)
            {
                double energy(pvars::energy(p));
                if((p.pid == 1 && energy > 70) || (p.pid != 1 && p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 40))
                    passes = true;
            }
            return passes;
        }
}
#endif // PARTICLE_CUTS_H