/**
 * @file cuts_muon2024.h
 * @brief Header file for definitions of analysis cuts specific to the muon2024
 * analysis.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the muon2024 analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_NUE2024_H
#define CUTS_NUE2024_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "include/utilities.h"

/**
 * @namespace cuts::muon2024
 * @brief Namespace for organizing cuts specific to the muon2024 analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the muon2024 analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::nue2024
{
    /**
     * @brief Apply a 1mu1p topological (final state) cut.
     * @details The interaction must have a topology matching 1mu1p as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1mu1p topology.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
        bool topological_1e1p_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[0] == 0 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] == 1;
        }

    /**
     * @brief Apply a 1muNp topological (final state) cut.
     * @details The interaction must have a topology matching 1muNp as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1muNp topology.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
        bool topological_1eNp_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[0] == 0 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] >= 1;
        }
    
    /**
     * @brief Apply a 1muX topological (final state) cut.
     * @details The interaction must have a topology matching 1muX as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1muX topology.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
        bool topological_1eX_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[1] == 1;
        }

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1mu1p
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1mu1p topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
        bool all_1e1p_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut_numi<T>(obj) && topological_1e1p_cut<T>(obj); }

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1muNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1muNp topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
        bool all_1eNp_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut_numi<T>(obj) && topological_1eNp_cut<T>(obj); }

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1muX
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1muX topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
        bool all_1eX_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut_numi<T>(obj) && topological_1eX_cut<T>(obj); }

    /**
     * @brief Apply a cut to select the 1mu1p signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1mu1p signal.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    bool signal_1e1p(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1e1p_cut(obj); }

    /**
     * @brief Apply a cut to select the 1mu1p non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1mu1p non-signal 
     * (1mu1p topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    bool nonsignal_1e1p(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1e1p_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muNp signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muNp signal.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    bool signal_1eNp(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1eNp_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muNp non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muNp non-signal
     * (1muNp topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    bool nonsignal_1eNp(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1eNp_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muX signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muX signal
     * definition.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    bool signal_1eX(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1eX_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muX non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muX non-signal
     * (1muX topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    bool nonsignal_1eX(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1eX_cut(obj); }
}
#endif // CUTS_MUON2024_H