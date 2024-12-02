/**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_H
#define CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "utilities.h"

/**
 * @namespace cuts
 * @brief Namespace for organizing generic cuts which act on interactions.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. The function should
 * be templated on the type of interaction object if the cut is intended to be
 * used on both true and reconstructed interactions.
 */
namespace cuts
{   
    /**
     * @brief Apply a cut on the validity of the flash match.
     * @details A "valid" flash match is defined as a flash-interaction
     * association with a flash time that is not NaN and a flash match
     * status of 1. The upstream flash matching algorithm (OpT0Finder) has a
     * flash filter that restricts candidate flashes to near the beam window,
     * which means that the majority of cosmogenic interactions are not
     * flash matched. If no flash match is found, the flash time is NaN. This
     * cut is intended to be applied as a preselection cut to reduce comparisons
     * to NaN values, which tend to be noisy on stderr.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
     */
    template<class T>
        bool valid_flashmatch(const T & obj) { return !std::isnan(obj.flash_time) && obj.fmatched == 1; }

    /**
     * @brief Apply no cut; all interactions passed.
     * @details This is a placeholder function for a cut which does not apply
     * any selection criteria. It is intended to be used in cases where a cut
     * function is required, but no selection is desired.
     * @tparam T the type of object (true or reco).
     * @param obj the interaction to select on.
     * @return true (always).
     */
    template<class T>
        bool no_cut(const T & obj) { return true; }

    /**
     * @brief Apply a cut to select neutrinos.
     * @details This function applies a cut to select neutrinos. This cut
     * makes use of the is_neutrino flag in the true interaction object and is
     * intended to be used to identify signal neutrinos.
     * @param obj the interaction to select on.
     * @return true if the interaction is a neutrino.
     * @note This cut is intended to be used for identifying neutrinos in
     * truth, which is useful for making signal definitions.
     */
    bool neutrino(const caf::SRInteractionTruthDLPProxy & obj) { return obj.is_neutrino; }

    /**
     * @brief Apply a cut to select cosmogenic interactions.
     * @details This function applies a cut to select cosmogenic interactions.
     * This cut makes use of the is_neutrino flag in the true interaction
     * object and is intended to be used to identify cosmogenic interactions.
     * @param obj the interaction to select on.
     * @return true if the interaction is a cosmogenic interaction.
     * @note This cut is intended to be used for identifying cosmogenic
     * interactions in truth, which is useful for making background definitions.
     */
    bool cosmic(const caf::SRInteractionTruthDLPProxy & obj) { return !obj.is_neutrino; }

    /**
     * @brief Apply a fiducial volume cut; the interaction vertex must be
     * reconstructed within the fiducial volume.
     * @details The fiducial volume cut is applied on the reconstructed
     * interaction vertex upstream in SPINE. The fiducial volume is defined
     * (in a SPINE post-processor) as a 25 cm border around the x and y
     * detector faces, a 50 cm border around the downstream (+) z face, and a
     * 30 cm border around the upstream (-) z face. The fiducial volume is
     * intended to reduce the impact of detector edge effects on the analysis.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
        bool fiducial_cut(const T & obj)
        {
            return obj.is_fiducial && !(obj.vertex[0] > 210.215 && obj.vertex[1] > 60 && (obj.vertex[2] > 290 && obj.vertex[2] < 390));
        }
    
    /**
     * @brief Apply a containment cut on the entire interaction.
     * @details The containment cut is applied on the entire interaction. The
     * interaction is considered contained if all particles and all spacepoints
     * are contained within 5cm of the detector edges (configured in a SPINE 
     * post-processor). Additionally, no spacepoints are allowed to be
     * reconstructed in a TPC that did not create it. This is an unphysical
     * condition that can occur when a cosmic muon is moved according to an
     * assumed t0 that is very out-of-time.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
        bool containment_cut(const T & obj) { return obj.is_contained; }

    /**
     * @brief Apply a flash time cut on the interaction.
     * @details The flash time cut is applied on the interaction. The flash time
     * is required to be within the beam window, which is defined as 0 to 1.6
     * microseconds. This cut is intended to reduce the impact of cosmogenic
     * interactions on analyses.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     * @note This cut is intended to be used for BNB analyses.
     * @note The cut window has been widened to reconcile the beam window as
     * observed in data and simulation.
     */
    template<class T>
        bool flash_cut_bnb(const T & obj)
        {
            if(!valid_flashmatch(obj))
                return false;
            else
                return (obj.flash_time >= -0.5) && (obj.flash_time <= 1.6);
        }
    
    /**
     * @brief Apply a flash time cut on the interaction.
     * @details The flash time cut is applied on the interaction. The flash time
     * is required to be within the beam window, which is defined as 0 to 9.6
     * microseconds. This cut is intended to reduce the impact of cosmogenic
     * interactions on analyses.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     * @note This cut is intended to be used for NuMI analyses.
     */
    template<class T>
        bool flash_cut_numi(const T & obj)
        {
            if(!valid_flashmatch(obj))
                return false;
            else
                return (obj.flash_time >= 0) && (obj.flash_time <= 9.6);
        }

    /**
     * @brief Apply a fiducial and containment cut (logical "and" of both).
     * @details This function applies a fiducial and containment cut on the
     * interaction using the logical "and" of both previously defined cuts.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
     */
    template<class T>
        bool fiducial_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj); }

    /**
     * @brief Apply a fiducial, containment, and flash time cut (BNB) (logical
     * "and" of each).
     * @details This function applies a fiducial, containment, and flash time cut
     * (BNB) on the interaction using the logical "and" of each previously
     * defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * flash time cut.
     * @note This cut is intended to be used for BNB analyses.
     */
    template<class T>
        bool fiducial_containment_flash_cut_bnb(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut_bnb<T>(obj); }

    /**
     * @brief Apply a fiducial, containment, and flash time cut (NuMI) (logical
     * "and" of each).
     * @details This function applies a fiducial, containment, and flash time cut
     * (NuMI) on the interaction using the logical "and" of each previously
     * defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * flash time cut.
     * @note This cut is intended to be used for NuMI analyses.
     */
    template<class T>
        bool fiducial_containment_flash_cut_numi(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut_numi<T>(obj); }

    /**
     * @brief Apply a fiducial and neutrino cut (logical "and" of each).
     * @details This function applies a fiducial and neutrino cut on the
     * interaction using the logical "and" of each previously defined cut.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial and neutrino cut.
     * @note This cut is intended to be used to select signal interactions
     * (neutrinos) that occur within the fiducial volume.
     */
    bool fiducial_neutrino_cut(const caf::SRInteractionTruthDLPProxy & obj) { return fiducial_cut(obj) && neutrino(obj); }

    /**
     * @brief Apply a fiducial, containment, and neutrino cut (logical "and" of
     * each).
     * @details This function applies a fiducial, containment, and neutrino cut
     * on the interaction using the logical "and" of each previously defined cut.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * neutrino cut.
     * @note This cut is intended to be used to select signal interactions
     * (neutrinos) that are fully contained within the detector and that occur
     * within the fiducial volume.
     */
    bool fiducial_containment_neutrino_cut(const caf::SRInteractionTruthDLPProxy & obj) { return fiducial_cut(obj) && containment_cut(obj) && neutrino(obj); }

    /**
     * @brief Apply a cut to select interactions with no primary charged pions.
     * @details This function applies a cut to select interactions with no
     * primary charged pions in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has zero charged pions.
     */
    template<class T>
        bool no_charged_pions(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[3] == 0;
        }

    /**
     * @brief Apply a cut to select interactions with no primary showers.
     * @details This function applies a cut to select interactions with no
     * primary showers in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has zero showers.
     */
    template<class T>
        bool no_showers(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[0] == 0 && c[1] == 0;
        }

    /**
     * @brief Apply a cut to select interactions with a single primary muon.
     * @details This function applies a cut to select interactions with a
     * single primary muon in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a single primary muon.
     */
    template<class T>
        bool has_single_muon(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[2] == 1;
        }

    /**
     * @brief Apply a cut to select interactions with a single primary proton.
     * @details This function applies a cut to select interactions with a
     * single primary proton in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a single primary proton.
     */
    template<class T>
        bool has_single_proton(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[4] == 1;
        }

    /**
     * @brief Apply a cut to select interactions with nonzero primary protons.
     * @details This function applies a cut to select interactions with
     * nonzero primary protons in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has nonzero primary protons.
     */
    template<class T>
        bool has_nonzero_protons(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return c[4] > 0;
        }
}
#endif
