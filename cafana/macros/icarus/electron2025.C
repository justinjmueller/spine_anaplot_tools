/**
 * @file electron2025.C
 * @brief The main analysis macro for the electron2025 benchmarking on ICARUS Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author jzettle@fnal.gov
*/

/**
 * @brief Block of preprocessor definitions for the analysis.
 * @details This block of preprocessor definitions is used to configure the
 * analysis. The definitions control the behavior of the analysis, such as
 * which beam is used, which cuts are applied, and which trees are created.
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PIDFUNC pvars::custom_pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/electron2025/variables_electron2025.h"
#include "include/cuts.h"
#include "include/electron2025/cuts_electron2025.h"
#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

void electron2025()
{
    ana::Analysis analysis("electron2025_rev1_icarus_testbenchmark");

    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/spineprod/mpv_boostedshower/mpv_boosted_ee.flat.root");
    analysis.AddLoader("mc", &mc, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */

    #define CUT cuts::electron2025::all_1shower_cut

    std::map<std::string, ana::SpillMultiVar> vars_selected_ee;
    vars_selected_ee.insert({"category", SpineVar<TTYPE, RTYPE>(&vars::electron2025::category, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"category_reco", SpineVar<RTYPE, RTYPE>(&vars::electron2025::category_templated, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::electron2025::visible_energy_ee, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"true_edep", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"nshowers", SpineVar<RTYPE,RTYPE>(&vars::electron2025::nshowers, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"nelectrons", SpineVar<RTYPE,RTYPE>(&vars::electron2025::nelectrons, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"nphotons", SpineVar<RTYPE,RTYPE>(&vars::electron2025::nphotons, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"ntracks", SpineVar<RTYPE,RTYPE>(&vars::electron2025::ntracks, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"opening_angle", SpineVar<TTYPE,RTYPE>(&vars::electron2025::opening_angle_ee, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"opening_angle_reco", SpineVar<RTYPE,RTYPE>(&vars::electron2025::opening_angle_ee, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"leading_shower_energy", SpineVar<TTYPE,RTYPE>(&vars::electron2025::leading_shower_energy, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"subleading_shower_energy", SpineVar<TTYPE,RTYPE>(&vars::electron2025::subleading_shower_energy, &CUT, &cuts::no_cut)}); 
    vars_selected_ee.insert({"leading_shower_energy_reco", SpineVar<RTYPE,RTYPE>(&vars::electron2025::leading_shower_energy, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"subleading_shower_energy_reco", SpineVar<RTYPE,RTYPE>(&vars::electron2025::subleading_shower_energy, &CUT, &cuts::no_cut)}); 
    vars_selected_ee.insert({"invariant_mass", SpineVar<RTYPE,RTYPE>(&vars::electron2025::invariant_mass, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"invariant_mass_true", SpineVar<TTYPE,RTYPE>(&vars::electron2025::invariant_mass, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"energy_asymmetry", SpineVar<RTYPE,RTYPE>(&vars::electron2025::energy_asymmetry, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"energy_asymmetry_true", SpineVar<TTYPE,RTYPE>(&vars::electron2025::energy_asymmetry, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"leading_electron_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_electron_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_electron_electron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::electron_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_electron_photon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::photon_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_px", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::px, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_py", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::py, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_pz", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pz, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_true_px", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::px, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_true_py", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::py, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_true_pz", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::pz, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"IoU", SpineVar<RTYPEP,RTYPE,RTYPE>(&vars::electron2025::iou, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_px_dir", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::px_dir, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_py_dir", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::py_dir, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_pz_dir", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pz_dir, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_true_px_dir", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::px_dir, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_true_py_dir", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::py_dir, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"leading_true_pz_dir", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::pz_dir, &CUT, &cuts::no_cut, &utilities::electron2025::leading_shower_index)});
    vars_selected_ee.insert({"subleading_electron_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_electron_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_electron_electron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::electron_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_electron_photon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::photon_softmax, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_px", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::px, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_py", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::py, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_pz", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pz, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_true_px", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::px, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_true_py", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::py, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_true_pz", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::pz, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_px_dir", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::px_dir, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_py_dir", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::py_dir, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_pz_dir", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pz_dir, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_true_px_dir", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::px_dir, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_true_py_dir", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::py_dir, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"subleading_true_pz_dir", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::pz_dir, &CUT, &cuts::no_cut, &utilities::electron2025::subleading_shower_index)});
    vars_selected_ee.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &cuts::no_cut)});
    vars_selected_ee.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &cuts::no_cut)});  

    analysis.AddTree("selectedEvents", vars_selected_ee, false);

    std::map<std::string, ana::SpillMultiVar> vars_pid_electron;
    auto pid = [](const RTYPEP & p) -> double { return p.pid; };
    auto true_pid = [](const TTYPEP & p) -> double { return p.pid; };
    auto primary_shower = [](const TTYPEP & p) { return (p.pid == 1 || p.pid == 0) && p.is_primary; };
    vars_pid_electron.insert({"pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(pid, primary_shower, &cuts::no_cut)});
    vars_pid_electron.insert({"true_pid", SpineVar<TTYPEP,TTYPEP,TTYPE>(true_pid, primary_shower, &cuts::no_cut)});
    vars_pid_electron.insert({"primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), primary_shower, &cuts::no_cut)});
    vars_pid_electron.insert({"IoU", SpineVar<RTYPEP,TTYPEP,TTYPE>(&vars::electron2025::iou, primary_shower, &cuts::no_cut)});
    vars_selected_ee.insert({"px_dir", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::px_dir, primary_shower, &cuts::no_cut)});
    vars_selected_ee.insert({"py_dir", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::py_dir, primary_shower, &cuts::no_cut)});
    vars_selected_ee.insert({"pz_dir", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::pz_dir, primary_shower, &cuts::no_cut)});
    vars_selected_ee.insert({"true_px_dir", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::px_dir, primary_shower, &cuts::no_cut)});
    vars_selected_ee.insert({"true_py_dir", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::py_dir, primary_shower, &cuts::no_cut)});
    vars_selected_ee.insert({"true_pz_dir", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::pz_dir, primary_shower, &cuts::no_cut)});
    analysis.AddTree("pid_electron", vars_pid_electron, false);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
     
    analysis.Go();

}
