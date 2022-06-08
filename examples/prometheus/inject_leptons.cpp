#include <phys-services/CrossSection.h>

#include <LeptonInjector/Controller.h>
#include <LeptonInjector/Particle.h>
#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Constants.h>
#include <LeptonInjector/Weighter.h>
#include <earthmodel-service/Geometry.h>
#include <earthmodel-service/EulerQuaternionConversions.h>
#include <earthmodel-service/Placement.h>

#include <string>
#include <iomanip>
#include <memory>
#include <chrono>
#include <ctime>
#include "argagg.hpp"
#include "date.h"

struct ExitStatus {
    ExitStatus(int status) : status(status) {}
    int status;
};

enum class InjectionMode {
    Ranged,
    Volume
};

enum class InteractionType {
    CC,
    NC,
    GR_E,
    GR_MU,
    GR_TAU,
    GR_HAD
};

template <class Precision>
std::string getISOCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    return date::format("%FT%TZ", date::floor<Precision>(now));
}

std::string interaction_to_string(InteractionType type) {
    std::map<InteractionType, std::string> const interactions = {
        {InteractionType::CC, "cc"},
        {InteractionType::NC, "nc"},
        {InteractionType::GR_E, "gr_e"},
        {InteractionType::GR_MU, "gr_mu"},
        {InteractionType::GR_TAU, "gr_tau"}
    };
    return interactions.at(type);
}

InteractionType string_to_interaction(std::string s) {
    std::map<std::string, InteractionType> const interactions = {
        {"cc", InteractionType::CC},
        {"nc", InteractionType::NC},
        {"gr_e", InteractionType::GR_E},
        {"gr_mu", InteractionType::GR_MU},
        {"gr_tau", InteractionType::GR_TAU}
    };
    return interactions.at(s);
}

std::string nu_to_string(LeptonInjector::Particle::ParticleType type) {
    std::map<LeptonInjector::Particle::ParticleType, std::string> const interactions = {
        {LeptonInjector::Particle::ParticleType::NuE, "nue"},
        {LeptonInjector::Particle::ParticleType::NuMu, "numu"},
        {LeptonInjector::Particle::ParticleType::NuTau, "nutau"},
        {LeptonInjector::Particle::ParticleType::NuEBar, "nuebar"},
        {LeptonInjector::Particle::ParticleType::NuMuBar, "numubar"},
        {LeptonInjector::Particle::ParticleType::NuTauBar, "nutaubar"}
    };
    return interactions.at(type);
}

LeptonInjector::Particle::ParticleType string_to_nu(std::string s) {
    std::map<std::string, LeptonInjector::Particle::ParticleType> const interactions = {
        {"nue", LeptonInjector::Particle::ParticleType::NuE},
        {"numu", LeptonInjector::Particle::ParticleType::NuMu},
        {"nutau", LeptonInjector::Particle::ParticleType::NuTau},
        {"nuebar", LeptonInjector::Particle::ParticleType::NuEBar},
        {"numubar", LeptonInjector::Particle::ParticleType::NuMuBar},
        {"nutaubar", LeptonInjector::Particle::ParticleType::NuTauBar}
    };
    return interactions.at(s);
}

std::vector<LeptonInjector::Particle::ParticleType> gen_nu_types() {
    std::vector<LeptonInjector::Particle::ParticleType> nu_types = {
        LeptonInjector::Particle::ParticleType::NuE,
        LeptonInjector::Particle::ParticleType::NuMu,
        LeptonInjector::Particle::ParticleType::NuTau,
        LeptonInjector::Particle::ParticleType::NuEBar,
        LeptonInjector::Particle::ParticleType::NuMuBar,
        LeptonInjector::Particle::ParticleType::NuTauBar
    };
    return nu_types;
}

std::string get_total_xs_id(InteractionType type, LeptonInjector::Particle::ParticleType nu) {
    std::string nu_string = nu_to_string(nu);
    std::string interaction_string = interaction_to_string(type);
    std::string id = nu_string + "_" + interaction_string + "_total_xs";
    return id;
}

std::string get_diff_xs_id(InteractionType type, LeptonInjector::Particle::ParticleType nu) {
    std::string nu_string = nu_to_string(nu);
    std::string interaction_string = interaction_to_string(type);
    std::string id = nu_string + "_" + interaction_string + "_diff_xs";
    return id;
}

std::string get_total_xs_fname(std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> cross_sections, InteractionType type, LeptonInjector::Particle::ParticleType nu) {
    return cross_sections[{type, nu}].second;
}

std::string get_diff_xs_fname(std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> cross_sections, InteractionType type, LeptonInjector::Particle::ParticleType nu) {
    return cross_sections[{type, nu}].first;
}


std::tuple<argagg::parser, argagg::parser_results> parse_arguments(int argc, char ** argv) {
    argagg::parser argparser = {{
        {
            "help", {"-h", "--help"},
            "Print help and exit", 0,
        },
        {
            "output", {"--output"},
            "Path and prefix for the output. Used for h5 and lic.", 1,
        },
        {
            "xs_path", {"--cross-section-path"},
            "Path to the cross sections.", 1,
        },
        {
            "earth_model_path", {"--earth-path", "--earth-model-path"},
            "Path to the earth model resources.", 1,
        },
        {
            "earth_model", {"--earth-model", "--earth-density", "--earth-density-model"},
            "Name of the earth density model.", 1
        },
        {
            "materials_model", {"--materials-model", "--materials-density", "--materials-density-model"},
            "Name of the materials density model.", 1,
        },
        {
            "minE", {"--min-energy", "--minE"},
            "The minimum injection energy in GeV", 1,
        },
        {
            "maxE", {"--max-energy", "--maxE"},
            "The maximum injection energy in GeV", 1,
        },
        {
            "gamma", {"--gamma"},
            "The injection spectral index", 1,
        },
        {
            "minZenith", {"--min-zenith", "--minZen"},
            "The minimum injection zenith angle in degrees", 1,
        },
        {
            "maxZenith", {"--max-zenith", "--maxZen"},
            "The maximum injection zenith angle in degrees", 1,
        },
        {
            "minAzimuth", {"--min-azimuth", "--minAzi"},
            "The minimum injection azimuth angle in degrees", 1,
        },
        {
            "maxAzimuth", {"--max-azimuth", "--maxAzi"},
            "The maximum injection azimuth angle in degrees", 1,
        },
        {
            "ranged_radius", {"--ranged-radius"},
            "The ranged mode radius of injection in meters", 1,
        },
        {
            "ranged_length", {"--ranged-length", "--end-cap", "--endcap", "--endcap-length", "--end-cap-length"},
            "The ranged mode end cap length in meters", 1,
        },
        {
            "volume_radius", {"--volume-radius", "--cylinder-radius"},
            "The volume mode injection cylinder radius in meters", 1,
        },
        {
            "volume_height", {"--volume-height", "--cylinder-height"},
            "The volume mode injection cylinder height in meters", 1,
        },
        {
            "n_events", {"--n-events", "--num-events", "--number-of-events"},
            "Number of ranged events per neutrino type per interaction", 1,
        },
        {
            "cc", {"--cc"},
            "Inject CC events?", 0,
        },
        {
            "nc", {"--nc"},
            "Inject NC events?", 0,
        },
        {
            "gr", {"--gr"},
            "Inject GR events?", 0,
        },
        {
            "gr_e", {"--gr-e"},
            "Inject GR e events?", 0,
        },
        {
            "gr_mu", {"--gr-mu"},
            "Inject GR mu events?", 0,
        },
        {
            "gr_tau", {"--gr-tau"},
            "Inject GR tau events?", 0,
        },
        {
            "gr_had", {"--gr-had"},
            "Inject GR had events?", 0,
        },
        {
            "nue", {"--nue"},
            "Inject nue events?", 0,
        },
        {
            "nuebar", {"--nuebar"},
            "Inject nuebar events?", 0,
        },
        {
            "numu", {"--numu"},
            "Inject numu events?", 0,
        },
        {
            "numubar", {"--numubar"},
            "Inject numubar events?", 0,
        },
        {
            "nutau", {"--nutau"},
            "Inject nutau events?", 0,
        },
        {
            "nutaubar", {"--nutaubar"},
            "Inject nutaubar events?", 0,
        },
        {
            "ranged_mode", {"--ranged"},
            "Use ranged mode", 0,
        },
        {
            "volume_mode", {"--volume"},
            "Use volume mode", 0,
        },
        {
            "ranged_and_volume_mode", {"--ranged-volume"},
            "Use ranged and volume mode", 0,
        },
        {
            "automatic_injection_mode", {"--automatic-injection-mode"},
            "Use automatically select ranged and/or volume mode", 0,
        },
        {
            "nu_cc_diff_xs", {"--nu-cc-dsdxdy"},
            "Path to the nu CC differential cross section file", 1,
        },
        {
            "nu_cc_total_xs", {"--nu-cc-sigma"},
            "Path to the nu CC total cross section file", 1,
        },
        {
            "nu_nc_diff_xs", {"--nu-nc-dsdxdy"},
            "Path to the nu NC differential cross section file", 1,
        },
        {
            "nu_nc_total_xs", {"--nu-nc-sigma"},
            "Path to the nu NC total cross section file", 1,
        },
        {
            "nubar_cc_diff_xs", {"--nubar-cc-dsdxdy"},
            "Path to the nubar CC differential cross section file", 1,
        },
        {
            "nubar_cc_total_xs", {"--nubar-cc-sigma"},
            "Path to the nubar CC total cross section file", 1,
        },
        {
            "nubar_nc_diff_xs", {"--nubar-nc-dsdxdy"},
            "Path to the nubar NC differential cross section file", 1,
        },
        {
            "nubar_nc_total_xs", {"--nubar-nc-sigma"},
            "Path to the nubar NC total cross section file", 1,
        },
        {
            "nue_cc_diff_xs", {"--nue-cc-dsdxdy"},
            "Path to the nue CC differential cross section file", 1,
        },
        {
            "nue_cc_total_xs", {"--nue-cc-sigma"},
            "Path to the nue CC total cross section file", 1,
        },
        {
            "nuebar_cc_diff_xs", {"--nuebar-cc-dsdxdy"},
            "Path to the nuebar CC differential cross section file", 1,
        },
        {
            "nuebar_cc_total_xs", {"--nuebar-cc-sigma"},
            "Path to the nuebar CC total cross section file", 1,
        },
        {
            "numu_cc_diff_xs", {"--numu-cc-dsdxdy"},
            "Path to the numu CC differential cross section file", 1,
        },
        {
            "numu_cc_total_xs", {"--numu-cc-sigma"},
            "Path to the numu CC total cross section file", 1,
        },
        {
            "numubar_cc_diff_xs", {"--numubar-cc-dsdxdy"},
            "Path to the numubar CC differential cross section file", 1,
        },
        {
            "numubar_cc_total_xs", {"--numubar-cc-sigma"},
            "Path to the numubar CC total cross section file", 1,
        },
        {
            "nutau_cc_diff_xs", {"--nutau-cc-dsdxdy"},
            "Path to the nutau CC differential cross section file", 1,
        },
        {
            "nutau_cc_total_xs", {"--nutau-cc-sigma"},
            "Path to the nutau CC total cross section file", 1,
        },
        {
            "nutaubar_cc_diff_xs", {"--nutaubar-cc-dsdxdy"},
            "Path to the nutaubar CC differential cross section file", 1,
        },
        {
            "nutaubar_cc_total_xs", {"--nutaubar-cc-sigma"},
            "Path to the nutaubar CC total cross section file", 1,
        },
        {
            "nue_nc_diff_xs", {"--nue-nc-dsdxdy"},
            "Path to the nue NC differential cross section file", 1,
        },
        {
            "nue_nc_total_xs", {"--nue-nc-sigma"},
            "Path to the nue NC total cross section file", 1,
        },
        {
            "nuebar_nc_diff_xs", {"--nuebar-nc-dsdxdy"},
            "Path to the nuebar NC differential cross section file", 1,
        },
        {
            "nuebar_nc_total_xs", {"--nuebar-nc-sigma"},
            "Path to the nuebar NC total cross section file", 1,
        },
        {
            "numu_nc_diff_xs", {"--numu-nc-dsdxdy"},
            "Path to the numu NC differential cross section file", 1,
        },
        {
            "numu_nc_total_xs", {"--numu-nc-sigma"},
            "Path to the numu NC total cross section file", 1,
        },
        {
            "numubar_nc_diff_xs", {"--numubar-nc-dsdxdy"},
            "Path to the numubar NC differential cross section file", 1,
        },
        {
            "numubar_nc_total_xs", {"--numubar-nc-sigma"},
            "Path to the numubar NC total cross section file", 1,
        },
        {
            "nutau_nc_diff_xs", {"--nutau-nc-dsdxdy"},
            "Path to the nutau NC differential cross section file", 1,
        },
        {
            "nutau_nc_total_xs", {"--nutau-nc-sigma"},
            "Path to the nutau NC total cross section file", 1,
        },
        {
            "nutaubar_nc_diff_xs", {"--nutaubar-nc-dsdxdy"},
            "Path to the nutaubar NC differential cross section file", 1,
        },
        {
            "nutaubar_nc_total_xs", {"--nutaubar-nc-sigma"},
            "Path to the nutaubar NC total cross section file", 1,
        },
        {
            "nuebar_gr_diff_xs", {"--nuebar-gr-dsdy"},
            "Path to the nuebar GR differential cross section file", 1,
        },
        {
            "nuebar_gr_total_xs", {"--nuebar-gr-sigma"},
            "Path to the nuebar GR total cross section file", 1,
        },
        {
            "nuebar_gr_e_diff_xs", {"--nuebar-gr-e-dsdy"},
            "Path to the nuebar GR e differential cross section file", 1,
        },
        {
            "nuebar_gr_e_total_xs", {"--nuebar-gr-e-sigma"},
            "Path to the nuebar GR e total cross section file", 1,
        },
        {
            "nuebar_gr_mu_diff_xs", {"--nuebar-gr-mu-dsdy"},
            "Path to the nuebar GR mu differential cross section file", 1,
        },
        {
            "nuebar_gr_mu_total_xs", {"--nuebar-gr-mu-sigma"},
            "Path to the nuebar GR mu total cross section file", 1,
        },
        {
            "nuebar_gr_tau_diff_xs", {"--nuebar-gr-tau-dsdy"},
            "Path to the nuebar GR tau differential cross section file", 1,
        },
        {
            "nuebar_gr_tau_total_xs", {"--nuebar-gr-tau-sigma"},
            "Path to the nuebar GR tau total cross section file", 1,
        },
        {
            "nuebar_gr_had_diff_xs", {"--nuebar-gr-had-dsdy"},
            "Path to the nuebar GR had differential cross section file", 1,
        },
        {
            "nuebar_gr_had_total_xs", {"--nuebar-gr-had-sigma"},
            "Path to the nuebar GR had total cross section file", 1,
        },
        {
            "seed", {"--seed"},
            "Seed for the injector", 1,
        },
    }};

    // Get the argagg usage string
    std::ostringstream usage;
    usage << argv[0] << std::endl;

    // Print usage information when we encounter a parsing error
    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception& e) {
        argagg::fmt_ostream fmt(std::cerr);
        fmt << usage.str() << argparser << std::endl
            << "Encountered exception while parsing arguments: " << e.what()
            << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }

    // Print the usage information when we help is requested
    if(args["help"]) {
        std::cerr << argparser;
        throw ExitStatus(EXIT_SUCCESS);
    }

    return std::make_tuple(argparser, args);
}

double get_seed(argagg::parser & argparser, argagg::parser_results & args) {
    // Grab the seed
    double seed = args["seed"].as<int>(-1);
    // Require seed to be set and positive
    if(seed < 0) {
        std::cerr << "--seed requires positive integer!" << std::endl << argparser;
        throw ExitStatus(EXIT_FAILURE);
    }
    return seed;
}

std::string get_output(argagg::parser & argparser, argagg::parser_results & args) {
    // Require output prefix
    if(not args["output"]) {
        std::cerr << "--output required!" << std::endl << argparser;
        throw ExitStatus(EXIT_FAILURE);
    }
    std::string output = args["output"].as<std::string>("./injected/output_DUNE");
    return output;
}

std::string get_path(argagg::parser_results & args) {
    // Require LI path
    std::string path;
    if(args["earth_model_path"]) {
        path = args["earth_model_path"].as<std::string>();
    }
    else { // Try to grab it from the environment variables instead
        if (const char* env_p = getenv("GOLEMSOURCEPATH")){
            path = std::string( env_p ) + "/LeptonInjector/";
            path += "resources/earthparams";
        }
        else {
            std::cerr << "WARNING no lepton injector path specified and GOLEMSOURCEPATH not set! Assuming earth model information is in ./resources/earthparams/" << std::endl;
            path = "./";
        }
    }
    if(path[path.size() - 1] != '/')
        path += "/";

    return path;
}

std::vector<InteractionType> get_allowed_interactions(argagg::parser_results const & args, LeptonInjector::Particle::ParticleType nu, std::vector<InteractionType> interactions_to_consider) {
    std::vector<InteractionType> allowed;
    std::string nu_string = nu_to_string(nu);
    if(args[nu_string]) {
        for(auto const & type : interactions_to_consider) {
            std::string type_string = interaction_to_string(type);
            bool is_glashow = type_string.rfind("gr", 0) == 0;
            if(is_glashow and nu_string != "nuebar")
                continue;
            if(args[type_string] or (is_glashow and args["gr"])) {
                allowed.push_back(type);
            }
        }
    }
    return allowed;
}

std::vector<std::string> gen_possible_interactions() {
    return {"cc", "nc", "gr_e", "gr_mu", "gr_tau", "gr_had"};
}

std::vector<std::string> get_interactions(argagg::parser_results & args) {
    std::vector<std::string> interactions;
    std::vector<std::string> possible_interactions = {"cc", "nc"};
    for(unsigned int i=0; i<possible_interactions.size(); ++i) {
        if(args[possible_interactions[i]]) {
            interactions.push_back(possible_interactions[i]);
        }
    }

    if(args["gr"]) {
        interactions.push_back("gr_e");
        interactions.push_back("gr_mu");
        interactions.push_back("gr_tau");
        interactions.push_back("gr_had");
    } else {
        std::vector<std::string> possible_interaction_suffixes = {"e", "mu", "tau", "had"};
        for(unsigned int i=0; i<possible_interaction_suffixes.size(); ++i) {
            if(args["gr_" + possible_interaction_suffixes[i]]) {
                interactions.push_back("gr_" + possible_interaction_suffixes[i]);
            }
        }
    }
    return interactions;
}

std::vector<std::string> gen_possible_neutrinos() {
    return {"nue", "nuebar", "numu", "numubar", "nutau", "nutaubar"};
}

std::vector<std::string> get_neutrinos(argagg::parser_results & args) {
    std::vector<std::string> neutrinos;
    std::vector<std::string> possible_neutrinos = gen_possible_neutrinos();
    for(unsigned int i=0; i<possible_neutrinos.size(); ++i) {
        if(args[possible_neutrinos[i]]) {
            neutrinos.push_back(possible_neutrinos[i]);
        }
    }
    return neutrinos;
}

std::map<std::string, std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>> gen_secondaries() {
    std::map<std::string, std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>> secondaries = {
        {"nue_cc", {LeptonInjector::Particle::EMinus, LeptonInjector::Particle::Hadrons}},
        {"nuebar_cc", {LeptonInjector::Particle::EPlus, LeptonInjector::Particle::Hadrons}},
        {"numu_cc", {LeptonInjector::Particle::MuMinus, LeptonInjector::Particle::Hadrons}},
        {"numubar_cc", {LeptonInjector::Particle::MuPlus, LeptonInjector::Particle::Hadrons}},
        {"nutau_cc", {LeptonInjector::Particle::TauMinus, LeptonInjector::Particle::Hadrons}},
        {"nutaubar_cc", {LeptonInjector::Particle::TauPlus, LeptonInjector::Particle::Hadrons}},
        {"nue_nc", {LeptonInjector::Particle::NuE, LeptonInjector::Particle::Hadrons}},
        {"nuebar_nc", {LeptonInjector::Particle::NuEBar, LeptonInjector::Particle::Hadrons}},
        {"numu_nc", {LeptonInjector::Particle::NuMu, LeptonInjector::Particle::Hadrons}},
        {"numubar_nc", {LeptonInjector::Particle::NuMuBar, LeptonInjector::Particle::Hadrons}},
        {"nutau_nc", {LeptonInjector::Particle::NuTau, LeptonInjector::Particle::Hadrons}},
        {"nutaubar_nc", {LeptonInjector::Particle::NuTauBar, LeptonInjector::Particle::Hadrons}},
        {"nuebar_gr_e", {LeptonInjector::Particle::NuEBar, LeptonInjector::Particle::EMinus}},
        {"nuebar_gr_mu", {LeptonInjector::Particle::NuMuBar, LeptonInjector::Particle::MuMinus}},
        {"nuebar_gr_tau", {LeptonInjector::Particle::NuTauBar, LeptonInjector::Particle::TauMinus}},
        {"nuebar_gr_had", {LeptonInjector::Particle::Hadrons, LeptonInjector::Particle::Hadrons}},
    };
    return secondaries;
}

std::string get_xs_base(argagg::parser_results & args) {
    std::string xs_base;
    if(args["xs_path"]) {
        xs_base = args["xs_path"].as<std::string>();
    }
    else {
        if (const char* env_p = getenv("GOLEMSOURCEPATH")){
            xs_base = std::string( env_p ) + "/DUNEAtmo/";
        }
        else {
            std::cerr << "WARNING no DUNEAtmo path specified and GOLEMSOURCEPATH not set! Assuming the default cross section information is in ../cross_sections/" << std::endl;
            xs_base = "../";
        }
    }
    if(xs_base[-1] != '/')
        xs_base += "/";

    return xs_base;
}

std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> gen_default_cross_sections(std::string xs_base) {

    std::string nu_cc_diff_xs = xs_base + "dsdxdy_nu_CC_iso.fits";
    std::string nu_cc_total_xs = xs_base + "sigma_nu_CC_iso.fits";
    std::string nubar_cc_diff_xs = xs_base + "dsdxdy_nubar_CC_iso.fits";
    std::string nubar_cc_total_xs = xs_base + "sigma_nubar_CC_iso.fits";
    std::string nu_nc_diff_xs = xs_base + "dsdxdy_nu_NC_iso.fits";
    std::string nu_nc_total_xs = xs_base + "sigma_nu_NC_iso.fits";
    std::string nubar_nc_diff_xs = xs_base + "dsdxdy_nubar_NC_iso.fits";
    std::string nubar_nc_total_xs = xs_base + "sigma_nubar_NC_iso.fits";

    std::string nuebar_gr_diff_xs = xs_base + "dsdy_nuebar_GR.fits";
    std::string nuebar_gr_total_xs = xs_base + "sigma_nuebar_GR.fits";

    std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> default_cross_sections = {
        {{InteractionType::CC, LeptonInjector::Particle::ParticleType::NuE}, {nu_cc_diff_xs, nu_cc_total_xs}},
        {{InteractionType::CC, LeptonInjector::Particle::ParticleType::NuEBar}, {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {{InteractionType::CC, LeptonInjector::Particle::ParticleType::NuMu}, {nu_cc_diff_xs, nu_cc_total_xs}},
        {{InteractionType::CC, LeptonInjector::Particle::ParticleType::NuMuBar}, {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {{InteractionType::CC, LeptonInjector::Particle::ParticleType::NuTau}, {nu_cc_diff_xs, nu_cc_total_xs}},
        {{InteractionType::CC, LeptonInjector::Particle::ParticleType::NuTauBar}, {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {{InteractionType::NC, LeptonInjector::Particle::ParticleType::NuE}, {nu_nc_diff_xs, nu_nc_total_xs}},
        {{InteractionType::NC, LeptonInjector::Particle::ParticleType::NuEBar}, {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {{InteractionType::NC, LeptonInjector::Particle::ParticleType::NuMu}, {nu_nc_diff_xs, nu_nc_total_xs}},
        {{InteractionType::NC, LeptonInjector::Particle::ParticleType::NuMuBar}, {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {{InteractionType::NC, LeptonInjector::Particle::ParticleType::NuTau}, {nu_nc_diff_xs, nu_nc_total_xs}},
        {{InteractionType::NC, LeptonInjector::Particle::ParticleType::NuTauBar}, {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {{InteractionType::GR_E, LeptonInjector::Particle::ParticleType::NuEBar}, {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {{InteractionType::GR_MU, LeptonInjector::Particle::ParticleType::NuEBar}, {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {{InteractionType::GR_TAU, LeptonInjector::Particle::ParticleType::NuEBar}, {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {{InteractionType::GR_HAD, LeptonInjector::Particle::ParticleType::NuEBar}, {nuebar_gr_diff_xs, nuebar_gr_total_xs}}
    };
    return default_cross_sections;
}

void update_cross_sections(argagg::parser_results & args, std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> & cross_sections, std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, bool> & replaced_cross_sections) {
    std::function<void(InteractionType, LeptonInjector::Particle::ParticleType, std::string)> replace_total_xs = [&] (InteractionType type, LeptonInjector::Particle::ParticleType nu, std::string name) -> void {
        if(name != cross_sections[{type, nu}].second) {
            replaced_cross_sections[{type, nu}] = true;
        }
        cross_sections[{type, nu}].second = name;
    };

    std::function<void(InteractionType, LeptonInjector::Particle::ParticleType, std::string)> replace_diff_xs = [&] (InteractionType type, LeptonInjector::Particle::ParticleType nu, std::string name) -> void {
        if(name != cross_sections[{type, nu}].first) {
            replaced_cross_sections[{type, nu}] = true;
        }
        cross_sections[{type, nu}].first = name;
    };

    std::vector<InteractionType> gr_types = {InteractionType::GR_E, InteractionType::GR_MU, InteractionType::GR_TAU, InteractionType::GR_HAD};
    // Check the GR arguments
    if(args["nuebar_gr_diff_xs"]) {
        std::string nuebar_gr_diff_xs = args["nuebar_gr_diff_xs"].as<std::string>();
        for(unsigned int i=0; i<gr_types.size(); ++i) {
            replace_diff_xs(gr_types[i], LeptonInjector::Particle::ParticleType::NuEBar, nuebar_gr_diff_xs);
        }
    }
    if(args["nuebar_gr_total_xs"]) {
        std::string nuebar_gr_total_xs = args["nuebar_gr_total_xs"].as<std::string>();
        for(unsigned int i=0; i<gr_types.size(); ++i) {
            replace_total_xs(gr_types[i], LeptonInjector::Particle::ParticleType::NuEBar, nuebar_gr_total_xs);
        }
    }

    for(unsigned int i=0; i<gr_types.size(); ++i) {
        std::string total_id = get_total_xs_id(gr_types[i], LeptonInjector::Particle::ParticleType::NuEBar);
        std::string diff_id = get_diff_xs_id(gr_types[i], LeptonInjector::Particle::ParticleType::NuEBar);
        if(args[diff_id]) {
            std::string diff = args[diff_id].as<std::string>();
            replace_diff_xs(gr_types[i], LeptonInjector::Particle::ParticleType::NuEBar, diff);
        }
        if(args[total_id]) {
            std::string total = args[total_id].as<std::string>();
            replace_total_xs(gr_types[i], LeptonInjector::Particle::ParticleType::NuEBar, total);
        }
    }


    std::vector<InteractionType> cc_nc = {InteractionType::CC, InteractionType::NC};
    std::vector<LeptonInjector::Particle::ParticleType> nus = gen_nu_types();
    std::vector<LeptonInjector::Particle::ParticleType> matter_nus = {
        LeptonInjector::Particle::ParticleType::NuE,
        LeptonInjector::Particle::ParticleType::NuMu,
        LeptonInjector::Particle::ParticleType::NuTau
    };
    std::vector<LeptonInjector::Particle::ParticleType> antimatter_nus = {
        LeptonInjector::Particle::ParticleType::NuE,
        LeptonInjector::Particle::ParticleType::NuMu,
        LeptonInjector::Particle::ParticleType::NuTau
    };

    // Check generic flavor xs arguments
    for(InteractionType interaction : cc_nc) {
        std::string diff_key = std::string("nu_") + interaction_to_string(interaction) + "diff_xs";
        std::string total_key = std::string("nu_") + interaction_to_string(interaction) + "total_xs";
        for(LeptonInjector::Particle::ParticleType nu : matter_nus) {
            if(args[diff_key]) {
                std::string diff = args[diff_key].as<std::string>(cross_sections[{interaction, nu}].first);
                replace_diff_xs(interaction, nu, diff);
            }
            if(args[total_key]) {
                std::string total = args[total_key].as<std::string>(cross_sections[{interaction, nu}].second);
                replace_total_xs(interaction, nu, total);
            }
        }
        diff_key = std::string("nubar_") + interaction_to_string(interaction) + "diff_xs";
        total_key = std::string("nubar_") + interaction_to_string(interaction) + "total_xs";
        for(LeptonInjector::Particle::ParticleType nu : antimatter_nus) {
            if(args[diff_key]) {
                std::string diff = args[diff_key].as<std::string>(cross_sections[{interaction, nu}].first);
                replace_diff_xs(interaction, nu, diff);
            }
            if(args[total_key]) {
                std::string total = args[total_key].as<std::string>(cross_sections[{interaction, nu}].second);
                replace_total_xs(interaction, nu, total);
            }
        }
    }


    // Check specific flavor xs arguments
    for(LeptonInjector::Particle::ParticleType nu : nus) {
        for(InteractionType interaction : cc_nc) {
            std::string diff_key = get_diff_xs_id(interaction, nu);
            std::string total_key = get_total_xs_id(interaction, nu);
            if(args[diff_key]) {
                std::string diff = args[diff_key].as<std::string>(cross_sections[{interaction, nu}].first);
                replace_diff_xs(interaction, nu, diff);
            }
            if(args[total_key]) {
                std::string total = args[total_key].as<std::string>(cross_sections[{interaction, nu}].second);
                replace_total_xs(interaction, nu, total);
            }
        }
    }
}

void warn_user_cross_sections(argagg::parser_results & args, std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> & cross_sections, std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, bool> & replaced_cross_sections) {
    std::function<bool(LeptonInjector::Particle::ParticleType, std::vector<InteractionType>)> any_replaced = [&] (LeptonInjector::Particle::ParticleType nu, std::vector<InteractionType> types) -> bool {
        bool b = false;
        for(auto const & t : types) {
            b |= replaced_cross_sections.at({t, nu});
        }
        return b;
    };
    std::function<bool(LeptonInjector::Particle::ParticleType, std::vector<InteractionType>)> all_replaced = [&] (LeptonInjector::Particle::ParticleType nu, std::vector<InteractionType> types) -> bool {
        bool b = true;
        for(auto const & t : types) {
            b &= replaced_cross_sections.at({t, nu});
        }
        return b;
    };

    std::vector<std::pair<LeptonInjector::Particle::ParticleType, std::vector<InteractionType>>> all_tags = {
        {LeptonInjector::Particle::ParticleType::NuE, {InteractionType::CC, InteractionType::NC}},
        {LeptonInjector::Particle::ParticleType::NuMu, {InteractionType::CC, InteractionType::NC}},
        {LeptonInjector::Particle::ParticleType::NuTau, {InteractionType::CC, InteractionType::NC}},
        {LeptonInjector::Particle::ParticleType::NuEBar, {InteractionType::CC, InteractionType::NC, InteractionType::GR_E, InteractionType::GR_MU, InteractionType::GR_TAU, InteractionType::GR_HAD}},
        {LeptonInjector::Particle::ParticleType::NuMuBar, {InteractionType::CC, InteractionType::NC}},
        {LeptonInjector::Particle::ParticleType::NuTauBar, {InteractionType::CC, InteractionType::NC}},
    };

    std::vector<LeptonInjector::Particle::ParticleType> bad_nus;

    for(auto const & tags : all_tags) {
        if(any_replaced(tags.first, tags.second) and not all_replaced(tags.first, tags.second)) {
            bad_nus.push_back(tags.first);
        }
    }
    if(bad_nus.size() > 0) {
        std::stringstream ss;
        for(unsigned int i=0; i<bad_nus.size()-1; ++i) {
            ss << nu_to_string(bad_nus[i]) << ",";
        }
        ss << nu_to_string(bad_nus[bad_nus.size()-1]);
        for(unsigned int i=0; i<5; ++i) {
            std::cerr << "WARNING!" << std::endl;
        }
        std::string flavor_list = ss.str();
        if(bad_nus.size() > 1) {
            flavor_list = "{" + flavor_list + "}";
        }
        for(unsigned int i=0; i<10; ++i) {
            std::cerr << "WARNING! At least one but not all " + flavor_list + " cross sections have been changed from the default! If you are modifying multiple interaction cross sections, be sure to provide the new cross sections to all jobs, even if they are not injecting with the modified cross section. Physically accurate weighting at the injection stage depends on knowning all cross sections apriori." << std::endl;
        }
        for(unsigned int i=0; i<5; ++i) {
            std::cerr << "WARNING!" << std::endl;
        }
    }
}

std::vector<std::tuple<InjectionMode, LeptonInjector::Particle::ParticleType, std::vector<InteractionType>, int>> get_injectors(argagg::parser_results & args) {
    bool ranged = bool(args["ranged_mode"]);
    bool volume = bool(args["volume_mode"]);
    bool ranged_and_volume = bool(args["ranged_and_volume_mode"]);
    bool automatic = bool(args["automatic_injection_mode"]);

    int n_true = int(ranged) + int(volume) + int(ranged_and_volume) + int(automatic);
    if(n_true == 0) {
        std::cerr << "WARNING! Defaulting to automatic injection mode selection!" << std::endl;
        automatic = true;
    } else if(n_true > 1) {
        std::cerr << "Fatal Error: More than one inejction mode option specified. You must choose a single injection mode option." << std::endl;
        throw std::runtime_error("Fatal Error: More than one inejction mode option specified. You must choose a single injection mode option.");
    }

    std::vector<std::tuple<InjectionMode, LeptonInjector::Particle::ParticleType, std::vector<InteractionType>, int>> injectors;
    if(automatic) {
        if(args["nue"]) {
            if(args["cc"] and args["nc"]) {
                std::cout << "Injecting NuE CC+NC simultaneously in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuE, {InteractionType::CC, InteractionType::NC}, int(args["n_events"].as<float>())});
            }
            else if(args["cc"]) {
                std::cout << "Injecting NuE CC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuE, {InteractionType::CC}, int(args["n_events"].as<float>())});
            }
            else if(args["nc"]) {
                std::cout << "Injecting NuE NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuE, {InteractionType::NC}, int(args["n_events"].as<float>())});
            }
        }
        if(args["numu"]) {
            if(args["cc"] and args["nc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuMu events with CC in Ranged mode, and half of NuMu events with CC+NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuMu, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuMu, {InteractionType::CC, InteractionType::NC}, volume_events});
            }
            else if(args["cc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuMu events with CC in Ranged mode, and half of NuMu events with CC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuMu, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuMu, {InteractionType::CC}, volume_events});
            }
            else if(args["nc"]) {
                std::cout << "Injecting NuMu NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuMu, {InteractionType::NC}, int(args["n_events"].as<float>())});
            }
        }
        if(args["nutau"]) {
            if(args["cc"] and args["nc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuTau events with CC in Ranged mode, and half of NuTau events with CC+NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuTau, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuTau, {InteractionType::CC, InteractionType::NC}, volume_events});
            }
            else if(args["cc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuTau events with CC in Ranged mode, and half of NuTau events with CC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuTau, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuTau, {InteractionType::CC}, volume_events});
            }
            else if(args["nc"]) {
                std::cout << "Injecting NuTau NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuTau, {InteractionType::NC}, int(args["n_events"].as<float>())});
            }
        }

        if(args["nuebar"]) {

            // gr or gr_mu or gr_tau --> ranged

            std::vector<InteractionType> nue_interaction_types = {InteractionType::CC, InteractionType::NC, InteractionType::GR_E, InteractionType::GR_MU, InteractionType::GR_TAU, InteractionType::GR_HAD};
            std::map<InteractionType, std::string> interaction_type_strings;
            std::vector<InteractionType> gr_interactions = {InteractionType::GR_E, InteractionType::GR_MU, InteractionType::GR_TAU, InteractionType::GR_HAD};
            std::vector<InteractionType> ranged_interactions = {InteractionType::GR_MU, InteractionType::GR_TAU};

            std::map<InteractionType, bool> enabled_interactions;
            std::vector<InteractionType> enabled_interaction_types;
            std::vector<InteractionType> enabled_ranged_types;
            for(unsigned int i=0; i<nue_interaction_types.size(); ++i) {
                InteractionType const & t = nue_interaction_types[i];
                std::string s = interaction_to_string(t);
                interaction_type_strings.emplace(t, s);
                enabled_interactions[t] = bool(args[s]);
                if(args[s]) {
                    enabled_interaction_types.push_back(t);
                    if(std::find(ranged_interactions.begin(), ranged_interactions.end(), t) != ranged_interactions.end()) {
                        enabled_ranged_types.push_back(t);
                    }
                }
            }

            bool gr = bool(args["gr"]);

            if(gr) {
                for(auto const & t : gr_interactions) {
                    enabled_interactions[t] = true;
                }
            }

            bool need_ranged = false;
            for(auto const & t : ranged_interactions) {
                need_ranged |= enabled_interactions[t];
            }

            std::stringstream ss;
            for(unsigned int i=0; i<enabled_interaction_types.size()-1; ++i) {
                std::string str = interaction_type_strings[enabled_interaction_types[i]];
                std::transform(str.begin(), str.end(),str.begin(), [](auto c) { return std::toupper(c);});
                ss << str << "+";
            }
            std::string str = interaction_type_strings[enabled_interaction_types[enabled_interaction_types.size()-1]];
            std::transform(str.begin(), str.end(),str.begin(), [](auto c) { return std::toupper(c);});
            ss << str;
            if(enabled_interaction_types.size() > 1) {
                ss << " simultaneously";
            }

            if(need_ranged) {
                std::stringstream ss_ranged;
                for(unsigned int i=0; i<enabled_ranged_types.size()-1; ++i) {
                    std::string str = interaction_type_strings[enabled_ranged_types[i]];
                    std::transform(str.begin(), str.end(),str.begin(), [](auto c) { return std::toupper(c);});
                    ss_ranged << str << "+";
                }
                str = interaction_type_strings[enabled_interaction_types[enabled_ranged_types.size()-1]];
                std::transform(str.begin(), str.end(),str.begin(), [](auto c) { return std::toupper(c);});
                ss_ranged << str;
                if(enabled_ranged_types.size() > 1) {
                    ss_ranged << " simultaneously";
                }
                std::cout << "Injecting NuEBar " + ss.str() + " in Volume mode, and NuEBar " + ss_ranged.str() + " in Ranged mode." << std::endl;

                int n_events = int(args["n_events"].as<float>());
                int n_ranged = int(0.107 * enabled_ranged_types.size() * n_events);
                int n_volume = n_events - n_ranged;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuEBar, enabled_ranged_types, n_ranged});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuEBar, enabled_interaction_types, n_volume});
            } else {
                std::cout << "Injecting NuEBar in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuEBar, enabled_interaction_types, int(args["n_events"].as<float>())});
            }
        }
        if(args["numubar"]) {
            if(args["cc"] and args["nc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuMuBar events with CC in Ranged mode, and half of NuMuBar events with CC+NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuMuBar, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuMuBar, {InteractionType::CC, InteractionType::NC}, volume_events});
            }
            else if(args["cc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuMuBar events with CC in Ranged mode, and half of NuMuBar events with CC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuMuBar, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuMuBar, {InteractionType::CC}, volume_events});
            }
            else if(args["nc"]) {
                std::cout << "Injecting NuMuBar NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuMuBar, {InteractionType::NC}, int(args["n_events"].as<float>())});
            }
        }
        if(args["nutaubar"]) {
            if(args["cc"] and args["nc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuTauBar events with CC in Ranged mode, and half of NuTauBar events with CC+NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuTauBar, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuTauBar, {InteractionType::CC, InteractionType::NC}, volume_events});
            }
            else if(args["cc"]) {
                int total_events = int(args["n_events"].as<float>());
                int ranged_events = total_events / 2;
                int volume_events = total_events - ranged_events;
                std::cout << "Injecting half of NuTauBar events with CC in Ranged mode, and half of NuTauBar events with CC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Ranged, LeptonInjector::Particle::ParticleType::NuTauBar, {InteractionType::CC}, ranged_events});
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuTauBar, {InteractionType::CC}, volume_events});
            }
            else if(args["nc"]) {
                std::cout << "Injecting NuTauBar NC in Volume mode." << std::endl;
                injectors.push_back({InjectionMode::Volume, LeptonInjector::Particle::ParticleType::NuTauBar, {InteractionType::NC}, int(args["n_events"].as<float>())});
            }
        }
    }
    else {
        int total_events = int(args["n_events"].as<float>());
        std::vector<int> events_to_inject;
        std::vector<InjectionMode> modes;
        std::map<InjectionMode, std::string> mode_to_str = {{InjectionMode::Ranged, "Ranged"}, {InjectionMode::Volume, "Volume"}};
        if(ranged) {
            modes.push_back(InjectionMode::Ranged);
            events_to_inject.push_back(total_events);
        } else if(volume) {
            modes.push_back(InjectionMode::Volume);
            events_to_inject.push_back(total_events);
        } else if(ranged_and_volume) {
            int n = total_events / 2;
            events_to_inject.push_back(n);
            events_to_inject.push_back(total_events - n);
            modes.push_back(InjectionMode::Ranged);
            modes.push_back(InjectionMode::Volume);
        }
        std::vector<LeptonInjector::Particle::ParticleType> nu_types = {
            LeptonInjector::Particle::ParticleType::NuE,
            LeptonInjector::Particle::ParticleType::NuMu,
            LeptonInjector::Particle::ParticleType::NuTau,
            LeptonInjector::Particle::ParticleType::NuEBar,
            LeptonInjector::Particle::ParticleType::NuMuBar,
            LeptonInjector::Particle::ParticleType::NuTauBar
        };
        std::vector<InteractionType> interaction_types = {InteractionType::CC, InteractionType::NC, InteractionType::GR_E, InteractionType::GR_MU, InteractionType::GR_TAU, InteractionType::GR_HAD};
        for(unsigned int m=0; m<modes.size(); ++m) {
            for(unsigned int n=0; n<nu_types.size(); ++n) {
                std::vector<InteractionType> interactions_to_inject = get_allowed_interactions(args, nu_types[n], interaction_types);
                if(interactions_to_inject.size() > 0 and events_to_inject[m] > 0) {
                    injectors.push_back({modes[m], nu_types[n], interactions_to_inject, events_to_inject[m]});
                }
            }
        }
    }
    return injectors;
}

std::map<std::pair<std::string, std::string>, std::shared_ptr<LeptonInjector::CrossSection>> get_xs_ptrs_by_fname(std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> & cross_sections) {
    std::map<std::pair<std::string, std::string>, std::vector<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>>> xs_fname_to_type;
    for(auto const & xs : cross_sections) {
		xs_fname_to_type[xs.second].push_back(xs.first);
    }
    std::map<std::pair<std::string, std::string>, std::shared_ptr<LeptonInjector::CrossSection>> xs_ptrs_by_fname;

    std::vector<LeptonInjector::Particle::ParticleType> dis_target_types = {LeptonInjector::Particle::ParticleType::Nucleon};
    for(auto const & xs : xs_fname_to_type) {
        std::shared_ptr<LeptonInjector::CrossSection> cross_section;
        std::vector<LeptonInjector::Particle::ParticleType> dis_primary_types;
        for(auto const & itype_ptype : xs.second) {
            dis_primary_types.push_back(itype_ptype.second);
        }
        std::shared_ptr<LeptonInjector::DISFromSpline> dis_cross_section = std::make_shared<LeptonInjector::DISFromSpline>(xs.first.first, xs.first.second, dis_primary_types, dis_target_types);
        xs_ptrs_by_fname[xs.first] = std::shared_ptr<LeptonInjector::CrossSection>(dis_cross_section);
    }
    return xs_ptrs_by_fname;
}
std::shared_ptr<earthmodel::EarthModel> get_earth_model(argagg::parser_results & args) {
    std::string path = get_path(args);
    std::string earth_file = args["earth_model"].as<std::string>("PREM_mmc");
    std::string materials_file = args["materials_model"].as<std::string>("Standard");
    std::shared_ptr<earthmodel::EarthModel> earth_model = std::make_shared<earthmodel::EarthModel>();
    earth_model->SetPath(path);
    earth_model->LoadMaterialModel(materials_file);
    earth_model->LoadEarthModel(earth_file);
    return earth_model;
}

void inject_events(std::shared_ptr<LeptonInjector::InjectorBase> injector) {// , std::shared_ptr<LeptonInjector::LeptonWeighter> weighter) {
    while(*injector) {
        LeptonInjector::InteractionRecord event = injector->GenerateEvent();
        //double weight = weighter->SimplifiedEventWeight(event);
    }
}

int main(int argc, char ** argv) {
try {
    std::tuple<argagg::parser, argagg::parser_results> argarg = parse_arguments(argc, argv);
    argagg::parser & argparser = std::get<0>(argarg);
    argagg::parser_results & args = std::get<1>(argarg);

    double seed = get_seed(argparser, args);
    std::string output = get_output(argparser, args);
    std::string path = get_path(args);

    int n_events = int(args["n_events"].as<float>(float(1e5)));

    double minE = args["minE"].as<double>(1e2)*LeptonInjector::Constants::GeV;
    double maxE = args["maxE"].as<double>(1e6)*LeptonInjector::Constants::GeV;
    double gamma = args["gamma"].as<double>(2.);
    double minZenith = args["minZenith"].as<double>(0.)*LeptonInjector::Constants::degrees;
    double maxZenith = args["maxZenith"].as<double>(180.)*LeptonInjector::Constants::degrees;
    double minAzimuth = args["minAzimuth"].as<double>(0.)*LeptonInjector::Constants::degrees;
    double maxAzimuth = args["maxAzimuth"].as<double>(360.)*LeptonInjector::Constants::degrees;
    double ranged_radius = args["ranged_radius"].as<double>(32)*LeptonInjector::Constants::m;
    double ranged_length = args["ranged_length"].as<double>(32)*LeptonInjector::Constants::m;
    double volume_radius = args["volume_radius"].as<double>(32)*LeptonInjector::Constants::m;
    double volume_height = args["volume_height"].as<double>(20)*LeptonInjector::Constants::m;

    std::vector<std::string> interactions = get_interactions(args);
    std::vector<std::string> neutrinos = get_neutrinos(args);
    std::map<std::string, std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>> secondaries = gen_secondaries();
    std::string xs_base = get_xs_base(args);
    std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, std::pair<std::string, std::string>> cross_sections = gen_default_cross_sections(xs_base);
    std::map<std::pair<InteractionType, LeptonInjector::Particle::ParticleType>, bool> replaced_cross_sections;
    for(auto const & default_xs : cross_sections) {
        replaced_cross_sections.emplace(default_xs.first, false);
    }

    update_cross_sections(args, cross_sections, replaced_cross_sections);
    warn_user_cross_sections(args, cross_sections, replaced_cross_sections);

    std::cout << getISOCurrentTimestamp<std::chrono::seconds>() << std::endl;

    std::vector<std::tuple<InjectionMode, LeptonInjector::Particle::ParticleType, std::vector<InteractionType>, int>> injectors = get_injectors(args);

    if(injectors.size() == 0) {
        std::cerr << "Must have at least one neutrino and one interaction type specified!" << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }

    std::map<std::pair<std::string, std::string>, std::shared_ptr<LeptonInjector::CrossSection>> xs_ptrs_by_fname = get_xs_ptrs_by_fname(cross_sections);
    std::vector<std::shared_ptr<LeptonInjector::CrossSection>> physical_xs_ptrs;
    for(auto xs : xs_ptrs_by_fname) {
        physical_xs_ptrs.push_back(xs.second);
    }

    std::shared_ptr<earthmodel::EarthModel> earth_model = get_earth_model(args);

    // Pick energy distribution
    std::shared_ptr<LeptonInjector::PrimaryEnergyDistribution> edist = std::make_shared<LeptonInjector::PowerLaw>(gamma, minE, maxE);

    // Choose injection direction
    std::shared_ptr<LeptonInjector::PrimaryDirectionDistribution> ddist = std::make_shared<LeptonInjector::Cone>(earthmodel::Vector3D{0.0, 0.0, 1.0}, minZenith, maxZenith, minAzimuth, maxAzimuth);

    // Targets should be stationary
    std::shared_ptr<LeptonInjector::TargetMomentumDistribution> target_momentum_distribution = std::make_shared<LeptonInjector::TargetAtRest>();

    // Let us inject according to the decay distribution
    std::shared_ptr<LeptonInjector::DepthFunction> depth_func = std::make_shared<LeptonInjector::LeptonDepthFunction>();

    // Helicity distribution
    std::shared_ptr<LeptonInjector::PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<LeptonInjector::PrimaryNeutrinoHelicityDistribution>();

    // Random number service
    std::shared_ptr<LeptonInjector::LI_random> random = std::make_shared<LeptonInjector::LI_random>();

    std::vector<std::shared_ptr<LeptonInjector::InjectorBase>> injector_pointers;

    for(unsigned inj_i=0; inj_i<injectors.size(); ++inj_i) {
        std::tuple<InjectionMode, LeptonInjector::Particle::ParticleType, std::vector<InteractionType>, int> injector_config = injectors[inj_i];
        InjectionMode injection_mode = std::get<0>(injector_config);
        LeptonInjector::Particle::ParticleType primary_type = std::get<1>(injector_config);
        std::vector<InteractionType> interaction_types = std::get<2>(injector_config);
        int n_events = std::get<3>(injector_config);
        std::vector<std::shared_ptr<LeptonInjector::CrossSection>> injector_cross_sections;
        std::map<std::pair<std::string, std::string>, std::shared_ptr<LeptonInjector::CrossSection>> injector_cross_sections_by_fname;
        for(auto const & interaction : interaction_types) {
            std::pair<InteractionType, LeptonInjector::Particle::ParticleType> xs_key = {interaction, primary_type};
            std::pair<std::string, std::string> xs_fname = cross_sections[xs_key];
            injector_cross_sections_by_fname[xs_fname] = xs_ptrs_by_fname[xs_fname];
        }
        for(auto & xs : xs_ptrs_by_fname) {
            injector_cross_sections.push_back(xs.second);
        }

        std::shared_ptr<LeptonInjector::PrimaryInjector> primary_injector = std::make_shared<LeptonInjector::PrimaryInjector>(primary_type);

        std::shared_ptr<LeptonInjector::InjectorBase> injector;
        if(injection_mode == InjectionMode::Ranged) {
            injector_pointers.push_back(std::make_shared<LeptonInjector::ColumnDepthLeptonInjector>(n_events, primary_injector, injector_cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, depth_func, ranged_radius, ranged_length, helicity_distribution));
        } else {
            // Placement of cylinder in detector coordinates
            earthmodel::Placement placement(earthmodel::Vector3D(0, 0, 0));
            earthmodel::Cylinder injection_cylinder(placement, volume_radius, 0.0, volume_height);
            injector_pointers.push_back(std::make_shared<LeptonInjector::VolumeLeptonInjector>(n_events, primary_injector, injector_cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, injection_cylinder, helicity_distribution));
        }
    }

    // Change the flux units from cm^-2 to m^-2
    std::shared_ptr<LeptonInjector::WeightableDistribution> flux_units = std::make_shared<LeptonInjector::NormalizationConstant>(1e4);

    std::shared_ptr<LeptonInjector::IsotropicDirection> physical_ddist = std::make_shared<LeptonInjector::IsotropicDirection>();

    std::vector<std::shared_ptr<LeptonInjector::WeightableDistribution>> physical_distributions = {
        std::shared_ptr<LeptonInjector::WeightableDistribution>(edist),
        std::shared_ptr<LeptonInjector::WeightableDistribution>(flux_units),
        std::shared_ptr<LeptonInjector::WeightableDistribution>(physical_ddist),
        std::shared_ptr<LeptonInjector::WeightableDistribution>(target_momentum_distribution),
        std::shared_ptr<LeptonInjector::WeightableDistribution>(helicity_distribution)
    };

    for(auto & injector : injector_pointers) {
        inject_events(injector);
    }

} catch(ExitStatus const & status) {
    return status.status;
}
}
