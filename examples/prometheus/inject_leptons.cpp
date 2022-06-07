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

template <class Precision>
std::string getISOCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    return date::format("%FT%TZ", date::floor<Precision>(now));
}

argagg::parser_results parse_arguments(int argc, char ** argv) {
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
            "dune_atmo_path", {"--dune-atmo-path"},
            "Path to the DUNEAtmo source directory. Used for the injection cross section.", 1,
        },
        {
            "lepton_injector_path", {"--li", "--li-path", "--lepton-injector-path"},
            "Path to the lepton injector source directory. Used for earth model resources.", 1,
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
            "ice_type", {"--ice", "--ice-type"},
            "The type of ice model.", 1,
        },
        {
            "ice_angle", {"--ice-angle", "--ice-cap-angle"},
            "Angle of the ice cap in degrees if it exists.", 1,
        },
        {
            "depth", {"--depth", "--detector-depth"},
            "Depth of the detector origin in meters.", 1,
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
        return EXIT_FAILURE;
    }

    // Print the usage information when we help is requested
    if(args["help"]) {
        std::cerr << argparser;
        return EXIT_SUCCESS;
    }

    return args;
}

double get_seed(argagg::parser_results & args) {
    // Grab the seed
    double seed = args["seed"].as<int>(-1);
    // Require seed to be set and positive
    if(seed < 0) {
        std::cerr << "--seed requires positive integer!" << std::endl << argparser;
        return EXIT_FAILURE;
    }
    return seed;
}

std::string get_output(argagg::parser_results & args) {
    // Require output prefix
    if(not args["output"]) {
        std::cerr << "--output required!" << std::endl << argparser;
        return EXIT_FAILURE;
    }
    std::string output = args["output"].as<std::string>("./injected/output_DUNE");
    return output;
}

std::string get_path(argagg::parser_results & args) {
    // Require LI path
    std::string path;
    if(args["lepton_injector_path"]) {
        path = args["lepton_injector_path"].as<std::string>();
    }
    else { // Try to grab it from the environment variables instead
        if (const char* env_p = getenv("GOLEMSOURCEPATH")){
            path = std::string( env_p ) + "/LeptonInjectorDUNE/";
        }
        else {
            std::cerr << "WARNING no lepton injector path specified and GOLEMSOURCEPATH not set! Assuming earth model information is in ./resources/earthparams/" << std::endl;
            path = "./";
        }
    }
    if(path[path.size() - 1] != '/')
        path += "/";
    path += "resources/earthparams";

    return path;
}

std::vector<std::string> gen_possible_interactions() {
    return {"cc", "nc", "gr_e", "gr_mu", "gr_tau", "gr_had"}
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

std::string get_xs_base(argagg::parser_result & args) {
    std::string xs_base;
    if(args["dune_atmo_path"]) {
        xs_base = args["dune_atmo_path"].as<std::string>();
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
    xs_base += "cross_sections/csms_differential_v1.0/";

    return xs_base;
}

std::map<std::string, std::pair<std::string, std::string>> gen_default_cross_sections(std::string xs_base) {

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

    std::map<std::string, std::pair<std::string, std::string>> default_cross_sections = {
        {"nue_cc", {nu_cc_diff_xs, nu_cc_total_xs}},
        {"nuebar_cc", {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {"numu_cc", {nu_cc_diff_xs, nu_cc_total_xs}},
        {"numubar_cc", {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {"nutau_cc", {nu_cc_diff_xs, nu_cc_total_xs}},
        {"nutaubar_cc", {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {"nue_nc", {nu_nc_diff_xs, nu_nc_total_xs}},
        {"nuebar_nc", {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {"numu_nc", {nu_nc_diff_xs, nu_nc_total_xs}},
        {"numubar_nc", {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {"nutau_nc", {nu_nc_diff_xs, nu_nc_total_xs}},
        {"nutaubar_nc", {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {"nuebar_gr_e", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {"nuebar_gr_mu", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {"nuebar_gr_tau", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {"nuebar_gr_had", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
    };
    return default_cross_sections;
}

void update_cross_sections(argagg::parser_result & args, std::map<std::string, std::pair<std::string, std::string>> & cross_sections, std::map<std::string, bool> & replaced_cross_sections) {
    std::function<void(std::string, std::string)> replace_total_xs = [&] (std::string key, std::string name) -> void {
        if(name != cross_sections[key].second) {
            replaced_cross_sections[key] = true;
        }
        cross_sections[key].second = name;
    };

    std::function<void(std::string, std::string)> replace_diff_xs = [&] (std::string key, std::string name) -> void {
        if(name != cross_sections[key].first) {
            replaced_cross_sections[key] = true;
        }
        cross_sections[key].first = name;
    };

    // Check the GR arguments
    if(args["nuebar_gr_diff_xs"]) {
        std::string nuebar_gr_diff_xs = args["nuebar_gr_diff_xs"].as<std::string>();
        replace_diff_xs("nuebar_gr", nuebar_gr_diff_xs);
    }
    if(args["nuebar_gr_total_xs"]) {
        std::string nuebar_gr_total_xs = args["nuebar_gr_total_xs"].as<std::string>();
        replace_total_xs("nuebar_gr", nuebar_gr_total_xs);
    }

    std::vector<std::string> possible_interaction_suffixes = {"e", "mu", "tau", "had"};
    for(unsigned int i=0; i<possible_interaction_suffixes.size(); ++i) {
        std::string diff = nuebar_gr_diff_xs;
        std::string total = nuebar_gr_total_xs;
        std::string id = "nuebar_gr_" + possible_interaction_suffixes[i];
        if(args[id + "_diff_xs"]) {
            diff = args[id + "_diff_xs"].as<std::string>();
            replace_diff_xs(id, diff);
        }
        if(args[id + "_total_xs"]) {
            total = args[id + "_total_xs"].as<std::string>();
            replace_total_xs(id, total);
        }
    }

    // Check generic flavor xs arguments
    std::vector<std::string> cc_nc = {"cc", "nc"};
    std::vector<std::string> flavors = {"e", "mu", "tau"};
    std::vector<std::string> prefixes = {"nu", "nubar"}
    for(std::string prefix : prefixes) {
        for(std::string & interaction : cc_nc) {
            std::string key = prefix + "_" + interaction;
            if(args[key + "_diff_xs"]) {
                std::string diff = args[key + "_diff_xs"].as<std::string>(cross_sections[key].first);
                for(std::string & flavor : flavors) {
                    std::string id = prefix + flavor + "_" + interaction;
                    replace_diff_xs(id, diff);
                }
            }
            if(args[key + "_total_xs"]) {
                std::string total = args[key + "_total_xs"].as<std::string>(cross_sections[key].second);
                for(std::string & flavor : flavors) {
                    std::string id = prefix + flavor + "_" + interaction;
                    replace_total_xs(id, total);
                }
            }
        }
    }

    // Check specific flavor xs arguments
    std::vector<std::string> interactions = gen_possible_interactions();
    std::vector<std::string> neutrinos = gen_possible_neutrinos();
    for(std::string & interaction : interactions) {
        bool is_glashow = interaction.rfind("gr", 0) == 0;
        for(std::string & neutrino : neutrinos) {
            if(not args[neutrino])
                continue;
            if(is_glashow and neutrino != "nuebar")
                continue;
            std::string id = neutrino + "_" + interaction;

            std::string diff_xs_id = id + "_diff_xs";
            std::string diff_xs =
                args[diff_xs_id].as<std::string>(cross_sections[id].first);
            replace_diff_xs(id, diff_xs);

            std::string total_xs_id = id + "_total_xs";
            std::string total_xs =
                args[total_xs_id].as<std::string>(cross_sections[id].second);
            replace_total_xs(id, total_xs);
        }
    }
}

void warn_user_cross_sections(argagg::parser_result & args, std::map<std::string, std::pair<std::string, std::string>> & cross_sections, std::map<std::string, bool> & replaced_cross_sections) {
    std::function<bool(std::vector<std::string>)> any_replaced = [] (std::vector<std::string> keys) -> bool {
        bool b = false;
        for(auto const & k : keys) {
            b |= replaced_cross_sections.at(k);
        }
        return b;
    };
    std::function<bool(std::vector<std::string>)> all_replaced = [] (std::vector<std::string> keys) -> bool {
        bool b = true;
        for(auto const & k : keys) {
            b &= replaced_cross_sections.at(k);
        }
        return b;
    };

    std::vector<std::pair<std::string, std::vector<std::string>>> all_tags = {
        {"nue", {"nue_cc", "nue_nc"}},
        {"numu", {"numu_cc", "numu_nc"}},
        {"nutau", {"nutau_cc", "nutau_nc"}},
        {"nubar", {"nuebar_cc", "nuebar_nc", "nuebar_gr_e", "nuebar_gr_mu", "nuebar_gr_tau", "nuebar_gr_had"}},
        {"numubar", {"numubar_cc", "numubar_nc"}},
        {"nutaubar", {"nutaubar_cc", "nutaubar_nc"}},
    };

    std::vector<std::string> bad_tags;

    for(auto const & tags : all_tags) {
        if(any_replaced(tags.second) and not all_replaced(tags.second)) {
            bad_tags.push_back(tags.first);
        }
    }
    if(bad_tags.size() > 0) {
        std::stringstream ss;
        for(unsigned int i=0; i<bad_tags.size()-1; ++i) {
            ss << bad_tags[i] << ",";
        }
        ss << bad_tags[bad_tags.size()-1];
        for(unsigned int i=0; i<5; ++i) {
            std::cerr << "WARNING!" << std::endl;
        }
        std::string flavor_list = ss.str();
        if(bad_tags.size() > 1) {
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

std::vector<std::tuple<bool, LeptonInjector::Particle::ParticleType, std::string, int>> get_injectors(argagg::parser_result & args) {
    std::vector<std::tuple<bool, LeptonInjector::Particle::ParticleType, std::string, int>> injectors;
    if(args["nue"]) {
        if(args["cc"] and args["nc"]) {
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "cc_nc", int(args["n_events"].as<float>())});
        }
        else if(args["cc"]) {
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "cc", int(args["n_events"].as<float>())});
        }
        else if(args["nc"]) {
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "nc", int(args["n_events"].as<float>())});
        }
    }
    if(args["numu"]) {
        if(args["cc"] and args["nc"]) {
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "cc_nc", int(args["n_events"].as<float>())});
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "cc_nc", int(args["n_events"].as<float>())});
        }
        else if(args["cc"]) {
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "cc", int(args["n_events"].as<float>())});
        }
        else if(args["nc"]) {
            injectors.push_back({false, LeptonInjector::Particle::ParticleType::NuE, "nc", int(args["n_events"].as<float>())});
        }
    }
}

int main(int argc, char ** argv) {

    argagg::parser_results args = parse_arguments(argc, argv);

    double seed = get_seed(args);
    std::string output = get_output(args);
    std::string path = get_path(args);

    std::string earth_model = args["earth_model"].as<std::string>("PREM_mmc");
    std::string materials_model = args["materials_model"].as<std::string>("Standard");
    std::string ice_type = args["ice_type"].as<std::string>("SimpleIceCap");
    double ice_angle = args["ice_angle"].as<double>(20.);
    double depth = args["depth"].as<double>(1480.);
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
    std::map<std::string, std::pair<std::string, std::string>> cross_sections = gen_default_cross_sectionst(xs_base);
    std::map<std::string, bool> replaced_cross_sections;
    for(auto const & default_xs : default_cross_sections) {
        replaced_cross_sections.emplace(default_xs.first, false);
    }

    update_cross_sections(args, cross_sections, replaced_cross_sections);
    warn_user_cross_sections(args, cross_sections, replaced_cross_sections);


    std::cout << getISOCurrentTimestamp<std::chrono::seconds>() << std::endl;


    std::vector<LeptonInjector::InjectorBase> injectors;
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    for(std::string & interaction : interactions) {
        bool is_glashow = interaction.rfind("gr", 0) == 0;
        if(args[interaction]) {
            for(std::string & neutrino : neutrinos) {
                if(args[neutrino]) {
                    if(is_glashow and neutrino != "nuebar") {
                        continue;
                    }
                    std::string id = neutrino + "_" + interaction;

                    std::string diff_xs_id = id + "_diff_xs";
                    std::string diff_xs =
                        args[diff_xs_id].as<std::string>(default_cross_sections[id].first);

                    std::string total_xs_id = id + "_total_xs";
                    std::string total_xs =
                        args[total_xs_id].as<std::string>(default_cross_sections[id].second);

                    LeptonInjector::Particle::ParticleType first_particle = secondaries[id].first;
                    LeptonInjector::Particle::ParticleType second_particle = secondaries[id].second;

                    if(n_ranged_events > 0) {
                        std::cout << "loading ranged injector\n";
                        LeptonInjector::RangedInjector injector(n_ranged_events,
                                first_particle,
                                //second_particle,
                                diff_xs,
                                total_xs,
                                true); // is_ranged
                        injectors.push_back(injector);
                    }
                    if(n_volume_events > 0) {
                        LeptonInjector::Injector injector(n_volume_events,
                                first_particle,
                                second_particle,
                                diff_xs,
                                total_xs,
                                false); // is_ranged
                        injectors.push_back(injector);
                    }

                }
                else {
                    continue;
                }
            }
        }
        else {
            continue;
        }
    }

    if(injectors.size() == 0) {
        std::cerr << "Must have at least one neutrino and one interaction type specified!" << std::endl;
        return EXIT_FAILURE;
    }

    // build the Controller object. This will facilitate the simulation itself
    // We need to pass the first injector while building this Controller
    // Dimensions of DUNE module is 58.2 x 3.5 x 12 x 4
    LeptonInjector::Controller cont(injectors[0],
            minE, maxE,
            gamma,
            minAzimuth, maxAzimuth,
            minZenith, maxZenith,
            ranged_radius, // injection radius
            ranged_length, // injection length
            volume_radius, // cylinder radius
            volume_height); // cylinder height


    for(unsigned int i=1; i<injectors.size(); ++i) {
        cont.AddInjector(injectors[i]);
    }


    cont.setSeed(seed);


    std::cout << "Mat model: " << materials_model << std::endl;
    try {
        earthmodel::EarthModel earthModel;
        earthModel.SetPath(path);
        earthModel.LoadMaterialModel(materials_model);
        //earthModel.LoadConcentricShellsFromLegacyFile(earth_model, depth*LeptonInjector::Constants::m, ice_angle*LeptonInjector::Constants::degrees);
        try{
					earthModel.LoadEarthModel(earth_model);
        }
        catch(const std::exception& err){
					std::cout << err.what();
        }

        std::vector<earthmodel::EarthSector> sectors = earthModel.GetSectors();
        for(unsigned int i=0; i<sectors.size(); ++i) {
            std::cerr << "Sector " << sectors[i].name << std::endl;
            std::cerr << "id: " << sectors[i].material_id << std::endl;
            std::cerr << "level: " << sectors[i].level << std::endl;
            std::shared_ptr<const earthmodel::Geometry> geo = sectors[i].geo;
            std::cerr << "Geo placement: " << geo->GetPlacement() << std::endl;
            std::cerr << "density: " << sectors[i].density->Evaluate(geo->GetPlacement().GetPosition()) << std::endl;
        }

				/*
        earthmodel::EarthModelService old_earthModel(
                "DUNE",
                path,
                std::vector<std::string>({earth_model}),
                std::vector<std::string>({materials_model}),
                ice_type,
                ice_angle*LeptonInjector::Constants::degrees,
                depth*LeptonInjector::Constants::m);

        */

        cont.SetEarthModel(std::shared_ptr<earthmodel::EarthModel>(&earthModel));
        //cont.SetEarthModel(std::shared_ptr<earthmodel::EarthModelService>(&old_earthModel));

        cont.NameOutfile(output + ".h5");
        cont.NameLicFile(output + ".lic");

        // Run the program.
        cont.Execute();
    } catch(char const * s) {
        std::cout << "Failure!" << std::endl;
        std::cout << s << std::endl;
    }
}
