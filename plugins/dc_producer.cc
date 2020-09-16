// Copyright [2020] Tyler Mitchell

#include <dirent.h>
#include <sys/types.h>

#include <ctime>
#include <fstream>
#include <iostream>
#include <locale>
#include <memory>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

#include "../include/CLParser.h"
#include "../include/json.hpp"
#include "TFile.h"
#include "TH2F.h"
#include "TStopwatch.h"
#include "TTree.h"

using std::string;
using std::vector;

void read_directory(const string &, vector<string> *, string match = "");
std::unordered_map<string, vector<string>> build_file_paths(string);
string format_output_name(string, bool, bool, string, int, int, string);
TH2F *build_histogram(string, vector<double>, vector<double>);
void process_file(string name, string channel, vector<TH2F *> histograms, string xvar_name, string yvar_name, string zvar_name, vector<double> edges,
                  string fake_weight_name, int DCP_idx);

int main(int argc, char *argv[]) {
    auto watch = TStopwatch();
    watch.Start();
    CLParser parser(argc, argv);
    bool is_ztt = parser.Flag("-z");
    bool do_syst = parser.Flag("-s");
    string dir = parser.Option("-d");
    string channel = parser.Option("-l");
    string config_name = parser.Option("-c");
    string year = parser.Option("-y");
    string suffix = parser.Option("-x");

    // get input file directory
    if (dir.empty()) {
        std::cerr << "You must give an input directory" << std::endl;
        return -1;
    }

    // get TTree name
    auto tree_name = channel + "_tree";

    // build file paths to process
    auto file_paths = build_file_paths(dir);

    std::cout << "Printing files for directory: " << dir << std::endl;
    for (auto &fp : file_paths) {
        std::cout << fp.first << std::endl;
        for (auto &f : fp.second) {
            std::cout << "\t" << f << std::endl;
        }
    }

    time_t now = time(0);
    tm *ltm = localtime(&now);
    auto month = 1 + ltm->tm_mon;
    auto day = ltm->tm_mday;

    auto output_file_name = format_output_name(channel, is_ztt, do_syst, year, month, day, suffix);
    std::cout << "creating output file " << output_file_name << std::endl;
    auto fout = std::make_shared<TFile>(output_file_name.c_str(), "recreate");

    // process json configs
    std::ifstream bp_file("configs/boilerplate.json");
    nlohmann::json bp_json;
    bp_file >> bp_json;

    std::ifstream binning_file("configs/binning.json");
    nlohmann::json binning_json;
    binning_file >> binning_json;

    // binnings
    auto tau_pt_bins = binning_json.at(config_name).at("tau_pt_bins");
    auto m_sv_bins_0jet = binning_json.at(config_name).at("m_sv_bins_0jet");
    auto higgs_pT_bins_boost = binning_json.at(config_name).at("higgs_pT_bins_boost");
    auto m_sv_bins_boost = binning_json.at(config_name).at("m_sv_bins_boost");
    auto vbf_cat_x_bins = binning_json.at(config_name).at("vbf_cat_x_bins");
    auto vbf_cat_y_bins = binning_json.at(config_name).at("vbf_cat_y_bins");
    auto vbf_cat_edges = binning_json.at(config_name).at("vbf_cat_edges");

    // boilerplate
    std::unordered_map<string, string> syst_name_map;
    bp_json.at("syst_name_map").get_to(syst_name_map);

    vector<string> output_file_directories, vbf_directories;
    bp_json.at("categories").get_to(output_file_directories);

    vector<string> fake_factor_systematics;
    bp_json.at("fake_factor_systematics").get_to(fake_factor_systematics);

    std::string vbf_sep_var = vbf_cat_x_bins.at(0);
    vector<string> vbf_cat_keys =
        vbf_sep_var.find("D0") == string::npos ? vector<string>{"vbf_sub_cats"} : vector<string>{"vbf_sub_cats_plus", "vbf_sub_cats_minus"};

    auto DCP_idx(0);
    auto n_edges(vbf_cat_edges.at(1).size() - 1);
    vector<string> temp_cats;
    for (auto &k : vbf_cat_keys) {
        auto tmp_idx(0);
        temp_cats.clear();
        bp_json.at(k).get_to(temp_cats);
        for (auto &t : temp_cats) {
            if (tmp_idx == n_edges) {
                break;
            }
            vbf_directories.push_back(t);
            ++tmp_idx;
        }
    }

    vector<string> all_output_directories;
    all_output_directories.reserve(output_file_directories.size() + vbf_directories.size());
    all_output_directories.insert(all_output_directories.end(), output_file_directories.begin(), output_file_directories.end());
    all_output_directories.insert(all_output_directories.end(), vbf_directories.begin(), vbf_directories.end());

    for (auto d : all_output_directories) {
        fout->cd();
        fout->mkdir((channel + "_" + d).c_str());
    }
    fout->cd();

    vector<TH2F *> histograms;
    for (auto &fp : file_paths) {
        // only process nominal unless user requested all systematics
        if (!do_syst && fp.first.find("nominal") == string::npos) {
            continue;
        }

        // make sure we know how to map this systematic
        std::string orig_syst_name("");
        if (syst_name_map.find(fp.first) != syst_name_map.end()) {
            orig_syst_name = syst_name_map.at(fp.first);
        } else if (fp.first != "nominal") {
            std::cerr << "\t \033[91m[INFO]  " << fp.first << " is unknown. Skipping...\033[0m" << std::endl;
            return -1;
        }

        orig_syst_name = std::regex_replace(orig_syst_name, std::regex("YEAR"), year);
        orig_syst_name = std::regex_replace(orig_syst_name, std::regex("LEP"), (channel == "et" ? "ele" : "mu"));
        orig_syst_name = std::regex_replace(orig_syst_name, std::regex("CHAN"), channel);

        std::cout << fp.first << " -> " << orig_syst_name << std::endl;
        for (auto &file : fp.second) {
            if (is_ztt && file.find("embed") != string::npos) {
                continue;
            } else if (!is_ztt && file.find("ZTT") != string::npos) {
                continue;
            }

            std::string syst_name = orig_syst_name;
            if (file.find("embed") != string::npos) {
                if (syst_name.find("CMS_tauideff") != string::npos) {
                    syst_name = std::regex_replace(syst_name, std::regex("tauideff"), "eff_t_embedded");
                } else if (syst_name.find("CMS_scale_e_") != string::npos) {
                    syst_name = std::regex_replace(syst_name, std::regex("scale_e_"), "scale_emb_e");
                } else if ((syst_name.find("CMS_single") != string::npos && syst_name.find("trg") != string::npos) ||
                           syst_name.find("tautrg_") != string::npos) {
                    syst_name = std::regex_replace(syst_name, std::regex("trg"), "trg_emb");
                } else if (syst_name.find("CMS_scale_t_") != string::npos) {
                    syst_name = std::regex_replace(syst_name, std::regex("scale_t_"), "scale_emb_t_");
                }
            }

            string name = std::regex_replace(file, std::regex(".root"), "") + syst_name;
            std::cout << name << std::endl;

            fout->cd((channel + "_0jet").c_str());
            histograms.push_back(build_histogram(name, tau_pt_bins, m_sv_bins_0jet));

            fout->cd((channel + "_boosted").c_str());
            histograms.push_back(build_histogram(name, higgs_pT_bins_boost, m_sv_bins_boost));

            fout->cd((channel + "_vbf").c_str());
            histograms.push_back(build_histogram(name, vbf_cat_x_bins.at(1), vbf_cat_y_bins.at(1)));

            for (auto &d : vbf_directories) {
                fout->cd((channel + "_" + d).c_str());
                histograms.push_back(build_histogram(name, vbf_cat_x_bins.at(1), vbf_cat_y_bins.at(1)));
            }

            // provide nominal fake weight for jetFakes
            auto fake_weight_name = name.find("jetFakes") != string::npos ? "fake_weight" : "None";
            process_file(dir + "/" + fp.first + "/" + file, channel, histograms, vbf_cat_x_bins.at(0), vbf_cat_y_bins.at(0), vbf_cat_edges.at(0),
                         vbf_cat_edges.at(1), fake_weight_name, DCP_idx);

            for (unsigned j = 0; j < histograms.size(); j++) {
                fout->cd((channel + "_" + all_output_directories.at(j)).c_str());
                histograms.at(j)->Write();
            }
            histograms.clear();

            if (fake_weight_name != "None" && do_syst) {
                for (auto &s : fake_factor_systematics) {
                    // make sure we know how to map this systematic
                    std::string syst_name("");
                    if (syst_name_map.find(s) != syst_name_map.end()) {
                        syst_name = syst_name_map.at(s);
                    } else if (fp.first != "nominal") {
                        std::cerr << "\t \033[91m[INFO]  " << fp.first << " is unknown. Skipping...\033[0m" << std::endl;
                        return -1;
                    }

                    syst_name = std::regex_replace(syst_name, std::regex("YEAR"), year);
                    syst_name = std::regex_replace(syst_name, std::regex("LEP"), (channel == "et" ? "ele" : "mu"));
                    syst_name = std::regex_replace(syst_name, std::regex("CHAN"), channel);

                    std::cout << "jetFakes " << syst_name << std::endl;

                    histograms.clear();
                    fout->cd((channel + "_0jet").c_str());
                    histograms.push_back(build_histogram("jetFakes_" + s, tau_pt_bins, m_sv_bins_0jet));

                    fout->cd((channel + "_boosted").c_str());
                    histograms.push_back(build_histogram("jetFakes_" + s, higgs_pT_bins_boost, m_sv_bins_boost));

                    fout->cd((channel + "_vbf").c_str());
                    histograms.push_back(build_histogram("jetFakes_" + s, vbf_cat_x_bins.at(1), vbf_cat_y_bins.at(1)));

                    for (auto &d : vbf_directories) {
                        fout->cd((channel + "_" + d).c_str());
                        histograms.push_back(build_histogram("jetFakes_" + s, vbf_cat_x_bins.at(1), vbf_cat_y_bins.at(1)));
                    }
                    process_file(dir + "/" + fp.first + "/" + file, channel, histograms, vbf_cat_x_bins.at(0), vbf_cat_y_bins.at(0),
                                 vbf_cat_edges.at(0), vbf_cat_edges.at(1), s, DCP_idx);

                    for (unsigned j = 0; j < histograms.size(); j++) {
                        fout->cd((channel + "_" + all_output_directories.at(j)).c_str());
                        histograms.at(j)->Write();
                    }
                    histograms.clear();
                }
            }
        }
    }
    fout->Close();
    std::cout << "Processing time: " << watch.RealTime() << std::endl;
}

void process_file(string name, string channel, vector<TH2F *> histograms, string xvar_name, string yvar_name, string zvar_name, vector<double> edges,
                  string fake_weight_name, int DCP_idx) {
    auto ifile = TFile::Open(name.c_str());
    auto tree = reinterpret_cast<TTree *>(ifile->Get((channel + "_tree").c_str()));
    Int_t isolation, contamination;
    Double_t fake_weight, NN_disc;
    Float_t evtwt, njets, mjj, t1_pt, m_sv, higgs_pt;

    std::unordered_map<string, Float_t> vbf_vars = {{"NN_disc", 0},  {"MELA_D2j", 0},   {"D0_ggH", 0},  {"D_a2_VBF", 0}, {"D0_VBF", 0},
                                                    {"D_l1_VBF", 0}, {"D_l1zg_VBF", 0}, {"DCP_ggH", 0}, {"DCP_VBF", 0}};

    bool is_jetFakes(false);
    string isolation_name("is_signal");
    if (name.find("jetFakes") != string::npos) {
        if (fake_weight_name == "None") {
            std::cerr << "\t \033[91m[INFO] can't process jetFakes without providing: fake_weight_name \033[0m" << std::endl;
        }
        is_jetFakes = true;
        isolation_name = "is_antiTauIso";
    }

    tree->SetBranchAddress("evtwt", &evtwt);
    tree->SetBranchAddress(isolation_name.c_str(), &isolation);
    tree->SetBranchAddress("contamination", &contamination);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("mjj", &mjj);
    tree->SetBranchAddress("t1_pt", &t1_pt);
    tree->SetBranchAddress("m_sv", &m_sv);
    tree->SetBranchAddress("higgs_pT", &higgs_pt);
    tree->SetBranchAddress("NN_disc", &NN_disc);
    tree->SetBranchAddress("MELA_D2j", &vbf_vars.at("MELA_D2j"));
    tree->SetBranchAddress("D0_ggH", &vbf_vars.at("D0_ggH"));
    tree->SetBranchAddress("D_a2_VBF", &vbf_vars.at("D_a2_VBF"));
    tree->SetBranchAddress("D0_VBF", &vbf_vars.at("D0_VBF"));
    tree->SetBranchAddress("D_l1_VBF", &vbf_vars.at("D_l1_VBF"));
    tree->SetBranchAddress("D_l1zg_VBF", &vbf_vars.at("D_l1zg_VBF"));
    tree->SetBranchAddress("DCP_ggH", &vbf_vars.at("DCP_ggH"));
    tree->SetBranchAddress("DCP_VBF", &vbf_vars.at("DCP_VBF"));

    if (fake_weight_name != "None") {
        tree->SetBranchAddress(fake_weight_name.c_str(), &fake_weight);
    }

    std::string dcp_name("None");
    if (xvar_name == "D0_ggH") {
        dcp_name = "DCP_ggH";
    } else if (xvar_name == "D0_VBF") {
        dcp_name = "DCP_VBF";
    }

    double final_evtwt(1.);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (isolation < 1 || contamination > 0) {
            continue;
        }

        vbf_vars["NN_disc"] = NN_disc;
        final_evtwt = is_jetFakes ? evtwt * fake_weight : evtwt;

        if (njets == 0) {
            histograms.at(0)->Fill(t1_pt, m_sv, final_evtwt);
        } else if (njets == 1 || (njets > 1 && mjj < 300)) {
            histograms.at(1)->Fill(higgs_pt, m_sv, final_evtwt);
        } else if (njets > 1 && mjj > 300) {
            histograms.at(2)->Fill(vbf_vars.at(xvar_name), vbf_vars.at(yvar_name), final_evtwt);
            if (edges.size() > 0) {
                for (auto j = 0; j < edges.size() - 1; j++) {
                    if (vbf_vars.at(zvar_name) < edges.at(j + 1)) {
                        auto curr_idx = 3 + j;
                        if (DCP_idx > 0 && vbf_vars.at(dcp_name) < 0) {
                            curr_idx += DCP_idx;
                        }
                        histograms.at(curr_idx)->Fill(vbf_vars.at(xvar_name), vbf_vars.at(yvar_name), final_evtwt);
                        break;
                    }
                }
            }
        }
    }
    ifile->Close();
}

TH2F *build_histogram(string name, vector<double> x_bins, vector<double> y_bins) {
    double *xbins = &x_bins[0];
    double *ybins = &y_bins[0];
    return new TH2F(name.c_str(), name.c_str(), x_bins.size() - 1, xbins, y_bins.size() - 1, ybins);
}

string format_output_name(string channel, bool is_ztt, bool is_syst, string year, int month, int day, string suffix) {
    auto ztt_name = is_ztt ? "ztt" : "emb";
    auto syst_suffix = is_syst ? "Sys" : "noSys";

    return "Output/templates/2D_htt_" + channel + "_" + ztt_name + "_" + syst_suffix + "_fa3_" + year + "_" + std::to_string(month) + "-" +
           std::to_string(day) + suffix + ".root";
}

std::unordered_map<string, vector<string>> build_file_paths(string dir) {
    // read all files from input directory
    vector<string> directories;
    read_directory(dir, &directories);

    vector<string> files;
    std::unordered_map<string, vector<string>> file_paths;
    for (string d : directories) {
        if (d == "." || d == "..") {
            continue;
        }
        files.clear();
        read_directory(dir + "/" + d, &files, ".root");
        file_paths[d] = files;
    }

    return file_paths;
}

// read all *.root files in the given directory and put them in the provided vector
void read_directory(const string &name, vector<string> *v, string match) {
    DIR *dirp = opendir(name.c_str());
    struct dirent *dp;
    while ((dp = readdir(dirp)) != 0) {
        if (match == "" || static_cast<string>(dp->d_name).find(match) != string::npos) {
            v->push_back(dp->d_name);
        }
    }
    closedir(dirp);
}
