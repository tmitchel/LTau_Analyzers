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
#include <utility>
#include <vector>

#include "../include/CLParser.h"
#include "../include/json.hpp"
#include "TFile.h"
#include "TH2F.h"
#include "TStopwatch.h"
#include "TTree.h"

using std::string;
using std::unordered_map;
using std::vector;

class file_processor {
   private:
    std::shared_ptr<TFile> fout;
    bool is_jetFakes;
    Int_t isolation, contamination;
    Double_t fake_weight, NN_disc;
    Float_t evtwt, njets, mjj, t1_pt, m_sv, higgs_pt;
    vector<string> vbf_cats;
    vector<double> tau_pt_bins, m_sv_bins_0jet, higgs_pT_bins_boost, m_sv_bins_boost, vbf_cat_x_bins, vbf_cat_y_bins;
    unordered_map<string, Float_t> vbf_vars;
    unordered_map<string, Double_t> other_vars;
    unordered_map<string, vector<TH2F *> *> all_histograms;

   public:
    string xvar_name, yvar_name, zvar_name, dcp_name, channel;
    vector<double> edges;

    file_processor(std::shared_ptr<TFile>, string, nlohmann::json, vector<string>);
    void register_branches(TTree *, bool);
    bool register_new_branch(TTree *, string);
    void create_histograms(string);
    void process_file(TTree *, string, int);
    void process_file_with_weights(TTree *, string, int, vector<std::pair<string, string>>, bool);
    void write(vector<string>);
};

void read_directory(const string &, vector<string> *, string match = "");
unordered_map<string, vector<string>> build_file_paths(string);
string format_output_name(string, bool, bool, string, int, int, string);
TH2F *build_histogram(string, vector<double>, vector<double>);

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

    auto vbf_cat_x_bins = binning_json.at(config_name).at("vbf_cat_x_bins");
    auto vbf_cat_y_bins = binning_json.at(config_name).at("vbf_cat_y_bins");
    auto vbf_cat_edges = binning_json.at(config_name).at("vbf_cat_edges");

    // boilerplate
    unordered_map<string, string> syst_name_map;
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

    auto p = new file_processor(fout, channel, binning_json.at(config_name), vbf_directories);

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

        bool is_jetFakes(false), is_embed(false);
        for (auto &file : fp.second) {
            if (file == "embed.root") {
                is_jetFakes = false;
                is_embed = true;
            } else if (file == "jetFakes.root") {
                is_jetFakes = true;
                is_embed = false;
            } else {
                is_jetFakes = false;
                is_embed = false;
            }

            if (is_ztt && is_embed) {
                continue;
            } else if (!is_ztt && file.find("ZTT") != string::npos) {
                continue;
            }

            std::string syst_name = orig_syst_name;
            if (is_embed) {
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

            p->create_histograms(name);

            auto fin = TFile::Open((dir + "/" + fp.first + "/" + file).c_str());
            auto tree = reinterpret_cast<TTree *>(fin->Get((channel + "_tree").c_str()));
            p->register_branches(tree, is_jetFakes);
            p->process_file(tree, name, DCP_idx);
            fin->Close();

            // handle jet systematics
            if (is_jetFakes && do_syst) {
                auto fin = TFile::Open((dir + "/" + fp.first + "/" + file).c_str());
                auto tree = reinterpret_cast<TTree *>(fin->Get((channel + "_tree").c_str()));
                p->register_branches(tree, is_jetFakes);
                vector<std::pair<string, string>> weights;
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

                    auto exists = p->register_new_branch(tree, s);
                    if (exists) {
                        p->create_histograms("jetFakes" + syst_name);
                        weights.push_back(std::make_pair(s, "jetFakes" + syst_name));
                    }
                }
                p->process_file_with_weights(tree, name, DCP_idx, weights, true);
                fin->Close();
            }
        }
    }

    p->write(all_output_directories);
    fout->Close();
    std::cout << "Processing time: " << watch.RealTime() << std::endl;
}

file_processor::file_processor(std::shared_ptr<TFile> _fout, string _channel, nlohmann::json json, vector<string> _vbf_cats)
    : fout(_fout), is_jetFakes(false), dcp_name("None"), vbf_cats(_vbf_cats), channel(_channel) {
    auto in_tau_pt_bins = json.at("tau_pt_bins");
    auto in_m_sv_bins_0jet = json.at("m_sv_bins_0jet");
    auto in_higgs_pT_bins_boost = json.at("higgs_pT_bins_boost");
    auto in_m_sv_bins_boost = json.at("m_sv_bins_boost");
    auto in_vbf_cat_x_bins = json.at("vbf_cat_x_bins");
    auto in_vbf_cat_y_bins = json.at("vbf_cat_y_bins");
    auto in_vbf_cat_edges = json.at("vbf_cat_edges");

    xvar_name = in_vbf_cat_x_bins.at(0);
    yvar_name = in_vbf_cat_y_bins.at(0);
    zvar_name = in_vbf_cat_edges.at(0);
    in_tau_pt_bins.get_to<std::vector<double>>(tau_pt_bins);
    in_m_sv_bins_0jet.get_to<std::vector<double>>(m_sv_bins_0jet);
    in_higgs_pT_bins_boost.get_to<std::vector<double>>(higgs_pT_bins_boost);
    in_m_sv_bins_boost.get_to<std::vector<double>>(m_sv_bins_boost);
    in_vbf_cat_x_bins.at(1).get_to<std::vector<double>>(vbf_cat_x_bins);
    in_vbf_cat_y_bins.at(1).get_to<std::vector<double>>(vbf_cat_y_bins);
    in_vbf_cat_edges.at(1).get_to<std::vector<double>>(edges);

    if (xvar_name == "D0_ggH") {
        dcp_name = "DCP_ggH";
    } else if (xvar_name == "D0_VBF") {
        dcp_name = "DCP_VBF";
    }

    vbf_vars = {{"NN_disc", 0},  {"MELA_D2j", 0},   {"D0_ggH", 0},  {"D_a2_VBF", 0}, {"D0_VBF", 0},
                {"D_l1_VBF", 0}, {"D_l1zg_VBF", 0}, {"DCP_ggH", 0}, {"DCP_VBF", 0}};
}

void file_processor::create_histograms(string name) {
    auto histograms = new vector<TH2F *>();
    fout->cd((channel + "_0jet").c_str());
    histograms->push_back(build_histogram(name, tau_pt_bins, m_sv_bins_0jet));

    fout->cd((channel + "_boosted").c_str());
    histograms->push_back(build_histogram(name, higgs_pT_bins_boost, m_sv_bins_boost));

    fout->cd((channel + "_vbf").c_str());
    histograms->push_back(build_histogram(name, vbf_cat_x_bins, vbf_cat_y_bins));

    for (auto &d : vbf_cats) {
        fout->cd((channel + "_" + d).c_str());
        histograms->push_back(build_histogram(name, vbf_cat_x_bins, vbf_cat_y_bins));
    }

    all_histograms.insert(std::make_pair(name, histograms));
}

void file_processor::register_branches(TTree *tree, bool _is_jetFakes = false) {
    is_jetFakes = _is_jetFakes;
    if (is_jetFakes) {
        tree->SetBranchAddress("is_antiTauIso", &isolation);
        tree->SetBranchAddress("fake_weight", &fake_weight);
    } else {
        tree->SetBranchAddress("is_signal", &isolation);
    }
    tree->SetBranchAddress("evtwt", &evtwt);
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
}

void file_processor::register_new_branch(TTree *tree, string vname) {
    if (tree->GetBranch(vname.c_str()) == nullptr) {
        return false
    }

    other_vars[vname] = 0;
    tree->SetBranchAddress(vname.c_str(), &other_vars.at(vname));
    return true
}

void file_processor::process_file(TTree *tree, string name, int DCP_idx) {
    double final_evtwt(1.);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (isolation < 1 || contamination > 0) {
            continue;
        }

        vbf_vars["NN_disc"] = NN_disc;

        // jetFake systematics include evtwt, others don't
        final_evtwt = include_evtwt ? evtwt : 1.;

        if (njets == 0) {
            all_histograms.at(name)->at(0)->Fill(t1_pt, m_sv, final_evtwt);
        } else if (njets == 1 || (njets > 1 && mjj < 300)) {
            all_histograms.at(name)->at(1)->Fill(higgs_pt, m_sv, final_evtwt);
        } else if (njets > 1 && mjj > 300) {
            all_histograms.at(name)->at(2)->Fill(vbf_vars.at(xvar_name), vbf_vars.at(yvar_name), final_evtwt);
            for (auto j = 0; j < edges.size() - 1; j++) {
                if (vbf_vars.at(zvar_name) < edges.at(j + 1)) {
                    auto curr_idx = 3 + j;
                    if (DCP_idx > 0 && vbf_vars.at(dcp_name) < 0) {
                        curr_idx += DCP_idx;
                    }
                    all_histograms.at(name)->at(curr_idx)->Fill(vbf_vars.at(xvar_name), vbf_vars.at(yvar_name), final_evtwt);
                    break;
                }
            }
        }
    }
}

void file_processor::process_file_with_weights(TTree *tree, string name, int DCP_idx, vector<std::pair<string, string>> weights,
                                               bool include_evtwt = false) {
    double final_evtwt(1.);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (isolation < 1 || contamination > 0) {
            continue;
        }

        vbf_vars["NN_disc"] = NN_disc;

        // jetFake systematics include evtwt, others don't
        final_evtwt = include_evtwt ? evtwt : 1.;

        if (njets == 0) {
            for (auto &w : weights) {
                all_histograms.at(w.second)->at(0)->Fill(t1_pt, m_sv, final_evtwt * other_vars.at(w.first));
            }
        } else if (njets == 1 || (njets > 1 && mjj < 300)) {
            for (auto &w : weights) {
                all_histograms.at(w.second)->at(1)->Fill(higgs_pt, m_sv, final_evtwt * other_vars.at(w.first));
            }
        } else if (njets > 1 && mjj > 300) {
            for (auto &w : weights) {
                all_histograms.at(w.second)->at(2)->Fill(vbf_vars.at(xvar_name), vbf_vars.at(yvar_name), final_evtwt * other_vars.at(w.first));
                for (auto j = 0; j < edges.size() - 1; j++) {
                    if (vbf_vars.at(zvar_name) < edges.at(j + 1)) {
                        auto curr_idx = 3 + j;
                        if (DCP_idx > 0 && vbf_vars.at(dcp_name) < 0) {
                            curr_idx += DCP_idx;
                        }
                        all_histograms.at(w.second)->at(curr_idx)->Fill(vbf_vars.at(xvar_name), vbf_vars.at(yvar_name),
                                                                        final_evtwt * other_vars.at(w.first));
                        break;
                    }
                }
            }
        }
    }
}

void file_processor::write(vector<string> all_output_directories) {
    fout->cd();
    for (auto &name : all_histograms) {
        for (unsigned i = 0; i < name.second->size(); i++) {
            fout->cd((channel + "_" + all_output_directories.at(i)).c_str());
            name.second->at(i)->Write();
        }
    }
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

unordered_map<string, vector<string>> build_file_paths(string dir) {
    // read all files from input directory
    vector<string> directories;
    read_directory(dir, &directories);

    vector<string> files;
    unordered_map<string, vector<string>> file_paths;
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
