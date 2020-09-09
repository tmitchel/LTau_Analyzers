// Copyright [2020] Tyler Mitchell

#include <fstream>
#include <iostream>

#include "../include/CLParser.h"
#include "../include/json.hpp"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

int main(int argc, char *argv[]) {
    CLParser parser(argc, argv);
    std::string input_name = parser.Option("-n");
    std::string tree_name = parser.Option("-t");
    std::string output_path = parser.Option("-o");

    // open file for processing
    auto fin = TFile::Open(input_name.c_str());
    auto tree = reinterpret_cast<TTree *>(fin->Get(tree_name.c_str()));
    auto counts = reinterpret_cast<TH1D *>(fin->Get("nevents"));

    // process json config
    std::ifstream config_file("configs/boilerplate.json");
    nlohmann::json config_json;
    config_file >> config_json;

    // use file name to get correct keys for json config
    std::string signal_map("None");
    if (input_name.find("madgraph") != std::string::npos) {
        signal_map = "mg_ac_reweighting_map";
    } else if (input_name.find("JHU") != std::string::npos) {
        signal_map = "jhu_ac_reweighting_map";
    }

    std::string signal("None");
    if (input_name.find("ggh125") != std::string::npos) {
        signal = "ggh";
    } else if (input_name.find("vbf125") != std::string::npos) {
        signal = "vbf";
    } else if (input_name.find("wh125") != std::string::npos) {
        signal = "wh";
    } else if (input_name.find("zh125") != std::string::npos) {
        signal = "zh";
    }

    // get the set of <weight, out_name> pairs to process
    std::vector<std::pair<std::string, std::string>> reweighting_map;
    config_json.at(signal_map).at(signal).get_to(reweighting_map);

    Float_t evtwt, new_evtwt, coupling_weight;

    // start reweighting the file to different coupling scenarios
    for (auto coupling : reweighting_map) {
        // create output file
        auto fout = new TFile((output_path + "/" + coupling.second + ".root").c_str(), "RECREATE");

        // copy old file EXCEPT evtwt branch
        tree->SetBranchStatus("evtwt", 0);
        auto new_tree = tree->CloneTree(-1, "fast");
        tree->SetBranchStatus("evtwt", 1);

        // create new evtwt in new tree
        auto new_evtwt_branch = new_tree->Branch("evtwt", &new_evtwt, "evtwt/F");

        // register variables needed to compute new branch
        tree->SetBranchAddress("evtwt", &evtwt);
        new_tree->SetBranchAddress(coupling.first.c_str(), &coupling_weight);

        // loop through entries to fill new branch
        Long64_t nentries = new_tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            new_tree->GetEntry(i);
            new_evtwt = evtwt * coupling_weight;  // new event weight for this coupling scenario
            new_evtwt_branch->Fill();
        }
        new_tree->Write();
        fout->Close();
    }

    // close the original file
    fin->Close();
}
