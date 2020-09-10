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
    std::string output_name = parser.Option("-o");
    std::string weight_name = parser.Option("-w");


    // open file for processing
    auto fin = TFile::Open(input_name.c_str());
    auto tree = reinterpret_cast<TTree *>(fin->Get(tree_name.c_str()));
    auto counts = reinterpret_cast<TH1D *>(fin->Get("nevents"));

    Float_t evtwt, new_evtwt, coupling_weight;

    // create output file
    auto fout = new TFile(output_name.c_str(), "RECREATE");

    // copy old file EXCEPT evtwt branch
    tree->SetBranchStatus("evtwt", 0);
    auto new_tree = tree->CloneTree(-1, "fast");
    tree->SetBranchStatus("evtwt", 1);

    // create new evtwt in new tree
    auto new_evtwt_branch = new_tree->Branch("evtwt", &new_evtwt, "evtwt/F");

    // register variables needed to compute new branch
    tree->SetBranchAddress("evtwt", &evtwt);
    new_tree->SetBranchAddress(weight_name.c_str(), &coupling_weight);

    // loop through entries to fill new branch
    Long64_t nentries = new_tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        new_tree->GetEntry(i);
        new_evtwt = evtwt * coupling_weight;  // new event weight for this coupling scenario
        new_evtwt_branch->Fill();
    }
    counts->Write();
    new_tree->Write();
    fout->Close();

    // close the original file
    fin->Close();
}
