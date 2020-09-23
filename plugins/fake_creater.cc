// Copyright [2020] Tyler Mitchell

#include <unordered_map>
#include <vector>

#include "../include/ApplyFF.h"
#include "../include/CLParser.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

std::vector<std::string> ff_syst = {"ff_qcd_0jet_unc1", "ff_qcd_0jet_unc2",    "ff_qcd_1jet_unc1",  "ff_qcd_1jet_unc2",   "ff_qcd_2jet_unc1",
                                    "ff_qcd_2jet_unc2", "ff_w_0jet_unc1",      "ff_w_0jet_unc2",    "ff_w_1jet_unc1",     "ff_w_1jet_unc2",
                                    "ff_w_2jet_unc1",   "ff_w_2jet_unc2",      "ff_tt_0jet_unc1",   "ff_tt_0jet_unc2",    "mtclosure_w_unc1",
                                    "mtclosure_w_unc2", "lptclosure_xtrg_qcd", "lptclosure_xtrg_w", "lptclosure_xtrg_tt", "lptclosure_qcd",
                                    "lptclosure_w",     "lptclosure_tt",       "osssclosure_qcd"};

class fake_map {
   private:
    std::vector<std::string> categories;
    std::unordered_map<std::string, std::vector<TH1F*>> fractions;

   public:
    explicit fake_map(TFile*);
    ~fake_map() {}

    std::vector<Float_t> get_fractions(std::string, Float_t, Float_t);
};

fake_map::fake_map(TFile* ff_file) : categories({"0jet", "boosted", "vbf"}) {
    for (auto cat : categories) {
        fractions[cat] = std::vector<TH1F*>{
            reinterpret_cast<TH1F*>(ff_file->Get((cat + "/frac_w").c_str())),
            reinterpret_cast<TH1F*>(ff_file->Get((cat + "/frac_tt").c_str())),
            reinterpret_cast<TH1F*>(ff_file->Get((cat + "/frac_qcd").c_str())),
            reinterpret_cast<TH1F*>(ff_file->Get((cat + "/frac_data").c_str())),
        };
    }
}

std::vector<Float_t> fake_map::get_fractions(std::string cat, Float_t x, Float_t y) {
    auto hist = fractions.at(cat).at(3);
    auto xbin = hist->GetXaxis()->FindBin(x);
    auto ybin = hist->GetYaxis()->FindBin(y);

    Float_t frac_w = fractions.at(cat).at(0)->GetBinContent(xbin, ybin);
    Float_t frac_tt = fractions.at(cat).at(1)->GetBinContent(xbin, ybin);
    Float_t frac_qcd = fractions.at(cat).at(2)->GetBinContent(xbin, ybin);

    return std::vector<Float_t>{frac_w, frac_tt, frac_qcd};
}

int main(int argc, char* argv[]) {
    CLParser parser(argc, argv);
    std::string input_path = parser.Option("-i");
    std::string fake_fraction_path = parser.Option("-p");
    std::string fake_factor_path = parser.Option("-f");
    std::string channel = parser.Option("-c");
    bool syst = parser.Flag("-s");

    apply_ff ffer(fake_factor_path, channel);

    std::string lpt_name("mu_pt");
    if (channel == "et") {
        lpt_name = "el_pt";
    }

    // read input file
    TFile* fin = TFile::Open((input_path + "/pre_jetFakes.root").c_str());
    TTree* tree = reinterpret_cast<TTree*>(fin->Get((channel + "_tree").c_str()));

    // read fake fractions
    TFile* ff_file = TFile::Open(fake_fraction_path.c_str());
    auto fractions = fake_map(ff_file);

    // create output file
    TFile* fout = new TFile((input_path + "/jetFakes.root").c_str(), "RECREATE");

    // Clone the original tree
    auto new_tree = tree->CloneTree(-1, "fast");

    // create new evtwt branch
    Float_t fake_weight;
    auto bfake_weight = new_tree->Branch("fake_weight", &fake_weight, "fake_weight/F");

    // systematic branches
    std::vector<Float_t> fake_weight_systs_up(23, 1), fake_weight_systs_dn(23, 1);
    std::vector<TBranch*> bfake_weight_systs_up, bfake_weight_systs_dn;
    if (syst) {
        for (auto i = 0; i < ff_syst.size(); i++) {
            bfake_weight_systs_up.push_back(new_tree->Branch((ff_syst.at(i) + "_up").c_str(), &fake_weight_systs_up[i]));
            bfake_weight_systs_dn.push_back(new_tree->Branch((ff_syst.at(i) + "_down").c_str(), &fake_weight_systs_dn[i]));
        }
    }

    // set addresses for inputs to apply_ff
    Float_t pt, mt, vis_mass, lpt, dr, met, njets, xtrg, mjj;
    tree->SetBranchAddress("t1_pt", &pt);
    tree->SetBranchAddress("mt", &mt);
    tree->SetBranchAddress("vis_mass", &vis_mass);
    tree->SetBranchAddress(lpt_name.c_str(), &lpt);
    tree->SetBranchAddress("lep_dr", &dr);
    tree->SetBranchAddress("met", &met);
    tree->SetBranchAddress("njets", &njets);
    tree->SetBranchAddress("cross_trigger", &xtrg);
    tree->SetBranchAddress("mjj", &mjj);

    // used to pass fractions to apply_ff
    Float_t frac_tt, frac_qcd, frac_w;
    std::string event_cat;
    std::vector<Float_t> event_fractions;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // categorize event to get correct fraction
        if (njets == 0) {
            event_cat = "0jet";
        } else if (njets == 1 || (njets > 1 && mjj <= 300)) {
            event_cat = "boosted";
        } else if (njets > 1 && mjj > 300) {
            event_cat = "vbf";
        }

        // get fractions
        event_fractions = fractions.get_fractions(event_cat, vis_mass, njets);
        frac_w = event_fractions.at(0);
        frac_tt = event_fractions.at(1);
        frac_qcd = event_fractions.at(2);

        // fill the weights
        fake_weight = ffer.get_ff(std::vector<Float_t>{pt, mt, vis_mass, lpt, dr, met, njets, xtrg, frac_tt, frac_qcd, frac_w});
        bfake_weight->Fill();

        if (syst) {
            for (auto j = 0; j < ff_syst.size(); j++) {
                fake_weight_systs_up[j] =
                    ffer.get_ff(std::vector<Float_t>{pt, mt, vis_mass, lpt, dr, met, njets, xtrg, frac_tt, frac_qcd, frac_w}, ff_syst.at(j), "up");
                bfake_weight_systs_up.at(j)->Fill();

                fake_weight_systs_dn[j] =
                    ffer.get_ff(std::vector<Float_t>{pt, mt, vis_mass, lpt, dr, met, njets, xtrg, frac_tt, frac_qcd, frac_w}, ff_syst.at(j), "down");
                bfake_weight_systs_dn.at(j)->Fill();
            }
        }
    }
    new_tree->Write();
    fout->Close();
    fin->Close();
    ff_file->Close();
}
