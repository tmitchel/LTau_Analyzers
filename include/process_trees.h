
// system includes
#include <dirent.h>
#include <sys/types.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"

// FF
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"
#include "HTTutilities/Jet2TauFakes/interface/IFunctionWrapper.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTFormula.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTGraph.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH2F.h"
#include "HTTutilities/Jet2TauFakes/interface/WrapperTH3D.h"

enum Categories { zeroJet,
                  boosted,
                  vbf,
                  vbf_ggHMELA_bin1_NN_bin1,
                  vbf_ggHMELA_bin2_NN_bin1,
                  vbf_ggHMELA_bin3_NN_bin1,
                  vbf_ggHMELA_bin4_NN_bin1,
                  vbf_ggHMELA_bin5_NN_bin1,
                  vbf_ggHMELA_bin6_NN_bin1,
                  vbf_ggHMELA_bin7_NN_bin1,
                  vbf_ggHMELA_bin8_NN_bin1,
                  vbf_ggHMELA_bin9_NN_bin1,
                  vbf_ggHMELA_bin10_NN_bin1,
                  vbf_ggHMELA_bin11_NN_bin1,
                  vbf_ggHMELA_bin12_NN_bin1};

// read all *.root files in the given directory and put them in the provided vector
void read_directory(const std::string &name, std::vector<std::string> *v) {
  DIR *dirp = opendir(name.c_str());
  struct dirent *dp;
  while ((dp = readdir(dirp)) != 0) {
    if (static_cast<std::string>(dp->d_name).find("root") != std::string::npos) {
      v->push_back(dp->d_name);
    }
  }
  closedir(dirp);
}

// class to hold the histograms until I'm ready to write them
class HistTool {
 public:
  HistTool(std::string, std::string, std::string, bool, bool);
  ~HistTool() { delete ff_weight; }
  void writeHistos();
  void writeTemplates();
  void initVectors2d(std::string);
  void initSystematics(std::string);
  void fillFraction(int, std::string, double, double, double);
  void convertDataToFake(Categories, std::string, double, double, double, double, double, double, double, double, double);  // 2d
  void histoLoop(std::vector<std::string>, std::string, std::string, std::string);
  void getJetFakes(std::vector<std::string>, std::string, std::string, bool);
  Categories getCategory(double, double);

  bool doNN, old_selection;
  TFile *fout;
  FakeFactor *ff_weight;
  std::string channel_prefix, var;
  std::vector<std::string> categories, systematics;
  std::vector<float> mvis_bins, njets_bins;
  std::map<std::string, std::string> acNameMap;
  std::map<std::string, std::vector<TH2F *>> hists_2d, FF_systs;
  std::vector<TH2F *> data, fakes_2d, frac_w, frac_tt, frac_real, frac_qcd;

  // binning
  std::vector<Float_t> bins_l2, bins_hpt, bins_vbf_var1, bins_lpt, bins_msv1, bins_vbf_var2, bins_hpt2;
};

// HistTool contructor to create the output file, the qcd histograms with the correct binning
// and the map from categories to vectors of TH2F*'s. Each TH2F* in the vector corresponds to
// one file that is being put into that categories directory in the output tempalte
HistTool::HistTool(std::string channel_prefix, std::string year, std::string suffix = "final", bool doNN = false, bool old = false)
    : fout(new TFile(("Output/templates/" + channel_prefix + year + "_" + suffix + ".root").c_str(), "recreate")),
      mvis_bins({0, 50, 80, 100, 110, 120, 130, 150, 170, 200, 250, 1000}),
      njets_bins({-0.5, 0.5, 1.5, 15}),
      // x-axis
      bins_l2{0, 1, 10, 11},
      bins_hpt{0, 100, 150, 200, 250, 300, 5000},
      // bins_vbf_var1{300, 500, 10000},  // real mjj
      bins_vbf_var1{0, 0.25, 0.5, 0.75, 1.},  // actually VBF MELA

      // y-axis
      bins_lpt{0, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 400},
      bins_msv1{0, 80, 90, 100, 110, 120, 130, 140, 150, 160, 300},
      // bins_vbf_var2{0, 80, 100, 115, 130, 150, 1000},
      // bins_vbf_var2{-5, -1.25, -0.75, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2., 3.},  // Fisher Disc
      // bins_vbf_var2{0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5,  0.6, 0.65, 0.7, 0.8},  // Perceptron
      // bins_vbf_var2{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1.},  // NN including m_sv et2016/mt2017
      // bins_vbf_var2{0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 0.95, 1.}, // mt2016
      bins_vbf_var2{0, 0.25, 0.5, 0.75, 1.},
      channel_prefix(channel_prefix),
      acNameMap{
          {"wt_ggH_a1", "JHU_GGH2Jets_sm_M125"},
          {"wt_ggH_a3", "JHU_GGH2Jets_pseudoscalar_M125"},
          {"wt_ggH_a3int", "JHU_GGH2Jets_pseudoscalar_Mf05ph0125"},
          {"wt_wh_a1", "reweighted_WH_htt_0PM125"},
          {"wt_wh_a2", "reweighted_WH_htt_0PH125"},
          {"wt_wh_a2int", "reweighted_WH_htt_0PHf05ph0125"},
          {"wt_wh_a3", "reweighted_WH_htt_0M125"},
          {"wt_wh_a3int", "reweighted_WH_htt_0Mf05ph0125"},
          {"wt_wh_L1", "reweighted_WH_htt_0L1125"},
          {"wt_wh_L1int", "reweighted_WH_htt_0L1f05ph0125"},
          {"wt_wh_L1Zg", "reweighted_WH_htt_0L1Zg125"},
          {"wt_wh_L1Zgint", "reweighted_WH_htt_0L1Zgf05ph0125"},
          {"wt_zh_a1", "reweighted_ZH_htt_0PM125"},
          {"wt_zh_a2", "reweighted_ZH_htt_0PH125"},
          {"wt_zh_a2int", "reweighted_ZH_htt_0PHf05ph0125"},
          {"wt_zh_a3", "reweighted_ZH_htt_0M125"},
          {"wt_zh_a3int", "reweighted_ZH_htt_0Mf05ph0125"},
          {"wt_zh_L1", "reweighted_ZH_htt_0L1125"},
          {"wt_zh_L1int", "reweighted_ZH_htt_0L1f05ph0125"},
          {"wt_zh_L1Zg", "reweighted_ZH_htt_0L1Zg125"},
          {"wt_zh_L1Zgint", "reweighted_ZH_htt_0L1Zgf05ph0125"},
          {"wt_a1", "reweighted_qqH_htt_0PM125"},
          {"wt_a2", "reweighted_qqH_htt_0PH125"},
          {"wt_a2int", "reweighted_qqH_htt_0PHf05ph0125"},
          {"wt_a3", "reweighted_qqH_htt_0M125"},
          {"wt_a3int", "reweighted_qqH_htt_0Mf05ph0125"},
          {"wt_L1", "reweighted_qqH_htt_0L1125"},
          {"wt_L1int", "reweighted_qqH_htt_0L1f05ph0125"},
          {"wt_L1Zg", "reweighted_qqH_htt_0L1Zg125"},
          {"wt_L1Zgint", "reweighted_qqH_htt_0L1Zgf05ph0125"},
      },
      categories{
          channel_prefix + "_0jet",
          channel_prefix + "_boosted",
          channel_prefix + "_vbf",
          channel_prefix + "_vbf_ggHMELA_bin1_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin2_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin3_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin4_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin5_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin6_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin7_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin8_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin9_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin10_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin11_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin12_NN_bin1",
          channel_prefix + "_vbf_ggHMELA_bin1_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin2_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin3_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin4_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin5_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin6_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin7_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin8_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin9_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin10_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin11_NN_bin2",
          channel_prefix + "_vbf_ggHMELA_bin12_NN_bin2",},
      systematics{
          "ff_qcd_syst_up", "ff_qcd_syst_down", "ff_qcd_dm0_njet0_stat_up",
          "ff_qcd_dm0_njet0_stat_down", "ff_qcd_dm0_njet1_stat_up", "ff_qcd_dm0_njet1_stat_down",
          "ff_qcd_dm1_njet0_stat_up", "ff_qcd_dm1_njet0_stat_down", "ff_qcd_dm1_njet1_stat_up",
          "ff_qcd_dm1_njet1_stat_down", "ff_w_syst_up", "ff_w_syst_down", "ff_w_dm0_njet0_stat_up",
          "ff_w_dm0_njet0_stat_down", "ff_w_dm0_njet1_stat_up", "ff_w_dm0_njet1_stat_down",
          "ff_w_dm1_njet0_stat_up", "ff_w_dm1_njet0_stat_down", "ff_w_dm1_njet1_stat_up",
          "ff_w_dm1_njet1_stat_down", "ff_tt_syst_up", "ff_tt_syst_down", "ff_tt_dm0_njet0_stat_up",
          "ff_tt_dm0_njet0_stat_down", "ff_tt_dm0_njet1_stat_up", "ff_tt_dm0_njet1_stat_down",
          "ff_tt_dm1_njet0_stat_up", "ff_tt_dm1_njet0_stat_down", "ff_tt_dm1_njet1_stat_up", "ff_tt_dm1_njet1_stat_down"} {

  // Create empty histograms for each category to fill later.
  for (auto cat : categories) {
    // make a 2d template
    hists_2d[cat.c_str()] = std::vector<TH2F *>();

    if (cat.find("0jet") != std::string::npos) {
      fakes_2d.push_back(new TH2F("fake_0jet", "fake_SS", bins_l2.size() - 1, &bins_l2[0], bins_lpt.size() - 1, &bins_lpt[0]));
    } else if (cat.find("boosted") != std::string::npos) {
      fakes_2d.push_back(new TH2F("fake_boosted", "fake_SS", bins_hpt.size() - 1, &bins_hpt[0], bins_msv1.size() - 1, &bins_msv1[0]));
    } else {
      fakes_2d.push_back(new TH2F(("fake_" + cat).c_str(), "fake_SS", bins_vbf_var1.size() - 1, &bins_vbf_var1[0], bins_vbf_var2.size() - 1, &bins_vbf_var2[0]));
    }

    // histograms for fake-factor are always 2d
    FF_systs[cat.c_str()] = std::vector<TH2F *>();

    data.push_back(new TH2F(("data_" + cat).c_str(), ("data_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_w.push_back(new TH2F(("frac_w_" + cat).c_str(), ("frac_w_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_tt.push_back(new TH2F(("frac_tt_" + cat).c_str(), ("frac_tt_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_real.push_back(new TH2F(("frac_real_" + cat).c_str(), ("frac_real_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_qcd.push_back(new TH2F(("frac_qcd_" + cat).c_str(), ("frac_qcd_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
  }

  // make all of the directories for templates
  for (auto it = hists_2d.begin(); it != hists_2d.end(); it++) {
    fout->cd();
    fout->mkdir((it->first).c_str());
    fout->cd();
  }

  // get FakeFactor workspace
  TFile *ff_file;
  if (year == "2017") {
    ff_file = new TFile(("${CMSSW_BASE}/src/HTTutilities/Jet2TauFakes/data2017/SM2017/tight/vloose/" + channel_prefix + "/fakeFactors.root").c_str(), "READ");
  } else if (year == "2016") {
    ff_file = new TFile(("${CMSSW_BASE}/src/HTTutilities/Jet2TauFakes/data2016/SM2016_ML/tight/" + channel_prefix + "/fakeFactors_tight.root").c_str(), "READ");
  } else {
    std::cerr << "Bad year" << std::endl;
  }
  ff_weight = reinterpret_cast<FakeFactor *>(ff_file->Get("ff_comb"));
  ff_file->Close();
}

// change to the correct output directory then create a new TH2F that will be filled for the current input file
void HistTool::initVectors2d(std::string name) {
  for (auto key : hists_2d) {
    fout->cd(key.first.c_str());
    if (name.find("Data") != std::string::npos) {
      name = "data_obs";
    }
    if (key.first == channel_prefix + "_0jet") {
      hists_2d.at(key.first.c_str()).push_back(new TH2F(name.c_str(), name.c_str(), bins_l2.size() - 1, &bins_l2[0], bins_lpt.size() - 1, &bins_lpt[0]));
    } else if (key.first == channel_prefix + "_boosted") {
      hists_2d.at(key.first.c_str()).push_back(new TH2F(name.c_str(), name.c_str(), bins_hpt.size() - 1, &bins_hpt[0], bins_msv1.size() - 1, &bins_msv1[0]));
    } else if (key.first.find("_vbf") != std::string::npos) {
      hists_2d.at(key.first.c_str()).push_back(new TH2F(name.c_str(), name.c_str(), bins_vbf_var1.size() - 1, &bins_vbf_var1[0], bins_vbf_var2.size() - 1, &bins_vbf_var2[0]));
    }
  }
}

// change to the correct output directory then create a new TH1F that will be filled for the current input file
void HistTool::initSystematics(std::string name) {
  for (auto key : FF_systs) {
    fout->cd(key.first.c_str());
    std::string name = "jetFakes_";
    for (auto syst : systematics) {
      if (key.first == channel_prefix + "_0jet") {
        FF_systs.at(key.first.c_str()).push_back(new TH2F((name + syst).c_str(), name.c_str(), bins_l2.size() - 1, &bins_l2[0], bins_lpt.size() - 1, &bins_lpt[0]));
      } else if (key.first == channel_prefix + "_boosted") {
        FF_systs.at(key.first.c_str()).push_back(new TH2F((name + syst).c_str(), name.c_str(), bins_hpt.size() - 1, &bins_hpt[0], bins_msv1.size() - 1, &bins_msv1[0]));
      } else if (key.first.find("_vbf") != std::string::npos) {
        FF_systs.at(key.first.c_str()).push_back(new TH2F((name + syst).c_str(), name.c_str(), bins_vbf_var1.size() - 1, &bins_vbf_var1[0], bins_vbf_var2.size() - 1, &bins_vbf_var2[0]));
      }
    }
  }
  std::cout << "initialized systematics" << std::endl;
}

void HistTool::fillFraction(int cat, std::string name, double var1, double var2, double weight) {
  TH2F *hist;
  if (name == "Data") {
    hist = data.at(cat);
  } else if (name == "W" || name == "ZJ" || name == "VVJ") {
    hist = frac_w.at(cat);
  } else if (name == "TTJ") {
    hist = frac_tt.at(cat);
  } else if (name == "embedded" || name == "TTT" || name == "VVT") {
    hist = frac_real.at(cat);
  }
  hist->Fill(var1, var2, weight);
}

void HistTool::convertDataToFake(Categories cat, std::string name, double var1, double var2, double vis_mass, double njets, double t1_pt, double t1_decayMode, double mt, double lep_iso, double weight) {
  auto bin_x = data.at(cat)->GetXaxis()->FindBin(vis_mass);
  auto bin_y = data.at(cat)->GetYaxis()->FindBin(njets);
  auto fakeweight = ff_weight->value({t1_pt, t1_decayMode, njets, vis_mass, mt, lep_iso,
                                      frac_w.at(cat)->GetBinContent(bin_x, bin_y),
                                      frac_tt.at(cat)->GetBinContent(bin_x, bin_y),
                                      frac_qcd.at(cat)->GetBinContent(bin_x, bin_y)});
  if (name.find("Data") != std::string::npos) {
    fakes_2d.at(cat)->Fill(var1, var2, weight * fakeweight);
  }
}

// write output histograms including the QCD histograms after scaling by OS/SS ratio
void HistTool::writeTemplates() {
  for (auto cat : hists_2d) {
    fout->cd(cat.first.c_str());
    for (auto hist : cat.second) {
      hist->Write();
    }
  }

  for (auto cat = 0; cat < fakes_2d.size(); cat++) {
    fout->cd(categories.at(cat).c_str());
    auto fake_hist = fakes_2d.at(cat);
    fake_hist->SetName("jetFakes");

    // if fake yield is negative, make it zero
    for (auto i = 0; i < fake_hist->GetNbinsX(); i++) {
      for (auto j = 0; j < fake_hist->GetNbinsY(); j++) {
        if (fake_hist->GetBinContent(i, j) < 0) {
          fake_hist->SetBinContent(i, j, 0);
        }
      }
    }
    fake_hist->Write();

    for (auto &hist : FF_systs.at(categories.at(cat))) {
      for (auto i = 0; i < hist->GetNbinsX(); i++) {
        for (auto j = 0; j < hist->GetNbinsY(); j++) {
          if (hist->GetBinContent(i, j) < 0) {
            hist->SetBinContent(i, j, 0);
          }
        }
      }
      hist->Write();
    }
  }
}

// basically a map from 2 inputs -> 1 Category
Categories HistTool::getCategory(double vbf_var3, double vbf_var4 = -1) {
  double edge = 1./6.;
  if (vbf_var3 > 0 && vbf_var3 <= 1.*edge) {
    return vbf_ggHMELA_bin1_NN_bin1;
  } else if (vbf_var3 <= 2.*edge) {
    return vbf_ggHMELA_bin2_NN_bin1;
  } else if (vbf_var3 <= 3.*edge) {
    return vbf_ggHMELA_bin3_NN_bin1;
  } else if (vbf_var3 <= 4.*edge) {
    return vbf_ggHMELA_bin4_NN_bin1;
  } else if (vbf_var3 <= 5.*edge) {
    return vbf_ggHMELA_bin5_NN_bin1;
  } else if (vbf_var3 <= 6.*edge) {
    return vbf_ggHMELA_bin6_NN_bin1;
  }

//  else if (D0_ggH <= 7.*edge) {
//    return vbf_ggHMELA_bin7_NN_bin1;
//  } else if (D0_ggH <= 8.*edge) {
//    return vbf_ggHMELA_bin8_NN_bin1;
//  } else if (D0_ggH <= 9.*edge) {
//    return vbf_ggHMELA_bin9_NN_bin1;
//  } else if (D0_ggH <= 10.*edge) {
//    return vbf_ggHMELA_bin10_NN_bin1;
//  } else if (D0_ggH <= 11.*edge) {
//    return vbf_ggHMELA_bin11_NN_bin1;
//  } else if (D0_ggH <= 12.*edge) {
//    return vbf_ggHMELA_bin12_NN_bin1;
//  }
}