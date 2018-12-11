
// system includes
#include <dirent.h>
#include <sys/types.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>

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

enum categories { zeroJet,
                  boosted,
                  vbf,
                  vbf_D0_0p0to0p2,
                  vbf_D0_0p2to0p4,
                  vbf_D0_0p4to0p8,
                  vbf_D0_0p8to1p0,
                  vbf_D0_0p0to0p2_DCPp,
                  vbf_D0_0p2to0p4_DCPp,
                  vbf_D0_0p4to0p8_DCPp,
                  vbf_D0_0p8to1p0_DCPp,
                  vbf_D0_0p0to0p2_DCPm,
                  vbf_D0_0p2to0p4_DCPm,
                  vbf_D0_0p4to0p8_DCPm,
                  vbf_D0_0p8to1p0_DCPm };

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
class histHolder {
 public:
  histHolder(std::string, std::string, std::string, bool, bool);
  ~histHolder() { delete ff_weight; }
  void writeHistos();
  void initVectors(std::string);
  void initSystematics(std::string);
  void fillFraction(int, std::string, double, double, double);
  void convertDataToFake(TH2F *, std::string, double, double, double);
  void histoLoop(std::vector<std::string>, std::string, std::string, std::string);
  void getJetFakes(std::vector<std::string>, std::string, std::string, bool);

  bool doNN, old_selection;
  TFile *fout;
  FakeFactor *ff_weight;
  std::string channel_prefix;
  std::vector<std::string> categories, systematics;
  std::vector<float> mvis_bins, njets_bins;
  std::map<std::string, std::string> acNameMap;
  std::map<std::string, std::vector<TH2F *>> hists, FF_systs;
  TH2F *fake_0jet, *fake_boosted, *fake_vbf;
  std::vector<TH2F *> data, fakes, frac_w, frac_tt, frac_real, frac_qcd;

  // binning
  std::vector<Float_t> bins_l2, bins_hpt, bins_mjj, bins_lpt, bins_msv1, bins_msv2;
};

// histHolder contructor to create the output file, the qcd histograms with the correct binning
// and the map from categories to vectors of TH2F*'s. Each TH2F* in the vector corresponds to
// one file that is being put into that categories directory in the output tempalte
histHolder::histHolder(std::string channel_prefix, std::string year, std::string suffix = "final", bool doNN = false, bool old = false) :
  fout( new TFile(("Output/templates/template_"+channel_prefix+year+"_"+suffix+".root").c_str(), "recreate") ),
  mvis_bins({0, 50, 80, 100, 110, 120, 130, 150, 170, 200, 250, 1000}),
  njets_bins({-0.5, 0.5, 1.5, 15}),
  // x-axis
  bins_l2 {0, 1, 10, 11},
  bins_hpt {0, 100, 150, 200, 250, 300, 5000},
  bins_mjj {300, 700, 1100, 1500, 10000},

  // y-axis
  bins_lpt {0, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 400},
  bins_msv1 {0, 80, 90, 100, 110, 120, 130, 140, 150, 160, 300},
  bins_msv2 {0, 95, 115, 135, 155, 400},
  channel_prefix(channel_prefix),
  doNN(doNN),
  old_selection(old),
  acNameMap {
    {"wt_ggH_a1", "reweighted_GGH2Jets_0PM125"},
    {"wt_ggH_a3", "reweighted_GGH2Jets_0M125"},
    {"wt_ggH_a3int", "reweighted_GGH2Jets_0Mf05ph0125"},
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
  categories {
    channel_prefix+"_0jet",
    channel_prefix+"_boosted",
    channel_prefix+"_vbf",
    channel_prefix+"_vbf_D0_0p0to0p2",
    channel_prefix+"_vbf_D0_0p2to0p4",
    channel_prefix+"_vbf_D0_0p4to0p8",
    channel_prefix+"_vbf_D0_0p8to1p0",
    channel_prefix+"_vbf_D0_0p0to0p2_DCPp",
    channel_prefix+"_vbf_D0_0p2to0p4_DCPp",
    channel_prefix+"_vbf_D0_0p4to0p8_DCPp",
    channel_prefix+"_vbf_D0_0p8to1p0_DCPp",
    channel_prefix+"_vbf_D0_0p0to0p2_DCPm",
    channel_prefix+"_vbf_D0_0p2to0p4_DCPm",
    channel_prefix+"_vbf_D0_0p4to0p8_DCPm",
    channel_prefix+"_vbf_D0_0p8to1p0_DCPm"   },
  systematics {
    "ff_qcd_syst_up"            , "ff_qcd_syst_down"           , "ff_qcd_dm0_njet0_stat_up"   ,
    "ff_qcd_dm0_njet0_stat_down", "ff_qcd_dm0_njet1_stat_up"   , "ff_qcd_dm0_njet1_stat_down" ,
    "ff_qcd_dm1_njet0_stat_up"  , "ff_qcd_dm1_njet0_stat_down" , "ff_qcd_dm1_njet1_stat_up"   ,
    "ff_qcd_dm1_njet1_stat_down", "ff_w_syst_up"               , "ff_w_syst_down"             , "ff_w_dm0_njet0_stat_up",
    "ff_w_dm0_njet0_stat_down"  , "ff_w_dm0_njet1_stat_up"     , "ff_w_dm0_njet1_stat_down"   ,
    "ff_w_dm1_njet0_stat_up"    , "ff_w_dm1_njet0_stat_down"   , "ff_w_dm1_njet1_stat_up"     ,
    "ff_w_dm1_njet1_stat_down"  , "ff_tt_syst_up"              , "ff_tt_syst_down"            , "ff_tt_dm0_njet0_stat_up",
    "ff_tt_dm0_njet0_stat_down" , "ff_tt_dm0_njet1_stat_up"    , "ff_tt_dm0_njet1_stat_down"  ,
    "ff_tt_dm1_njet0_stat_up"   , "ff_tt_dm1_njet0_stat_down"  , "ff_tt_dm1_njet1_stat_up"    ,  "ff_tt_dm1_njet1_stat_down"
  }
{
  if (doNN) {
    bins_mjj = {0., 0.1, 0.5, 0.9, 1.};
  }

  // Create empty histograms for each category to fill later.
  for (auto cat : categories) {
    hists[cat.c_str()] = std::vector<TH2F *>();
    FF_systs[cat.c_str()] = std::vector<TH2F *>();

    data.push_back(new TH2F(("data_" + cat).c_str(), ("data_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_w.push_back(new TH2F(("frac_w_" + cat).c_str(), ("frac_w_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_tt.push_back(new TH2F(("frac_tt_" + cat).c_str(), ("frac_tt_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_real.push_back(new TH2F(("frac_real_" + cat).c_str(), ("frac_real_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));
    frac_qcd.push_back(new TH2F(("frac_qcd_" + cat).c_str(), ("frac_qcd_" + cat).c_str(), mvis_bins.size() - 1, &mvis_bins[0], njets_bins.size() - 1, &njets_bins[0]));

    if (cat.find("0jet") != std::string::npos) {
      fakes.push_back(new TH2F("fake_0jet", "fake_SS", bins_l2.size() - 1, &bins_l2[0], bins_lpt.size() - 1, &bins_lpt[0]));
    } else if (cat.find("boosted") != std::string::npos) {
      fakes.push_back(new TH2F("fake_boosted", "fake_SS", bins_hpt.size() - 1, &bins_hpt[0], bins_msv1.size() - 1, &bins_msv1[0]));
    } else {
      fakes.push_back(new TH2F(("fake_" + cat).c_str(), "fake_SS", bins_mjj.size() - 1, &bins_mjj[0], bins_msv2.size() - 1, &bins_msv2[0]));
    }
  }

  for (auto it = hists.begin(); it != hists.end(); it++) {
    fout->cd();
    fout->mkdir((it->first).c_str());
    fout->cd();
  }

  // get FakeFactor workspace
  TFile *ff_file;
  if (year == "2017") {
    ff_file = new TFile(("${CMSSW_BASE}/src/HTTutilities/Jet2TauFakes/data/SM2017/tight/vloose/" + channel_prefix + "/fakeFactors.root").c_str(), "READ");
  } else if (year == "2016") {
    ff_file = new TFile(("${CMSSW_BASE}/src/HTTutilities/Jet2TauFakes/data/SM2016_ML/tight/" + channel_prefix + "/fakeFactors_20180831_tight.root").c_str(), "READ");
  } else {
    std::cerr << "Bad year" << std::endl;
  }
  ff_weight = reinterpret_cast<FakeFactor *>(ff_file->Get("ff_comb"));
  ff_file->Close();
}

// change to the correct output directory then create a new TH1F that will be filled for the current input file
void histHolder::initVectors(std::string name) {
  for (auto key : hists) {
    fout->cd(key.first.c_str());
    if (name.find("Data") != std::string::npos) {
      name = "data_obs";
    }
    if (key.first == channel_prefix + "_0jet") {
      hists.at(key.first.c_str()).push_back(new TH2F(name.c_str(), name.c_str(), bins_l2.size() - 1, &bins_l2[0], bins_lpt.size() - 1, &bins_lpt[0]));
    } else if (key.first == channel_prefix + "_boosted") {
      hists.at(key.first.c_str()).push_back(new TH2F(name.c_str(), name.c_str(), bins_hpt.size() - 1, &bins_hpt[0], bins_msv1.size() - 1, &bins_msv1[0]));
    } else if (key.first.find("_vbf") != std::string::npos) {
      hists.at(key.first.c_str()).push_back(new TH2F(name.c_str(), name.c_str(), bins_mjj.size() - 1, &bins_mjj[0], bins_msv2.size() - 1, &bins_msv2[0]));
    }
  }
}

// change to the correct output directory then create a new TH1F that will be filled for the current input file
void histHolder::initSystematics(std::string name) {
  for (auto key : FF_systs) {
    fout->cd(key.first.c_str());
    std::string name = "jetFakes_";
    for (auto syst : systematics) {
      if (key.first == channel_prefix + "_0jet") {
        FF_systs.at(key.first.c_str()).push_back(new TH2F((name + syst).c_str(), name.c_str(), bins_l2.size() - 1, &bins_l2[0], bins_lpt.size() - 1, &bins_lpt[0]));
      } else if (key.first == channel_prefix + "_boosted") {
        FF_systs.at(key.first.c_str()).push_back(new TH2F((name + syst).c_str(), name.c_str(), bins_hpt.size() - 1, &bins_hpt[0], bins_msv1.size() - 1, &bins_msv1[0]));
      } else if (key.first.find("_vbf") != std::string::npos) {
        FF_systs.at(key.first.c_str()).push_back(new TH2F((name + syst).c_str(), name.c_str(), bins_mjj.size() - 1, &bins_mjj[0], bins_msv2.size() - 1, &bins_msv2[0]));
      }
    }
  }
  std::cout << "initialized systematics" << std::endl;
}

void histHolder::fillFraction(int cat, std::string name, double var1, double var2, double weight) {
  TH2F *hist;
  if (name == "Data") {
    hist = data.at(cat);
  } else if (name == "W" || name == "ZJ" || name == "VVJ") {
    hist = frac_w.at(cat);
  } else if (name == "TTJ") {
    hist = frac_tt.at(cat);
  } else if (name == "ZTT" || name == "TTT" || name == "VVT") {
    hist = frac_real.at(cat);
  }
  hist->Fill(var1, var2, weight);
}

void histHolder::convertDataToFake(TH2F *hist, std::string name, double var1, double var2, double weight) {
  if (name.find("Data") != std::string::npos) {
    hist->Fill(var1, var2, weight);
  }
}

// write output histograms including the QCD histograms after scaling by OS/SS ratio
void histHolder::writeHistos() {
  for (auto cat : hists) {
    fout->cd(cat.first.c_str());
    for (auto hist : cat.second) {
      std::cout << hist->GetName() << std::endl;
      hist->Write();
    }
  }

  for (auto cat = 0; cat < fakes.size(); cat++) {
    fout->cd(categories.at(cat));
    auto fake_hist = fakes.at(i);
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
  }

  // fout->cd((channel_prefix + "_0jet").c_str());
  // fake_0jet->SetName("jetFakes");
  // for (auto i = 0; i < fake_0jet->GetNbinsX(); i++) {
  //   for (auto j = 0; j < fake_0jet->GetNbinsY(); j++) {
  //     if (fake_0jet->GetBinContent(i, j) < 0) {
  //       fake_0jet->SetBinContent(i, j, 0);
  //     }
  //   }
  // }
  // fake_0jet->Write();

  // for (auto &hist : FF_systs.at(channel_prefix + "_0jet")) {
  //   for (auto i = 0; i < hist->GetNbinsX(); i++) {
  //     for (auto j = 0; j < hist->GetNbinsY(); j++) {
  //       if (hist->GetBinContent(i, j) < 0) {
  //         hist->SetBinContent(i, j, 0);
  //       }
  //     }
  //   }
  //   hist->Write();
  // }

  // fout->cd((channel_prefix + "_boosted").c_str());
  // fake_boosted->SetName("jetFakes");
  // for (auto i = 0; i < fake_boosted->GetNbinsX(); i++) {
  //   for (auto j = 0; j < fake_boosted->GetNbinsY(); j++) {
  //     if (fake_boosted->GetBinContent(i, j) < 0) {
  //       fake_boosted->SetBinContent(i, j, 0);
  //     }
  //   }
  // }
  // fake_boosted->Write();

  // for (auto &hist : FF_systs.at(channel_prefix + "_boosted")) {
  //   for (auto i = 0; i < hist->GetNbinsX(); i++) {
  //     for (auto j = 0; j < hist->GetNbinsY(); j++) {
  //       if (hist->GetBinContent(i, j) < 0) {
  //         hist->SetBinContent(i, j, 0);
  //       }
  //     }
  //   }
  //   hist->Write();
  // }

  // fout->cd((channel_prefix + "_vbf").c_str());
  // fake_vbf->SetName("jetFakes");
  // for (auto i = 0; i < fake_vbf->GetNbinsX(); i++) {
  //   for (auto j = 0; j < fake_vbf->GetNbinsY(); j++) {
  //     if (fake_vbf->GetBinContent(i, j) < 0) {
  //       fake_vbf->SetBinContent(i, j, 0);
  //     }
  //   }
  // }
  // fake_vbf->Write();

  // for (auto &hist : FF_systs.at(channel_prefix + "_vbf")) {
  //   for (auto i = 0; i < hist->GetNbinsX(); i++) {
  //     for (auto j = 0; j < hist->GetNbinsY(); j++) {
  //       if (hist->GetBinContent(i, j) < 0) {
  //         hist->SetBinContent(i, j, 0);
  //       }
  //     }
  //   }
  //   hist->Write();
  // }

  fout->Close();
}
