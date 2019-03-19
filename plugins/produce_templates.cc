// Copyright 2019 Tyler Mitchell

#include <memory>

// user includes
#include "../include/CLParser.h"
#include "../include/sample_template.h"
#include "TMath.h"
#include "TStopwatch.h"

using std::string;
using std::vector;

int main(int argc, char *argv[]) {
  auto watch = TStopwatch();
  // get CLI arguments
  CLParser parser(argc, argv);
  bool doSyst = parser.Flag("-s");
  string dir = parser.Option("-d");
  string year = parser.Option("-y");
  string suffix = parser.Option("--suf");
  string tree_name = parser.Option("-t");
  string ff_name = parser.Option("-f");

  // get input file directory
  if (dir.empty()) {
    std::cerr << "You must give an input directory" << std::endl;
    return -1;
  }

  // get channel info
  string channel_prefix, lep_charge;
  if (tree_name.find("etau_tree") != string::npos) {
    channel_prefix = "et";
  } else if (tree_name.find("mutau_tree") != string::npos) {
    channel_prefix = "mt";
  } else if (tree_name.find("tautau_tree") != string::npos) {
    channel_prefix = "tt";
  } else {
    std::cerr << "Um. I don't know that tree. Sorry...";
    return -1;
  }

  // read all files from input directory
  vector<string> files;
  read_directory(dir, &files);

  // make output file and TemplateTool containing useful information
  auto fout = std::make_shared<TFile>(("Output/templates/" + channel_prefix + year + "_" + suffix + ".root").c_str(), "recreate");
  auto info = std::make_unique<TemplateTool>(channel_prefix);
  info->make_extension_map(channel_prefix);

  // make all of the TDirectoryFiles we need
  for (auto cat : info->get_categories()) {
    fout->cd();
    fout->mkdir(cat.c_str());
  }
  fout->cd();

  // loop through all files in the directory
  for (auto ifile : files) {
    // get the sample name
    auto name = ifile.substr(0, ifile.find("."));
    std::cout << name << std::endl;

    // open the file
    auto fin = std::unique_ptr<TFile>(TFile::Open((dir + "/" + ifile).c_str()));

    // run for nominal case first
    auto tree = std::shared_ptr<TTree>(reinterpret_cast<TTree *>(fin->Get(tree_name.c_str())));  // open TTree
    auto sample = std::make_unique<Sample_Template>(channel_prefix, year, name, suffix, fout);   // create Sample_Template
    sample->load_fake_fractions(ff_name);                                                        // load fractions from input file
    sample->fill_histograms(tree, doSyst);                                                       // do event loop and fill histos
    sample->write_histograms(doSyst);                                                            // write all to the file
    sample->Close();                                                                             // delete the ff_weight pointer

    // run all systematics stored in tree
    if (doSyst) {
      for (auto key : (*fin->GetListOfKeys())) {
        std::string syst_tree_name = key->GetName();
        if (syst_tree_name.find(tree_name) != std::string::npos && syst_tree_name != tree_name) {
          std::cout << syst_tree_name << std::endl;
          auto ext = info->get_extension(syst_tree_name);
          auto tree = std::shared_ptr<TTree>(reinterpret_cast<TTree *>(fin->Get(syst_tree_name.c_str())));
          auto sample = std::make_unique<Sample_Template>(channel_prefix, year, name, suffix, fout, ext);
          sample->load_fake_fractions(ff_name);
          sample->fill_histograms(tree, doSyst);
          sample->write_histograms(doSyst);
          sample->Close();
        }
      }
    }

    // AC reweighting for JHU samples only
    if (ifile.find("_inc.root") != std::string::npos) {
      auto ac_weights = info->get_AC_weights(ifile);
      for (auto ac_weight : ac_weights) {
        auto fin = std::unique_ptr<TFile>(TFile::Open((dir + "/" + ifile).c_str()));
        auto tree = std::shared_ptr<TTree>(reinterpret_cast<TTree *>(fin->Get(tree_name.c_str())));
        auto jhu_sample = std::make_unique<Sample_Template>(channel_prefix, year, ac_weight.second, suffix, fout);
        jhu_sample->load_fake_fractions(ff_name);
        jhu_sample->fill_histograms(tree, doSyst, ac_weight.first);
        jhu_sample->write_histograms(doSyst, ac_weight.second);
        jhu_sample->Close();
      }
    }
  }
}