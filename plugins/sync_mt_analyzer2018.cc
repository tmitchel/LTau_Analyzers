// Copyright [2019] Tyler Mitchell

// system includes
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

// ROOT includes
#include "RooFunctor.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"

// user includes
#include "../include/ACWeighter.h"
#include "../include/CLParser.h"
#include "../include/ComputeWG1Unc.h"
#include "../include/LumiReweightingStandAlone.h"
#include "../include/event_info.h"
#include "../include/jet_factory.h"
#include "../include/met_factory.h"
#include "../include/muon_factory.h"
#include "../include/slim_tree.h"
#include "../include/swiss_army_class.h"
#include "../include/tau_factory.h"

typedef std::vector<double> NumV;

int main(int argc, char *argv[]) {
    ////////////////////////////////////////////////
    // Initial setup:                             //
    // Get file names, normalization, paths, etc. //
    ////////////////////////////////////////////////

    CLParser parser(argc, argv);
    std::string name = parser.Option("-n");
    std::string path = parser.Option("-p");
    std::string syst = parser.Option("-u");
    std::string sample = parser.Option("-s");
    std::string output_dir = parser.Option("-d");
    std::string signal_type = parser.Option("--stype");
    std::string fname = path + sample + ".root";
    bool isData = sample.find("data") != std::string::npos;
    bool isEmbed = sample.find("embed") != std::string::npos || name.find("embed") != std::string::npos;
    bool isMG = sample.find("madgraph") != std::string::npos;
    bool doAC = signal_type != "None";

    // get systematic shift name
    std::string systname = "NOMINAL";
    if (!syst.empty()) {
        systname = "SYST_" + syst;
    }

    // create output path
    auto suffix = "_output.root";
    auto prefix = "Output/trees/" + output_dir + "/" + systname + "/";
    std::string filename, logname;
    if (name == sample) {
        filename = prefix + name + systname + suffix;
        logname = "Output/trees/" + output_dir + "/logs/" + name + systname + ".txt";
    } else {
        filename = prefix + sample + std::string("_") + name + "_" + systname + suffix;
        logname = "Output/trees/" + output_dir + "/logs/" + sample + std::string("_") + name + "_" + systname + ".txt";
    }

    // create the log file
    std::ofstream logfile;
    logfile.open(logname, std::ios::out | std::ios::trunc);

    // open log file and log some things
    logfile << "Opening file... " << sample << std::endl;
    logfile << "With name...... " << name << std::endl;
    logfile << "And running systematic " << systname << std::endl;
    logfile << "Using options: " << std::endl;
    logfile << "\t name: " << name << std::endl;
    logfile << "\t path: " << path << std::endl;
    logfile << "\t syst: " << syst << std::endl;
    logfile << "\t sample: " << sample << std::endl;
    logfile << "\t output_dir: " << output_dir << std::endl;
    logfile << "\t signal_type: " << signal_type << std::endl;
    logfile << "\t isData: " << isData << " isEmbed: " << isEmbed << " doAC: " << doAC << std::endl;

    auto fin = TFile::Open(fname.c_str());
    auto ntuple = reinterpret_cast<TTree *>(fin->Get("mutau_tree"));

    // get number of generated events
    auto counts = reinterpret_cast<TH1D *>(fin->Get("nevents"));
    auto gen_number = counts->GetBinContent(2);

    auto fout = new TFile(filename.c_str(), "RECREATE");
    counts->Write();
    fout->mkdir("grabbag");
    fout->cd("grabbag");

    // initialize Helper class
    Helper *helper = new Helper(fout, name, syst);

    // cd to root of output file and create tree
    fout->cd();
    slim_tree *st = new slim_tree("mt_tree", doAC);

    std::string original = sample;
    if (name == "VBF125") {
        sample = "vbf125";
    } else if (name == "ggH125" && signal_type != "madgraph") {
        sample = "ggh125";
    } else if (name == "WH125") {
        sample = "wh125";
    } else if (name == "WHsigned125") {
        sample = sample.find("plus") == std::string::npos ? "wplus125" : "wminus125";
    } else if (name == "ZH125") {
        sample = "zh125";
    }

    if (signal_type == "JHU" && (sample == "ggh125" || sample == "vbf125")) {
        gen_number = 1.;
    }

    // reweighter for anomolous coupling samples
    ACWeighter ac_weights = ACWeighter(original, sample, signal_type, "2018");
    ac_weights.fillWeightMap();

    // get normalization (lumi & xs are in util.h)
    double norm(1.);
    if (!isData && !isEmbed) {
        norm = helper->getLuminosity2018() * helper->getCrossSection(sample) / gen_number;
    }

    ///////////////////////////////////////////////
    // Scale Factors:                            //
    // Read weights, hists, graphs, etc. for SFs //
    ///////////////////////////////////////////////

    auto lumi_weights =
        new reweight::LumiReWeighting("/hdfs/store/user/tmitchel/HTT_ScaleFactors/pu_distributions_mc_2018.root",
                                      "/hdfs/store/user/tmitchel/HTT_ScaleFactors/pu_distributions_data_2018.root", "pileup", "pileup");

    // legacy sf's
    TFile htt_sf_file("/hdfs/store/user/tmitchel/HTT_ScaleFactors/htt_scalefactors_legacy_2018.root");
    RooWorkspace *htt_sf = reinterpret_cast<RooWorkspace *>(htt_sf_file.Get("w"));
    htt_sf_file.Close();

    // MadGraph Higgs pT file
    RooWorkspace *mg_sf;
    if (signal_type == "madgraph") {
        TFile mg_sf_file("/hdfs/store/user/tmitchel/HTT_ScaleFactors/htt_scalefactors_2017_MGggh.root");
        mg_sf = reinterpret_cast<RooWorkspace *>(mg_sf_file.Get("w"));
        mg_sf_file.Close();
    }

    TFile *f_NNLOPS = new TFile("/hdfs/store/user/tmitchel/HTT_ScaleFactors/NNLOPS_reweight.root");
    TGraph *g_NNLOPS_0jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_powheg_0jet"));
    TGraph *g_NNLOPS_1jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_powheg_1jet"));
    TGraph *g_NNLOPS_2jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_powheg_2jet"));
    TGraph *g_NNLOPS_3jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_powheg_3jet"));

    //////////////////////////////////////
    // Final setup:                     //
    // Declare histograms and factories //
    //////////////////////////////////////

    // declare histograms (histogram initializer functions in util.h)
    fout->cd("grabbag");
    auto histos = helper->getHistos1D();

    // construct factories
    event_info event(ntuple, lepton::MUON, 2018, isMG, syst);
    muon_factory muons(ntuple, 2018, syst);
    tau_factory taus(ntuple, 2018, syst);
    jet_factory jets(ntuple, 2018, syst);
    met_factory met(ntuple, 2018, syst);

    if (sample == "ggh125" && signal_type == "powheg") {
        event.setRivets(ntuple);
    }

    // begin the event loop
    Int_t nevts = ntuple->GetEntries();
    int progress(0), fraction((nevts - 1) / 10);
    for (Int_t i = 0; i < nevts; i++) {
        ntuple->GetEntry(i);
        if (i == progress * fraction) {
            logfile << "LOG: Processing: " << progress * 10 << "% complete." << std::endl;
            progress++;
        }

        // find the event weight (not lumi*xs if looking at W or Drell-Yan)
        Float_t evtwt(norm), corrections(1.), sf_trig(1.), sf_id(1.), sf_iso(1.), sf_reco(1.);
        if (name == "W") {
            if (event.getNumGenJets() == 1) {
                evtwt = 9.679;
            } else if (event.getNumGenJets() == 2) {
                evtwt = 4.808;
            } else if (event.getNumGenJets() == 3) {
                evtwt = 3.290;
            } else if (event.getNumGenJets() == 4) {
                evtwt = 3.435;
            } else {
                evtwt = 55.160;
            }
        }

        if (name == "ZTT" || name == "ZLL" || name == "ZL" || name == "ZJ") {
            if (event.getNumGenJets() == 1) {
                evtwt = 0.671;
            } else if (event.getNumGenJets() == 2) {
                evtwt = 0.588;
            } else if (event.getNumGenJets() == 3) {
                evtwt = 0.640;
            } else if (event.getNumGenJets() == 4) {
                evtwt = 1.172;
            } else {
                evtwt = 3.865;
            }
        }
        histos->at("cutflow")->Fill(1., 1.);

        // run factories
        auto muon = muons.run_factory();
        auto tau = taus.run_factory();
        jets.run_factory();

        // event flags
        if (event.getPassFlags(isData)) {
            histos->at("cutflow")->Fill(2., 1.);
        } else {
            continue;
        }

        // Separate processes
        if ((name == "ZL" || name == "TTL" || name == "VVL" || name == "STL") && tau.getGenMatch() > 4) {
            continue;
        } else if ((name == "ZTT" || name == "TTT" || name == "VVT" || name == "STT") && tau.getGenMatch() != 5) {
            continue;
        } else if ((name == "ZJ" || name == "TTJ" || name == "VVJ" || name == "STJ") && tau.getGenMatch() != 6) {
            continue;
        } else {
            histos->at("cutflow")->Fill(3., 1.);
        }

        // only opposite-sign
        int evt_charge = tau.getCharge() + muon.getCharge();
        if (evt_charge == 0) {
            histos->at("cutflow")->Fill(4., 1.);
        } else {
            continue;
        }

        // build Higgs
        TLorentzVector Higgs = muon.getP4() + tau.getP4() + met.getP4();

        // calculate mt
        double met_x = met.getMet() * cos(met.getMetPhi());
        double met_y = met.getMet() * sin(met.getMetPhi());
        double met_pt = sqrt(pow(met_x, 2) + pow(met_y, 2));
        double mt = sqrt(pow(muon.getPt() + met_pt, 2) - pow(muon.getP4().Px() + met_x, 2) - pow(muon.getP4().Py() + met_y, 2));

        // now do mt selection
        if (mt < 50) {
            histos->at("cutflow")->Fill(5., 1.);
        } else {
            continue;
        }

        if (tau.getAgainstVVLooseElectronDeep() > 0) {
            histos->at("cutflow")->Fill(5., 1.);
        } else {
            continue;
        }

        // b-jet veto
        if (jets.getNbtagLoose() < 2 && jets.getNbtagMedium() < 1) {
            histos->at("cutflow")->Fill(6., 1.);
        } else {
            continue;
        }

        // create regions
        bool signalRegion = (tau.getMediumIsoDeep() && muon.getIso() < 0.15);
        bool antiTauIsoRegion = (tau.getMediumIsoDeep() == 0 && tau.getVVVLooseIsoDeep() > 0 && muon.getIso() < 0.15);

        if (signal_type != "None") {
            antiTauIsoRegion = false;  // don't need anti-tau iso region in signal
        }

        // only keep the regions we need
        if (signalRegion || antiTauIsoRegion) {
            histos->at("cutflow")->Fill(7., 1.);
        } else {
            continue;
        }

        // apply all scale factors/corrections/etc.
        if (!isData && !isEmbed) {
            // pileup reweighting
            evtwt *= lumi_weights->weight(event.getNPU());

            // generator weights
            evtwt *= event.getGenWeight();

            // b-tagging scale factor goes here
            evtwt *= jets.getBWeight();

            // set workspace variables
            htt_sf->var("m_pt")->setVal(muon.getPt());
            htt_sf->var("m_eta")->setVal(muon.getEta());
            htt_sf->var("t_pt")->setVal(tau.getPt());
            htt_sf->var("t_eta")->setVal(tau.getEta());
            htt_sf->var("t_phi")->setVal(tau.getPhi());
            htt_sf->var("t_dm")->setVal(tau.getDecayMode());
            htt_sf->var("z_gen_mass")->setVal(event.getGenM());
            htt_sf->var("z_gen_pt")->setVal(event.getGenPt());

            // start applying weights from workspace
            evtwt *= htt_sf->function("m_trk_ratio")->getVal();
            evtwt *= htt_sf->function("m_idiso_ic_ratio")->getVal();

            // tau ID efficiency SF and systematics
            std::string id_name = "t_deeptauid_pt_medium";  // nominal
            if (syst.find("tau_id_") != std::string::npos) {
                id_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
            }
            evtwt *= htt_sf->function(id_name.c_str())->getVal();

            // muon fake rate SF
            if (tau.getDecayMode() == 2 || tau.getDecayMode() == 4) {
                evtwt *= htt_sf->function("t_id_vs_mu_eta_tight")->getVal();
            }

            // trigger scale factors
            if (muon.getPt() < 25) {  // cross-trigger
                evtwt *= htt_sf->function("m_trg_20_ic_ratio")->getVal();
                if (syst == "trigger_up") {
                    evtwt *= htt_sf->function("t_trg_pog_deeptau_medium_mutau_ratio_up")->getVal();
                } else if (syst == "trigger_down") {
                    evtwt *= htt_sf->function("t_trg_pog_deeptau_medium_mutau_ratio_down")->getVal();
                } else {
                    evtwt *= htt_sf->function("t_trg_pog_deeptau_medium_mutau_ratio")->getVal();
                }
            } else {  // muon trigger
                evtwt *= htt_sf->function("m_trg_ic_ratio")->getVal();
            }

            // Z-pT Reweighting
            if (name == "EWKZ2l" || name == "EWKZ2nu" || name == "ZTT" || name == "ZLL" || name == "ZL" || name == "ZJ") {
                auto nom_zpt_weight = htt_sf->function("zptmass_weight_nom")->getVal();
                if (syst == "dyShape_Up") {
                    nom_zpt_weight = nom_zpt_weight + ((nom_zpt_weight - 1) * 0.1);
                } else if (syst == "dyShape_Down") {
                    nom_zpt_weight = nom_zpt_weight - ((nom_zpt_weight - 1) * 0.1);
                }
                evtwt *= nom_zpt_weight;
            }

            // top-pT Reweighting
            if (name == "TTT" || name == "TTJ" || name == "TTL" || name == "STT" || name == "STJ" || name == "STL") {
                float pt_top1 = std::min(static_cast<float>(400.), jets.getTopPt1());
                float pt_top2 = std::min(static_cast<float>(400.), jets.getTopPt2());
                if (syst == "ttbarShape_Up") {
                    evtwt *= (2 * sqrt(exp(0.0615 - 0.0005 * pt_top1) * exp(0.0615 - 0.0005 * pt_top2)) - 1);  // 2*√[e^(..)*e^(..)] - 1
                } else if (syst == "ttbarShape_Up") {
                    // no weight for shift down
                } else {
                    evtwt *= sqrt(exp(0.0615 - 0.0005 * pt_top1) * exp(0.0615 - 0.0005 * pt_top2));  // √[e^(..)*e^(..)]
                }
            }

            // NNLOPS ggH reweighting
            if (sample == "ggh125" && signal_type == "powheg") {
                if (event.getNjetsRivet() == 0) evtwt *= g_NNLOPS_0jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(125.0)));
                if (event.getNjetsRivet() == 1) evtwt *= g_NNLOPS_1jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(625.0)));
                if (event.getNjetsRivet() == 2) evtwt *= g_NNLOPS_2jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(800.0)));
                if (event.getNjetsRivet() >= 3) evtwt *= g_NNLOPS_3jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(925.0)));
                NumV WG1unc = qcd_ggF_uncert_2017(event.getNjetsRivet(), event.getHiggsPtRivet(), event.getJetPtRivet());
                if (syst.find("Rivet") != std::string::npos) {
                    evtwt *= (1 + event.getRivetUnc(WG1unc, syst));
                }
            }

            // MadGraph Higgs pT correction
            if (signal_type == "madgraph") {
                mg_sf->var("HpT")->setVal(Higgs.Pt());
                evtwt *= mg_sf->function("ggH_quarkmass_corr")->getVal();
            }

            // begin systematics

            // jet to tau fake rate systematics
            if (tau.getGenMatch() == 6 && name == "TTJ" || name == "ZJ" || name == "W" || name == "VVJ") {
                auto temp_tau_pt = std::min(200., static_cast<double>(tau.getPt()));
                if (syst == "jetToTauFake_Up") {
                    evtwt *= (1 - (0.2 * temp_tau_pt / 100));
                } else if (syst == "jetToTauFake_Down") {
                    evtwt *= (1 + (0.2 * temp_tau_pt / 100));
                }
            }
        } else if (!isData && isEmbed) {
            event.setEmbed();

            // embedded generator weights
            auto genweight(event.getGenWeight());
            if (genweight > 1 || genweight < 0) {
                genweight = 0;
            }
            evtwt *= genweight;

            // tracking sf
            evtwt *= helper->embed_tracking(tau.getDecayMode());

            // set workspace variables
            htt_sf->var("m_pt")->setVal(muon.getPt());
            htt_sf->var("m_eta")->setVal(muon.getEta());
            htt_sf->var("t_pt")->setVal(tau.getPt());
            htt_sf->var("t_eta")->setVal(tau.getEta());
            htt_sf->var("t_phi")->setVal(tau.getPhi());
            htt_sf->var("t_dm")->setVal(tau.getDecayMode());
            htt_sf->var("gt1_pt")->setVal(muon.getGenPt());
            htt_sf->var("gt1_eta")->setVal(muon.getGenEta());
            htt_sf->var("gt2_pt")->setVal(tau.getGenPt());
            htt_sf->var("gt2_eta")->setVal(tau.getGenEta());

            // start applying weights from workspace
            evtwt *= htt_sf->function("m_trk_ratio")->getVal();
            evtwt *= htt_sf->function("m_idiso_ic_embed_ratio")->getVal();

            // tau ID eff SF
            std::string embed_id_name = "t_deeptauid_pt_embed_medium";
            if (syst.find("tau_id_") != std::string::npos) {
                embed_id_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
            }
            evtwt *= htt_sf->function(embed_id_name.c_str())->getVal();

            // trigger scale factors
            if (muon.getPt() < 25) {  // cross-trigger
                // muon-leg
                evtwt *= htt_sf->function("m_trg_20_ic_embed_ratio")->getVal();
                // tau-leg
                std::string tau_leg_name("t_trg_mediumDeepTau_mutau_embed_ratio");
                if (syst.find("trigger") != std::string::npos) {
                    tau_leg_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
                }
                evtwt *= htt_sf->function(tau_leg_name.c_str())->getVal();
            } else {  // muon trigger
                evtwt *= htt_sf->function("m_trg_ic_embed_ratio")->getVal();
            }

            // muon fake rate SF
            if (tau.getDecayMode() == 2 || tau.getDecayMode() == 4) {
                evtwt *= htt_sf->function("t_id_vs_mu_eta_tight")->getVal();
            }

            // double muon trigger eff in selection
            evtwt *= htt_sf->function("m_sel_trg_ratio")->getVal();

            // muon ID eff in selection (leg 1)
            htt_sf->var("gt_pt")->setVal(muon.getGenPt());
            htt_sf->var("gt_eta")->setVal(muon.getGenEta());
            evtwt *= htt_sf->function("m_sel_id_ic_ratio")->getVal();

            // muon ID eff in selection (leg 2)
            htt_sf->var("gt_pt")->setVal(tau.getGenPt());
            htt_sf->var("gt_eta")->setVal(tau.getGenEta());
            evtwt *= htt_sf->function("m_sel_id_ic_ratio")->getVal();
        }
        fout->cd();

        std::vector<std::string> tree_cat;

        // regions
        if (signalRegion) {
            tree_cat.push_back("signal");
        } else if (antiTauIsoRegion) {
            tree_cat.push_back("antiTauIso");
        }

        // event charge
        if (evt_charge == 0) {
            tree_cat.push_back("OS");
        }

        std::shared_ptr<std::vector<double>> weights(nullptr);
        Long64_t currentEventID = event.getLumi();
        currentEventID = currentEventID * 1000000 + event.getEvt();
        if (doAC) {
            weights = std::make_shared<std::vector<double>>(ac_weights.getWeights(currentEventID));
        }

        // fill the tree
        st->fillTree(tree_cat, &muon, &tau, &jets, &met, &event, mt, evtwt, weights, name);
    }  // close event loop

    fin->Close();
    fout->cd();
    fout->Write();
    fout->Close();
    logfile << "Finished processing " << sample << std::endl;
    logfile.close();
    return 0;
}
