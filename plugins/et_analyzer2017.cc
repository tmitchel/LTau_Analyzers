// Copyright [2018] Tyler Mitchell

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
#include "../include/electron_factory.h"
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
    bool condor = parser.Flag("--condor");
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
    bool doAC = signal_type != "None" && signal_type != "powheg";

    std::string systname = "NOMINAL";
    if (!syst.empty()) {
        systname = "SYST_" + syst;
    }

    // create output path
    auto suffix = "_output.root";
    auto prefix = "Output/trees/" + output_dir;
    std::string filename, logname;
    filename = prefix + "/" + systname + "/" + sample + std::string("_") + name + "_" + systname + suffix;
    logname = prefix + "/logs/" + sample + std::string("_") + name + "_" + systname + ".txt";

    if (condor) {
        filename = sample + std::string("_") + name + "_" + systname + suffix;
    }

    // create the log file
    std::ofstream logfile;
    if (!condor) {
        logfile.open(logname, std::ios::out | std::ios::trunc);
    }

    std::ostream &running_log = (condor ? std::cout : logfile);

    // open log file and log some things
    running_log << "Opening file... " << sample << std::endl;
    running_log << "With name...... " << name << std::endl;
    running_log << "And running systematic " << systname << std::endl;
    running_log << "Using options: " << std::endl;
    running_log << "\t name: " << name << std::endl;
    running_log << "\t path: " << path << std::endl;
    running_log << "\t syst: " << syst << std::endl;
    running_log << "\t sample: " << sample << std::endl;
    running_log << "\t output_dir: " << output_dir << std::endl;
    running_log << "\t signal_type: " << signal_type << std::endl;
    running_log << "\t isData: " << isData << " isEmbed: " << isEmbed << " doAC: " << doAC << std::endl;

    auto fin = TFile::Open(fname.c_str());
    auto ntuple = reinterpret_cast<TTree *>(fin->Get("etau_tree"));

    // get number of generated events
    auto counts = reinterpret_cast<TH1D *>(fin->Get("nevents"));
    auto gen_number = counts->GetBinContent(2);

    // create output file
    auto fout = new TFile(filename.c_str(), "RECREATE");
    counts->Write();
    fout->mkdir("grabbag");
    fout->cd("grabbag");

    // initialize Helper class
    Helper *helper = new Helper(fout, name, syst);

    // cd to root of output file and create tree
    fout->cd();
    slim_tree *st = new slim_tree("et_tree", doAC);

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
    ACWeighter ac_weights = ACWeighter(original, sample, signal_type, "2017");
    ac_weights.fillWeightMap();

    // get normalization (lumi & xs are in util.h)
    double norm(1.);
    if (!isData && !isEmbed) {
        norm = helper->getLuminosity2017() * helper->getCrossSection(sample) / gen_number;
    }

    ///////////////////////////////////////////////
    // Scale Factors:                            //
    // Read weights, hists, graphs, etc. for SFs //
    ///////////////////////////////////////////////

    reweight::LumiReWeighting *lumi_weights;
    // read inputs for lumi reweighting
    if (!isData && !isEmbed && !doAC && !isMG) {
        TH1F *dbsName = reinterpret_cast<TH1F *>(fin->Get("MiniAOD_name"));
        std::string datasetName = dbsName->GetTitle();
        if (datasetName.find("Not Found") != std::string::npos && !isEmbed && !isData) {
            fin->Close();
            fout->Close();
            return 2;
        }
        std::replace(datasetName.begin(), datasetName.end(), '/', '#');
        lumi_weights = new reweight::LumiReWeighting("/hdfs/store/user/tmitchel/HTT_ScaleFactors/pu_distributions_mc_2017.root",
                                                     "/hdfs/store/user/tmitchel/HTT_ScaleFactors/pu_distributions_data_2017.root",
                                                     ("pua/#" + datasetName).c_str(), "pileup");
        running_log << "using PU dataset name: " << datasetName << std::endl;
    }

    // legacy sf's
    TFile htt_sf_file("/hdfs/store/user/tmitchel/HTT_ScaleFactors/htt_scalefactors_legacy_2017.root");
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
    TGraph *g_mcatnlo_NNLOPS_0jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_mcatnlo_0jet"));
    TGraph *g_mcatnlo_NNLOPS_1jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_mcatnlo_1jet"));
    TGraph *g_mcatnlo_NNLOPS_2jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_mcatnlo_2jet"));
    TGraph *g_mcatnlo_NNLOPS_3jet = reinterpret_cast<TGraph *>(f_NNLOPS->Get("gr_NNLOPSratio_pt_mcatnlo_3jet"));

    //////////////////////////////////////
    // Final setup:                     //
    // Declare histograms and factories //
    //////////////////////////////////////

    // declare histograms (histogram initializer functions in util.h)
    fout->cd("grabbag");
    auto histos = helper->getHistos1D();

    // construct factories
    event_info event(ntuple, lepton::ELECTRON, 2017, isMG, syst);
    electron_factory electrons(ntuple, 2017, syst);
    tau_factory taus(ntuple, 2017, syst);
    jet_factory jets(ntuple, 2017, syst);
    met_factory met(ntuple, 2017, syst);

    if (signal_type == "powheg" || signal_type == "madgraph") {
        event.setRivets(ntuple);
    }

    // begin the event loop
    Int_t nevts = ntuple->GetEntries();
    int progress(0), fraction((nevts - 1) / 10);
    for (Int_t i = 0; i < nevts; i++) {
        ntuple->GetEntry(i);
        if (i == progress * fraction) {
            running_log << "LOG: Processing: " << progress * 10 << "% complete." << std::endl;
            progress++;
        }

        // find the event weight (not lumi*xs if looking at W or Drell-Yan)
        Float_t evtwt(norm), corrections(1.), sf_trig(1.), sf_id(1.), sf_iso(1.), sf_reco(1.);
        if (name == "W") {
            if (event.getNumGenJets() == 1) {
                evtwt = 3.656;
            } else if (event.getNumGenJets() == 2) {
                evtwt = 3.383;
            } else if (event.getNumGenJets() == 3) {
                evtwt = 2.145;
            } else if (event.getNumGenJets() == 4) {
                evtwt = 1.954;
            } else {
                evtwt = 25.609;
            }
        }

        if (name == "ZTT" || name == "ZLL" || name == "ZL" || name == "ZJ") {
            if (event.getNumGenJets() == 1) {
                evtwt = 0.710;
            } else if (event.getNumGenJets() == 2) {
                evtwt = 0.921;
            } else if (event.getNumGenJets() == 3) {
                evtwt = 1.651;
            } else if (event.getNumGenJets() == 4) {
                evtwt = 0.220;
            } else {
                evtwt = 2.581;
            }
        }
        histos->at("cutflow")->Fill(1., 1.);

        // run factories
        auto electron = electrons.run_factory();
        auto tau = taus.run_factory();
        jets.run_factory();

        // pass event flags
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
        int evt_charge = tau.getCharge() + electron.getCharge();
        if (evt_charge == 0) {
            histos->at("cutflow")->Fill(4., 1.);
        } else {
            continue;
        }

        // build Higgs
        TLorentzVector Higgs = electron.getP4() + tau.getP4() + met.getP4();

        // calculate mt
        double met_x = met.getMet() * cos(met.getMetPhi());
        double met_y = met.getMet() * sin(met.getMetPhi());
        double met_pt = sqrt(pow(met_x, 2) + pow(met_y, 2));
        double mt = sqrt(pow(electron.getPt() + met_pt, 2) - pow(electron.getP4().Px() + met_x, 2) - pow(electron.getP4().Py() + met_y, 2));

        // now do mt selection
        if (mt < 50) {
            histos->at("cutflow")->Fill(5., 1.);
        } else {
            continue;
        }

        // b-jet veto
        if (jets.getNbtagLoose() < 2 && jets.getNbtagMedium() < 1) {
            histos->at("cutflow")->Fill(7., 1.);
        } else {
            continue;
        }

        // create regions
        bool signalRegion = (tau.getMediumIsoDeep() && electron.getIso() < 0.15);
        bool antiTauIsoRegion = (tau.getMediumIsoDeep() == 0 && tau.getVVVLooseIsoDeep() > 0 && electron.getIso() < 0.15);
        if (signal_type != "None") {
            antiTauIsoRegion = false;  // don't need anti-tau iso region in signal
        }

        // only keep the regions we need
        if (signalRegion || antiTauIsoRegion) {
            histos->at("cutflow")->Fill(8., 1.);
        } else {
            continue;
        }

        // apply all scale factors/corrections/etc.
        if (!isData && !isEmbed) {
            // pileup reweighting
            if (!doAC && !isMG) {
                evtwt *= lumi_weights->weight(event.getNPU());
            }

            // generator weights
            evtwt *= event.getGenWeight();

            // prefiring weight (systematics are taken care of already)
            evtwt *= event.getPrefiringWeight();

            // b-tagging scale factor goes here
            evtwt *= jets.getBWeight();

            // Z-Vtx HLT Correction
            evtwt *= 0.991;

            // set workspace variables
            htt_sf->var("e_pt")->setVal(electron.getPt());
            htt_sf->var("e_eta")->setVal(electron.getEta());
            htt_sf->var("t_pt")->setVal(tau.getPt());
            htt_sf->var("t_eta")->setVal(tau.getEta());
            htt_sf->var("t_phi")->setVal(tau.getPhi());
            htt_sf->var("t_dm")->setVal(tau.getDecayMode());
            htt_sf->var("z_gen_mass")->setVal(event.getGenM());
            htt_sf->var("z_gen_pt")->setVal(event.getGenPt());

            // start applying weights from workspace
            evtwt *= htt_sf->function("e_trk_ratio")->getVal();
            evtwt *= htt_sf->function("e_idiso_ic_ratio")->getVal();

            // tau ID efficiency SF and systematics
            std::string id_name = "t_deeptauid_pt_medium";  // nominal
            if (syst.find("tau_id_") != std::string::npos) {
                if ((syst.find("30to35") != std::string::npos && tau.getPt() >= 30 && tau.getPt() < 35) ||
                    (syst.find("35to40") != std::string::npos && tau.getPt() >= 35 && tau.getPt() < 40) ||
                    (syst.find("ptgt40") != std::string::npos && tau.getPt() >= 40)) {
                    id_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
                }
            }
            if (tau.getGenMatch() == 5) {
                evtwt *= htt_sf->function(id_name.c_str())->getVal();
            }

            // electron fake rate SF
            std::string e_fake_id_name = "t_id_vs_e_eta_tight";
            if (syst.find("tau_id_el_disc") != std::string::npos) {
                if ((syst.find("DM0_barrel") != std::string::npos && tau.getDecayMode() == 0 && fabs(tau.getEta()) < 1.479) ||
                    (syst.find("DM0_endcap") != std::string::npos && tau.getDecayMode() == 0 && fabs(tau.getEta()) >= 1.479) ||
                    (syst.find("DM1_barrel") != std::string::npos && tau.getDecayMode() == 1 && fabs(tau.getEta()) < 1.479) ||
                    (syst.find("DM1_endcap") != std::string::npos && tau.getDecayMode() == 1 && fabs(tau.getEta()) >= 1.479)) {
                    e_fake_id_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
                }
            }
            if (tau.getGenMatch() == 1 || tau.getGenMatch() == 3) {
                evtwt *= htt_sf->function(e_fake_id_name.c_str())->getVal();
            }

            // trigger scale factors
            if (electron.getPt() < 33) {
                // electron leg with systematics
                evtwt *= htt_sf->function("e_trg_24_ic_ratio")->getVal();
                if (syst == "mc_cross_trigger_up") {
                    evtwt *= 1.02;  // 2% per light lepton leg
                } else if (syst == "mc_cross_trigger_down") {
                    evtwt *= 0.98;
                }

                // tau leg with systematics
                if (syst == "mc_cross_trigger_up") {
                    evtwt *= htt_sf->function("t_trg_pog_deeptau_medium_etau_ratio_up")->getVal();
                } else if (syst == "mc_cross_trigger_down") {
                    evtwt *= htt_sf->function("t_trg_pog_deeptau_medium_etau_ratio_down")->getVal();
                } else {
                    evtwt *= htt_sf->function("t_trg_pog_deeptau_medium_etau_ratio")->getVal();
                }
            } else {
                evtwt *= htt_sf->function("e_trg_ic_ratio")->getVal();
                if (syst == "mc_single_trigger_up") {
                    evtwt *= 1.02;  // 2% per light lepton leg
                } else if (syst == "mc_single_trigger_down") {
                    evtwt *= 0.98;
                }
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
                float pt_top1 = std::min(static_cast<float>(470.), jets.getTopPt1());
                float pt_top2 = std::min(static_cast<float>(470.), jets.getTopPt2());
                auto top_pt_weight = sqrt(exp(0.088 - 0.00087 * pt_top1 + 0.00000092 * pt_top1 * pt_top1) *
                                          exp(0.088 - 0.00087 * pt_top2 + 0.00000092 * pt_top2 * pt_top2));
                if (syst == "ttbarShape_Up") {
                    top_pt_weight = 2 * top_pt_weight - 1;
                } else if (syst == "ttbarShape_Up") {
                    top_pt_weight = 1.;
                }
                evtwt *= top_pt_weight;
            }

            // ggH theory uncertainty
            if (sample == "ggh125" && signal_type == "powheg") {
                if (event.getNjetsRivet() == 0) evtwt *= g_NNLOPS_0jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(125.0)));
                if (event.getNjetsRivet() == 1) evtwt *= g_NNLOPS_1jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(625.0)));
                if (event.getNjetsRivet() == 2) evtwt *= g_NNLOPS_2jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(800.0)));
                if (event.getNjetsRivet() >= 3) evtwt *= g_NNLOPS_3jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(925.0)));
                NumV WG1unc = qcd_ggF_uncert_2017(event.getNjetsRivet(), event.getHiggsPtRivet(), event.getJetPtRivet());
                if (syst.find("ggH_Rivet") != std::string::npos) {
                    evtwt *= (1 + event.getRivetUnc(WG1unc, syst));
                }
            } else if (signal_type == "madgraph") {
                if (event.getNjetsRivet() == 0) evtwt *= g_mcatnlo_NNLOPS_0jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(125.0)));
                if (event.getNjetsRivet() == 1) evtwt *= g_mcatnlo_NNLOPS_1jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(625.0)));
                if (event.getNjetsRivet() == 2) evtwt *= g_mcatnlo_NNLOPS_2jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(800.0)));
                if (event.getNjetsRivet() >= 3) evtwt *= g_mcatnlo_NNLOPS_3jet->Eval(std::min(event.getHiggsPtRivet(), static_cast<float>(925.0)));
                NumV WG1unc = qcd_ggF_uncert_2017(event.getNjetsRivet(), event.getHiggsPtRivet(), event.getJetPtRivet());
                if (syst.find("ggH_Rivet") != std::string::npos) {
                    evtwt *= (1 + event.getRivetUnc(WG1unc, syst));
                }
            }

            // VBF theory uncertainty
            if (sample == "vbf125" && signal_type == "powheg" && syst.find("VBF_Rivet") != std::string::npos) {
                evtwt *= event.getVBFTheoryUnc(syst);
            }

            auto efake_pt_shift(1.);
            if (syst.find("efaket_norm_ptgt50") != std::string::npos && tau.getPt() > 50) {
                efake_pt_shift = (syst == "efaket_norm_ptgt50_Up" ? 1.1 : 0.9);
            } else if (syst.find("efaket_norm_pt40to50") != std::string::npos && tau.getPt() > 40) {
                efake_pt_shift = (syst == "efaket_norm_pt40to50_Up" ? 1.1 : 0.9);
            } else if (syst.find("efaket_norm_pt30to40") != std::string::npos && tau.getPt() > 30) {
                efake_pt_shift = (syst == "efaket_norm_pt30to40_Up" ? 1.1 : 0.9);
            }
            evtwt *= efake_pt_shift;
        } else if (!isData && isEmbed) {
            event.setEmbed();
            // embedded cross-triggers not applied in skimmer
            if (electron.getPt() < 28 && !event.getPassElEmbedCross()) {
                continue;
            }

            // embedded generator weights
            auto genweight(event.getGenWeight());
            if (genweight > 1 || genweight < 0) {
                genweight = 0;
            }
            evtwt *= genweight;

            // tracking sf
            evtwt *= helper->embed_tracking(tau.getDecayMode(), syst);

            // set workspace variables
            htt_sf->var("e_pt")->setVal(electron.getPt());
            htt_sf->var("e_eta")->setVal(electron.getEta());
            htt_sf->var("t_pt")->setVal(tau.getPt());
            htt_sf->var("t_eta")->setVal(tau.getEta());
            htt_sf->var("t_phi")->setVal(tau.getPhi());
            htt_sf->var("t_dm")->setVal(tau.getDecayMode());
            htt_sf->var("gt1_pt")->setVal(electron.getGenPt());
            htt_sf->var("gt1_eta")->setVal(electron.getGenEta());
            htt_sf->var("gt2_pt")->setVal(tau.getGenPt());
            htt_sf->var("gt2_eta")->setVal(tau.getGenEta());

            evtwt *= htt_sf->function("e_trk_embed_ratio")->getVal();
            evtwt *= htt_sf->function("e_idiso_ic_embed_ratio")->getVal();

            // tau ID eff SF
            std::string id_name = "t_deeptauid_pt_tightvse_embed_medium";
            if (syst.find("tau_id_") != std::string::npos) {
                if ((syst.find("30to35") != std::string::npos && tau.getPt() >= 30 && tau.getPt() < 35) ||
                    (syst.find("35to40") != std::string::npos && tau.getPt() >= 35 && tau.getPt() < 40) ||
                    (syst.find("ptgt40") != std::string::npos && tau.getPt() >= 40)) {
                    id_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
                }
            }
            if (tau.getGenMatch() == 5) {
                evtwt *= htt_sf->function(id_name.c_str())->getVal();
            }

            // electron fake rate SF
            std::string e_fake_id_name = "t_id_vs_e_eta_tight";
            if (syst.find("tau_id_el_disc") != std::string::npos) {
                if ((syst.find("DM0_barrel") != std::string::npos && tau.getDecayMode() == 0 && fabs(tau.getEta()) < 1.479) ||
                    (syst.find("DM0_endcap") != std::string::npos && tau.getDecayMode() == 0 && fabs(tau.getEta()) >= 1.479) ||
                    (syst.find("DM1_barrel") != std::string::npos && tau.getDecayMode() == 1 && fabs(tau.getEta()) < 1.479) ||
                    (syst.find("DM1_endcap") != std::string::npos && tau.getDecayMode() == 1 && fabs(tau.getEta()) >= 1.479)) {
                    e_fake_id_name += syst.find("Up") != std::string::npos ? "_up" : "_down";
                }
            }
            if (tau.getGenMatch() == 1 || tau.getGenMatch() == 3) {
                evtwt *= htt_sf->function(e_fake_id_name.c_str())->getVal();
            }

            // trigger scale factors
            bool fireSingle = electron.getPt() > 28;
            bool fireCross = electron.getPt() < 28;
            std::string single_eff_name = fabs(electron.getEta()) < 1.479 ? "e_trg_ic_embed_ratio" : "e_trg_ic_data";
            std::string el_leg_eff_name = fabs(electron.getEta()) < 1.479 ? "e_trg_24_ic_embed_ratio" : "e_trg_24_ic_data";
            std::string tau_leg_eff_name = fabs(electron.getEta()) < 1.479 ? "t_trg_mediumDeepTau_etau_embed_ratio" : "t_trg_mediumDeepTau_etau_data";
            if (syst == "embed_cross_trigger_up") {
                tau_leg_eff_name += "_up";
            } else if (syst == "embed_cross_trigger_down") {
                tau_leg_eff_name += "_down";
            }

            auto single_eff = htt_sf->function(single_eff_name.c_str())->getVal();
            if (syst == "embed_single_trigger_up") {
                single_eff *= 1.02;  // 2% per light lepton leg
            } else if (syst == "embed_single_trigger_down") {
                single_eff *= 0.98;
            }

            auto el_leg_eff = htt_sf->function(el_leg_eff_name.c_str())->getVal();
            if (syst == "embed_cross_trigger_up") {
                el_leg_eff *= 1.02;  // 2% per light lepton leg
            } else if (syst == "embed_cross_trigger_down") {
                el_leg_eff *= 0.98;
            }

            auto tau_leg_eff = htt_sf->function(tau_leg_eff_name.c_str())->getVal();
            evtwt *= (single_eff * fireSingle + el_leg_eff * tau_leg_eff * fireCross);

            // double muon trigger eff in selection
            evtwt *= htt_sf->function("m_sel_trg_ratio")->getVal();

            // muon ID eff in selection (leg 1)
            htt_sf->var("gt_pt")->setVal(electron.getGenPt());
            htt_sf->var("gt_eta")->setVal(electron.getGenEta());
            evtwt *= htt_sf->function("m_sel_id_ic_ratio")->getVal();

            // muon ID eff in selection (leg 1)
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
        st->fillTree(tree_cat, &electron, &tau, &jets, &met, &event, mt, evtwt, weights, name);
    }  // close event loop

    fin->Close();
    fout->cd();
    fout->Write();
    fout->Close();
    running_log << "Finished processing " << sample << std::endl;
    if (!condor) {
        logfile.close();
    }
    return 0;
}
