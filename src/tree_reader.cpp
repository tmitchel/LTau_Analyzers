#include "../include/tree_reader.h"
#include "TMath.h"

tree_reader::tree_reader(std::map<std::string, std::vector<float>> vars) : acWeightVal(1.), variables(vars) {}

void tree_reader::setBranches(TTree* tree, std::string acName) {
  tree->SetBranchAddress("evtwt", &evtwt);

  tree->SetBranchAddress("el_pt", &el_pt);
  tree->SetBranchAddress("el_eta", &el_eta);
  tree->SetBranchAddress("el_phi", &el_phi);
  tree->SetBranchAddress("el_mass", &el_mass);
  tree->SetBranchAddress("el_charge", &el_charge);
  tree->SetBranchAddress("el_iso", &el_iso);
  tree->SetBranchAddress("mu_pt", &mu_pt);
  tree->SetBranchAddress("mu_eta", &mu_eta);
  tree->SetBranchAddress("mu_phi", &mu_phi);
  tree->SetBranchAddress("mu_mass", &mu_mass);
  tree->SetBranchAddress("mu_charge", &mu_charge);
  tree->SetBranchAddress("mu_iso", &mu_iso);
  tree->SetBranchAddress("t1_pt", &t1_pt);
  tree->SetBranchAddress("t1_eta", &t1_eta);
  tree->SetBranchAddress("t1_phi", &t1_phi);
  tree->SetBranchAddress("t1_mass", &t1_mass);
  tree->SetBranchAddress("t1_charge", &t1_charge);
  tree->SetBranchAddress("t1_iso", &t1_iso);
  tree->SetBranchAddress("t1_iso_VL", &t1_iso_VL);
  tree->SetBranchAddress("t1_iso_L", &t1_iso_L);
  tree->SetBranchAddress("t1_iso_M", &t1_iso_M);
  tree->SetBranchAddress("t1_iso_T", &t1_iso_T);
  tree->SetBranchAddress("t1_iso_VT", &t1_iso_VT);
  tree->SetBranchAddress("t1_iso_VVT", &t1_iso_VVT);
  tree->SetBranchAddress("t1_decayMode", &t1_decayMode);
  tree->SetBranchAddress("t1_dmf", &dmf);
  tree->SetBranchAddress("t1_dmf_new", &dmf_new);
  tree->SetBranchAddress("t1_genMatch", &t1_genMatch);
  tree->SetBranchAddress("t2_pt", &t2_pt);
  tree->SetBranchAddress("t2_eta", &t2_eta);
  tree->SetBranchAddress("t2_phi", &t2_phi);
  tree->SetBranchAddress("t2_mass", &t2_mass);
  tree->SetBranchAddress("t2_charge", &t2_charge);
  tree->SetBranchAddress("t2_iso", &t2_iso);
  tree->SetBranchAddress("t2_iso_VL", &t2_iso_VL);
  tree->SetBranchAddress("t2_iso_L", &t2_iso_L);
  tree->SetBranchAddress("t2_iso_M", &t2_iso_M);
  tree->SetBranchAddress("t2_iso_T", &t2_iso_T);
  tree->SetBranchAddress("t2_iso_VT", &t2_iso_VT);
  tree->SetBranchAddress("t2_iso_VVT", &t2_iso_VVT);
  tree->SetBranchAddress("t2_decayMode", &t2_decayMode);
  tree->SetBranchAddress("t2_dmf", &dmf);
  tree->SetBranchAddress("t2_dmf_new", &dmf_new);
  tree->SetBranchAddress("t2_genMatch", &t2_genMatch);

  tree->SetBranchAddress("njets", &njets);
  tree->SetBranchAddress("nbjets", &nbjets);
  tree->SetBranchAddress("j1_pt", &j1_pt);
  tree->SetBranchAddress("j1_eta", &j1_eta);
  tree->SetBranchAddress("j1_phi", &j1_phi);
  tree->SetBranchAddress("j2_pt", &j2_pt);
  tree->SetBranchAddress("j2_eta", &j2_eta);
  tree->SetBranchAddress("j2_phi", &j2_phi);
  tree->SetBranchAddress("b1_pt", &b1_pt);
  tree->SetBranchAddress("b1_eta", &b1_eta);
  tree->SetBranchAddress("b1_phi", &b1_phi);
  tree->SetBranchAddress("b2_pt", &b2_pt);
  tree->SetBranchAddress("b2_eta", &b2_eta);
  tree->SetBranchAddress("b2_phi", &b2_phi);

  tree->SetBranchAddress("met", &met);
  tree->SetBranchAddress("metphi", &metphi);
  tree->SetBranchAddress("mjj", &mjj);
  tree->SetBranchAddress("mt", &mt);

  tree->SetBranchAddress("numGenJets", &numGenJets);

  tree->SetBranchAddress("pt_sv", &pt_sv);
  tree->SetBranchAddress("m_sv", &m_sv);
  tree->SetBranchAddress("Dbkg_VBF", &Dbkg_VBF);
  tree->SetBranchAddress("Dbkg_ggH", &Dbkg_ggH);
  tree->SetBranchAddress("D0_VBF", &D0_VBF);
  tree->SetBranchAddress("DCP_VBF", &DCP_VBF);
  tree->SetBranchAddress("D0_ggH", &D0_ggH);
  tree->SetBranchAddress("DCP_ggH", &DCP_ggH);

  tree->SetBranchAddress("Phi", &Phi);
  tree->SetBranchAddress("Phi1", &Phi1);
  tree->SetBranchAddress("costheta1", &costheta1);
  tree->SetBranchAddress("costheta2", &costheta2);
  tree->SetBranchAddress("costhetastar", &costhetastar);
  tree->SetBranchAddress("Q2V1", &Q2V1);
  tree->SetBranchAddress("Q2V2", &Q2V2);
  tree->SetBranchAddress("ME_sm_VBF", &ME_sm_VBF);
  tree->SetBranchAddress("ME_sm_ggH", &ME_sm_ggH);
  tree->SetBranchAddress("ME_sm_WH", &ME_sm_WH);
  tree->SetBranchAddress("ME_sm_ZH", &ME_sm_ZH);
  tree->SetBranchAddress("ME_bkg", &ME_bkg);
  tree->SetBranchAddress("ME_bkg1", &ME_bkg1);
  tree->SetBranchAddress("ME_bkg2", &ME_bkg2);
  tree->SetBranchAddress("VBF_MELA", &VBF_MELA);

  tree->SetBranchAddress("higgs_pT", &higgs_pT);
  tree->SetBranchAddress("higgs_m", &higgs_m);
  tree->SetBranchAddress("hjj_pT", &hjj_pT);
  tree->SetBranchAddress("hjj_m", &hjj_m);
  tree->SetBranchAddress("vis_mass", &vis_mass);
  tree->SetBranchAddress("dEtajj", &dEtajj);
  tree->SetBranchAddress("dPhijj", &dPhijj);
  tree->SetBranchAddress("NN_disc", &NN_disc);
  tree->SetBranchAddress("NN_disc_vbf", &NN_disc_vbf);
  tree->SetBranchAddress("NN_disc_boost", &NN_disc_boost);

  tree->SetBranchAddress("lt_dphi", &lt_dphi);
  tree->SetBranchAddress("MT_lepMET", &MT_lepMET);
  tree->SetBranchAddress("MT_HiggsMET", &MT_HiggsMET);
  tree->SetBranchAddress("hj_dphi", &hj_dphi);
  tree->SetBranchAddress("jmet_dphi", &jmet_dphi);
  tree->SetBranchAddress("MT_t2MET", &MT_t2MET);
  tree->SetBranchAddress("hj_deta", &hj_deta);
  tree->SetBranchAddress("hmet_dphi", &hmet_dphi);
  tree->SetBranchAddress("hj_dr", &hj_dr);

  tree->SetBranchAddress("is_signal", &is_signal);
  tree->SetBranchAddress("is_antiLepIso", &is_antiLepIso);
  tree->SetBranchAddress("is_antiTauIso", &is_antiTauIso);
  tree->SetBranchAddress("is_qcd", &is_qcd);
  tree->SetBranchAddress("is_looseIso", &is_looseIso);
  tree->SetBranchAddress("cat_0jet", &cat_0jet);
  tree->SetBranchAddress("cat_boosted", &cat_boosted);
  tree->SetBranchAddress("cat_vbf", &cat_vbf);
  tree->SetBranchAddress("cat_VH", &cat_VH);
  tree->SetBranchAddress("OS", &OS);
  tree->SetBranchAddress("SS", &SS);

  // include weights for anomolous coupling
  if (acName != "None") {
    tree->SetBranchAddress(acName.c_str(), &acWeightVal);

    tree->SetBranchAddress("wt_a1", &wt_a1);
    tree->SetBranchAddress("wt_a2", &wt_a2);
    tree->SetBranchAddress("wt_a3", &wt_a3);
    tree->SetBranchAddress("wt_L1", &wt_L1);
    tree->SetBranchAddress("wt_L1Zg", &wt_L1Zg);
    tree->SetBranchAddress("wt_a2int", &wt_a2int);
    tree->SetBranchAddress("wt_a3int", &wt_a3int);
    tree->SetBranchAddress("wt_L1int", &wt_L1int);
    tree->SetBranchAddress("wt_L1Zgint", &wt_L1Zgint);

    tree->SetBranchAddress("wt_ggH_a1", &wt_ggH_a1);
    tree->SetBranchAddress("wt_ggH_a3", &wt_ggH_a3);
    tree->SetBranchAddress("wt_ggH_a3int", &wt_ggH_a3int);

    tree->SetBranchAddress("wt_wh_a1", &wt_wh_a1);
    tree->SetBranchAddress("wt_wh_a2", &wt_wh_a2);
    tree->SetBranchAddress("wt_wh_a3", &wt_wh_a3);
    tree->SetBranchAddress("wt_wh_L1", &wt_wh_L1);
    tree->SetBranchAddress("wt_wh_L1Zg", &wt_wh_L1Zg);
    tree->SetBranchAddress("wt_wh_a2int", &wt_wh_a2int);
    tree->SetBranchAddress("wt_wh_a3int", &wt_wh_a3int);
    tree->SetBranchAddress("wt_wh_L1int", &wt_wh_L1int);
    tree->SetBranchAddress("wt_wh_L1Zgint", &wt_wh_L1Zgint);

    tree->SetBranchAddress("wt_zh_a1", &wt_zh_a1);
    tree->SetBranchAddress("wt_zh_a2", &wt_zh_a2);
    tree->SetBranchAddress("wt_zh_a3", &wt_zh_a3);
    tree->SetBranchAddress("wt_zh_L1", &wt_zh_L1);
    tree->SetBranchAddress("wt_zh_L1Zg", &wt_zh_L1Zg);
    tree->SetBranchAddress("wt_zh_a2int", &wt_zh_a2int);
    tree->SetBranchAddress("wt_zh_a3int", &wt_zh_a3int);
    tree->SetBranchAddress("wt_zh_L1int", &wt_zh_L1int);
    tree->SetBranchAddress("wt_zh_L1Zgint", &wt_zh_L1Zgint);
  }
}

// map names to variables in the tree
// new variables can be created here
Float_t tree_reader::getVar(std::string var) {
  if (var == "NN_disc_vbf") {
    return this->NN_disc_vbf;
  } else if (var == "NN_disc_boost") {
    return this->NN_disc_boost;
  } else if (var == "NN_disc") {
    return this->NN_disc;
  } else if (var == "el_pt") {
    return this->el_pt;
  } else if (var == "mu_pt") {
    return this->mu_pt;
  } else if (var == "t1_pt") {
    return this->t1_pt;
  } else if (var == "j1_pt") {
    return this->j1_pt;
  } else if (var == "j2_pt") {
    return this->j2_pt;
  } else if (var == "VBF_MELA") {
    return this->VBF_MELA;
  } else if (var == "D0_ggH") {
    return this->D0_ggH;
  } else if (var == "dPhijj") {
    return TMath::ACos(TMath::Cos(this->j1_phi - this->j2_phi));
  } else if (var == "m_sv") {
    return this->m_sv;
  } else if (var == "lt_dphi") {
    return this->lt_dphi;
  } else if (var == "MT_lepMET") {
    return this->MT_lepMET;
  } else if (var == "MT_HiggsMET") {
    return this->MT_HiggsMET;
  } else if (var == "hj_dphi") {
    return this->hj_dphi;
  } else if (var == "jmet_dphi") {
    return this->jmet_dphi;
  } else if (var == "MT_t2MET") {
    return this->MT_t2MET;
  } else if (var == "hj_deta") {
    return this->hj_deta;
  } else if (var == "hmet_dphi") {
    return this->hmet_dphi;
  } else if (var == "hj_dr") {
    return this->hj_dr;
  } else if (var == "met") {
    return this->met;
  }
}

