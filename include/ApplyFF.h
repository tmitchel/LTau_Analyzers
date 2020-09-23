// Copyright [2020] Tyler Mitchell

#ifndef INCLUDE_APPLYFF_H_
#define INCLUDE_APPLYFF_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"

class apply_ff {
   private:
    std::string channel, year;

    // input files
    TFile *raw_file, *vis_mass_file, *osss_file, *tpt_file;

    // TF1 from files
    TF1 *ff_qcd_0jet, *ff_qcd_0jet_unc1_up, *ff_qcd_0jet_unc1_down, *ff_qcd_0jet_unc2_up, *ff_qcd_0jet_unc2_down, *ff_qcd_1jet, *ff_qcd_1jet_unc1_up,
        *ff_qcd_1jet_unc1_down, *ff_qcd_1jet_unc2_up, *ff_qcd_1jet_unc2_down, *ff_qcd_2jet, *ff_qcd_2jet_unc1_up, *ff_qcd_2jet_unc1_down,
        *ff_qcd_2jet_unc2_up, *ff_qcd_2jet_unc2_down, *ff_w_0jet, *ff_w_0jet_unc1_up, *ff_w_0jet_unc1_down, *ff_w_0jet_unc2_up, *ff_w_0jet_unc2_down,
        *ff_w_1jet, *ff_w_1jet_unc1_up, *ff_w_1jet_unc1_down, *ff_w_1jet_unc2_up, *ff_w_1jet_unc2_down, *ff_w_2jet, *ff_w_2jet_unc1_up,
        *ff_w_2jet_unc1_down, *ff_w_2jet_unc2_up, *ff_w_2jet_unc2_down, *ff_tt_0jet, *ff_tt_0jet_unc1_up, *ff_tt_0jet_unc1_down, *ff_tt_0jet_unc2_up,
        *ff_tt_0jet_unc2_down;

    TF1 *mVisClosure_QCD_0jet, *mVisClosure_QCD_1jet, *mVisClosure_QCD_2jet, *mVisClosure_W_0jet, *mVisClosure_W_1jet, *mVisClosure_W_2jet,
        *mVisClosure_TT;

    TF1 *lptClosure_W_taupt30to50, *lptClosure_W_taupt50to70, *lptClosure_W_tauptgt70, *lptClosure_QCD_taupt30to50, *lptClosure_QCD_taupt50to70,
        *lptClosure_QCD_tauptgt70, *lptClosure_TT_taupt30to50, *lptClosure_TT_taupt50to70, *lptClosure_TT_tauptgt70, *lptClosure_W_xtrg_taupt30to50,
        *lptClosure_W_xtrg_taupt50to70, *lptClosure_W_xtrg_tauptgt70, *lptClosure_QCD_xtrg_taupt30to50, *lptClosure_QCD_xtrg_taupt50to70,
        *lptClosure_QCD_xtrg_tauptgt70, *lptClosure_TT_xtrg_taupt30to50, *lptClosure_TT_xtrg_taupt50to70, *lptClosure_TT_xtrg_tauptgt70;

    TF1 *lptClosure_W_taupt30to40, *lptClosure_W_taupt40to50, *lptClosure_W_tauptgt50, *lptClosure_QCD_taupt30to40, *lptClosure_QCD_taupt40to50,
        *lptClosure_QCD_tauptgt50, *lptClosure_TT_taupt30to40, *lptClosure_TT_taupt40to50, *lptClosure_TT_tauptgt50, *lptClosure_W_xtrg_taupt30to40,
        *lptClosure_W_xtrg_taupt40to50, *lptClosure_W_xtrg_tauptgt50, *lptClosure_QCD_xtrg_taupt30to40, *lptClosure_QCD_xtrg_taupt40to50,
        *lptClosure_QCD_xtrg_tauptgt50, *lptClosure_TT_xtrg_taupt30to40, *lptClosure_TT_xtrg_taupt40to50, *lptClosure_TT_xtrg_tauptgt50;

    TF1 *OSSSClosure_QCD, *MTClosure_W, *MTClosure_W_unc1_up, *MTClosure_W_unc1_down, *MTClosure_W_unc2_up, *MTClosure_W_unc2_down,
        *tauPtCorrection_qcd, *tauPtCorrection_w;

    Float_t raw_qcd(Float_t, Float_t, std::string, std::string);
    Float_t raw_w(Float_t, Float_t, std::string, std::string);
    Float_t raw_tt(Float_t, std::string, std::string);

    Float_t lpt_cls_corr_w(Float_t, Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_tt(Float_t, Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_qcd(Float_t, Float_t, Float_t, std::string, std::string);

    Float_t lpt_cls_corr_xtrg_mt_w(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_xtrg_mt_qcd(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_xtrg_mt_tt(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_mt_w(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_mt_qcd(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_mt_tt(Float_t, Float_t, std::string, std::string);

    Float_t lpt_cls_corr_xtrg_et_w(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_xtrg_et_qcd(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_xtrg_et_tt(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_et_w(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_et_qcd(Float_t, Float_t, std::string, std::string);
    Float_t lpt_cls_corr_et_tt(Float_t, Float_t, std::string, std::string);

    Float_t mt_closure_corr(Float_t, std::string, std::string);
    Float_t osss_closure_corr(Float_t, std::string, std::string);

   public:
    apply_ff(std::string, std::string);
    ~apply_ff() {}

    Float_t get_ff(std::vector<Float_t>, std::string, std::string);
};

apply_ff::apply_ff(std::string path, std::string _channel) : channel(_channel) {
    raw_file = TFile::Open((path + "uncorrected_fakefactors_" + channel + ".root").c_str());
    vis_mass_file = TFile::Open((path + "FF_corrections_1.root").c_str());
    osss_file = TFile::Open((path + "FF_QCDcorrectionOSSS.root").c_str());
    tpt_file = TFile::Open((path + "tauptcorrection_" + channel + ".root").c_str());

    // make sure there were no problems opening the file
    if (raw_file->IsZombie() || vis_mass_file->IsZombie() || osss_file->IsZombie() || tpt_file->IsZombie()) {
        std::cerr << "cannot open input FF files" << std::endl;
    }

    // read all TF1's
    ff_qcd_0jet = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_0jet").c_str()));
    ff_qcd_0jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_0jet_unc1_up").c_str()));
    ff_qcd_0jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_0jet_unc1_down").c_str()));
    ff_qcd_0jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_0jet_unc2_up").c_str()));
    ff_qcd_0jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_0jet_unc2_down").c_str()));

    ff_qcd_1jet = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_1jet").c_str()));
    ff_qcd_1jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_1jet_unc1_up").c_str()));
    ff_qcd_1jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_1jet_unc1_down").c_str()));
    ff_qcd_1jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_1jet_unc2_up").c_str()));
    ff_qcd_1jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_1jet_unc2_down").c_str()));

    ff_qcd_2jet = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_2jet").c_str()));
    ff_qcd_2jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_2jet_unc1_up").c_str()));
    ff_qcd_2jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_2jet_unc1_down").c_str()));
    ff_qcd_2jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_2jet_unc2_up").c_str()));
    ff_qcd_2jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_qcd_2jet_unc2_down").c_str()));

    ff_w_0jet = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_0jet").c_str()));
    ff_w_0jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_0jet_unc1_up").c_str()));
    ff_w_0jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_0jet_unc1_down").c_str()));
    ff_w_0jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_0jet_unc2_up").c_str()));
    ff_w_0jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_0jet_unc2_down").c_str()));

    ff_w_1jet = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_1jet").c_str()));
    ff_w_1jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_1jet_unc1_up").c_str()));
    ff_w_1jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_1jet_unc1_down").c_str()));
    ff_w_1jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_1jet_unc2_up").c_str()));
    ff_w_1jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_1jet_unc2_down").c_str()));

    ff_w_2jet = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_2jet").c_str()));
    ff_w_2jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_2jet_unc1_up").c_str()));
    ff_w_2jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_2jet_unc1_down").c_str()));
    ff_w_2jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_2jet_unc2_up").c_str()));
    ff_w_2jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("rawFF_" + channel + "_w_2jet_unc2_down").c_str()));

    ff_tt_0jet = reinterpret_cast<TF1 *>(raw_file->Get(("mc_rawFF_" + channel + "_tt").c_str()));
    ff_tt_0jet_unc1_up = reinterpret_cast<TF1 *>(raw_file->Get(("mc_rawFF_" + channel + "_tt_unc1_up").c_str()));
    ff_tt_0jet_unc1_down = reinterpret_cast<TF1 *>(raw_file->Get(("mc_rawFF_" + channel + "_tt_unc1_down").c_str()));
    ff_tt_0jet_unc2_up = reinterpret_cast<TF1 *>(raw_file->Get(("mc_rawFF_" + channel + "_tt_unc2_up").c_str()));
    ff_tt_0jet_unc2_down = reinterpret_cast<TF1 *>(raw_file->Get(("mc_rawFF_" + channel + "_tt_unc2_down").c_str()));

    mVisClosure_QCD_0jet = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_0jet_qcd").c_str()));
    mVisClosure_QCD_1jet = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_1jet_qcd").c_str()));
    mVisClosure_QCD_2jet = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_2jet_qcd").c_str()));
    mVisClosure_W_0jet = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_0jet_w").c_str()));
    mVisClosure_W_1jet = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_1jet_w").c_str()));
    mVisClosure_W_2jet = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_2jet_w").c_str()));
    mVisClosure_TT = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_mvis_" + channel + "_ttmc").c_str()));

    if (channel == "mt") {
        lptClosure_W_taupt30to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to50_" + channel + "_w").c_str()));
        lptClosure_W_taupt50to70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt50to70_" + channel + "_w").c_str()));
        lptClosure_W_tauptgt70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt70_" + channel + "_w").c_str()));
        lptClosure_QCD_taupt30to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to50_" + channel + "_qcd").c_str()));
        lptClosure_QCD_taupt50to70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt50to70_" + channel + "_qcd").c_str()));
        lptClosure_QCD_tauptgt70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt70_" + channel + "_qcd").c_str()));
        lptClosure_TT_taupt30to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to50_" + channel + "_ttmc").c_str()));
        lptClosure_TT_taupt50to70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt50to70_" + channel + "_ttmc").c_str()));
        lptClosure_TT_tauptgt70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt70_" + channel + "_ttmc").c_str()));

        lptClosure_W_xtrg_taupt30to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to50_xtrg_" + channel + "_w").c_str()));
        lptClosure_W_xtrg_taupt50to70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt50to70_xtrg_" + channel + "_w").c_str()));
        lptClosure_W_xtrg_tauptgt70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt70_xtrg_" + channel + "_w").c_str()));
        lptClosure_QCD_xtrg_taupt30to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to50_xtrg_" + channel + "_qcd").c_str()));
        lptClosure_QCD_xtrg_taupt50to70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt50to70_xtrg_" + channel + "_qcd").c_str()));
        lptClosure_QCD_xtrg_tauptgt70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt70_xtrg_" + channel + "_qcd").c_str()));
        lptClosure_TT_xtrg_taupt30to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to50_xtrg_" + channel + "_ttmc").c_str()));
        lptClosure_TT_xtrg_taupt50to70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt50to70_xtrg_" + channel + "_ttmc").c_str()));
        lptClosure_TT_xtrg_tauptgt70 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt70_xtrg_" + channel + "_ttmc").c_str()));
    } else {
        lptClosure_W_taupt30to40 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to40_" + channel + "_w").c_str()));
        lptClosure_W_taupt40to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt40to50_" + channel + "_w").c_str()));
        lptClosure_W_tauptgt50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt50_" + channel + "_w").c_str()));
        lptClosure_QCD_taupt30to40 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to40_" + channel + "_qcd").c_str()));
        lptClosure_QCD_taupt40to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt40to50_" + channel + "_qcd").c_str()));
        lptClosure_QCD_tauptgt50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt50_" + channel + "_qcd").c_str()));
        lptClosure_TT_taupt30to40 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to40_" + channel + "_ttmc").c_str()));
        lptClosure_TT_taupt40to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt40to50_" + channel + "_ttmc").c_str()));
        lptClosure_TT_tauptgt50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt50_" + channel + "_ttmc").c_str()));

        lptClosure_W_xtrg_taupt30to40 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to40_xtrg_" + channel + "_w").c_str()));
        lptClosure_W_xtrg_taupt40to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt40to50_xtrg_" + channel + "_w").c_str()));
        lptClosure_W_xtrg_tauptgt50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt50_xtrg_" + channel + "_w").c_str()));
        lptClosure_QCD_xtrg_taupt30to40 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to40_xtrg_" + channel + "_qcd").c_str()));
        lptClosure_QCD_xtrg_taupt40to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt40to50_xtrg_" + channel + "_qcd").c_str()));
        lptClosure_QCD_xtrg_tauptgt50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt50_xtrg_" + channel + "_qcd").c_str()));
        lptClosure_TT_xtrg_taupt30to40 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt30to40_xtrg_" + channel + "_ttmc").c_str()));
        lptClosure_TT_xtrg_taupt40to50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_taupt40to50_xtrg_" + channel + "_ttmc").c_str()));
        lptClosure_TT_xtrg_tauptgt50 = reinterpret_cast<TF1 *>(vis_mass_file->Get(("closure_lpt_tauptgt50_xtrg_" + channel + "_ttmc").c_str()));
    }

    OSSSClosure_QCD = reinterpret_cast<TF1 *>(osss_file->Get(("closure_OSSS_dr_flat_" + channel + "_qcd").c_str()));

    MTClosure_W = reinterpret_cast<TF1 *>(osss_file->Get(("closure_mt_" + channel + "_w").c_str()));
    MTClosure_W_unc1_up = reinterpret_cast<TF1 *>(osss_file->Get(("closure_mt_" + channel + "_w_unc1_up").c_str()));
    MTClosure_W_unc1_down = reinterpret_cast<TF1 *>(osss_file->Get(("closure_mt_" + channel + "_w_unc1_down").c_str()));
    MTClosure_W_unc2_up = reinterpret_cast<TF1 *>(osss_file->Get(("closure_mt_" + channel + "_w_unc2_up").c_str()));
    MTClosure_W_unc2_down = reinterpret_cast<TF1 *>(osss_file->Get(("closure_mt_" + channel + "_w_unc2_down").c_str()));

    tauPtCorrection_qcd = reinterpret_cast<TF1 *>(tpt_file->Get("mt_0jet_qcd_taupt_iso"));
    tauPtCorrection_w = reinterpret_cast<TF1 *>(tpt_file->Get("mt_0jet_w_taupt_iso"));
}

Float_t apply_ff::raw_w(Float_t pt, Float_t njets, std::string unc, std::string dir) {
    if (njets == 0) {
        if (unc == "ff_w_0jet_unc1") {
            return dir == "up" ? ff_w_0jet_unc1_up->Eval(pt) : ff_w_0jet_unc1_down->Eval(pt);
        } else if (unc == "ff_w_0jet_unc2") {
            return dir == "up" ? ff_w_0jet_unc2_up->Eval(pt) : ff_w_0jet_unc2_down->Eval(pt);
        }
        return ff_w_0jet->Eval(pt);
    } else if (njets == 1) {
        if (unc == "ff_w_1jet_unc1") {
            return dir == "up" ? ff_w_1jet_unc1_up->Eval(pt) : ff_w_1jet_unc1_down->Eval(pt);
        } else if (unc == "ff_w_1jet_unc2") {
            return dir == "up" ? ff_w_1jet_unc2_up->Eval(pt) : ff_w_1jet_unc2_down->Eval(pt);
        }
        return ff_w_1jet->Eval(pt);
    }

    // must be 2+ jets
    if (unc == "ff_w_2jet_unc1") {
        return dir == "up" ? ff_w_2jet_unc1_up->Eval(pt) : ff_w_2jet_unc1_down->Eval(pt);
    } else if (unc == "ff_w_2jet_unc2") {
        return dir == "up" ? ff_w_2jet_unc2_up->Eval(pt) : ff_w_2jet_unc2_down->Eval(pt);
    }
    return ff_w_2jet->Eval(pt);
}

Float_t apply_ff::raw_qcd(Float_t pt, Float_t njets, std::string unc, std::string dir) {
    if (njets == 0) {
        if (unc == "ff_qcd_0jet_unc1") {
            return dir == "up" ? ff_qcd_0jet_unc1_up->Eval(pt) : ff_qcd_0jet_unc1_down->Eval(pt);
        } else if (unc == "ff_qcd_0jet_unc2") {
            return dir == "up" ? ff_qcd_0jet_unc2_up->Eval(pt) : ff_qcd_0jet_unc2_down->Eval(pt);
        }
        return ff_qcd_0jet->Eval(pt);
    } else if (njets == 1) {
        if (unc == "ff_qcd_1jet_unc1") {
            return dir == "up" ? ff_qcd_1jet_unc1_up->Eval(pt) : ff_qcd_1jet_unc1_down->Eval(pt);
        } else if (unc == "ff_qcd_1jet_unc2") {
            return dir == "up" ? ff_qcd_1jet_unc2_up->Eval(pt) : ff_qcd_1jet_unc2_down->Eval(pt);
        }

        return ff_qcd_1jet->Eval(pt);
    }

    // must be 2+ jets
    if (unc == "ff_qcd_2jet_unc1") {
        return dir == "up" ? ff_qcd_2jet_unc1_up->Eval(pt) : ff_qcd_2jet_unc1_down->Eval(pt);
    } else if (unc == "ff_qcd_2jet_unc2") {
        return dir == "up" ? ff_qcd_2jet_unc2_up->Eval(pt) : ff_qcd_2jet_unc2_down->Eval(pt);
    }
    return ff_qcd_2jet->Eval(pt);
}

Float_t apply_ff::raw_tt(Float_t pt, std::string unc, std::string dir) {
    if (unc == "ff_tt_0jet_unc1") {
        return dir == "up" ? ff_tt_0jet_unc1_up->Eval(pt) : ff_tt_0jet_unc1_down->Eval(pt);
    } else if (unc == "ff_tt_0jet_unc2") {
        return dir == "up" ? ff_tt_0jet_unc2_up->Eval(pt) : ff_tt_0jet_unc2_down->Eval(pt);
    }
    return ff_tt_0jet->Eval(pt);
}

Float_t apply_ff::lpt_cls_corr_w(Float_t pt, Float_t lpt, Float_t xtrg, std::string unc, std::string dir) {
    // call sub-functions so that if mt/et synchronize at some point we can continue to just call
    // this function to handle things
    if (channel == "mt") {
        return xtrg ? lpt_cls_corr_xtrg_mt_w(pt, lpt, unc, dir) : lpt_cls_corr_mt_w(pt, lpt, unc, dir);
    } else {
        return xtrg ? lpt_cls_corr_xtrg_et_w(pt, lpt, unc, dir) : lpt_cls_corr_et_w(pt, lpt, unc, dir);
    }
}

Float_t apply_ff::lpt_cls_corr_tt(Float_t pt, Float_t lpt, Float_t xtrg, std::string unc, std::string dir) {
    // call sub-functions so that if mt/et synchronize at some point we can continue to just call
    // this function to handle things
    if (channel == "mt") {
        return xtrg ? lpt_cls_corr_xtrg_mt_tt(pt, lpt, unc, dir) : lpt_cls_corr_mt_tt(pt, lpt, unc, dir);
    } else {
        return xtrg ? lpt_cls_corr_xtrg_et_tt(pt, lpt, unc, dir) : lpt_cls_corr_et_tt(pt, lpt, unc, dir);
    }
}

Float_t apply_ff::lpt_cls_corr_qcd(Float_t pt, Float_t lpt, Float_t xtrg, std::string unc, std::string dir) {
    // call sub-functions so that if mt/et synchronize at some point we can continue to just call
    // this function to handle things
    if (channel == "mt") {
        return xtrg ? lpt_cls_corr_xtrg_mt_qcd(pt, lpt, unc, dir) : lpt_cls_corr_mt_qcd(pt, lpt, unc, dir);
    } else {
        return xtrg ? lpt_cls_corr_xtrg_et_qcd(pt, lpt, unc, dir) : lpt_cls_corr_et_qcd(pt, lpt, unc, dir);
    }
}

Float_t apply_ff::lpt_cls_corr_mt_qcd(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_qcd") {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_QCD_taupt30to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_QCD_taupt50to70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 70.) {
            return lptClosure_QCD_tauptgt70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_QCD_taupt30to50->Eval(lpt);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_QCD_taupt50to70->Eval(lpt);
        } else if (pt > 70.) {
            return lptClosure_QCD_tauptgt70->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_mt_w(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_w") {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_W_taupt30to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_W_taupt50to70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 70.) {
            return lptClosure_W_tauptgt70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_W_taupt30to50->Eval(lpt);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_W_taupt50to70->Eval(lpt);
        } else if (pt > 70.) {
            return lptClosure_W_tauptgt70->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_mt_tt(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_tt") {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_TT_taupt30to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_TT_taupt50to70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 70.) {
            return lptClosure_TT_tauptgt70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_TT_taupt30to50->Eval(lpt);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_TT_taupt50to70->Eval(lpt);
        } else if (pt > 70.) {
            return lptClosure_TT_tauptgt70->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_xtrg_mt_qcd(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_qcd") {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_QCD_xtrg_taupt30to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_QCD_xtrg_taupt50to70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 70.) {
            return lptClosure_QCD_xtrg_tauptgt70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_QCD_xtrg_taupt30to50->Eval(lpt);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_QCD_xtrg_taupt50to70->Eval(lpt);
        } else if (pt > 70.) {
            return lptClosure_QCD_xtrg_tauptgt70->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_xtrg_mt_w(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_w") {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_W_xtrg_taupt30to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_W_xtrg_taupt50to70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 70.) {
            return lptClosure_W_xtrg_tauptgt70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_W_xtrg_taupt30to50->Eval(lpt);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_W_xtrg_taupt50to70->Eval(lpt);
        } else if (pt > 70.) {
            return lptClosure_W_xtrg_tauptgt70->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_xtrg_mt_tt(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_tt") {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_TT_xtrg_taupt30to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_TT_xtrg_taupt50to70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 70.) {
            return lptClosure_TT_xtrg_tauptgt70->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 50.) {
            return lptClosure_TT_xtrg_taupt30to50->Eval(lpt);
        } else if (pt > 50. && pt <= 70.) {
            return lptClosure_TT_xtrg_taupt50to70->Eval(lpt);
        } else if (pt > 70.) {
            return lptClosure_TT_xtrg_tauptgt70->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_et_qcd(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_qcd") {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_QCD_taupt30to40->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_QCD_taupt40to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50.) {
            return lptClosure_QCD_tauptgt50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_QCD_taupt30to40->Eval(lpt);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_QCD_taupt40to50->Eval(lpt);
        } else if (pt > 50.) {
            return lptClosure_QCD_tauptgt50->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_et_w(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_w") {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_W_taupt30to40->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_W_taupt40to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50.) {
            return lptClosure_W_tauptgt50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_W_taupt30to40->Eval(lpt);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_W_taupt40to50->Eval(lpt);
        } else if (pt > 50.) {
            return lptClosure_W_tauptgt50->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_et_tt(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_tt") {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_TT_taupt30to40->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_TT_taupt40to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50.) {
            return lptClosure_TT_tauptgt50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_TT_taupt30to40->Eval(lpt);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_TT_taupt40to50->Eval(lpt);
        } else if (pt > 50.) {
            return lptClosure_TT_tauptgt50->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_xtrg_et_qcd(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_qcd") {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_QCD_xtrg_taupt30to40->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_QCD_xtrg_taupt40to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50.) {
            return lptClosure_QCD_xtrg_tauptgt50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_QCD_xtrg_taupt30to40->Eval(lpt);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_QCD_xtrg_taupt40to50->Eval(lpt);
        } else if (pt > 50.) {
            return lptClosure_QCD_xtrg_tauptgt50->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_xtrg_et_w(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_w") {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_W_xtrg_taupt30to40->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_W_xtrg_taupt40to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50.) {
            return lptClosure_W_xtrg_tauptgt50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_W_xtrg_taupt30to40->Eval(lpt);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_W_xtrg_taupt40to50->Eval(lpt);
        } else if (pt > 50.) {
            return lptClosure_W_xtrg_tauptgt50->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::lpt_cls_corr_xtrg_et_tt(Float_t pt, Float_t lpt, std::string unc, std::string dir) {
    if (unc == "lptclosure_tt") {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_TT_xtrg_taupt30to40->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_TT_xtrg_taupt40to50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        } else if (pt > 50.) {
            return lptClosure_TT_xtrg_tauptgt50->Eval(lpt) * (dir == "up" ? 1.1 : 0.9);
        }
    } else {
        if (pt > 30. && pt <= 40.) {
            return lptClosure_TT_xtrg_taupt30to40->Eval(lpt);
        } else if (pt > 40. && pt <= 50.) {
            return lptClosure_TT_xtrg_taupt40to50->Eval(lpt);
        } else if (pt > 50.) {
            return lptClosure_TT_xtrg_tauptgt50->Eval(lpt);
        }
    }

    return 1.;
}

Float_t apply_ff::mt_closure_corr(Float_t mt, std::string unc, std::string dir) {
    if (unc == "mtclosure_w_unc1") {
        return dir == "up" ? MTClosure_W_unc1_up->Eval(mt) : MTClosure_W_unc1_down->Eval(mt);
    } else if (unc == "mtclosure_w_unc2") {
        return dir == "up" ? MTClosure_W_unc2_up->Eval(mt) : MTClosure_W_unc2_down->Eval(mt);
    }

    return MTClosure_W->Eval(mt);
}

Float_t apply_ff::osss_closure_corr(Float_t dr, std::string unc, std::string dir) {
    if (unc == "osssclosure_qcd") {
        return OSSSClosure_QCD->Eval(dr) * (dir == "up" ? 1.1 : 0.9);
    }

    return OSSSClosure_QCD->Eval(dr);
}

Float_t apply_ff::get_ff(std::vector<Float_t> kin, std::string unc = "", std::string dir = "") {
    // kin = pt(0), mt(1), mvis(2), lpt(3), dr(4), met(5), njets(6), xtrg(7), frac_tt(8), frac_qcd(9), frac_w(10)
    Float_t eff_pt(std::min(static_cast<Float_t>(100), kin.at(0)));
    Float_t eff_mvis(std::min(static_cast<Float_t>(250), kin.at(2)));
    Float_t eff_lpt(std::min(static_cast<Float_t>(150), kin.at(3)));
    Float_t eff_met(std::max(static_cast<Float_t>(0), kin.at(5)));

    // get raw fake factors
    Float_t ff_qcd(1.), ff_w(1.), ff_tt(1.);
    ff_qcd = raw_qcd(eff_pt, kin.at(6), unc, dir);
    ff_w = raw_w(eff_pt, kin.at(6), unc, dir);
    ff_tt = raw_tt(eff_pt, unc, dir);

    // lepton pT closure (depends on channel and cross trigger)
    ff_w *= lpt_cls_corr_w(eff_pt, eff_lpt, kin.at(7), unc, dir);
    ff_tt *= lpt_cls_corr_tt(eff_pt, eff_lpt, kin.at(7), unc, dir);
    ff_qcd *= lpt_cls_corr_qcd(eff_pt, eff_lpt, kin.at(7), unc, dir);

    // other closure corrections
    ff_w *= mt_closure_corr(kin.at(1), unc, dir);
    ff_qcd *= osss_closure_corr(kin.at(4), unc, dir);
    
    return kin.at(8) * ff_tt + kin.at(9) * ff_qcd + kin.at(10) * ff_w;
}

#endif  // INCLUDE_APPLYFF_H_
