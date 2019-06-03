import ROOT as r
from glob import glob
r.gROOT.SetBatch(True)

def grab(fin, category, var, name, rebin=None):
    ihist = fin.Get('plots/{}/{}/{}'.format(var, category, name)).Clone()
    if rebin != None:
        ihist.Rebin(rebin)
    return ihist

def get_histos(fout, fname, category, variable, rebin=None):
    """Read all histograms from the input directory and pack
    them into a dictionary"""

    hists = {}
    ss_iso, ss_anti = {}, {}
    fin = r.TFile(fname, 'READ')
    for sname in fin.Get('plots/{}/{}'.format(variable, category)).GetListOfKeys():
        fout.cd()
        hists[sname.GetName()] = grab(fin, category, variable, sname.GetName(), rebin)
        ss_iso['SS_iso_{}'.format(sname.GetName())] = grab(fin, 'SS_iso_{}'.format(category), variable, 'SS_iso_{}'.format(sname.GetName()), rebin)
        ss_anti['SS_anti_{}'.format(sname.GetName())] = grab(fin, 'SS_anti_{}'.format(category), variable, 'SS_anti_{}'.format(sname.GetName()), rebin)

    fin.Close()
    return hists, {
        'ss_iso': ss_iso,
        'ss_anti': ss_anti,
    }

def build_qcd(cat, osss_histos):
    # copy the data histograms
    ss_iso = osss_histos['ss_iso']['SS_iso_data_obs'].Clone()
    ss_anti = osss_histos['ss_anti']['SS_anti_data_obs'].Clone()

    # do data - bkgs
    bkgs = ['ZJ', 'ZL', 'TTT', 'TTJ', 'VVT', 'VVJ', 'W', 'embedded']
    for bkg in bkgs:
        ss_iso.Add(osss_histos['ss_iso']['SS_iso_' + bkg], -1)
        ss_anti.Add(osss_histos['ss_anti']['SS_anti_' + bkg], -1)

    qcd = ss_anti.Clone()
    qcd.SetName('QCD')
    if '0jet' in cat:
        qcd.Scale(1.0 * ss_iso.Integral() / ss_anti.Integral())
    elif 'boosted' in cat:
        qcd.Scale(1.28 * ss_iso.Integral() / ss_anti.Integral())
    elif 'vbf' in cat:
        qcd.Scale(1.0 * ss_iso.Integral() / ss_anti.Integral())

    for ibin in range(0, qcd.GetNbinsX() + 1):
        if qcd.GetBinContent(ibin) < 0:
            qcd.SetBinContent(ibin, 0)

    return qcd

def format_data(data, legend):
    """Make the data look pretty"""

    data.SetLineColor(r.kBlack)
    data.SetFillColor(0)
    data.SetMarkerStyle(8)
    legend.AddEntry(data, 'Data', 'lep')
    return data

def format_backgrounds(histos, legend):
    """Make all the backgrounds look pretty then build the
    stack and the stat uncertainty histogram"""

    histos['embedded'].SetFillColor(r.TColor.GetColor("#f9cd66"))
    histos['ZJ'].SetFillColor(r.TColor.GetColor("#fcc894"))
    histos['ZL'].SetFillColor(r.TColor.GetColor("#de5a6a"))
    histos['TTT'].SetFillColor(r.TColor.GetColor("#cfe87f"))
    histos['TTJ'].SetFillColor(r.TColor.GetColor("#a0abff"))
    histos['VVT'].SetFillColor(r.TColor.GetColor("#9feff2"))
    histos['VVJ'].SetFillColor(r.TColor.GetColor("#9feff2"))
    histos['W'].SetFillColor(r.TColor.GetColor("#d1c7be"))
    histos['QCD'].SetFillColor(r.TColor.GetColor("#ffccff"))

    stack = r.THStack()
    stat = histos['embedded'].Clone()
    stat.Reset()
    stat.SetName('stat')
    sorted_histos = sorted(histos.iteritems(), key=lambda kv: kv[1].Integral())
    for name, ihist in sorted_histos:
        if 'ggh' in name.lower() or 'vbf' in name.lower() or 'wh' in name.lower() or 'zh' in name.lower() or 'qqh' in name.lower():
          continue
        stack.Add(ihist)
        stat.Add(ihist)
        legend.AddEntry(ihist, name, 'f')

    return stack, stat

def formatStat(stat):
    """Make the stat uncertainty histogram pretty"""

    stat.SetMarkerStyle(0)
    stat.SetLineWidth(2)
    stat.SetLineColor(0)
    stat.SetFillStyle(3004)
    stat.SetFillColor(r.kBlack)


def createCanvas():
    """Create the canvas with the right margins, etc."""

    can = r.TCanvas('can', 'can', 432, 451)
    can.Divide(2, 1)

    pad1 = can.cd(1)
    pad1.SetLeftMargin(.12)
    pad1.cd()
    pad1.SetPad(0, .3, 1, 1)
    pad1.SetTopMargin(.1)
    pad1.SetBottomMargin(0.02)
    pad1.SetLogy()
    pad1.SetTickx(1)
    pad1.SetTicky(1)

    pad2 = can.cd(2)
    pad2.SetLeftMargin(.12)
    pad2.SetPad(0, 0, 1, .3)
    pad2.SetTopMargin(0.06)
    pad2.SetBottomMargin(0.35)
    pad2.SetTickx(1)
    pad2.SetTicky(1)

    can.cd(1)

    return can

def format_labels(can):
    """Draw LaTeX things on the canvas"""

    can.cd(1)
    ll = r.TLatex()
    ll.SetNDC(r.kTRUE)
    ll.SetTextSize(0.06)
    ll.SetTextFont(42)
    lepLabel = "#mu#tau_{#mu}"
    lumi = "41.5 fb^{-1}"
    ll.DrawLatex(0.42, 0.92, "{} {}, {} (13 TeV)".format(
        lepLabel, '2017', lumi))

    cms = r.TLatex()
    cms.SetNDC(r.kTRUE)
    cms.SetTextFont(61)
    cms.SetTextSize(0.09)
    cms.DrawLatex(0.16, 0.8, "CMS")

    prel = r.TLatex()
    prel.SetNDC(r.kTRUE)
    prel.SetTextFont(52)
    prel.SetTextSize(0.06)
    prel.DrawLatex(0.16, 0.74, "Preliminary")


def formatRatio(data, stat):
    """Calculate the data/MC ratio then plot it nicely"""

    ratio = data.Clone()
    ratio.Divide(stat)
    ratio.SetTitle('')
    ratio.SetMaximum(2)
    ratio.SetMinimum(0)
    ratio.GetXaxis().SetTitle()
    ratio.SetMarkerStyle(21)
    ratio.GetXaxis().SetTitleSize(0.18)
    ratio.GetXaxis().SetTitleOffset(0.8)
    ratio.GetXaxis().SetTitleFont(42)
    ratio.GetXaxis().SetLabelFont(42)
    ratio.GetXaxis().SetLabelSize(.111)
    ratio.GetXaxis().SetNdivisions(505)

    ratio.GetYaxis().SetTitle('Data / MC')
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetTitleFont(42)
    ratio.GetYaxis().SetTitleOffset(.475)
    ratio.GetYaxis().SetLabelSize(.12)
    ratio.GetYaxis().SetNdivisions(204)

    ratio_unc = ratio.Clone()
    for ibin in range(1, ratio_unc.GetNbinsX()+1):
        ratio_unc.SetBinContent(ibin, 1)
        ratio_unc.SetBinError(ibin, ratio.GetBinError(ibin))
    ratio_unc.SetMarkerSize(0)
    ratio_unc.SetMarkerStyle(8)
    ratio_unc.SetFillColor(r.kGray)

    return ratio, ratio_unc


def sigmaLines(data):
    """Draw lines on ratio plot"""

    low = data.GetBinLowEdge(1)
    high = data.GetBinLowEdge(data.GetNbinsX()) + \
        data.GetBinWidth(data.GetNbinsX())

    # high line
    line1 = r.TLine(low, 1.5, high, 1.5)
    line1.SetLineWidth(1)
    line1.SetLineStyle(3)
    line1.SetLineColor(r.kBlack)

    # low line
    line2 = r.TLine(low, 0.5, high, 0.5)
    line2.SetLineWidth(1)
    line2.SetLineStyle(3)
    line2.SetLineColor(r.kBlack)

    return line1, line2

def main(args):
    # get the input files
    fout = r.TFile('out.root', 'RECREATE')
    histos, osss_histos = get_histos(fout, args.input_dir, args.cat, args.var)
    histos['QCD'] = build_qcd(args.cat, osss_histos)
    # create the histogram legend
    legend = r.TLegend(0.7, 0.65, 0.85, 0.85)
    legend.SetLineColor(0)
    legend.SetFillColor(0)
    legend.SetTextFont(61)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)

    data = format_data(histos.pop('data_obs'), legend)
    stack, stat = format_backgrounds(histos, legend)
    formatStat(stat)

    can = createCanvas()
    r.gStyle.SetOptStat(0)
    high = max(data.GetMaximum(), stat.GetMaximum()) * args.scale
    stack.SetMaximum(high)

    stack.Draw('hist')
    stat.Draw('same e2')
    data.Draw('lep same')
    legend.Draw('same')

    format_labels(can)

    can.cd(2)

    ratio, ratio_unc = formatRatio(data, stat)
    ratio_unc.Draw('same e2')
    ratio.Draw('same lep')

    line1, line2 = sigmaLines(data)
    line1.Draw()
    line2.Draw()

    can.SaveAs('Output/plots/{}_{}_{}_{}.pdf'.format(args.prefix, args.var, args.cat, '2017'))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input-dir', '-i', action='store',
                        dest='input_dir', help='path to input files')
    parser.add_argument('--var', '-v', action='store',
                        dest='var', help='variable to plot')
    parser.add_argument('--cat', '-c', action='store',
                        dest='cat', help='category to plot')
    parser.add_argument('--scale', '-s', action='store', default=1., type=float,
                        dest='scale', help='category to plot')
    parser.add_argument('--prefix', '-p', action='store',
                        dest='prefix', help='prefix for file name')
    main(parser.parse_args())