import ROOT
from pprint import pprint
from collections import namedtuple
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)


# short-hand
GetColor = ROOT.TColor.GetColor
black = ROOT.kBlack
no_color = 0
signal_scaling = 50.

style_map_tuple = namedtuple('style_map_tuple', [
    'fill_color', 'line_color', 'line_style', 'line_width', 'marker_style'
])
style_map = {
    "data_obs": style_map_tuple(no_color, black, 1, 1, 8),
    "backgrounds": {
        "embedded": style_map_tuple(GetColor("#f9cd66"), black, 1, 1, 1),
        "ZTT": style_map_tuple(GetColor("#f9cd66"), black, 1, 1, 1),
        "jetFakes": style_map_tuple(GetColor("#ffccff"), black, 1, 1, 1),
        "ZL": style_map_tuple(GetColor("#de5a6a"), black, 1, 1, 1),
        "TTT": style_map_tuple(GetColor("#cfe87f"), black, 1, 1, 1),
        "TTL": style_map_tuple(GetColor("#cfe87f"), black, 1, 1, 1),
        "VVT": style_map_tuple(GetColor("#9feff2"), black, 1, 1, 1),
        "VVL": style_map_tuple(GetColor("#9feff2"), black, 1, 1, 1),
        "STT": style_map_tuple(GetColor("#9feff2"), black, 1, 1, 1),
        "STL": style_map_tuple(GetColor("#9feff2"), black, 1, 1, 1),
    },
    "signals": {
        "ggh125_powheg": style_map_tuple(no_color, no_color, 0, 0, 1),  # don't show powheg
        "reweighted_ggH_htt_0PM125": style_map_tuple(no_color, GetColor("#0000FF"), 1, 3, 1),
        "reweighted_ggH_htt_0Mf05ph0125": style_map_tuple(no_color, GetColor("#AA00FF"), 1, 3, 1),
        "reweighted_ggH_htt_0M125": style_map_tuple(no_color, GetColor("#00AAFF"), 1, 3, 1),

        "vbf125_powheg": style_map_tuple(no_color, no_color, 0, 0, 1),  # don't show powheg
        "reweighted_qqH_htt_0PM125": style_map_tuple(no_color, GetColor("#FF0000"), 1, 3, 1),
        "reweighted_qqH_htt_0Mf05ph0125": style_map_tuple(no_color, GetColor("#FF00AA"), 1, 3, 1),
        "reweighted_qqH_htt_0M125": style_map_tuple(no_color, GetColor("#ff5e00"), 1, 3, 1),
        "reweighted_qqH_htt_0PH125": style_map_tuple(no_color, GetColor("#ff5e00"), 1, 3, 1),
        "reweighted_qqH_htt_0L1125": style_map_tuple(no_color, GetColor("#ff5e00"), 1, 3, 1),
        "reweighted_qqH_htt_0L1Zg125": style_map_tuple(no_color, GetColor("#ff5e00"), 1, 3, 1),
    }
}


def clean_samples(input_histos):
    merged = {
        'ZTT': input_histos['ZTT'].Clone(),
        'ZL': input_histos['ZL'].Clone(),
        'jetFakes': input_histos['jetFakes'].Clone(),
        'tt': input_histos['TTT'].Clone(),
        'Others': input_histos['STT'].Clone()
    }

    merged['tt'].Add(input_histos['TTL'])
    for name in ['STL', 'VVT', 'VVL']:
        merged['Others'].Add(input_histos[name])

    return merged


def ApplyStyle(ihist, styles):
    """Apply styling to histogram."""
    ihist.SetFillColor(styles.fill_color)
    ihist.SetLineColor(styles.line_color)
    ihist.SetLineStyle(styles.line_style)
    ihist.SetLineWidth(styles.line_width)
    ihist.SetMarkerStyle(styles.marker_style)
    return ihist


def createCanvas():
    """Build the TCanvas and split into pads."""
    can = ROOT.TCanvas('can', 'can', 432, 451)
    can.Divide(2, 1)

    pad1 = can.cd(1)
    pad1.SetLeftMargin(.12)
    pad1.cd()
    pad1.SetPad(0, .3, 1, 1)
    pad1.SetTopMargin(.1)
    pad1.SetBottomMargin(0.02)
#    pad1.SetLogy()
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


def formatStat(stat):
    """Format the statistical uncertainty histogram"""
    stat.SetMarkerStyle(0)
    stat.SetLineWidth(2)
    stat.SetLineColor(0)
    stat.SetFillStyle(3004)
    stat.SetFillColor(ROOT.kBlack)
    return stat


def formatStack(stack):
    """Format the stacked background histograms."""
    stack.GetXaxis().SetLabelSize(0)
    stack.GetYaxis().SetTitle('Events / Bin')
    stack.GetYaxis().SetTitleFont(42)
    stack.GetYaxis().SetTitleSize(.05)
    stack.GetYaxis().SetTitleOffset(1.2)
    stack.SetTitle('')
    stack.GetXaxis().SetNdivisions(505)


def fillLegend(data, backgrounds, signals, stat, signal_plots):
    """Fill the legend with appropriate processes."""
    leg = ROOT.TLegend(0.5, 0.45, 0.85, 0.85)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetTextFont(61)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)

    # data
    leg.AddEntry(data, 'Data', 'lep')

    # signals
    leg.AddEntry(signals[signal_plots[0][0]], signal_plots[0][1], 'l')
    leg.AddEntry(signals[signal_plots[1][0]], signal_plots[1][1], 'l')

    # backgrounds
    leg.AddEntry(backgrounds['ZTT'], 'ZTT', 'f')
    leg.AddEntry(backgrounds['ZL'], 'ZL', 'f')
    leg.AddEntry(backgrounds['jetFakes'], 'Jet Mis-ID', 'f')
    leg.AddEntry(backgrounds['tt'], 'tt', 'f')
    leg.AddEntry(backgrounds['Others'], 'Others', 'f')

    # stat. uncertainty
    leg.AddEntry(stat, 'Uncertainty', 'f')

    return leg


def formatPull(pull, title):
    """Format the pull (or ratio) histogram in the lower pad."""
    pull.SetTitle('')
    pull.SetMaximum(1.3)
    pull.SetMinimum(0.7)
    pull.GetXaxis().SetTitle(title)
    pull.SetMarkerStyle(21)
    pull.GetXaxis().SetTitleSize(0.18)
    pull.GetXaxis().SetTitleOffset(0.8)
    pull.GetXaxis().SetTitleFont(42)
    pull.GetXaxis().SetLabelFont(42)
    pull.GetXaxis().SetLabelSize(.111)
    pull.GetXaxis().SetNdivisions(505)
    # pull.GetXaxis().SetLabelSize(0)
    # pull.GetXaxis().SetTitleSize(0)

    pull.GetYaxis().SetTitle('Obs. / Exp.')
    pull.GetYaxis().SetTitleSize(0.12)
    pull.GetYaxis().SetTitleFont(42)
    pull.GetYaxis().SetTitleOffset(.475)
    pull.GetYaxis().SetLabelSize(.12)
    pull.GetYaxis().SetNdivisions(204)
    return pull


def sigmaLines(data):
    """Draw lines on pull (or ratio) plot."""
    low = data.GetBinLowEdge(1)
    high = data.GetBinLowEdge(data.GetNbinsX()) + \
        data.GetBinWidth(data.GetNbinsX())

    # high line
    line1 = ROOT.TLine(low, 1.2, high, 1.2)
    line1.SetLineWidth(1)
    line1.SetLineStyle(3)
    line1.SetLineColor(ROOT.kBlack)

    # low line
    line2 = ROOT.TLine(low, 0.8, high, 0.8)
    line2.SetLineWidth(1)
    line2.SetLineStyle(3)
    line2.SetLineColor(ROOT.kBlack)

    # nominal line
    line3 = ROOT.TLine(low, 1., high, 1.)
    line3.SetLineWidth(1)
    line3.SetLineStyle(3)
    line3.SetLineColor(ROOT.kBlack)

    return line1, line2, line3


def blindData(data, signal, background):
    """Apply blinding procedure to data."""
    for ibin in range(data.GetNbinsX()+1):
        sig = signal.GetBinContent(ibin) / signal_scaling
        bkg = background.GetBinContent(ibin)
        if bkg > 0 and sig / ROOT.TMath.Sqrt(bkg + pow(0.09*bkg, 2)) >= .3:
            err = data.GetBinError(ibin)
            data.SetBinContent(ibin, -1)
            data.SetBinError(ibin, err)

    return data


def BuildPlot(args):
    """
    Build the stacked plot with everything included and formatted then save as PDF.

    Variables (inside args):
    input       -- input TFile full of histograms
    category    -- which TDirectory to read
    variable    -- which variable to plot
    scale       -- value to scale the top of the plot (keep histograms from being cutoff)
    year        -- which era is this?
    label       -- LaTeX label for variable on x-axis
    prefix      -- name to attach to output file
    """
    ifile = ROOT.TFile(args.input)
    category = ifile.Get(args.category)
    variable = category.Get(args.variable)

    # start getting histograms
    data_hist = variable.Get('data_obs').Clone()
    signals = {}
    backgrounds = {}

    # loop through histograms to read and store to dict
    for hkey in variable.GetListOfKeys():
        hname = hkey.GetName()
        ihist = variable.Get(hname).Clone()
        if hname in style_map['backgrounds']:
            ihist = ApplyStyle(ihist, style_map['backgrounds'][hname])
            backgrounds[hname] = ihist
        elif hname in style_map['signals']:
            ihist = ApplyStyle(ihist, style_map['signals'][hname])
            signals[hname] = ihist

    # merge backgrounds
    backgrounds = clean_samples(backgrounds)

    # now get stat and stack filled
    stat = data_hist.Clone()  # sum of all backgrounds
    stat.Reset()
    stack = ROOT.THStack()  # stack of all backgrounds
    for bkg in sorted(backgrounds.itervalues(), key=lambda hist: hist.Integral()):
        stat.Add(bkg)
        stack.Add(bkg)

    sig_yields = [ihist.GetMaximum() for ihist in [signals['ggh125_powheg'], signals['vbf125_powheg']]] + [data_hist.GetMaximum(), stat.GetMaximum()]
    stack.SetMaximum(max(sig_yields) * args.scale)

    # format the plots
    can = createCanvas()
    data_hist = ApplyStyle(data_hist, style_map['data_obs'])
    stat = formatStat(stat)
    stack.Draw('hist')
    formatStack(stack)

    for sig_name, sig_hist in signals.iteritems():
        if 'ggH' in sig_name:
            sig_hist.Scale(signal_scaling*signals['ggh125_powheg'].Integral()/sig_hist.Integral())
        if 'qqH' in sig_name:
            sig_hist.Scale(signal_scaling*signals['vbf125_powheg'].Integral()/sig_hist.Integral())

    signal_plots = []
    if args.variable == 'D0_ggH':
        signal_plots = [('reweighted_ggH_htt_0PM125', 'ggH SM Higgs(125)x50'), ('reweighted_ggH_htt_0M125', 'ggH PS Higgs(125)x50')]
    elif args.variable == 'D0_VBF':
        signal_plots = [('reweighted_qqH_htt_0PM125', 'VBF SM Higgs(125)x50'), ('reweighted_qqH_htt_0M125', 'VBF PS Higgs(125)x50')]
    elif args.variable == 'D_a2_VBF':
        signal_plots = [('reweighted_qqH_htt_0PM125', 'VBF SM Higgs(125)x50'), ('reweighted_qqH_htt_0PH125', 'VBF a2 Higgs(125)x50')]
    elif args.variable == 'D_l1_VBF':
        signal_plots = [('reweighted_qqH_htt_0PM125', 'VBF SM Higgs(125)x50'), ('reweighted_qqH_htt_0L1125', 'VBF #lambda_{1} Higgs(125)x50')]
    elif args.variable == 'D_l1zg_VBF':
        signal_plots = [('reweighted_qqH_htt_0PM125', 'VBF SM Higgs(125)x50'), ('reweighted_qqH_htt_0L1Zg125', 'VBF #lambda_{1zg} Higgs(125)x50')]
    elif args.variable == 'DCP_ggH':
        signal_plots = [('reweighted_ggH_htt_0PM125', 'ggH SM Higgs(125)x50'), ('reweighted_ggH_htt_0Mf05ph0125', 'ggH MM Higgs(125)x50')]
    elif args.variable == 'DCP_VBF':
        signal_plots = [('reweighted_qqH_htt_0PM125', 'VBF SM Higgs(125)x50'), ('reweighted_qqH_htt_0Mf05ph0125', 'VBF MM Higgs(125)x50')]
    else:
        signal_plots = [('reweighted_ggH_htt_0PM125', 'ggH SM Higgs(125)x50'), ('reweighted_qqH_htt_0PM125', 'VBF SM Higgs(125)x50')]

    combo_signal = signals[signal_plots[0][0]].Clone()
    combo_signal.Add(signals[signal_plots[1][0]])
    data_hist = blindData(data_hist, combo_signal, stat)

    # draw the plots
    data_hist.Draw('same lep')
    stat.Draw('same e2')

    signals[signal_plots[0][0]].Draw('same hist')
    signals[signal_plots[1][0]].Draw('same hist')

    legend = fillLegend(data_hist, backgrounds, signals, stat, signal_plots)
    legend.Draw('same')

    # do some printing on the canvas
    ll = ROOT.TLatex()
    ll.SetNDC(ROOT.kTRUE)
    ll.SetTextSize(0.06)
    ll.SetTextFont(42)
    if 'et_' in args.category:
        lepLabel = "#tau_{e}#tau_{h}"
    elif 'mt_' in args.category:
        lepLabel = "#tau_{#mu}#tau_{h}"
    if args.year == '2016':
        lumi = "35.9 fb^{-1}"
    elif args.year == '2017':
        lumi = "41.5 fb^{-1}"
    elif args.year == '2018':
        lumi = "59.7 fb^{-1}"
    elif args.year == 'all':
        lumi = "137.0 fb^{-1}"

    year_label = args.year
    if args.year == 'all':
      year_label = ''
    ll.DrawLatex(0.42, 0.94, "{} {}, {} (13 TeV)".format(lepLabel, year_label, lumi))

    cms = ROOT.TLatex()
    cms.SetNDC(ROOT.kTRUE)
    cms.SetTextFont(61)
    cms.SetTextSize(0.09)
    cms.DrawLatex(0.16, 0.8, "CMS")

    prel = ROOT.TLatex()
    prel.SetNDC(ROOT.kTRUE)
    prel.SetTextFont(52)
    prel.SetTextSize(0.06)
    prel.DrawLatex(0.16, 0.74, "Preliminary")

    if args.category == 'et_inclusive' or args.category == 'mt_inclusive':
        catName = 'Inclusive'
    elif args.category == 'et_0jet' or args.category == 'mt_0jet':
        catName = '0-Jet'
    elif args.category == 'et_boosted' or args.category == 'mt_boosted':
        catName = 'Boosted'
    elif args.category == 'et_vbf' or args.category == 'mt_vbf':
        catName = 'VBF Category'
    else:
        catName = ''

    lcat = ROOT.TLatex()
    lcat.SetNDC(ROOT.kTRUE)
    lcat.SetTextFont(42)
    lcat.SetTextSize(0.06)
    lcat.DrawLatex(0.16, 0.68, catName)

    # now work on ratio plot
    can.cd(2)
    ratio = data_hist.Clone()
    ratio.Divide(stat)
    ratio = formatPull(ratio, args.label)
    rat_unc = ratio.Clone()
    for ibin in range(1, rat_unc.GetNbinsX()+1):
        rat_unc.SetBinContent(ibin, 1)
        rat_unc.SetBinError(ibin, ratio.GetBinError(ibin))
    rat_unc = formatStat(rat_unc)
    # rat_unc.SetMarkerSize(0)
    # rat_unc.SetMarkerStyle(8)

    # rat_unc.SetFillColor(ROOT.kGray)
    rat_unc.Draw('same e2')
    ratio.Draw('same lep')

    line1, line2, line3 = sigmaLines(data_hist)
    line1.Draw()
    line2.Draw()
    line3.Draw()

    # save the pdf
    can.SaveAs('Output/plots/{}_{}_{}_{}.pdf'.format(args.prefix, args.variable, args.category, args.year))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', required=True, action='store', help='input file')
    parser.add_argument('--year', '-y', required=True, action='store', help='year to plot')
    parser.add_argument('--category', '-c', required=True, action='store', help='category to plot')
    parser.add_argument('--variable', '-v', required=True, action='store', help='variable to plot')
    parser.add_argument('--label', '-l', required=True, action='store', help='label for plot')
    parser.add_argument('--prefix', '-p', action='store', help='prefix for output name')
    parser.add_argument('--scale', '-s', default=1.2, type=float, action='store', help='scale max by x')
    BuildPlot(parser.parse_args())
