import ROOT
import time
import pandas
import numpy
import uproot
import os
from array import array
from subprocess import call

mvis_bins = [0, 50, 80, 100, 110, 120, 130, 150, 170, 200, 250, 1000]
njets_bins = [-0.5, 0.5, 1.5, 15]


def get_categories(channel):
    """Return list of categories with the lepton prefix added."""
    return ['inclusive', '0jet', 'boosted', 'vbf']


def build_histogram(name):
    """Build TH2F to fill with fake fraction."""
    return ROOT.TH2F(name, name, len(mvis_bins) - 1, array('d', mvis_bins), len(njets_bins) - 1, array('d', njets_bins))


def fill_fraction(df, fraction):
    """Use visible Higgs mass and njets to fill this histogram."""
    vis_mass = df['vis_mass'].values
    njets = df['njets'].values
    evtwt = df['evtwt'].values
    for i in xrange(len(df.index)):
        fraction.Fill(vis_mass[i], njets[i], evtwt[i])


def parse_tree_name(keys):
    """Take list of keys in the file and search for our TTree"""
    if 'et_tree;1' in keys:
        return 'et_tree'
    elif 'mt_tree;1' in keys:
        return 'mt_tree'
    else:
        raise Exception('Can\t find et_tree or mt_tree in keys: {}'.format(keys))


def main(args):
    start = time.time()

    # read info from data file
    open_file = uproot.open('{}/data_obs.root'.format(args.input))
    keys = open_file.keys()
    tree_name = parse_tree_name(keys)
    oldtree = open_file[tree_name].arrays(['*'])
    treedict = {ikey: oldtree[ikey].dtype for ikey in oldtree.keys()}
    pre_jet_fakes = open_file[tree_name].arrays('*', outputtype=pandas.DataFrame)

    channel_prefix = tree_name[:2]
    fake_file_path = '/hdfs/store/user/tmitchel/HTT_FakeFactors/ff_files_{}_{}/'.format(channel_prefix, args.year)
    fake_fraction_output_name = 'Output/fake_fractions/{}{}_{}.root'.format(channel_prefix, args.year, args.suffix)

    fout = ROOT.TFile(fake_fraction_output_name, 'recreate')
    categories = get_categories(channel_prefix)
    for cat in categories:
        fout.cd()
        fout.mkdir(cat)
        fout.cd()

    inputs = {
        'frac_w': ['W', 'ZJ', 'VVJ', 'STJ'],
        'frac_tt': ['TTJ'],
        'frac_data': ['data_obs'],
        'frac_real': ['STL', 'VVL', 'TTL', 'ZL', 'STT', 'VVT', 'TTT', 'embed'],
    }

    fractions = {
        'frac_w': {cat: build_histogram('frac_w_{}'.format(cat)) for cat in categories},
        'frac_tt': {cat: build_histogram('frac_tt_{}'.format(cat)) for cat in categories},
        'frac_qcd': {cat: build_histogram('frac_qcd_{}'.format(cat)) for cat in categories},
        'frac_data': {cat: build_histogram('frac_data_{}'.format(cat)) for cat in categories},
        'frac_real': {cat: build_histogram('frac_real_{}'.format(cat)) for cat in categories},
    }

    for frac, samples in inputs.iteritems():
        for sample in samples:
            print sample
            events = uproot.open('{}/{}.root'.format(args.input, sample)
                                 )[tree_name].arrays('*', outputtype=pandas.DataFrame)

            anti_iso_events = events[
                (events['is_antiTauIso'] > 0) & (events['contamination'] == 0)
            ]

            zero_jet_events = anti_iso_events[anti_iso_events['njets'] == 0]
            boosted_events = anti_iso_events[
                (anti_iso_events['njets'] == 1) |
                ((anti_iso_events['njets'] > 1) & anti_iso_events['mjj'] < 300)
            ]
            vbf_events = anti_iso_events[(anti_iso_events['njets'] > 1) & (anti_iso_events['mjj'] > 300)]

            # inclusive region
            fill_fraction(anti_iso_events, fractions[frac]['inclusive'])
            fill_fraction(zero_jet_events, fractions[frac]['0jet'])
            fill_fraction(boosted_events, fractions[frac]['boosted'])
            fill_fraction(vbf_events, fractions[frac]['vbf'])

            # fill pre_jet_fakes with anti-isolated events
            if sample in inputs['frac_real']:
                events['evtwt'] = -1 * events['evtwt']
                pre_jet_fakes = pandas.concat([pre_jet_fakes, events])

    pre_jet_fakes = pre_jet_fakes[(pre_jet_fakes['is_antiTauIso'] > 0)]

    for cat in categories:
        fractions['frac_qcd'][cat] = fractions['frac_data'][cat].Clone()
        fractions['frac_qcd'][cat].Add(fractions['frac_w'][cat], -1)
        fractions['frac_qcd'][cat].Add(fractions['frac_tt'][cat], -1)
        fractions['frac_qcd'][cat].Add(fractions['frac_real'][cat], -1)

        # handle bins that go negative
        for xbin in range(0, fractions['frac_qcd'][cat].GetNbinsX() + 1):
            for ybin in range(0, fractions['frac_qcd'][cat].GetNbinsY() + 1):
                if fractions['frac_qcd'][cat].GetBinContent(xbin, ybin) < 0:
                    fractions['frac_qcd'][cat].SetBinContent(xbin, ybin, 0.)

        denom = fractions['frac_qcd'][cat].Clone()
        denom.Add(fractions['frac_w'][cat])
        denom.Add(fractions['frac_tt'][cat])

        print 'Category: {}'.format(cat)
        print '\tfrac_w: {}'.format(fractions['frac_w'][cat].Integral() / denom.Integral())
        print '\tfrac_tt: {}'.format(fractions['frac_tt'][cat].Integral() / denom.Integral())
        print '\tfrac_qcd: {}'.format(fractions['frac_qcd'][cat].Integral() / denom.Integral())
        print '\tfrac_real: {}'.format(fractions['frac_real'][cat].Integral() / denom.Integral())

        fractions['frac_w'][cat].Divide(denom)
        fractions['frac_tt'][cat].Divide(denom)
        fractions['frac_qcd'][cat].Divide(denom)

    # write the fractions to disk
    for frac_name, categories in fractions.iteritems():
        for cat_name, ihist in categories.iteritems():
            fout.cd(cat_name)
            ihist.Write(frac_name)

    fout.Close()

    # write the prefakes file to disk
    tmp_path = 'tmp/{}/'.format(args.suffix)
    tmp_file_path = tmp_path + 'pre_jetFakes.root'
    call('mkdir -p {}'.format(tmp_path), shell=True)
    with uproot.recreate(tmp_file_path) as f:
        f[tree_name] = uproot.newtree(treedict)
        f[tree_name].extend(pre_jet_fakes.to_dict('list'))


    # call C++ binary to fill weights
    callstring = './bin/create-fakes -i {} -p {} -f {} -c {}'.format(
        tmp_path, fake_fraction_output_name, fake_file_path, channel_prefix)
    if args.syst:
        callstring += ' -s '
    print callstring
    call(callstring, shell=True)
    os.system('mv -v {} {}'.format(tmp_file_path.replace('pre_jetFakes', 'jetFakes'), args.input),)
    print 'Finished in {} seconds'.format(time.time() - start)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='path to input files')
    parser.add_argument('--suffix', '-s', required=True, help='string to append to output file name')
    parser.add_argument('--year', '-y', required=True, help='year being processed')
    parser.add_argument('--syst', action='store_true', help='process systematics too')
    main(parser.parse_args())
