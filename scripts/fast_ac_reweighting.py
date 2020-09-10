from glob import glob
from tqdm import tqdm
from pprint import pprint
from subprocess import call

def to_reweight(ifile):
    """List of signal samples. Only processes these files."""
    for name in [
        'vbf125_JHU.root', 'wh125_JHU.root', 'zh125_JHU.root'
    ]:
        if name in ifile:
            return True
    return False


def main(args):
    input_directories = [idir for idir in glob('{}/*'.format(args.input)) if not 'logs' in idir]
    input_files = {
        idir: [ifile for ifile in glob('{}/merged/*.root'.format(idir)) if to_reweight(ifile)]
        for idir in input_directories
    }
    print 'Directory structure to process'
    pprint(input_files, width=150)

    temp_name = args.input
    if '/hdfs' in temp_name:
        temp_name = 'tmp/{}'.format(args.input.split('/')[-1])
        call('mkdir -p {}'.format(temp_name), shell=True)

    i = 0
    ndir = len(input_files.keys())
    for idir, files in tqdm(input_files.iteritems()):
        i += 1
        for ifile in tqdm(files, leave=False):
            if '/hdfs' in args.input:
                call('bin/ac-reweight -n {} -t {} -o {}/'.format(ifile, args.tree_name, temp_name), shell=True)
                call('mv {}/*.root {}'.format(temp_name, idir), shell=True)
            else:
                call('bin/ac-reweight -n {} -t {} -o {}/merged'.format(ifile, args.tree_name, idir), shell=True)
            fname = ifile.split('/')[-1]
            call('mv {} {}'.format(ifile, ifile.replace(fname, '/merged')), shell=True) # move from "merged" to parent directory


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='path to input files')
    parser.add_argument('--tree-name', '-t', required=True, help='name of TTree')
    main(parser.parse_args())
