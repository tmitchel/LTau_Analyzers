from os import system
from glob import glob
from pprint import pprint


def main(args):
    to_migrate = {
        idir.split('/')[-1]: [
            ifile.split('/')[-1] for ifile in glob('{}/*.root'.format(idir))
        ] for idir in glob('{}/*'.format(args.input))
    }

    for idir, files in to_migrate.iteritems():
        for ifile in files:
            system('cp -v {0}/{1}/{2} {3}/{1}/{2}'.format(args.input, idir, ifile, args.target))

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='input directory')
    parser.add_argument('--target', '-t', required=True, help='directory to merge into')
    main(parser.parse_args())
