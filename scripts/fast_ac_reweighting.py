import multiprocessing
from glob import glob
from tqdm import tqdm
from pprint import pprint
from subprocess import call, Popen


def to_reweight(ifile):
    """List of signal samples. Only processes these files."""
    for name in [
        'vbf125_JHU.root', 'wh125_JHU.root', 'zh125_JHU.root'
    ]:
        if name in ifile:
            return True
    return False


def call_cmd(cmd):
    call(cmd, shell=True)
    return None


def main(args):
    input_directories = [idir for idir in glob('{}/*'.format(args.input)) if not 'logs' in idir]
    input_files = {
        idir: [ifile for ifile in glob('{}/merged/*.root'.format(idir)) if to_reweight(ifile)]
        for idir in input_directories
    }
    input_files = dict(filter(lambda x: len(x[1]) > 0, input_files.items()))

    print 'Directory structure to process'
    pprint(input_files, width=150)

    temp_name = 'tmp/{}'.format(args.input.split('/')[-1])
    call('mkdir -p {}'.format(temp_name), shell=True)

    n_processes = min(8, multiprocessing.cpu_count() / 2)
    pool = multiprocessing.Pool(processes=n_processes)

    ndir = len(input_files.keys())
    pbar = tqdm(input_files.items())
    for idir, files in pbar:
        pbar.set_description('Processing: {}'.format(idir.split('/')[-1]))
        # start reweighting things and wait to complete
        jobs = [
            pool.apply_async(call_cmd, ('bin/ac-reweight -n {} -t {} -o {}/'.format(ifile, args.tree_name, temp_name))) for ifile in files
        ]
        [j.wait() for j in jobs]

        # move output files
        call('mv {}/*.root {}/merged'.format(temp_name, idir), shell=True)
        procs = [Popen('mv {} {}'.format(ifile, ifile.replace('merged', '')), shell=True)
                 for ifile in files]  # move from "merged" to parent directory
        for p in procs:
            p.wait()

    # make sure everything in the pool is finished
    pool.close()
    pool.join()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='path to input files')
    parser.add_argument('--tree-name', '-t', required=True, help='name of TTree')
    main(parser.parse_args())
