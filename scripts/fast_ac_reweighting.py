import json
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


def recognize_signal(ifile):
    """Pick the correct keys for this sample."""
    process = ifile.split('/')[-1].split('125')[0]
    key = 'jhu_ac_reweighting_map'
    return key, process


def call_cmd(ifile, tree_name, temp_name, syst_dir, weight, name):
    call('bin/ac-reweight -n {} -t {} -o {}/{}/{}.root -w {}'.format(ifile, tree_name, temp_name, syst_dir, name, weight), shell=True)
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

    boilerplate = {}
    with open('configs/boilerplate.json', 'r') as config_file:
        boilerplate = json.load(config_file)

    temp_name = 'tmp/{}'.format(args.input.split('/')[-1])
    call('mkdir -p {}'.format(temp_name), shell=True)

    n_processes = min(8, multiprocessing.cpu_count() / 2)
    pool = multiprocessing.Pool(processes=n_processes)

    ndir = len(input_files.keys())
    pbar = tqdm(input_files.items())
    jobs = []
    for idir, files in pbar:
        pbar.set_description('Starting: {}'.format(idir.split('/')[-1]))
        syst_dir = idir.split('/')[-1]
        call('mkdir -p {}/{}'.format(temp_name, syst_dir), shell=True)

        for ifile in files:
            key, process = recognize_signal(ifile)
            weight_names = boilerplate[key][process]

            # start reweighting things and wait to complete
            jobs.append([
                pool.apply_async(call_cmd, (ifile, args.tree_name, temp_name, syst_dir, weight, name)) for weight, name in weight_names
            ])

    # make sure everything in the pool is finished
    [j.get() for job in jobs for j in job]
    print 'Waiting for processes to finish...'
    pool.close()
    pool.join()

    # move output files
    out_dirs = [idir for idir in glob('{}/'.format(temp_name))]
    pbar = tqdm(out_dirs)
    for idir in pbar:
        pbar.set_description('Moving: {}'.format(idir.split('/')[-1]))
        to_name = idir.replace('tmp/', '/hdfs/store/user/tmitchel/')
        call('mv {}/*.root {}/merged/'.format(idir, to_name), shell=True)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='path to input files')
    parser.add_argument('--tree-name', '-t', required=True, help='name of TTree')
    main(parser.parse_args())
