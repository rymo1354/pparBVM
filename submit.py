import os
import argparse


def write(nodes, cores, time, out, alloc, script, read_se, read_p, write_p, algo, kwargs, options):
    writelines = '#!/bin/bash' + '\n'
    writelines += '#SBATCH -J ' + out + '\n'
    writelines += '#SBATCH --time=' + str(time) + ':00:00' + '\n'
    writelines += '#SBATCH -N ' + str(nodes) + '\n'
    writelines += '#SBATCH --tasks ' + str(cores) + '\n'
    writelines += '#SBATCH -o ' + out + '-%j.out' + '\n'
    writelines += '#SBATCH -e ' + out + '-%j.err' + '\n'
    writelines += '#SBATCH --account=' + alloc + '\n'
    if time == 1:
        writelines += '#SBATCH --partition=debug\n'
    elif time <= 4:
        writelines += '#SBATCH --partition=short\n'
    elif time >= 48:
        writelines += '#SBATCH --partition=long\n'
    else:
        writelines += '#SBATCH --partition=standard\n'

    writelines += 'srun'+' -n '+str(nodes*cores)+' python '+script+' -rse '+read_se+' -rp '+read_p+' -wp '+write_p+' -algo '+algo+' -kw '+"\'{0}\'".format(kwargs)+' -opt '+"\'{0}\'".format(options)+'\n'
    writelines +='exit 0'+'\n'

    with open('submit.sh', 'w') as f:
        f.write(writelines)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    ### General SC run arguments ###
    parser.add_argument('-n', '--nodes', help = 'Number of nodes',
                        type = int, default = 1)
    parser.add_argument('-c', '--cores', help = 'Number of cores',
                        type = int, default = 36)
    parser.add_argument('-t', '--time', help = 'Time limit',
                        type = int, default = 4)
    parser.add_argument('-o', '--outfile', help = 'Outfile name',
                        type = str, required = True)
    parser.add_argument('-a', '--allocation', help = 'Allocation',
                        type = str, default = 'custws')
    parser.add_argument('-s', '--script', help = 'Python script to submit.',
                        type = str, required = True)

    ### run_parameterization.py arguments  ###
    parser.add_argument(
        '-rse', '--read_structures_energies', help='path to .json file with structures and energies', type=str, required=True)
    parser.add_argument(
        '-rp', '--read_parameters', help='path to .json file with BVM parameters', type=str, required=True)
    parser.add_argument(
        '-algo', '--algorithm', help='pyOpt algorithm to use', type=str, required=True)
    parser.add_argument(
        '-kw', '--optimizer_kwargs', help='.json convertible str of pyOpt optimizer kwargs, form \'{"key": "value"}\'', type=str, required=True)
    parser.add_argument(
        '-opt', '--optimizer_options', help='.json convertible str of pyOpt optimizer options, form \'{"key": "value"}\'', type=str, required=True)
    parser.add_argument(
        '-wp', '--write_parameters', help='path to .json file of parameterized bond valence parameters', type=str, required=False)
    args = parser.parse_args()

    write(args.nodes, args.cores, args.time, args.outfile, args.allocation,
          args.script, args.read_structures_energies, args.read_parameters, args.write_parameters, args.algorithm, args.optimizer_kwargs, args.optimizer_options)
    os.system('sbatch submit.sh')
