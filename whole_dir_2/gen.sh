#!/bin/sh

# SLURM options:

#SBATCH --job-name=generator_         # Job name
#SBATCH --output=save_logs/generator_%j.log     # Standard output and error log

#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=6                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=18000                   # Memory in MB per default
#SBATCH --time=1-00:00:00             # 7 days by default on htc partition

#SBATCH --mail-user=weymann@ijclab.in2p3.fr   # Where to send mail
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

# Commands to be submitted:

module load python

#python generator.py -infl_type='lle' -bp=1 -ln1010As=3.047 -ns=0.9665 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='lle' -bp=1 -ln1010As=3.047 -ns=0.9703 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='lle' -bp=1 -ln1010As=3.047 -ns=0.9627 -lnRrad=0 -xaxis='A6'

#python generator.py -infl_type='udd' -bp=1 -ln1010As=3.047 -ns=0.9665 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='udd' -bp=1 -ln1010As=3.047 -ns=0.9703 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='udd' -bp=1 -ln1010As=3.047 -ns=0.9627 -lnRrad=0 -xaxis='A6'

#python generator.py -infl_type='lle' -bp=2 -ln1010As=3.047 -ns=0.9665 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='lle' -bp=2 -ln1010As=3.047 -ns=0.9703 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='lle' -bp=2 -ln1010As=3.047 -ns=0.9627 -lnRrad=0 -xaxis='A6'

#python generator.py -infl_type='udd' -bp=2 -ln1010As=3.047 -ns=0.9665 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='udd' -bp=2 -ln1010As=3.047 -ns=0.9703 -lnRrad=0 -xaxis='A6'
#python generator.py -infl_type='udd' -bp=2 -ln1010As=3.047 -ns=0.9627 -lnRrad=0 -xaxis='A6'

python generator.py -infl_type='tree' -ln1010As=3.047 -ns=0.9665 -lnRrad=0 -xaxis='A6'
python generator.py -infl_type='tree' -ln1010As=3.047 -ns=0.9703 -lnRrad=0 -xaxis='A6'
python generator.py -infl_type='tree' -ln1010As=3.047 -ns=0.9627 -lnRrad=0 -xaxis='A6'
