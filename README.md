# strepto2
Evolution of mutation-driven division of labour in Streptomyces

This code accompanies the article:
Evolution of genome fragility enables microbial division of labor.
Authors: E.S. Colizzi, B. van Dijk, R.M.H. Merks, D.E. Rozen, R.M.A. Vroomans
Which is currently under revision. 
Preprint: https://www.biorxiv.org/content/10.1101/2021.06.04.447040v3

The source code is written in c, and uses a custom version of the cash library 
(which you can find here: https://github.com/escolizzi/cash2-s.BINF2018). 
It runs with Ubuntu (and, I expect, with any reasonable OS built on Linux) but I have not tested it with other OS. 

NOTICE: 
- Depending on your machine, a single simulation with default parameters will take anywhere between 3 days and 2 weeks. 
- Simulations with non-default parameters can take up to 3 months.
- The code produces A LOT of data, from a few GB up to 100GB per run. 

To run this code do the following:
- Install cash from: https://github.com/escolizzi/cash2-s.BINF2018
- Download strepto.c
- In the terminal, type: 
    make strepto
- Once the file has compiled, you can run it with: 
    ./strepto [options]
The default parameter values run the evolutionary simulation described in Fig. 2. 
For the various options you can open strepto.c and look for the function Initial(), where the command line options are screened. 

The output is organised as follows:
- a data file with data about the population (saved periodically)
- a directory with images that give a rough impression about the field. 

To reproduce the figures in the manuscript, you'll need the plotting scripts, which are in the directory "scripts".
They are written in Python3. 
To make the time plots (e.g. inf Fig. 2a), type: 
./plot_strepto_intime_abdiversity_v2.py your_data_file.txt

To make the analysis, you'll need to re-run one growth cycle with a slightly different version of the code that outputs more fine grained data.
For this first you'll need to extract the time point you are interested in from your data file

  ( which you can grep from the data file as follows: 
    say you want the growth cycle at time 3000000 from a simulation that lasted until 5000000 time steps, the type:
    grep 3000000 your_data_file.txt > your_data_file_time3000000.txt
  ) 

Then, you'll need to start a new simulation with: 
  ./strepto_oneseason_lotsadata -input your_data_file_time3000000.txt -change_season_at_initialisation_from_input 1 -maxtime 2500
  (you can add additional options)
This produces an additional log file which details every single replication and antibiotic production event.
With this file, you can run:
./correlate_F_and_ABorR_intime.py your_log_file.txt 
which produces a few figures similar to Fig. 2c, 2d, 3c and 3d. 

To compare one feature of multiple runs (e.g. see Fig. SF5), take all the population at one time step (see above) from all the files you want to compare. 
Say that you have a bunch of files like: data_1.txt data_2.txt data_3.txt etc 
then use the script: 
./plot_genomelendistr_v2.py [FEATURE POSITION IN FILENAME] "data_\*.txt"
Where [FEATURE POSITION IN FILENAME] in this case would be 1, as it is read from the file name, split with the character "_".

It is my experience that these descriptions are rarely complete, and for this code in particular, some previous experience with the cash libraries might be useful. For any question please contact me at: enricosandro [dot] colizzi [at] gmail [dot] com
Also, please let me know if you find bugs :)
