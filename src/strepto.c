/*
In this version antibiotic diversity can evolve. 
Antibiotic are 16 bits elements, resistance depends on dist. 
Antib are actually coded as integers. :P

strepto.c, a simulation system for streptomyces evolution, based on cash.
Bacteria have a genome made of three types of genes: 'F', 'A', 'B'. 
F genes are needed for growth, A genes for antibiotic production, B are break points.
In some cases there is a fourth gene type 'H'. H genes are for homeostasis.
Bacteria live on a grid, which is TYPE2** world.
Antibiotics are on a second TYPE2** plane, called antib.

The variable TYPE2 is defined in cash2.h, here it takes different meaning depending on whetehr it is used for the bacteria place (world), or the antibiotic plane antib

BACTERIA:

typedef struct __type2{
  int val;  //presence/absence of bacterium
  int val2; // 0 = no bacterium, > 0 indicates strain nr. (inherited from spore)
  int val3; auxiliary variable for nr. of F genes, calculated at birth from the genome (so that I don't have to recalculate everything every time step), see function Regulation0()
  int val4; same as above, for nr. of A
  int val5; same as above, nr. of H (if used)
  double fval; //unused
  double fval2;//unused
  double fval3;// auxiliary variable for replication rate, calculated from the genome, see function Regulation0
  double fval4;//antibiotic production rate
  double fval5; //unused
  char seq[1024]; // genome, a string composed of many F, A, B, H 
  char str[1024]; //unused
  int crow; //unused
  int ccol; //unused
  int valarray[1024]; // corresponding to the positions where there is an antibiotic gene in seq, there is a bitstring (represented as integer) here that specifies the antiobiotic type
  
  int valarray2[1024]; //unused
  int len_valarray2; //unused
  int ghash; //unused
//  REPLICATOR replicator[256]; // //unused
} TYPE2;

ANTIBIOTICS PLANE:

typedef struct __type2{
  int val; // used but for nothing really
  int val2; // nr of antibiotics at this location
  int val3; //unused
  int val4; //unused
  int val5; //unused
  double fval; //unused
  double fval2; //unused
  double fval3; //unused
  double fval4; //unused
  double fval5; //unused
  char seq[1024]; //unused
  char str[1024]; //unused
  int crow; //unused
  int ccol; //unused
  
  int valarray[1024]; //antibiotic types in the plane antib
  
  int valarray2[1024];  //unused
  int len_valarray2; //unused
  int ghash; //unused
//  REPLICATOR replicator[256]; //unused
} TYPE2;



*/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <grace_np.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <dirent.h>
#include "cash2003.h"
#include "cash2.h"
#include "mersenne.h"
#include "cash2-s.h"

#define MAXSIZE 1024 // same as cash2.h, max size of genome

struct point *ab_poslist; //2D array used to efficiently randomize antibiotics placement
int len_ab_poslist; //size of ab_poslist
void InitialiseABPosList(struct point **p_ab_poslist, int *p_len_ab_poslist, int MAXRADIUS); //initializes ab_poslist
void InitialiseFromInput(const char* par_fileinput_name,TYPE2 **world,TYPE2 **antib);
void InitialiseFromSingleGenome(const char* init_genome, char* init_ab_gen, TYPE2 **world,TYPE2 **antib);
void InitialiseFromScratch(TYPE2 **world,TYPE2 **antib);

//// Biological function declarations
void Mutate(TYPE2 **world, int row, int col);      // Mutate the genome at this position
void UpdateABproduction(int row, int col);        // Places new antibiotics in the field
void ChangeSeasonMix(TYPE2 **world);                 // "Sporulates" bacteria, shuffle spore position and restarts the field
void ChangeSeasonNoMix(TYPE2 **world);                 // "Sporulates" bacteria and restarts the field, NO spore shuffling

//Printing, plotting, painting.
void ColourPlanes(TYPE2 **world, TYPE2 **G, TYPE2 **A, TYPE2 **R); 
int Genome2genenumber(const char *seq, char gene);  // Convert genome to int of number of genes equal to char gene
int PrintPopFull(TYPE2 **world,TYPE2 **antib);
int tomovie_nrow;
int tomovie_ncol;
int ToMovie(TYPE2 **world, TYPE2 **antib, TYPE2** G, TYPE2** A, TYPE2** R);
void Pause();
void Exit(int exit_number); //just exits, but prints the program and arguments while doing so - why haven't I though of this before?

void BreakPoint_Recombination_LeftToRight_SemiHomog(TYPE2 *icel);
void BreakPointDeletion_RightToLeft(TYPE2 *icel);
void BreakPointDeletion_LeftToRight(TYPE2* icel);
void ScrambleGenomeBtwnSeasons(TYPE2* icel);
//takes care of regulation, and also set val3 val4 fval3 and fval4, 
// It is accessed by function pointers in the code - this is going to be fun
void Regulation0(TYPE2 *icel); //Old version of regulation - this works well

// Static global data structures and parameters
static TYPE2** world; // the world where bacteria live
static TYPE2** antib; // the antibiotic plane

static TYPE2** G; //auxiliary planes for visualising some features of bacteria, this colors the number of growth genes F
static TYPE2** A; //nr. of antibiotic genes 
static TYPE2** R; //nr. of break points

char par_movie_directory_name[MAXSIZE]="movie_strepto"; //genome alphabet
char par_fileoutput_name[MAXSIZE] = "data_strepto.txt";
char par_name[MAXSIZE] = "\0"; //name - it will create a par_fileoutput_name: data_[name].txt and a par_movie_directory_name: movie_[name];
int par_movie_period = 20;
int par_outputdata_period = 100;
char init_genome[MAXSIZE]="\0"; // initial genome, for specific experiments
char init_ab_gen[MAXSIZE]="\0"; // initial antibiotic bitstring,
int initialise_from_singlegenome=0;
char par_fileinput_name[MAXSIZE] = "\0"; //reads a one timestep snip of data_strepto.txt file
int initialise_from_input=0;
int antib_with_bitstring=1; // 1: Antibiotic with bistring; 0: antib without bitstrings - keep this to 1
int init_genome_size = 20; // initial genome size

static char* AZ="HFAB"; //genome alphabet
int MAXRADIUS = 10; //max distance from bacterium at which antibiotics are placed
double p_movement = 0.01; // probability of one step of random walk (if empty space available) per bacterium per unit time
int par_season_duration=2500; // nr. of time steps for one growth cycle

double ddrate=0.001; //per-gene duplication and deletion probability
double prob_new_brpoint = 0.01; //inflow of one new randomly placed breakpoints, per replication
double breakprob=0.01;//0.005; // probability of activating a break point
double prob_mut_antibtype_perbit = 0.05; //per bit probability of antibiotic type mutation

double spore_fraction=0.001; // fraction of nrow*ncol that sporulates -- i.e. probability that an individual is sampled to become a spore

int antib_bitstring_length = 16; // length of the antibiotic bitstring
double max_ab_prod_per_unit_time = -1.; // set either by command line argument, or as 1/len_ab_pos
double beta_antib_tradeoff = 1.; // tradeoff strength parameter, higher means stronger tradeoff -> more easy division of labor
int par_all_vs_all_competition = 1; // only set if we are not using bitstrings
double prob_noABspores_fromouterspace = 0.; //add spores from the outside of the grid
double tmp_prob_noABspores_fromouterspace=-1.; //auxiliart variable
int burn_in=0; // initial number of seasons where don't collect data - relaxation time for ad-hoc initial conditions
int par_burn_in_time=10; //this is going to be multiplied to season length

int mix_between_seasons = 0; //are spores mixed between seasons 0 = no, 1 = yes.
char breakpoint_mut_type = 'C'; // S: semi homolog recombination, T: telomeric deletion, c: centromeric towards telomeric (strepto like)
double par_beta_birthrate=0.3; //scaling factor for the birthrate as a function of the cumulative distance between antib and resistance
int const_tot_ab_mut=1;                 // if 1, the per AB gene mut rate is set - if 0 the per antibiotic gene bit mutrate is set: 
double prob_mut_antibtype_tot = 0.005;    // prob_mut_antibtype_tot is used instead of prob_mut_antibtype_perbit
                                        // to be precise: prob_mut_antibtype_perbit is set to prob_mut_antibtype_tot/antib_bitstring_length
int nr_H_genes_to_stay_alive=0;         // in experiments with Homeostatic genes - this is the min number of H genes to stay alive  
double h_growth=10.; // nr. of growth genes for half max growth rate
double h_antib_act=3.; // nr. of antibiotic genes for half max ab production rate

double constABprod = 0.; //constitutive antibiotic production, default = 0.
double scaling_factor_max_ab_prod_per_unit_time = 1.; // make this smaller to reduce AB prod, DO NOT MAKE THIS LARGER THAN 1. (default = 1.)
double max_repl_prob_per_unit_time=0.1; //scales the probability of replication per unit time

int which_regulation=0; //flag to tell the program which regulation function
void (*pt_Regulation)(TYPE2 *icel); // pointer to the regulation function (in case you want different ones)

int scramble_genome_btwn_seasons = 0; // if 1 it scrambles the gene order on the chromosome at the beginning of each season
int perfectmix=0; // if 1, mixes the world plane every time step

TYPE2 TYPE2_empty = { 0,0,0,0,0,0.,0.,0.,0.,0.,"\0","\0",0,0,{0} }; // zero the all values in a TYPE2 variable

void Initial(void)
{
	// readout parameters
	int myseed = 5431;//time(NULL);
	int myncol = 300; //nr. of rows in the field
	int mynrow = 300; //nr. of columns in the field
  char* readOut;
  display=0;
  MaxTime = 2147483647; // how long does the simulation lasts - you'll probably want to interrupt it by hand
  nrow = mynrow; /* # of row (default=100)*/
  ncol = myncol; /* # of column (default=100)*/
  init_genome[0]='\0';
  
  // COMMAND LINE ARGUMENTS
	for(int i = 1; i < (int)argc_g; i++)
	{
	  readOut = (char*)argv_g[i];
		if(strcmp(readOut, "-seed") == 0) myseed = atoi(argv_g[i+1]);
		else if(strcmp(readOut, "-c") == 0) ncol = atoi(argv_g[i+1]);
		else if(strcmp(readOut, "-r") == 0) nrow = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-moviedir") == 0) strcpy( par_movie_directory_name , argv_g[i+1] );
    else if(strcmp(readOut, "-datafile") == 0) strcpy( par_fileoutput_name , argv_g[i+1] );
    else if(strcmp(readOut, "-name") == 0) strcpy( par_name , argv_g[i+1] ); // gives a common name to all output files, which are then distinguished by extension
    else if(strcmp(readOut, "-movie_period") == 0) par_movie_period = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-data_period") == 0) par_outputdata_period = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-display") == 0) display = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-season") == 0) par_season_duration = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-breakpoint_inflow") == 0) prob_new_brpoint = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-spore_fraction") == 0) spore_fraction = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-maxtime") == 0) MaxTime = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-init_genome") == 0) {strcpy( init_genome , argv_g[i+1] ); 
                                                   strcpy( init_ab_gen , argv_g[i+2] ); i++; 
                                                   init_genome_size=strlen(init_genome);}
    else if(strcmp(readOut, "-max_ab_prod_per_unit_time") == 0) max_ab_prod_per_unit_time = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-beta_antib_tradeoff") == 0) beta_antib_tradeoff = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-antib_with_bitstring") == 0) antib_with_bitstring = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-antib_bitstring_length") == 0) antib_bitstring_length = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-prob_mut_antibtype_perbit") == 0) prob_mut_antibtype_perbit = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-par_all_vs_all_competition") == 0) par_all_vs_all_competition = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-prob_noABspores_fromouterspace") == 0) tmp_prob_noABspores_fromouterspace = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-mix_between_seasons") == 0) mix_between_seasons = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-breakpoint_mut_type") == 0) breakpoint_mut_type = *(argv_g[i+1]);
    else if(strcmp(readOut, "-beta_birthrate") == 0) par_beta_birthrate = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-const_tot_ab_mut") == 0) const_tot_ab_mut = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-prob_mut_antibtype_tot") == 0) prob_mut_antibtype_tot = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-nr_Hgenes_to_stay_alive") == 0) nr_H_genes_to_stay_alive = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-which_regulation") == 0) which_regulation = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-input") == 0) strcpy( par_fileinput_name , argv_g[i+1] );
    else if(strcmp(readOut, "-ddrate") == 0) ddrate = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-scramble_genome_btwn_seasons") == 0) scramble_genome_btwn_seasons = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-perfectmix") == 0) perfectmix = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-breakprob") == 0) breakprob = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-constABprod") == 0) constABprod = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-scaling_factor_max_ab_prod_per_unit_time") == 0) scaling_factor_max_ab_prod_per_unit_time = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-h_growth") == 0) h_growth= atof(argv_g[i+1]);
    else {fprintf(stderr,"Parameter number %d was not recognized, simulation not starting\n",i);
          fprintf(stderr,"It might help that parameter number %d was %s\n",i-1, argv_g[i-1]);
          Exit(1);}
    i++;
	}
  //makes output file
  if(strlen(par_name)!=0){
    strcpy(par_fileoutput_name,"data_");
    strcat(par_fileoutput_name,par_name);
    strcat(par_fileoutput_name,".txt");
    strcpy(par_movie_directory_name,"movie_");
    strcat(par_movie_directory_name,par_name);
    fprintf(stderr,"Output file name: %s\n Output movie dir name: %s\n",par_fileoutput_name,par_movie_directory_name);
  }
  //check if par_movie_directory_name and par_fileoutput_name already exist,
  // simulation not starting if that is the case
  DIR *dir = opendir(par_movie_directory_name);
  if(dir){ fprintf(stderr, "Initial(): Error. Directory %s already exists, simulation not starting\n",par_movie_directory_name); Exit(1);}
  FILE *fp = fopen(par_fileoutput_name,"r"); 
  if(fp){ fprintf(stderr, "Initial(): Error. File %s already exists, simulation not starting\n",par_fileoutput_name); Exit(1);}
  //File for input - if any
  if(strlen(par_fileinput_name)){
    FILE *fp = fopen(par_fileinput_name,"r"); 
    if(!fp){ fprintf(stderr, "Initial(): Error. Input file %s does not exist, simulation not starting\n",par_fileoutput_name); Exit(1);}
    else {initialise_from_input=1;fclose(fp);}
  }

  nplane = 5; /* # of planes (default=0)*/
  scale = 1; /* size of the window (default=2)*/
  margin=10;
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		                  Margolus diffusion is not supported for ECHO. */

  ulseedG =(myseed>0)?myseed:time(NULL);// time(NULL); /* random seed ... if you don't know the time */
  fprintf(stderr,"Seeding with: %ld\n", ulseedG);
  
  // Sets mutation rate
  if(const_tot_ab_mut){
    prob_mut_antibtype_perbit = prob_mut_antibtype_tot/(double)(antib_bitstring_length);
    fprintf(stderr,"Constant per AB mut rate used, prob_mut_antibtype_perbit reset to %f\n", prob_mut_antibtype_perbit);
  }
  
  // checks if simulation will be started from single genome
  if( strlen(init_genome) ){
    if(initialise_from_input){
      fprintf(stderr, "Initial(): Error. You cannot specify input file AND initial genome, simulation not starting\n"); Exit(1);
    }
    initialise_from_singlegenome=1; //this will start a sim with one genome in the center of the grid
  }
  
  if(!antib_with_bitstring){
    antib_bitstring_length=1; // one antibiotic
    prob_mut_antibtype_perbit = 0.; //which does not mutate, but should be re-assigned at the beg. of each season
    fprintf(stderr,"No bitstrings for antibiotic, all vs all competition set to: %d\n",par_all_vs_all_competition);
  }
  if(tmp_prob_noABspores_fromouterspace>0.){
    burn_in=1;
    fprintf(stderr,"Probability of 'noAB' spores from outer space set to: %f\n",tmp_prob_noABspores_fromouterspace);
    par_burn_in_time *= par_season_duration;
    fprintf(stderr,"par_burn_in_time set to: %d\n",par_burn_in_time);
  }
  // sets regulation pointer
  if(which_regulation==0) pt_Regulation = &Regulation0;
  else{
    fprintf(stderr,"Initial(): Error. which_regulation got unrecognised value: %d\n",which_regulation);
    Exit(1);
  }
  if(perfectmix){
    p_movement=0.;
    fprintf(stderr,"perfectmix = %d, mixing every step (also setting p_movement = 0)\n",perfectmix);
  }
  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.};
  
  tomovie_nrow=nrow;
  tomovie_ncol=nplane*ncol;
  OpenPNG(par_movie_directory_name, tomovie_nrow, tomovie_ncol); //recall that for some weird reason you must allocate a data structure
                                                                 // of size tomovie_nrow+1 and tomovie_ncol+1, if OpenPNG gets arguments
                                                                 // tomovie_nrow and tomovie_ncol 
}

void InitialPlane(void){
  MakePlane(&world,&antib,&G,&A,&R); //Allocates the planes

  InitialiseABPosList(&ab_poslist, &len_ab_poslist, MAXRADIUS); //initialises an auxiliary variable for antibiotic positions in the ab plane
  if(max_ab_prod_per_unit_time<0.) max_ab_prod_per_unit_time = 1/(double)len_ab_poslist; //scales the antibiotic production probability
  
  // Initialisation of the grid with a bunch of cells - or one, depending on the flags set in Initial()
  if(initialise_from_input) InitialiseFromInput(par_fileinput_name,world,antib);
  else if(initialise_from_singlegenome) InitialiseFromSingleGenome(init_genome, init_ab_gen, world, antib);
  else InitialiseFromScratch(world,antib);
  
  fprintf(stderr,"\n\nworld is ready. Let's go!\n\n");
  Boundaries2(world);
  
  //gradient takes c values from 100 to 255
  // colors cyan, blue, magenta, red i.e. 155/4 ~ 35 values per segment, i.e. 255/35 = 7
  for(int c=1;c<255;c++)
  {
    if(c<100) continue;
    else{
      if(c<135) ColorRGB(c, 0, 255 - 7*(c-100) , 255);
      else if(c<170) ColorRGB(c, 7*(c-135) , 0, 255);
      else if(c<205) ColorRGB(c, 255 , 0, 255- 7*(c-170));
      else ColorRGB(c, 255 , 7*(c-205),0);
    }
  }
}

// Calculates the Hamming distance between two bitstrings - represented as integers
int HammingDistance(int a, int b){
  int h=0;
  int xor = a ^ b;
  while (xor>0){
    h += xor & 1;
    xor >>= 1;
  }
  return h;
}

//Birthrate() checks, for every antib at the pos of an individual, the antib with the smallest Hamm Dist
// then sums all these hamm dists into a final score cumhammdist
// end returns the birthrate - a decreasing funct of cumhammdist
double BirthRate(TYPE2 *icel, TYPE2 *ab)
{
  int i, k, h;
  int hammdist, cumhammdist=0;
  float birthrate;
  if(ab->val2){
    for (i=0; i<ab->val2; i++){
      //check for each ab in the field whether individual has some resistance.
      hammdist=999;
      h=0;
      for(k=0; icel->seq[k]!='\0'; k++){
        if(icel->seq[k]=='A'){
          h = HammingDistance( icel->valarray[k] , ab->valarray[i] );  //finds HD
          if(hammdist>h) hammdist=h; //check if h is the smallest
        }
      }
      cumhammdist+=hammdist; // add as cumulative distance that particular distance
    }
  }
  if(!antib_with_bitstring) birthrate = (cumhammdist==0)?1.:0.; 
  else birthrate = exp(-par_beta_birthrate*cumhammdist*cumhammdist); 
  return birthrate;
}

void NextState(int row,int col)
{
  //int sum1,sum0;
  int k;
  int counter=0;

  TYPE2 *nei;
  // return;
  if(world[row][col].val2==0){
    int howmany_emptynei = CountMoore8(world,0,row,col); //counts how many neigh are empty
    if( howmany_emptynei==8) return; // if the neighbourhood is empty no replication is going to happen
    
    int dirarray[8];
    double grarray[8];
        
    double totgrowth=0;
    
    //we find who these neighbors are:
    for(k=1;k<9;k++){
      nei = GetNeighborP(world,row,col,k);
      
      //check if neighbour cell is inhibited 
      if(nei->val!=0){
                        
        double repprob= ( nei->val5 < nr_H_genes_to_stay_alive )? 0.: BirthRate(nei, &antib[row][col]); //calculates replication probability
        //if cell has no fitness genes in genome it cannot reproduce
        if(nei->val3==0 || repprob<=0.000000000001) continue;
        //save direction
        dirarray[counter]=k;

        double growth = repprob * nei->fval3; // fval3 contains the previously calculated prob of replication
        
        grarray[counter]=growth; //save its replication rate
        totgrowth+=growth;
        counter++; //counts how many neighbourhs we have
      }
      
    }
    //we filled up the neigh array, so the one with higher repl rate has higher chances of replicating
    if(counter){
      double rn = genrand_real2()*8.;
      int counter2=0;
      double cumprob=0.0;
      if(rn < totgrowth){
        cumprob += grarray[counter2];
        while(cumprob < rn ){
          counter2++;
          cumprob +=grarray[counter2];
        }
        //now we replicate the individual at direction dirarray[counter] into the focal pixel
        int direction=dirarray[counter2];

        TYPE2 winner = GetNeighborS(world,row,col,direction);
        world[row][col] = winner;
        
        Mutate(world,row,col);
        
        // UPON BIRTH WE SET A BUNCH OF PARAMETERS: val3 val4 val5 fval3 fval4
        (*pt_Regulation)(&world[row][col]);
        
      }
    }
    
  }
  // movement + continuous antib + death if you are on it prod maybe helps?
  else if(world[row][col].val2!=0){
    
    // DEATH
    double deathprob;
    deathprob = ( world[row][col].val5 < nr_H_genes_to_stay_alive )? 1. : 1.-BirthRate(&world[row][col], &antib[row][col]); //no H genes you die
    // deathprob = 1. - BirthRate( &world[row][col] , &antib[row][col] ) ; // in this version nr_H_genes_to_stay_alive controls only birth rate
                                                                           // this is the case if H are as biosynth genes but do not cause death when missing
    if(genrand_real2()<deathprob){
      world[row][col]=TYPE2_empty; // kill - i.e. we remove the bacterium
    }else{
      //MOVEMENT
      UpdateABproduction(row, col);
      //now we will try to move
      if(p_movement>0.){
        int neidir=1+(int)(8.*genrand_real2()); //picks random direction
        int neirow, neicol;
        GetNeighborC(world, row, col, neidir, &neirow, &neicol); //get the neighbour in that direction
        nei=&world[neirow][neicol];
        if(nei->val2==0){ //we can only move if neighbour is empty
          if(genrand_real2() < p_movement){
            *nei = world[row][col]; // move - i.e. copy content
            world[row][col]=TYPE2_empty; //don't forget to delete previous location
          }
        }  
      }
    } 
  }
}

void Update(void){
  static int how_long_no_antib=0;
  int n_antib;
  if(perfectmix) PerfectMix(world);
  
  Asynchronous(); 
  
  if(!burn_in){
    if(Time%par_movie_period==0){
      ColourPlanes(world,G,A,R); // Pretty big speed cost here, can be fixed later
      if(display) Display(world,antib,G,A,R);
      else ToMovie(world,antib,G,A,R);
    }
    if(Time%par_outputdata_period==0){
      if(display){
        fprintf(stderr,"Update(): Error. display = 1 is deprecated - just don't use me\nBye\n");
        Exit(1);
      }else{
        n_antib=PrintPopFull(world,antib);
        if(n_antib==0){
          how_long_no_antib++;
          if(how_long_no_antib==20){
            fprintf(stderr,"No antibiotic genes in the field for 20 saved data points, the simulation will terminate\n");
            Exit(1);
          }
        }
      }
    }
  }else{ // else if burn_in == 1
    if(Time>=par_burn_in_time){
      fprintf(stderr, "Update(): message. End of burn-in\n");
      burn_in=0;
      Time=0;
      prob_noABspores_fromouterspace = tmp_prob_noABspores_fromouterspace;   /////// OK BUT DOES THIS WORK ALSO WHEN WE DON'T WANT BURN IN?
    }
  }

  if(Time%par_season_duration==0 && Time>0){
    if(mix_between_seasons) ChangeSeasonMix(world);
    else ChangeSeasonNoMix(world);
  }

  //while( Mouse()==0) {}; // you can run the program continuously by commenting out this statement

}



////////////////////////
////////////////////////
////////////////////////
////////////////////////
////////////////////////
//BELOW IS OUR OWN STUFF//
////////////////////////
////////////////////////
////////////////////////
////////////////////////
////////////////////////
void Pause()
{
  printf("Paused. Where's the any key?\n");
  getchar();
}


int Genome2genenumber(const char *seq, char gene)
{
  int i,genecount=0;
  int seqlen=strlen(seq);
  for(i=0;i<seqlen;i++) if(seq[i]==gene) genecount++;
  return genecount;
}

void ScrambleGenomeBtwnSeasons(TYPE2* icel){
  int genlen = strlen(icel->seq);
  int i;
  printf("Got genome: %s\n",icel->seq);
  for(i=0;i<genlen;i++) if(icel->seq[i] == 'A') printf("%d,",icel->valarray[i]);
  printf("\n");
  for(i=genlen-1;i>=0;i--){
    int rpos = (i+1)*genrand_real2(); // 0 <= rpos <= i, with  0 <= i <= genlen-1
    char tmp_gene = icel->seq[rpos];
    int tmp_ab = icel->valarray[rpos];
    
    icel->seq[rpos] = icel->seq[i];
    icel->valarray[rpos] = icel->valarray[i];

    icel->seq[i] = tmp_gene;
    icel->valarray[i] = tmp_ab;
  }
  printf("--Scrabled: %s\n",icel->seq);
  for(i=0;i<genlen;i++) if(icel->seq[i] == 'A') printf("%d,",icel->valarray[i]);
  printf("\n");
}

//function to reseed plane with "spores", mixing the field in the process
void ChangeSeasonMix(TYPE2 **world)
{
  TYPE2 spores[MAXSIZE];
  int sporenr=0;
  int i,j,k,row,col;
  
  //attempt 2, pick 100 spores at random
  sporenr=spore_fraction*nrow*ncol;
  if(sporenr==0) {fprintf(stderr,"ChangeSeasonMix(): Error.No spores for next generation?\n"); Exit(1);}
  else if(sporenr>MAXSIZE) sporenr=MAXSIZE-1;
  int actual_sporenr=0;
  for(i=0;i<sporenr;i++){
    int attempt=0;
    while(attempt<20){
      row=1+(nrow*genrand_real2());
      col=1+(ncol*genrand_real2());
      if(world[row][col].val2){
        spores[actual_sporenr]=world[row][col];
        //if(1+i==250) printf("make spore 250, val2 was %d, genome %s\n",spores[i].val2,spores[i].seq );
        //for(k=0;k<MAXSIZE;k++) spores[i].seq[k]=world[row][col].seq[k];
        world[row][col].val = 0;
        world[row][col].val2 = 0;
        actual_sporenr++;
        break;
      }
      else attempt++;
    
    }
  }
  if(actual_sporenr==0){fprintf(stderr,"ChangeSeasonMix(): No spores were found for next generation? System extinct\n"); Exit(1);}
  //erase the plane
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
  {
    //empty the plane
    world[i][j].val = 0;
    world[i][j].val2 = 0;
    antib[i][j].val=0;
    antib[i][j].val2=0;
    for(k=0;k<MAXSIZE;k++) {
      antib[i][j].valarray[k]=0;
      world[i][j].seq[k]='\0';
      world[i][j].valarray[k]=-1;
    }
  }

  
  int ipos, jpos, zpos, attempt;
  //place spores
  zpos=(int)(genrand_real2()*nrow*ncol);
  ipos=zpos/ncol+1;
  jpos=zpos%ncol+1;
  for(i=0; i<actual_sporenr; i++){
    attempt=0;
    while (world[ipos][jpos].val && attempt<100){
      zpos=(int)(genrand_real2()*nrow*ncol);
      ipos=zpos/ncol+1;
      jpos=zpos%ncol+1;
      attempt++;
    }
    if(attempt==100){
      fprintf(stderr,"ChangeSeasonMix(): Error. Cannot place new spore %d\n",i);
      Exit(1);
    }else{
      world[ipos][jpos]=spores[i];
      if(scramble_genome_btwn_seasons) ScrambleGenomeBtwnSeasons(&world[ipos][jpos]); //scrambles genome and AB genome accordingly, does not alter anything other than order of genes
      //if(1+i==250) printf("val2 was %d, genome %s\n",spores[i].val2,spores[i].seq );
      if(genrand_real2()<prob_noABspores_fromouterspace){
        //instead of placing our own spores, we place an outsider, with no AB prod
        int howmanyFfromouterspace=100;
        for(k=0;k<howmanyFfromouterspace;k++) {
          world[ipos][jpos].seq[k]='F';
          world[ipos][jpos].valarray[k]=-1;
        }
        world[ipos][jpos].seq[howmanyFfromouterspace]='\0';
      }
      
      if(!antib_with_bitstring){ 
        if(par_all_vs_all_competition){
          for(k=0;k<MAXSIZE;k++) if(world[ipos][jpos].seq[k]=='A')world[ipos][jpos].valarray[k] = i;
        }
      }
      world[ipos][jpos].val=1+i%10;
      world[ipos][jpos].val2=1+i;
      
    }
  }
}

//function to reseed plane with "spores"
// without reshiffling their positions between seasons
void ChangeSeasonNoMix(TYPE2 **world)
{
  int i,j,k,bactnumber=0;
  
  //Go through the plane
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    // delete the antib plane
    antib[i][j].val=0;
    antib[i][j].val2=0;
    for(k=0;k<MAXSIZE;k++){
      antib[i][j].valarray[k]=0;
    }
    //empty the plane - with small prob, we do not wipe this guys out, because it sporulates
    // in other words, with high prob we wipe people out
    if(genrand_real2() < (1. - spore_fraction) ){
      world[i][j].val = 0;
      world[i][j].val2 = 0;
      for(k=0;k<MAXSIZE;k++) {
        world[i][j].seq[k]='\0';
        world[i][j].valarray[k]=-1;
      }
    }
    else{
      bactnumber++;
      //with small prob we did not kill the guy, then: 
      
      //instead of placing our own spores, we place an outsider, with no AB prod
      if(genrand_real2()<prob_noABspores_fromouterspace){
        int howmanyFfromouterspace =100;
        for(k=0;k<howmanyFfromouterspace;k++){
          world[i][j].seq[k]='F';
          world[i][j].valarray[k]=-1;
        }
        world[i][j].seq[howmanyFfromouterspace]='\0';
      }
      
      if(scramble_genome_btwn_seasons) ScrambleGenomeBtwnSeasons(&world[i][j]); //scrambles genome and AB genome accordingly, does not alter anything other than order of genes

      if(!antib_with_bitstring){
        if(par_all_vs_all_competition){
          for(k=0;k<MAXSIZE;k++) if(world[i][j].seq[k]=='A')world[i][j].valarray[k] = bactnumber;
        }
      }
      world[i][j].val = 1+bactnumber%10;
      world[i][j].val2 =1+bactnumber;
      
    }
  }

}



// What kind of mutations are possible?
// dupl, del, translocations and breaking mut. A-> a, G->g, R->r
// mutations are more frequet at the end of the chrom
// according to function: 1 - exp(-0.0001*genome_size^2)
// because this goes to zero in zero -> epsilon + (1-epsilon)(1 - exp(-0.0001*genome_size^2))
void Mutate(TYPE2** world, int row, int col)
{
  TYPE2 *icell;
  icell=&world[row][col];
  int i,j,genome_size;//, which_mutations;
  
  genome_size = strlen(icell->seq);
  
  // ANTIBIOTIC TYPE MUTATIONS
  if(antib_with_bitstring){
    for(int i=0;icell->seq[i]!='\0';i++){
      if(icell->seq[i]=='A'){
        if(genrand_real2() < prob_mut_antibtype_perbit) 
          icell->valarray[i]=icell->valarray[i] ^ (1<<(int)(genrand_real2()*antib_bitstring_length)); //this flips one random bit between 0 and 12( ^ is xor mask)
      }
    }
  }

  // DUPLICATION AND DELETIONS
  // fills up an array that tells us which genes we want to duplicate and delete
  int nrmuts=bnldev(2*ddrate,genome_size); //binomial sample of genes that will be duplicated or deleted
  int mutposarr[MAXSIZE];
  int imut;
  for(imut=0;imut<genome_size;imut++) mutposarr[imut]=imut;
  for(imut=0;imut<nrmuts;imut++) {
    int rpos = imut+ genrand_real2() * (genome_size - imut);
    int tmp = mutposarr[imut];
    mutposarr[imut] = mutposarr[rpos];
    mutposarr[rpos] = tmp;
  }
  for(imut=0; imut<nrmuts; imut++){
    int mutpos = mutposarr[imut];
    if(genrand_real2() < 0.5){ 
      // Deletions
      for(i=mutpos;i<genome_size;i++){ 
        icell->seq[i]=icell->seq[i+1];
        icell->valarray[i]=icell->valarray[i+1];
      }
      genome_size--;
      icell->seq[genome_size]='\0';
      icell->valarray[genome_size]=-1;
      //update mutpos array, decrease positions that are larger than mutpos
      for(j=imut;j<nrmuts;j++) if(mutposarr[j]>mutpos) mutposarr[j]--;
    }else{                             
      //Duplication
      if(genome_size>=MAXSIZE-2) continue;
      int duppos=genrand_real2()*genome_size;
      char insertgene=icell->seq[mutpos];
      int insertval = icell->valarray[mutpos];
      for(i=genome_size;i>duppos;i--){
        icell->seq[i] = icell->seq[i-1];
        icell->valarray[i]=icell->valarray[i-1];
      }
      icell->seq[duppos]=insertgene;
      icell->valarray[duppos]=insertval; 
      
      genome_size++;
      icell->seq[genome_size]='\0';
      for(j=imut;j<nrmuts;j++) if(mutposarr[j]>duppos) mutposarr[j]++;
    }
  }
  
  //goes from left to right, breaks are recombination mediated
  if(breakpoint_mut_type=='S') BreakPoint_Recombination_LeftToRight_SemiHomog(icell); //how it was  <---previous attempt that doesn't work super well
  else if(breakpoint_mut_type=='T') BreakPointDeletion_RightToLeft(icell); //no recomb, only 3'->5' instability  <---previous attempts that doesn't work super well
  else if(breakpoint_mut_type=='C') BreakPointDeletion_LeftToRight(icell);  //5'->3' instability    <--------  Default, used in all simulations
  else{ 
    fprintf(stderr,"Mutate(): Error. Unrecognised option for the type of breakpoint mutation\n");
    Exit(1);
  }
  
  //Random break point insertion
  if(genome_size<MAXSIZE-2 && genrand_real2()< prob_new_brpoint){ 
    int new_brpos=genrand_real2()*genome_size;
    char insertgene='B';
    for(i=genome_size;i>new_brpos;i--){ 
      icell->seq[i] = icell->seq[i-1];
      icell->valarray[i]=icell->valarray[i-1];
    }
    icell->seq[new_brpos]=insertgene;
    icell->valarray[new_brpos]=-1;
    genome_size++;
    icell->seq[genome_size]='\0';
  }  
  
}

void BreakPoint_Recombination_LeftToRight_SemiHomog(TYPE2* icel){

  char *seq=icel->seq;
  //break points
  int breakarray[MAXSIZE];
  int breaknr=0;
  int genome_size = strlen(seq);
  for (int ipos=0; ipos<genome_size; ipos++){
    if(seq[ipos]=='B'){
      breakarray[breaknr]=ipos;
      breaknr++;
    }
  }
  int b, match;//,lcount;
  // notice that like this breaks happen more frequently closer to 5' than to 3'
  //Also, one break happens, period.
  if(breaknr>1){
    for(b=0; b<breaknr-1; b++){
      if(genrand_real2()<breakprob){
        match=b+genrand_real2()*(breaknr-b);
        for(int i=breakarray[match],lcount=0 ; i<genome_size ; i++, lcount++){
          seq[ breakarray[b]+lcount ] = seq[i];
          icel->valarray[ breakarray[b]+lcount ] = icel->valarray[i];
        }
        genome_size-= (breakarray[match] - breakarray[b]) ;
        seq[genome_size] = '\0';
        break;
      }
    }
  }
}


void BreakPointDeletion_RightToLeft(TYPE2* icel){
  char *seq=icel->seq;

  //break points
  int genome_size = strlen(seq);
  for (int ipos=genome_size-1; ipos>=0; ipos--){
    if(seq[ipos]=='B' && genrand_real2()<breakprob){
      for (int k=ipos; k<MAXSIZE; k++){
        icel->valarray[k]=-1;
      }
      seq[ipos]='\0';
      break;
    }
  }
}

void BreakPointDeletion_LeftToRight(TYPE2* icel){
  char *seq=icel->seq;

  //break points
  int genome_size = strlen(seq);
  for(int ipos=0; ipos<genome_size; ipos++){
    if(seq[ipos]=='B' && genrand_real2()<breakprob){
      for(int k=ipos; k<MAXSIZE; k++){
        icel->valarray[k]=-1;
      }
      seq[ipos]='\0';
      break;
    }
  }
}

//takes care of "regulation", i.e. it sets fval3 and fval4, accessed by function pointers in the code
void Regulation0(TYPE2 *icel){
  // SET GROWTH AND AB PRODUCTION PARAMETERS
  icel->val3=Genome2genenumber(icel->seq,'F');
  icel->val4=Genome2genenumber(icel->seq,'A');
  icel->val5=Genome2genenumber(icel->seq,'H');

  double fg = icel->val3; 
  double ag = icel->val4; 
  
  icel->fval3 = max_repl_prob_per_unit_time * fg/(fg+h_growth);
  
  double constABprod_if_ag;
  constABprod_if_ag = (ag>0)?constABprod:0.; //check that there are antibiotics - otherwise it might make them out of thin air?
  
  icel->fval4 = max_ab_prod_per_unit_time * scaling_factor_max_ab_prod_per_unit_time *constABprod_if_ag + (1.-constABprod_if_ag)*max_ab_prod_per_unit_time * scaling_factor_max_ab_prod_per_unit_time * ag/(ag+h_antib_act) * (exp(-beta_antib_tradeoff*fg));
} 

void UpdateABproduction(int row, int col){
  TYPE2 *icell=&world[row][col];

  int i,k;
  if( icell->val4==0 || icell->fval4<0.000000000001) return;//if you don't have antib genes, surely no ab are placed

  int howmany_pos_get_ab = bnldev(icell->fval4,len_ab_poslist);
  int which_ab[MAXSIZE];
  int nrtypes=0;

  //which ab can this individual produce?
  int genome_length = strlen(icell->seq);
  for (i=0; i<genome_length; i++){
    if(icell->seq[i]=='A'){
      which_ab[nrtypes]=icell->valarray[i];
      nrtypes++;
    }
  }

  //ab_poslist
  struct point tmp;
  int pos_ab_poslist;
  //get a number of random positions from ab_poslist equal to howmany_pos_get_ab
  for(i=0; i< howmany_pos_get_ab; i++){
    pos_ab_poslist= i+(int)((len_ab_poslist-i)*genrand_real2());
    tmp = ab_poslist[pos_ab_poslist];
    ab_poslist[pos_ab_poslist] = ab_poslist[i];
    ab_poslist[i] = tmp;
  }
  //now the first "howmany_pos_get_ab" positions in "ab_poslist" are random
  for(i=0;i<howmany_pos_get_ab;i++){
    int ii = row+ab_poslist[i].row;
    int jj = col+ab_poslist[i].col;
    if(ii<1) ii = nrow+ii;    //check for boundaries
    else if(ii>nrow) ii = ii%nrow;
    if(jj<1) jj = ncol+jj; 
    else if(jj>ncol) jj = jj%ncol;
    
    //pick a random ab to place
    int ab_toplace=which_ab[(int)(nrtypes*genrand_real2())];
    int found=0;
    for (k=0; k<antib[ii][jj].val2;k++){
      if(antib[ii][jj].valarray[k]== ab_toplace){  //icell->val2){
        found=1;
        break;
      }
    }
    if(!found){
      antib[ii][jj].valarray[antib[ii][jj].val2]=ab_toplace;   //icell->val2;
      antib[ii][jj].val2++;
      antib[ii][jj].val=antib[ii][jj].val2%10;
    }
  }
  
}


char Num2Char(int c)
{
	switch(c){
		case 0: return 'F';
		case 1: return 'A';
		case 2: return 'B';
		default: {
      fprintf(stderr, "Num2Char(): Error, character not recognised, got number: %d\n",c); 
      Exit(1);
      exit(1);
    }
	}
}


// print to file cell by cell in order, with format:
// [Time] [tot_nr_antib] [antib,array,if,antibs,else,0,] (if val2)[[val2] [genome] [abgenes,if,ab,genes]](else)[0 n 0,]
//returns tot_antib_genes
int PrintPopFull(TYPE2 **world,TYPE2 **antib){
  FILE *fp;
  fp= fopen( par_fileoutput_name, "a" );
  int tot_antib_genes=0; //if this number == 0 only fitness genes... kinda pointless to go on with the sim

  for(int i=1;i<=nrow;i++)for(int j=1;j<=ncol;j++){
    fprintf(fp, "%d %d ",Time, antib[i][j].val2 );
    if(antib[i][j].val2){
      for(int k=0; k<antib[i][j].val2; k++){
         fprintf(fp,"%d," , antib[i][j].valarray[k]);
      }
    }else fprintf(fp,"0,"); 

    if(world[i][j].val2){
      tot_antib_genes+=world[i][j].val4;
      fprintf(fp, " %d %s ",world[i][j].val2, world[i][j].seq);
      if(world[i][j].val4){
        for (int k=0; k<strlen(world[i][j].seq); k++){
          if(world[i][j].seq[k]=='A') fprintf(fp, "%d,",world[i][j].valarray[k]);
        }
      }else{
        fprintf(fp,","); // this is so that you know there is nothing (!=0 - which is a valid AB ), and hopefully every piece of python breaks on this? 
      }
    }else{
      fprintf(fp," 0 n 0,"); 
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  return tot_antib_genes;
}

// Assigs a color (an RGB tuple) to a number, which means that you can use it to map, e.g. G[row][col].val = nr_of_growth_genes, to a colour in a color palette
// the colour palette is defined in InitialPlane()
//I am pretty sure this is kinda buggy - but it should be easy to fix
void ColourPlanes(TYPE2 **world, TYPE2 **G, TYPE2 **A, TYPE2 **R)
{
	int i,j;
  int colordisplacement = 100;
	for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
		if(world[i][j].val==0)
    {
		  G[i][j].val = 0;
      A[i][j].val = 0;
      R[i][j].val = 0;
    }
    else
    {
      int Gcolour = Genome2genenumber(world[i][j].seq,'F')+colordisplacement;//+200;
      if(Gcolour==colordisplacement) G[i][j].val = 1;
      else if (Gcolour<255) G[i][j].val = Gcolour;
      else G[i][j].val = 255;

      int Acolour = Genome2genenumber(world[i][j].seq,'A')+colordisplacement;//+200;
      if(Acolour==colordisplacement) A[i][j].val = 1;
      else if (Acolour<255) A[i][j].val = Acolour;
      else A[i][j].val = 255;

      int Rcolour = Genome2genenumber(world[i][j].seq,'B')+colordisplacement;//+200;
      if(Rcolour==colordisplacement) R[i][j].val = 1;
      else if (Rcolour<255) R[i][j].val = Rcolour;
      else R[i][j].val = 255;
    }
}

//Prepares the arrays for making pngs
int ToMovie(TYPE2 **world, TYPE2 **antib, TYPE2** G, TYPE2** A, TYPE2** R)
{
  int popcount=0;
  
  int i,j;
  int **tomovie;
  
  tomovie=(int **)malloc( (tomovie_nrow+1)*sizeof(int*));
  for(i=0; i<=nrow;i++) tomovie[i]=(int *)malloc( (tomovie_ncol+1)*sizeof(int));
    
  for(i=0;i<= tomovie_nrow ;i++)for(j=0;j<=tomovie_ncol;j++){
    tomovie[i][j]=0;
  }
  
  //here you'll save data
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++){
    if(world[i][j].val!=0) popcount++;
    tomovie[i-1][j-1 +0*(ncol)] = world[i][j].val;
    tomovie[i-1][j-1 +1*(ncol)] = antib[i][j].val;
    tomovie[i-1][j-1 +2*(ncol)] = G[i][j].val;
    tomovie[i-1][j-1 +3*(ncol)] = A[i][j].val;
    tomovie[i-1][j-1 +4*(ncol)] = R[i][j].val;
  }
  
  PlanePNG(tomovie,0);

  for(i=0; i<tomovie_nrow;i++) free(tomovie[i]);
  free(tomovie);
  
  return popcount;
}

void InitialiseABPosList(struct point **p_ab_poslist, int *p_len_ab_poslist, int MAXRADIUS){
  int i,j, howmany=0;
  for (i=-MAXRADIUS; i<MAXRADIUS+1; i++){
    for (j=-MAXRADIUS; j<MAXRADIUS+1; j++){
      if( sqrt((double)( (i*i)+(j*j) ))<=MAXRADIUS){
        howmany++; 
      }
    }
  }
  *p_len_ab_poslist = howmany; //set the len_ab_poslist variable
  *p_ab_poslist = malloc(howmany * sizeof(struct point)); //alllllocate the array
  int array_pos=0;
  //fill arrrrray
  for (i=-MAXRADIUS; i<MAXRADIUS+1; i++){
    for (j=-MAXRADIUS; j<MAXRADIUS+1; j++){
      if( sqrt((double)( (i*i)+(j*j) ))<=MAXRADIUS){
        (*p_ab_poslist)[array_pos].row = i;
        (*p_ab_poslist)[array_pos].col = j;
        array_pos++;
      }
    }
  }
  
  return;
}

void InitialiseFromInput(const char* par_fileinput_name, TYPE2 **world,TYPE2 **bact){
  
  FILE *fp = fopen(par_fileinput_name,"r"); 
  fprintf(stderr, "InitialiseFromInput(): Warning, for now only reads genomes, not AB already in the plane\n");
  fprintf(stderr, "InitialiseFromInput(): Warning, changing season within this function\n");
  
  char fline[5241], *token, *token2,fab_string[10241];
  char fline_cp[5241];
  const char sep[2]=" ";
  const char sepab[2]=",";
  int i,j,k;
  int sig_break = 0;
  int line_nr=0;
  for(i=1;i<=nrow;i++){
    for(j=1;j<=ncol;j++){
      world[i][j].val = 0;world[i][j].val2 = 0;
      antib[i][j].val = 0;antib[i][j].val2 = 0;
      for(k=0;k<MAXSIZE;k++) antib[i][j].valarray[k]=-1;
      for(k=0;k<MAXSIZE;k++) world[i][j].seq[k]='\0';
      for(k=0;k<MAXSIZE;k++) world[i][j].valarray[k]='\0';
      line_nr++;
      // fprintf(stderr,"Reading next line:\n");
      if(  fgets(fline,5240,fp) != NULL ){
        strcpy(fline_cp,fline); //for errors
        // fprintf(stderr,"The line is: %s",fline); // contains new line
        token = strtok(fline,sep);    // Time
        // fprintf(stderr,"First token %s\n",token);
        // fprintf(stderr,"First token %d\n",atoi(token));
        // exit(1);
        token = strtok(NULL, sep);      // nr of antib in field at this pos
        
        token = strtok(NULL, sep);      // antibs in field at this pos
        // fprintf(stderr,"Got token %s\n",token);
        
        token = strtok(NULL, sep);      // identifier of this guy
        // fprintf(stderr,"Got identifier %s\n",token);
        if(token==NULL) continue;
        int val2 = atoi(token);
        token = strtok(NULL, sep);      // genome of this guy
        if(token==NULL) continue;
        if(token[0]!='n' || token[0]!=','){
          world[i][j].val = 1+val2%10;
          world[i][j].val2 = val2;
          
          //world[i][j].fval5=(double)(global_tag++); //to be used for log file
          
          strcpy(world[i][j].seq,token);
          int genomelen = strlen(world[i][j].seq); world[i][j].seq[genomelen]='\0';
          if(Genome2genenumber(world[i][j].seq, 'A') ){
            token = strtok(NULL, sep);      //AB genes of this guy
            strcpy(fab_string, token);
            token2 = strtok(fab_string,sepab);
            
            for(k=0;k<MAXSIZE;k++){
              if( world[i][j].seq[k]=='\0' )break;
              if( world[i][j].seq[k]!='A' ) world[i][j].valarray[k]=-1;
              if( world[i][j].seq[k]=='A' ){
                // if(first_time)
                if(token2 != NULL) world[i][j].valarray[k] = atoi(token2); //finally sets AB
                else{
                  fprintf(stderr,"Boia, vecio!\n");
                  exit(1);
                }
                token2 = strtok(NULL, sepab);
              }
              
            }
          }
          
          (*pt_Regulation)(&world[i][j]); // sets regulation parameters

        }else{
          //empty genome
          continue;
        }
        
      }else{
        fprintf(stderr,"Got NULL gets or reached end of file, hope these i,j satisfy you: %d %d\n",i,j);
        fprintf(stderr,"%s\n",fline);
        fprintf(stderr,"At line = %d",line_nr);
        sig_break=1;
        break;
      }
    }
    if(sig_break) break;
  }
  // fprintf(stderr,"InitialiseFromInput(): Error. I don't exist yet!\n");
  // Exit(1);
  if(mix_between_seasons) ChangeSeasonMix(world);
  else ChangeSeasonNoMix(world);
}

void InitialiseFromSingleGenome(const char* init_genome, char* init_ab_gen, TYPE2 **world,TYPE2 **antib){
  int i,j,k;
  const char sepab[2]=",";
  char *token2;
   
  //Initialise the field
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++) {
    world[i][j].val=0;//only for colour
    world[i][j].val2=0;//strain indication
    antib[i][j].val2=0;//how many different AB strains
    antib[i][j].val=0;//only for colour
    for(k=0;k<MAXSIZE;k++) antib[i][j].valarray[k]=0;
    for(k=0;k<MAXSIZE;k++) world[i][j].seq[k]='\0';
  }
  
  // place sequence in the middle
  i=nrow/2; j=ncol/2;  
  world[i][j].val=1;
  world[i][j].val2=2;
  strcpy( world[i][i].seq, init_genome);
  world[i][i].seq[ strlen(init_genome) ]='\0'; //It may well be that this is not needed...
  
  //world[i][j].fval5=(double)(global_tag++);

  if(Genome2genenumber(world[i][j].seq, 'A') ){
    token2 = strtok(init_ab_gen,sepab);
    
    for(k=0;k<MAXSIZE;k++){
      if( world[i][j].seq[k]=='\0' )break;
      if( world[i][j].seq[k]!='A' ) world[i][j].valarray[k]=-1;
      if( world[i][j].seq[k]=='A' ){
        // if(first_time)
        if(token2 != NULL) world[i][j].valarray[k] = atoi(token2); //finally sets AB
        else{
          fprintf(stderr,"Boia, vecio!\n");
          exit(1);
        }
        token2 = strtok(NULL, sepab);
      }
      
    }
  }

  (*pt_Regulation)(&world[i][j]);

}

void InitialiseFromScratch(TYPE2 **world,TYPE2 **bact){
  int i,j,k;

  int count=1;
  // Initialise the grid with a bunch of genomes
  int antib_counter = 1;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++) {
    world[i][j].val=0;//only for colour
    world[i][j].val2=0;//strain indication
    antib[i][j].val2=0;//how many different AB strains
    antib[i][j].val=0;//only for colour
    for(k=0;k<MAXSIZE;k++) antib[i][j].valarray[k]=0;
    for(k=0;k<MAXSIZE;k++) world[i][j].seq[k]='\0';
    if( genrand_real1()<0.01) // i==nrow/2 && j==nrow/2)
    {
      world[i][j].val=1+count%10;
      world[i][j].val2=1+count;
      antib_counter += 17;
      
      for(k=0;k<init_genome_size+nr_H_genes_to_stay_alive;k++){
        // if(k<init_genome_size/2) world[i][j].seq[k]='F';
        // else world[i][j].seq[k]='A';
        // world[i][j].seq[0]='B';
        // world[i][j].seq[3+init_genome_size/2]='B';
      
        // world[i][j].seq[k]=AZ[(int)(2*genrand_real2())]; 
        if(k<nr_H_genes_to_stay_alive) world[i][j].seq[k] = AZ[0];
        else{
          world[i][j].seq[k] = AZ[1]; //if no homeost genes needed -> we put none
          if(genrand_real2() < 0.2) world[i][j].seq[k] = AZ[2]; // antib
          if(genrand_real2() < 0.2) world[i][j].seq[k] = AZ[3]; // give break points only once in a while  
        }
        
        if( world[i][j].seq[k]=='A') {
          if(antib_with_bitstring) world[i][j].valarray[k]=(int)( /*antib_counter*/ pow(2,antib_bitstring_length) * genrand_real2() ); //gives random antib type
          else world[i][j].valarray[k]=antib_counter; //same within the same genome
        }
        else world[i][j].valarray[k]=-1;
      }

      (*pt_Regulation)(&world[i][j]);  //sets regulation parameters - tested, works fine :)
      
      world[i][j].seq[init_genome_size+nr_H_genes_to_stay_alive]='\0';
      
      count++;
      
    }
  }
}

void Exit(int exit_number){
  
  fprintf(stderr,"Exit(): exit. The program:\n");
  
  for(int i = 0; i < (int)argc_g; i++){
    fprintf(stderr,"%s ",argv_g[i]);
  }
  fprintf(stderr,"\nwill exit now, with exit_number %d\n",exit_number);
  
  exit(exit_number);

}
