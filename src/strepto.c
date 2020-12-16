/*
In this version antibiotic diversity can evolve. 
Antibiotic are 12 bits elements, resistance depends on dist. 
Antib are actually coded as integers. :P

strepto.c, a simulation system for streptomyces evolution, based on cash.
Bacteria have a genome made of three types of genes: 'F', 'A', 'B'.
F genes are needed for growth, A genes for antibiotic production, B are break points.
Bacteria live on a grid, which is TYPE2** world.
Antibiotics are on int** antib.
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

//// Biological function declarations
int Mutate(TYPE2 **world, int row, int col);      // Mutate the genome at this position
void UpdateABproduction(int row, int col);        // Places new antibiotics in the field
void ChangeSeasonMix(TYPE2 **world);                 // "Sporulates" bacteria, shuffle spore position and restarts the field
void ChangeSeasonNoMix(TYPE2 **world);                 // "Sporulates" bacteria and restarts the field, NO spore shuffling

//Printing, plotting, painting.
void ColourPlanes(TYPE2 **world, TYPE2 **G, TYPE2 **A, TYPE2 **R); 
int Genome2genenumber(const char *seq, char gene);  // Convert genome to int of number of genes equal to char gene
void PrintPopStats(TYPE2 **world, TYPE2 **antib);
int PrintPopFull(TYPE2 **world,TYPE2 **antib);
int tomovie_nrow;
int tomovie_ncol;
int ToMovie(TYPE2 **world, TYPE2 **antib, TYPE2** G, TYPE2** A, TYPE2** R);
void Pause();

void BreakPoint_Recombination_LeftToRight_SemiHomog(TYPE2 *icel);
void BreakPoint_Recombination_Homog(TYPE2 *icel);
void BreakPointDeletion_RightToLeft(TYPE2 *icel);
void BreakPointDeletion_LeftToRight(TYPE2* icel);

// Static global data structures and parameters
static TYPE2** world;
static TYPE2** antib;
static TYPE2** G;
static TYPE2** A;
static TYPE2** R;
//static TYPE2 empty = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.,{'\0'},{'\0'}};

static char* AZ="FAB"; //genome alphabet
int MAXRADIUS = 10; //max distance from bacterium at which antibiotics are placed
int init_genome_size = 10; 
double rscale=10.; 
double p_movement = 0.5;
int par_season_duration=1000;
double ddrate=0.001; //per-gene dupdel prob
double prob_new_brpoint = 0.01; //inflow of new randomly placed breakpoints, per replication
double breakprob=0.01;//0.005; // probability of activating a break point
double spore_fraction=0.001; // fraction of nrow*ncol that sporulates
char par_movie_directory_name[MAXSIZE]="movie_strepto"; //genome alphabet
char par_fileoutput_name[MAXSIZE] = "data_strepto.txt";
char par_name[MAXSIZE] = "\0"; //name - it will create a par_fileoutput_name: data_[name].txt and a par_movie_directory_name: movie_[name];
int par_movie_period = 20;
int par_outputdata_period = 100;
char init_genome[MAXSIZE]; // initial genome, for specific experiments
int antib_with_bitstring=1; // 1: Antibiotic with bistring; 0: antib without bitstrings
int antib_bitstring_length = 6;
double prob_mut_antibtype_perbit = 0.05; //per bit probability of antibiotic type mutation
double h_ag=2.;
double max_ab_prod_per_unit_time = -1.; // set either by command line argument, or as 1/len_ab_pos
double beta_antib_tradeoff = 3.;
int par_all_vs_all_competition = 1; // only set if we are not using bitstrings
double prob_noABspores_fromouterspace = 0.;
double tmp_prob_noABspores_fromouterspace=-1.;
int burn_in=0;
int par_burn_in_time=10; //this is going to be multiplied to season length
int mix_between_seasons = 1;
char breakpoint_mut_type = 'C'; // S: semi homolog recombination, T: telomeric deletion, c: centromeric towards telomeric (strepto like)
double par_beta_birthrate=0.3;
int const_tot_ab_mut=0;                 // if 1, the per AB mut rate is constant - rather than the per bit mutrate: 
double prob_mut_antibtype_tot = 0.05;    // prob_mut_antibtype_tot is used instead of prob_mut_antibtype_perbit
                                        // to be precise: prob_mut_antibtype_perbit is set to prob_mut_antibtype_tot/antib_bitstring_length

void Initial(void)
{
	// readout parameters
	int myseed = 5431;//time(NULL);
	int myncol = 300;
	int mynrow = 300;
  char* readOut;
  display=0;
  MaxTime = 2147483647; /* default=2147483647 */
  nrow = mynrow; /* # of row (default=100)*/
  ncol = myncol; /* # of column (default=100)*/
  init_genome[0]='\0';
  
	for(int i = 1; i < (int)argc_g; i++)
	{
	  readOut = (char*)argv_g[i];
		if(strcmp(readOut, "-seed") == 0) myseed = atoi(argv_g[i+1]);
		else if(strcmp(readOut, "-c") == 0) ncol = atoi(argv_g[i+1]);
		else if(strcmp(readOut, "-r") == 0) nrow = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-moviedir") == 0) strcpy( par_movie_directory_name , argv_g[i+1] );
    else if(strcmp(readOut, "-datafile") == 0) strcpy( par_fileoutput_name , argv_g[i+1] );
    else if(strcmp(readOut, "-name") == 0) strcpy( par_name , argv_g[i+1] );
    else if(strcmp(readOut, "-movie_period") == 0) par_movie_period = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-data_period") == 0) par_outputdata_period = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-display") == 0) display = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-season") == 0) par_season_duration = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-breakpoint_inflow") == 0) prob_new_brpoint = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-spore_fraction") == 0) spore_fraction = atof(argv_g[i+1]);
    else if(strcmp(readOut, "-maxtime") == 0) MaxTime = atoi(argv_g[i+1]);
    else if(strcmp(readOut, "-init_genome") == 0) {strcpy( init_genome , argv_g[i+1] );init_genome_size=strlen(init_genome);}
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
    else {fprintf(stderr,"Parameter number %d was not recognized, simulation not starting\n",i);exit(1);}
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
  if(dir){ fprintf(stderr, "Directory %s already exists, simulation not starting\n",par_movie_directory_name); exit(1);}
  FILE *fp = fopen(par_fileoutput_name,"r"); 
  if(fp){ fprintf(stderr, "File %s already exists, simulation not starting\n",par_fileoutput_name); exit(1);}

  nplane = 5; /* # of planes (default=0)*/
  scale = 1; /* size of the window (default=2)*/
  margin=10;
  boundary = WRAP; /* the type of boundary: FIXED, WRAP, ECHO (default=WRAP). Note that
		      Margolus diffusion is not supported for ECHO. */

  ulseedG =(myseed>0)?myseed:time(NULL);// time(NULL); /* random seed ... if you don't know the time */
  fprintf(stderr,"Seeding with: %ld\n", ulseedG);
  
  if(const_tot_ab_mut){
    prob_mut_antibtype_perbit = prob_mut_antibtype_tot/(double)(antib_bitstring_length);
    fprintf(stderr,"Constant per AB mut rate used, prob_mut_antibtype_perbit reset to %f\n", prob_mut_antibtype_perbit);
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

  /* useally, one does not have to change the followings */
  /* the value of boundary (default=(TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.})*/
  boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.};
  
  tomovie_nrow=nrow;
  tomovie_ncol=nplane*ncol;
  OpenPNG(par_movie_directory_name, tomovie_nrow, tomovie_ncol); //recall that for some weird reason you must allocate a data structure
                                                                 // of size tomovie_nrow+1 and tomovie_ncol+1, if OpenPNG gets arguments
                                                                 // tomovie_nrow and tomovie_ncol 
}

void InitialPlane(void)
{
  MakePlane(&world,&antib,&G,&A,&R);

  InitialiseABPosList(&ab_poslist, &len_ab_poslist, MAXRADIUS);
  if(max_ab_prod_per_unit_time<0.) max_ab_prod_per_unit_time = 1/(double)len_ab_poslist;

  // for(int i=0; i<len_ab_poslist;i++) printf("%d %d\n",ab_poslist[i].row,ab_poslist[i].col);
  // exit(1);
  /* InitialSet(1,2,3,4,5)
    1: name of plane
    2: number of the state other than the background state
    3: background state or empty state
    4: state I want to put
    5: fraction of cells that get S1 state (0 to 1)
  */
  int i,j,k;
  int count=1;
  // Initialise the grid with a bunch of genomes
  int antib_counter = 1;
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
  {
    world[i][j].val=0;//only for colour
    world[i][j].val2=0;//strain indication
    antib[i][j].val2=0;//how many different AB strains
    antib[i][j].val=0;//only for colour
    for(k=0;k<MAXSIZE;k++) antib[i][j].valarray[k]=0;
    for(k=0;k<MAXSIZE;k++)
    	world[i][j].seq[k]='\0';
    if( genrand_real1()<0.01) // i==nrow/2 && j==nrow/2)
    {
      world[i][j].val=1+count%10;
      world[i][j].val2=1+count;
      antib_counter += 17;
      if(init_genome[0] == '\0'){
        for(k=0;k<init_genome_size;k++){
          // if(k<init_genome_size/2) world[i][j].seq[k]='F';
          // else world[i][j].seq[k]='A';
          // world[i][j].seq[0]='B';
          // world[i][j].seq[3+init_genome_size/2]='B';
        
        // world[i][j].seq[k]=AZ[(int)(2*genrand_real2())]; 
        world[i][j].seq[k] = AZ[0];
        if(genrand_real2() < 0.2) world[i][j].seq[k] = AZ[1];
        if(genrand_real2() < 0.2) world[i][j].seq[k] = AZ[2]; //give break points only once in a while  

        if( world[i][j].seq[k]=='A' ) {
          if(antib_with_bitstring) world[i][j].valarray[k]=(int)( /*antib_counter*/ pow(2,antib_bitstring_length) * genrand_real2() ); //gives random antib type
          else world[i][j].valarray[k]=antib_counter; //same within the same genome
        }
        else world[i][j].valarray[k]=-1;
        }
        world[i][j].val3=Genome2genenumber(world[i][j].seq,'F');
        world[i][j].val4=Genome2genenumber(world[i][j].seq,'A');
        world[i][j].val5=Genome2genenumber(world[i][j].seq,'B');
      }else{
        fprintf(stderr,"Don't use me\n");
        exit(1);
        strcpy(world[i][j].seq,init_genome);
      }
      world[i][j].seq[init_genome_size]='\0';
      // printf("Genome: %s\n", world[i][j].seq);
      count++;
      //world[i][j].val2 = Genome2genenumber(world[i][j].seq,'G');
       UpdateABproduction(i, j);

    }
  }
  //InitialSet(world,1,0,1,0.001);
  //ReadSavedData("glidergun.sav",1,world);
  fprintf(stderr,"\n\nworld is ready. Let's go!\n\n");
  Boundaries2(world);


  // Initialise the grid with a bunch of cells
  

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
  // ColorRGB(200,0,0,245);
  // ColorRGB(200,0,0,245);
}

void print_binary(unsigned int number){
    if(number >> 1){
      print_binary(number >> 1);
    }
    putc((number & 1) ? '1' : '0', stdout);
}

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
  //return (cumhammdist==0)?1.:0.;
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
    if( /*howmany_emptynei<7 ||*/ howmany_emptynei==8) return; // if number is less than this -> not enough resources, returns

    // int howmanynei = 8 - howmany_emptynei; //so this is number of competing cells (1 or 2)
    // printf("howmanynei: %d\n", howmanynei);
    double dirarray[8];
    double grarray[8];

    double totgrowth=0;
    //printf("Hello 0\n");
    //we find who these neighbors are:
    for(k=1;k<9;k++){
      nei = GetNeighborP(world,row,col,k);
      //printf("Getting fitness of genome:\t%s\n",nei->seq);

      //check if this cell is inhibited by any antibiotic present in focal point
      if(nei->val!=0){
        int i=0;
        int flag=0;
        double repprob=BirthRate(nei, &antib[row][col]);
        
        //check number of growth genes
        double fg=nei->val3; //Genome2genenumber(nei->seq,'F');
        double ag=nei->val4; //Genome2genenumber(nei->seq,'A');
        //cell has no fitness genes in genome: cannot reproduce
        if (fg==0) continue;
      
        //save direction
        dirarray[counter]=k;

        //double ratio=fg/(fg+ag);
        double fgscale=3.; //with 1 it was doing interesting things, with 5 ab production never happened
        double regulation_growth = fg/(fg+fgscale);
        double regulation_antib  = ag/(ag+fgscale) - regulation_growth; regulation_antib = (regulation_antib>0.)?regulation_antib:0.;
        
        double growth = repprob*0.1*regulation_growth/(regulation_growth+regulation_antib+0.1); // was -> =repprob*0.1*fg/(fg+fgscale);// was -> /(ratio+rscale));

        //double numgrowthgenes = Genome2genenumber(nei->seq,'G')*growthperG - Genome2genenumber(nei->seq,'p')*growthperG*1.5 - Genome2genenumber(nei->seq,'P')*growthperG*1.5 - costperR*Genome2genenumber(nei->seq,'R') -prodperA*Genome2genenumber(nei->seq,'A'); //- Genome2genenumber(nei->seq,'A')- Genome2genenumber(nei->seq,'R');
        //if(numgrowthgenes<0) numgrowthgenes=0;
        // printf("Fitness: \t%f, fit wo p: %f\n", numgrowthgenes, Genome2genenumber(nei->seq,'G')*growthperG - costperR*Genome2genenumber(nei->seq,'R') -prodperA*Genome2genenumber(nei->seq,'A') );
        grarray[counter]=growth;
        totgrowth+=growth;
        counter++;
      }
      // if(counter==howmanynei) break;
    }
    // printf("Hello 1\n");
    if(counter){
      //double totalgrowth= 1. - exp(-beta*totgrowthgenes);  //total prob of replication
      double rn = genrand_real2()*8.;//(double)counter;
      int counter2=0;
      double cumprob=0.0;
      if(rn < totgrowth /*totalgrowth*/){
        cumprob += grarray[counter2];
        while(cumprob < rn /*&& counter2<counter*/ ){
          counter2++;
          cumprob +=grarray[counter2];
        }
        //now we replicate the individual at direction dirarray[counter] into the focal pixel
        //Replicate( row,col,dirarray[counter]);
        int direction=dirarray[counter2];

        // printf("counter: %d, direction: %d\n",counter, direction);
        TYPE2 winner = GetNeighborS(world,row,col,direction);
        world[row][col] = winner;
        // printf("Hello 2\n");
        if(Mutate(world,row,col)==1){
        //   //printf("Genome after del at pos 0: %s\n",world[row][col].seq );
        //   //printf("Which would then have fitness: %f\n", Genome2genenumber(world[row][col].seq,'G')*growthperG - Genome2genenumber(world[row][col].seq,'p')*growthperG - costperR*Genome2genenumber(world[row][col].seq,'R') -prodperA*Genome2genenumber(nei->seq,'A'));
        //   //printf("From parent with genome: %s\n",winner.seq );
        //   //printf("While the parent had fitness:  %f\n", Genome2genenumber(winner.seq,'G')*growthperG - Genome2genenumber(winner.seq,'p')*growthperG - costperR*Genome2genenumber(winner.seq,'R') -prodperA*Genome2genenumber(winner.seq,'A'));
        //   //Pause();
        //   ;
        //   //printf("Hello\n");
        }
        world[row][col].val3=Genome2genenumber(nei->seq,'F');
        world[row][col].val4=Genome2genenumber(nei->seq,'A');
        world[row][col].val5=Genome2genenumber(nei->seq,'B');
        // printf("Hello 3\n");
        UpdateABproduction(row, col);
        // printf("Hello 4\n");

      }
    }
    
  }
  // movement + continuous antib + death if you are on it prod maybe helps?
  else if(world[row][col].val2!=0){
    //int flag=0;
    double deathprob=1.-BirthRate(&world[row][col], &antib[row][col]);
    if(genrand_real2()<deathprob){//death
      world[row][col].val=0;
      world[row][col].val2=0;
      world[row][col].seq[0]='\0';
      // fprintf(stderr,"Death\n");
    }else{
      
      UpdateABproduction(row, col);
      //now we will try to move
      if(p_movement>0.){
        int neidir=1+(int)(8.*genrand_real2());
        
        int neirow, neicol;
        GetNeighborC(world, row, col, neidir, &neirow, &neicol);
        nei=&world[neirow][neicol];
        if(nei->val2==0){ //we can only move if neighbour is empty
          //flag=0;
          // if(antib[neirow][neicol].valarray[0]!=0){ //if there are foreign ABs, we cannot move into the neighbouring pixel
          //   int i=0;
          //   while(antib[neirow][neicol].valarray[i]!=0 && i<antib[neirow][neicol].val2){
          //    if(world[row][col].val2!=antib[neirow][neicol].valarray[i]){
          //      flag=1;
          //      break;
          //    }
          //     i++;
          //   }
          // }
          if(/*!flag &&*/ genrand_real2() < p_movement){
            *nei = world[row][col];
            world[row][col].val=0;
            world[row][col].val2=0;
            world[row][col].seq[0]='\0';
            // fprintf(stderr,"movement from %d %d to %d %d\n", row,col, neirow, neicol);
          }
        }

  
      }
    } 
  }

}

int n_antib;
void Update(void){
  // PerfectMix(world);
  Asynchronous(); 
  
  if(!burn_in){  
    if(Time%par_movie_period==0){
      ColourPlanes(world,G,A,R); // Pretty big speed cost here, can be fixed later
      if(display) Display(world,antib,G,A,R);
      else ToMovie(world,antib,G,A,R);
    }
    if(Time%par_outputdata_period==0) {
      if(display){
        fprintf(stderr,"Time = %d\n",Time);
        PrintPopStats(world,antib);
      }else{
        n_antib=PrintPopFull(world,antib);
        if(n_antib==0){
          fprintf(stderr,"No antibiotic genes in the field, the simulation will terminate\n");
          exit(1);
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


//function to reseed plane with "spores"
void ChangeSeasonMix(TYPE2 **world)
{
  TYPE2 spores[MAXSIZE];
  int sporenr=0;
  int i,j,k,row,col;
  
  //attempt 2, pick 100 spores at random
  sporenr=spore_fraction*nrow*ncol;
  if(sporenr==0) {fprintf(stderr,"ChangeSeasonMix(): Error.No spores for next generation?\n"); exit(1);}
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
  if(actual_sporenr==0){fprintf(stderr,"ChangeSeasonMix(): No spores were found for next generation? System extinct\n"); exit(1);}
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

  // printf("I got %d spores\n", actual_sporenr);
  
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
      exit(1);
    }else{
      world[ipos][jpos]=spores[i];
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
      UpdateABproduction(ipos, jpos);
      //if(world[ipos][jpos].val2 == 250) printf("val2=250, genome %s\n",world[ipos][jpos].seq);
    
    }
  }
}

//function to reseed plane with "spores"
// without reshiffling their positions between seasons
void ChangeSeasonNoMix(TYPE2 **world)
{
  int i,j,k,row,col,bactnumber=0;
  
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
      
      if(!antib_with_bitstring){
        if(par_all_vs_all_competition){
          for(k=0;k<MAXSIZE;k++) if(world[i][j].seq[k]=='A')world[i][j].valarray[k] = bactnumber;
        }
      }
      world[i][j].val = 1+bactnumber%10;
      world[i][j].val2 =1+bactnumber;
      UpdateABproduction(i,j);
    }
  }

  // printf("I got %d spores\n", bactnumber);
  
}



// What kind of mutations are possible?
// dupl, del, translocations and breaking mut. A-> a, G->g, R->r
// mutations are more frequet at the end of the chrom
// according to function: 1 - exp(-0.0001*genome_size^2)
// because this goes to zero in zero -> epsilon + (1-epsilon)(1 - exp(-0.0001*genome_size^2))
int Mutate(TYPE2** world, int row, int col)
{
  //double epsilon=0.001;
  //int output=0;
  TYPE2 *icell;
  icell=&world[row][col];
  //char new_seq[MAXSIZE];
  int i,j,genome_size;//, which_mutations;
  //double tot_mutrate;
  //int swappos;

  genome_size = strlen(icell->seq);
  //strcpy(new_seq,icell->seq);

  //int number_pos_seen = 0;
  // int gsize_before = genome_size;

  // ANTIBIOTIC TYPE MUTATIONS
  if(antib_with_bitstring){
    for(int i=0;icell->seq[i]!='\0';i++){
      if(icell->seq[i]=='A'){
        if(genrand_real2() < prob_mut_antibtype_perbit) 
          icell->valarray[i]=icell->valarray[i] ^ (1<<(int)(genrand_real2()*antib_bitstring_length)); //this flips one random bit between 0 and 12( ^ is xor mask)
      }
    }
  }

  //DUPLICATION AND DELETIONS
  int nrmuts=bnldev(2*ddrate,genome_size);
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
    if(genrand_real2() < 0.5){          //Deletion
      for(i=mutpos;i<genome_size;i++){ 
        icell->seq[i]=icell->seq[i+1];
        icell->valarray[i]=icell->valarray[i+1];
      }
      genome_size--;
      icell->seq[genome_size]='\0';
      icell->valarray[genome_size]=-1;
      //update mutpos array, decrease positions that are larger than mutpos
      for(j=imut;j<nrmuts;j++) if(mutposarr[j]>mutpos) mutposarr[j]--;
    }else{                             //Duplication
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
      // if(insertgene=='A'){
      //   if(genrand_real2()<prob_mut_antibtype){
      //     // fprintf(stderr,"Time: %d, before,after %d, ",Time, icell->valarray[duppos]);
      //     icell->valarray[duppos]=icell->valarray[duppos] ^ (1<<(int)(genrand_real2()*antib_bitstring_length)); //this flips one random bit between 0 and 12( ^ is xor mask)
      //     // fprintf(stderr,"%d\n",icell->valarray[duppos]);
      //   }
      // }
      genome_size++;
      icell->seq[genome_size]='\0';
      for(j=imut;j<nrmuts;j++) if(mutposarr[j]>duppos) mutposarr[j]++;
    }
    
  }
  
  //goes from left to right, breaks are recombination mediated
  if(breakpoint_mut_type=='S') BreakPoint_Recombination_LeftToRight_SemiHomog(icell); //how it was
  else if(breakpoint_mut_type=='T') BreakPointDeletion_RightToLeft(icell); //no recomb, only 3'->5' instability
  else if(breakpoint_mut_type=='C') BreakPointDeletion_LeftToRight(icell);  //5'->3' instability
  else{ 
    fprintf(stderr,"Mutate(): Error. Unrecognised option for the type of breakpoint mutation\n");
    exit(1);
    if(0) BreakPoint_Recombination_Homog(icell);
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
  
  
  return 1;
}

void BreakPoint_Recombination_Homog(TYPE2* icel){
  
  //     UNDER CONSTRUCTION --- DO NOT USE ME --- //
  fprintf(stderr,"BreakPoint_Recombination_Homog(): Error. Function under construction");
  exit(1);

  // The whole thing is still to be constructed
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
        // printf("\n val2 = %d; Break genome \n%s at pos breakarray[b] = %d\n", icel->val2, icel->seq, breakarray[b]);
        // printf("Old antib gen:\n");
        // for(int bla=0; bla<genome_size;bla++) printf("%d ",icel->valarray[bla]);
        // printf("\n");

        match=b+genrand_real2()*(breaknr-b);
        // printf("match pos = %d\n", breakarray[match]);
        
        // int lcount=0;
        for(int i=breakarray[match],lcount=0 ; i<genome_size ; i++, lcount++){
          seq[ breakarray[b]+lcount ] = seq[i];
          icel->valarray[ breakarray[b]+lcount ] = icel->valarray[i];
          // lcount++;
        }
        genome_size-= (breakarray[match] - breakarray[b]) ;
        seq[genome_size] = '\0';
        // printf("New genome is \n%s\n", icel->seq);
        // printf("New antib gen:\n");
        // for(int bla=0; bla<genome_size;bla++) printf("%d ",icel->valarray[bla]);
        // printf("\n");
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

void UpdateABproduction(int row, int col){
  TYPE2 *icell=&world[row][col];

  int i,k;
  //check number of growth genes
  double fg=icell->val3; // Genome2genenumber(icell->seq,'F');
  double ag=icell->val4; // Genome2genenumber(icell->seq,'A');
  //cell has no genome: cannot reproduce
  if (ag==0) return;

  //int MAXRADIUS=10;
  // double ratio = ag/(fg+ag);
  
  double fgscale=3.; //with 1 it was doing interesting things, with 5 ab production never happened
  double regulation_growth = fg/(fg+fgscale);
  double regulation_antib  = ag/(ag+fgscale) - regulation_growth; regulation_antib = (regulation_antib>0.)?regulation_antib:0.;

  double agprod= max_ab_prod_per_unit_time*regulation_antib/(regulation_growth+regulation_antib+0.1);//was -> //max_ab_prod_per_unit_time*ag/(ag+h_ag)*(exp(-beta_antib_tradeoff*fg));
  
  double howmany_pos_get_ab = bnldev(agprod,len_ab_poslist);
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

int Char2Num(char c)
{
	switch(c){
		case 'F': return 0;
		case 'A': return 1;
		case 'B': return 2;
		default: {
      fprintf(stderr, "Char2Num(): Error, character not recognised, got %c\n",c); 
      return -1;
    }
	}
}
char Num2Char(int c)
{
	switch(c){
		case 0: return 'F';
		case 1: return 'A';
		case 2: return 'B';
		default: fprintf(stderr, "Num2Char(): Error, character not recognised, got number: %d\n",c), exit(1);
	}
}


void PrintPopStats(TYPE2 **world,TYPE2 **antib)
{
  int i,j,k,l;//,max;
  int countA = 0;
  int countB =0;
  int countF=0;
  double sumabs = 0.0;
  // char agenome[MAXSIZE];
	int av_genome[MAXSIZE][8] ={}; //prints consensus genome, rather than a single one
  int cellcount = 0;

  // printf("Print nothing and return for now\n");
  // return;

	for(k=0;k<MAXSIZE;k++)for(l=0;l<8;l++) av_genome[k][l]=0;
	
  for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
  {
    if(world[i][j].val2>0)
    {
      cellcount++;
      // countG += Genome2genenumber(world[i][j].seq,'G');
      // countR += Genome2genenumber(world[i][j].seq,'R');
      countF += world[i][j].val3; //Genome2genenumber(world[i][j].seq,'F');
      countA += world[i][j].val4; //Genome2genenumber(world[i][j].seq,'A');
      countB += world[i][j].val5; //Genome2genenumber(world[i][j].seq,'B');
      
      // countP += Genome2genenumber(world[i][j].seq,'p');
      // countP += Genome2genenumber(world[i][j].seq,'P');
      //strcpy(agenome,world[i][j].seq);
			for(k=0; world[i][j].seq[k] != '\0';k++){
        int char2num = Char2Num(world[i][j].seq[k]);
        if(char2num==-1){
          fprintf(stderr,"Pop stats, val2 =%d,  got garbage from genome. Time = %d, k (garbage pos)=%d, seq = %s\n" ,world[i][j].val2, Time, k, world[i][j].seq );
          exit(1);
        }
				av_genome[k][ char2num ]++; // clearly the better way would be an alignment, because this mismatches positions after dup/dels
			}
    }
    sumabs+=antib[i][j].fval;
  }
  printf("Population stats at T=%d\n",Time);
  printf("POPSIZE:\t%d\n",cellcount);
  printf("AVG F:  \t%f\n",(float)countF/cellcount);
  printf("AVG A:  \t%f\n",(float)countA/cellcount);
  printf("AVG B:  \t%f\n",(float)countB/cellcount);
  // printf("AVG P:  \t%f\n",(float)countP/cellcount);
  // printf("AVG AB: \t%f\n",sumabs/(nrow*ncol));
  // printf("A genome:\t%s\n",agenome);
	// Not very good
  // for(k=0;k<MAXSIZE;k++){
	// 	max=-1;
	// 	for(l=0;l<strlen(AZ);l++){
  //     // printf("seq pos %d, gene=%c, got howmany = %d\n",k,Num2Char( l ), av_genome[k][l]);
  //     if(av_genome[k][l] > max) {
  //       max = av_genome[k][l];
  //       if(max>0) agenome[k] = Num2Char( l ) ;
  //       else agenome[k] = '\0';
	// 	}
  //   }
  //   if(max<=0){
  //     agenome[k] = '\0';
  //     break;
  //   }
	// }
  // printf("A genome:\t%s\n",agenome);
  // printf("\n");

  int tot_per_pos[MAXSIZE]={0};
  for( k=0; k<MAXSIZE; k++ ){
		int tot=0;
    for(l=0; l<strlen(AZ); l++){
      tot+= av_genome[k][l];
    }
    tot_per_pos[k]=tot;
  }
  
  for(l=0;l<strlen(AZ);l++){
    for( k=0; k<MAXSIZE; k++ ){
      if(tot_per_pos[k]){
        printf("%d",(int)(9.999* av_genome[k][l]/((double)tot_per_pos[k]) )  );
      }else{
        printf("\n");
        break;
      }
    }
  }
  
  
  
  double array[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  array[0] = (float)countF/cellcount;
  array[1] = (float)countA/cellcount;
  array[2] = (float)countB/cellcount;
  // array[3] = (float)countP/cellcount;
  PlotArray(array);
}

//returns tot_antib_genes
int PrintPopFull(TYPE2 **world,TYPE2 **antib){
  FILE *fp;
  fp= fopen( par_fileoutput_name, "a" );
  int tot_antib_genes=0; //if this number == 0 only fitness genes... kinda pointless to go on with the sim

  for(int i=1;i<=nrow;i++)for(int j=1;j<=ncol;j++){
    fprintf(fp, "%d %d ",Time, antib[i][j].val2 );
    if(antib[i][j].val2){
      for (int k=0; k<antib[i][j].val2; k++) {
         fprintf(fp,"%d," , antib[i][j].valarray[k]);
      }
    }else fprintf(fp,"0,"); 

    if(world[i][j].val2){
      tot_antib_genes+=world[i][j].val4;
      fprintf(fp, " %d %s ",world[i][j].val2, world[i][j].seq);
      for (int k=0; k<strlen(world[i][j].seq); k++){
        if(world[i][j].seq[k]=='A') fprintf(fp, "%d,",world[i][j].valarray[k]);
      }
    }else{
      fprintf(fp," 0 n 0,"); 
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  return tot_antib_genes;
}

int GetGenomeColour(const char * seq){
    // Hash functie die doet positie^index+1 waarbij R=1, P=2, q=3 (dus die niet)
    //int gsize = strlen(seq);
	char course[2];
	course[0] = (Genome2genenumber(seq,'A')>0)?'A':'_';
	course[1] = (Genome2genenumber(seq,'R')>0)?'R':'_';

    unsigned int hash = 0;
    for(int i = 0; i < 2; i++)
    {
        hash += (int)course[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    hash = ((int)(hash%99))+1;
		if(0==hash) fprintf(stderr, "GetGenomeColour(): Error. hash = 0\n");
    //cout << hash << endl;
    return hash;
}


void ColourPlanes(TYPE2 **world, TYPE2 **G, TYPE2 **A, TYPE2 **R)
{
	int i,j;
  int colordisplacement = 100;
	for(i=1;i<=nrow;i++)for(j=1;j<=ncol;j++)
		if(world[i][j].val==0)
    {
      //world[i][j].val2 = 0;
		  G[i][j].val = 0;
      A[i][j].val = 0;
      R[i][j].val = 0;
    }
    else
    {
      //world[i][j].val = 1;
      //world[i][j].val = GetGenomeColour(world[i][j].seq);
      //if(world[i][j].val < 1) fprintf(stderr,"FUck.\n"), exit(0);
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
  
  // printf("Hello, howmany = %d\n", *p_len_ab_poslist);
  // printf("The last position looks like: row = %d, col = %d \n", (*p_ab_poslist)[len_ab_poslist-1].row,(*p_ab_poslist)[len_ab_poslist-1].col);
  
  // exit(1);
  return;
}

