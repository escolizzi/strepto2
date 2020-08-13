/*
strepto.c, a simulation system for streptomyces evolution, based on cash.
Bacteria have a genome made of two types of genes, 'g' and 'a'.
g genes are needed for growth, a genes for antibiotic production.
Bacteria live on a grid, which is TYPE2** world.
Antibiotics are on int** antib.
So far nothinf happen, except bacterial growth.
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
#include <cash2003.h>
#include <cash2.h>
#include <mersenne.h>
#include <cash2-s.h>

#define MAXSIZE 256 // same as cash2.h, max size of genome

struct point *ab_poslist; //2D array used to efficiently randomize antibiotics placement
int len_ab_poslist; //size of ab_poslist
void InitialiseABPosList(struct point **p_ab_poslist, int *p_len_ab_poslist, int MAXRADIUS); //initializes ab_poslist

//// Biological function declarations
int Mutate(TYPE2 **world, int row, int col);      // Mutate the genome at this position
void UpdateABproduction(int row, int col);        // Places new antibiotics in the field
void ChangeSeason(TYPE2 **world);                 // "Sporulates" bacteria and restarts the field

//Printing, plotting, painting.
void ColourPlanes(TYPE2 **world, TYPE2 **G, TYPE2 **A, TYPE2 **R); 
int Genome2genenumber(const char *seq, char gene);  // Convert genome to int of number of genes equal to char gene
void PrintPopStats(TYPE2 **world, TYPE2 **antib);
void PrintPopFull(TYPE2 **world,TYPE2 **antib);
int tomovie_nrow;
int tomovie_ncol;
int ToMovie(TYPE2 **world, TYPE2 **antib, TYPE2** G, TYPE2** A, TYPE2** R);
void Pause();


// Static global data structures and parameters
static TYPE2** world;
static TYPE2** antib;
static TYPE2** G;
static TYPE2** A;
static TYPE2** R;
//static TYPE2 empty = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.,{'\0'},{'\0'}};

static char* AZ="FAB"; //genome alphabet
int MAXRADIUS = 10; //max distance from bacterium at which antibiotics are placed
int init_genome_size = 15; 
double rscale=10.; 
double p_movement = 0.5;
int par_season_duration=1000;
double ddrate=0.001; //per-gene dupdel prob
double prob_new_brpoint = 0.01; //inflow of new randomly placed breakpoints, per replication
double breakprob=0.01;//0.005; // probability of activating a break point
double spore_fraction=0.001; // fraction of nrow*ncol that sporulates
char par_movie_directory_name[MAXSIZE]="movie_strepto"; //genome alphabet
char par_fileoutput_name[MAXSIZE] = "data_strepto.txt";
int par_movie_period = 20;
int par_outputdata_period = 100;

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
  
	for(int i = 0; i < (int)argc_g; i++)
	{
	  readOut = (char*)argv_g[i];
		if(strcmp(readOut, "-seed") == 0) myseed = atoi(argv_g[i+1]);
		if(strcmp(readOut, "-c") == 0) ncol = atoi(argv_g[i+1]);
		if(strcmp(readOut, "-r") == 0) nrow = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-moviedir") == 0) strcpy( par_movie_directory_name , argv_g[i+1] );
    if(strcmp(readOut, "-datafile") == 0) strcpy( par_fileoutput_name , argv_g[i+1] );
    if(strcmp(readOut, "-movie_period") == 0) par_movie_period = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-data_period") == 0) par_outputdata_period = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-display") == 0) display = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-season") == 0) par_season_duration = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-breakpoint_inflow") == 0) prob_new_brpoint = atof(argv_g[i+1]);
    if(strcmp(readOut, "-spore_fraction") == 0) spore_fraction = atof(argv_g[i+1]);
    if(strcmp(readOut, "-maxtime") == 0) MaxTime = atoi(argv_g[i+1]);
    
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

  ulseedG =myseed;// time(NULL); /* random seed ... if you don't know the time */
  fprintf(stderr,"Seeding with time: %ld\n", ulseedG);

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
      for(k=0;k<init_genome_size;k++){
        // if(k<init_genome_size/2) world[i][j].seq[k]='F';
        // else world[i][j].seq[k]='A';
        // world[i][j].seq[0]='B';
        // world[i][j].seq[3+init_genome_size/2]='B';
       world[i][j].seq[k]=AZ[(int)(2*genrand_real2())]; 
       if(genrand_real2() < 0.2) world[i][j].seq[k] = AZ[2]; //give break points only once in a while  
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
        while(antib[row][col].valarray[i]!=0 && i<antib[row][col].val2){
          if(nei->val2!=antib[row][col].valarray[i]){
            flag=1;
            break;
          }
          i++;
        }
        if(flag) continue;
        
        //check number of growth genes
        double fg=Genome2genenumber(nei->seq,'F');
        double ag=Genome2genenumber(nei->seq,'A');
        //cell has no genome: cannot reproduce
        if (fg+ag==0) continue;
      
        //save direction
        dirarray[counter]=k;

        //double ratio=fg/(fg+ag);
        double fgscale=3.; //with 1 it was doing interesting things, with 5 ab production never happened
        double growth=0.1*fg/(fg+fgscale);// was -> /(ratio+rscale));

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
        // printf("Hello 3\n");
        UpdateABproduction(row, col);
        // printf("Hello 4\n");

      }
    }
    
  }
  // movement + continuous antib + death if you are on it prod maybe helps?
  else if(world[row][col].val2!=0){
    int flag=0;
    if(antib[row][col].valarray[0]!=0){ //if there are antibiotic at all, we check if those produced by focal cells are there
      int i=0;
      while(antib[row][col].valarray[i]!=0 && i<antib[row][col].val2){
        if(world[row][col].val2!=antib[row][col].valarray[i]){
          flag=1;
          break;
        }
        i++;
      }
    }
    if(flag){//death
      world[row][col].val=0;
      world[row][col].val2=0;
      world[row][col].seq[0]='\0';
    }else{
      
      UpdateABproduction(row, col);
      //now we will try to move
      if(p_movement>0.){
        int neidir=1+(int)(8.*genrand_real2());
        
        int neirow, neicol;
        GetNeighborC(world, row, col, neidir, &neirow, &neicol);
        nei=&world[neirow][neicol];
        if(nei->val2==0){ //we can only move if neighbour is empty
          flag=0;
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
          if(!flag && genrand_real2() < p_movement){
            *nei = world[row][col];
            world[row][col].val=0;
            world[row][col].val2=0;
            world[row][col].seq[0]='\0';
          }
        }

  
      }
    } 
  }

}

void Update(void)
{
  
  Asynchronous(); 

  if(Time%par_movie_period==0) {
    ColourPlanes(world,G,A,R); // Pretty big speed cost here, can be fixed later
    if(display) Display(world,antib,G,A,R);
    else ToMovie(world,antib,G,A,R);
  }
  if(Time%par_outputdata_period==0) {
    if(display){
      fprintf(stderr,"Time = %d\n",Time);
      PrintPopStats(world,antib);
    }else{
      PrintPopFull(world,antib);
    }
  }
  
  if(Time%par_season_duration==0 && Time>0){
    ChangeSeason(world);
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
void ChangeSeason(TYPE2 **world)
{
  TYPE2 spores[MAXSIZE];
  int sporenr=0;
  int i,j,k,row,col;
  
  //attempt 2, pick 100 spores at random
  sporenr=spore_fraction*nrow*ncol;
  if(sporenr==0) {fprintf(stderr,"ChangeSeason(): Error.No spores for next generation?\n"); exit(1);}
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
  if(actual_sporenr==0){fprintf(stderr,"ChangeSeason(): No spores were found for next generation? System extinct\n"); exit(1);}
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
      fprintf(stderr,"Error: cannot place new spore %d\n",i);
      exit(1);
    }else{
      world[ipos][jpos]=spores[i];
      //if(1+i==250) printf("val2 was %d, genome %s\n",spores[i].val2,spores[i].seq );
      //for(k=0;k<MAXSIZE;k++) world[ipos][ipos].seq[k]=spores[i].seq[k];
      world[ipos][jpos].val=1+i%10;
      world[ipos][jpos].val2=1+i;
      UpdateABproduction(ipos, jpos);
      //if(world[ipos][jpos].val2 == 250) printf("val2=250, genome %s\n",world[ipos][jpos].seq);
    }

  }
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
  int i,j,ipos,genome_size;//, which_mutations;
  //double tot_mutrate;
  //int swappos;

  genome_size = strlen(icell->seq);
  //strcpy(new_seq,icell->seq);

  //int number_pos_seen = 0;
  // int gsize_before = genome_size;


  //Duplications and Deletions
  
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
    if(genrand_real2() < 0.5) {
      //Deletion
      for(i=mutpos;i<genome_size;i++) icell->seq[i]=icell->seq[i+1];
      genome_size--;
      icell->seq[genome_size]='\0';
      //update mutpos array, decrease positions that are larger than mutpos
      for(j=imut;j<nrmuts;j++) if(mutposarr[j]>mutpos) mutposarr[j]--;
    }else{
      //Duplication
      int duppos=genrand_real2()*genome_size;
      char insertgene=icell->seq[mutpos];
      for(i=genome_size;i>duppos;i--) icell->seq[i] = icell->seq[i-1];
      icell->seq[duppos]=insertgene;
      genome_size++;
      icell->seq[genome_size]='\0';
      for(j=imut;j<nrmuts;j++) if(mutposarr[j]>duppos) mutposarr[j]++;
    }
    
  }

  //break points
  int breakarray[MAXSIZE];
  int breaknr=0;
  for (ipos=0; ipos<genome_size; ipos++){
    if (icell->seq[ipos]=='B'){
      breakarray[breaknr]=ipos;
      breaknr++;
    }
  }
  int b, match,lcount;
  // notice that like this breaks happen more frequently closer to 5' than to 3'
  //Also, one break happens, period.
  if(breaknr>1){
    for(b=0; b<breaknr-1; b++){
      if(genrand_real2()<breakprob){
        // printf("val2 = %d; Break genome \n%s at pos breakarray[b] = %d\n", icell->val2, icell->seq, breakarray[b]);
        match=b+genrand_real2()*(breaknr-b);
        // printf("match pos = %d\n", breakarray[match]);
        
        // int lcount=0;
        for(i=breakarray[match],lcount=0 ; i<genome_size ; i++, lcount++){
          icell->seq[ breakarray[b]+lcount ] = icell->seq[i];
          // lcount++;
        }
        genome_size-= (breakarray[match] - breakarray[b]) ;
        icell->seq[genome_size] = '\0';
        // printf("New genome is \n%s\n", icell->seq);
        
        break;
      }
    }
  }
  //Random break point insertion
  if(genrand_real2()< prob_new_brpoint){ 
    int new_brpos=genrand_real2()*genome_size;
    char insertgene='B';
    for(i=genome_size;i>new_brpos;i--) icell->seq[i] = icell->seq[i-1];
    icell->seq[new_brpos]=insertgene;
    genome_size++;
    icell->seq[genome_size]='\0';
  }  
  
  
  return 1;
}

void UpdateABproduction(int row, int col){
  TYPE2 *icell=&world[row][col];

  int i,k;
  //check number of growth genes
  double fg=Genome2genenumber(icell->seq,'F');
  double ag=Genome2genenumber(icell->seq,'A');
  //cell has no genome: cannot reproduce
  if (ag==0) return;

  //int MAXRADIUS=10;
  double agscale=2.;
  // double ratio = ag/(fg+ag);
  double agprod=0.01*ag/(ag+agscale)*(exp(-3.*fg)); // was -> /(ratio+rscale));
  double howmany_pos_get_ab = bnldev(agprod,len_ab_poslist);
  //ab_poslist
  struct point tmp;
  int pos_ab_poslist;
  for(i=0; i< howmany_pos_get_ab; i++){
    pos_ab_poslist= i+(int)((len_ab_poslist-i)*genrand_real2());
    tmp = ab_poslist[pos_ab_poslist];
    ab_poslist[pos_ab_poslist] = ab_poslist[i];
    ab_poslist[i] = tmp;
  }
  for(i=0;i<howmany_pos_get_ab;i++){
    int ii = row+ab_poslist[i].row;
    int jj = col+ab_poslist[i].col;
    if(ii<1) ii = nrow+ii; 
    if(ii>nrow) ii = ii%nrow;
    if(jj<1) jj = ncol+jj; 
    if(jj>ncol) jj = jj%ncol;
    int found=0;
    for (k=0; k<antib[ii][jj].val2;k++){
      if(antib[ii][jj].valarray[k]==icell->val2){
        found=1;
        break;
      }
    }
    if(!found){
      antib[ii][jj].valarray[antib[ii][jj].val2]=icell->val2;
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
      countF += Genome2genenumber(world[i][j].seq,'F');
      countA += Genome2genenumber(world[i][j].seq,'A');
      countB += Genome2genenumber(world[i][j].seq,'B');
      
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

void PrintPopFull(TYPE2 **world,TYPE2 **antib){
  FILE *fp;
  fp= fopen( par_fileoutput_name, "a" );
  
  for(int i=1;i<=nrow;i++)for(int j=1;j<=ncol;j++){
    fprintf(fp, "%d %d %d ",Time, world[i][j].val2, antib[i][j].val2 );
    if(world[i][j].val2 != 0) fprintf(fp,"%s " , world[i][j].seq);
    else fprintf(fp,"n ");
    if(antib[i][j].val2 != 0) for(int k=0;k<antib[i][j].val2;k++) fprintf(fp,"%d,",antib[i][j].valarray[k] );
    else fprintf(fp,"0,");
    fprintf(fp,"\n");
  }
  fclose(fp);
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

