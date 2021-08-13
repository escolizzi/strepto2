#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "parameters.h"
#include "hash_table.h"

// Creates a pointer to a new hash table item
HT_ITEM* CreateItem(const int* genome, int len_genome, int ab, double birthrate) {
    //printf("CreateItem(): going to create element\n");
    //reserve memory for the whole item
    HT_ITEM* item = (HT_ITEM*) malloc (sizeof(HT_ITEM));
    //copy genome
    item->genome = (int*) malloc ( len_genome * sizeof(int) );
    memcpy(item->genome, genome, len_genome * sizeof(int));
    item->len_genome = len_genome;
    
    //reserve mem for ab,birthrate pairs
    item-> ab_minhd = (AB_BIRTHRATE*) malloc ( sizeof(AB_BIRTHRATE) );
    item-> ab_minhd->ab = ab;
    item-> ab_minhd->minhd = birthrate;
    item-> ab_minhd->next = NULL;
    
    item->left=NULL;
    item->right=NULL;
    
    //printf("CreateItem(): done. item->genome[0] = %d, ab = %d, minhd = %d\n",item->genome[0],item-> ab_minhd->ab,item-> ab_minhd->minhd);
    
    return item;
}

// int prime_powers[5]={923521,29791,961,31,1};
int prime_powers[5]={1,1,1,1,1};
int HashGenome(const int* genome, int len_genome){
    int i,ghash=0;
    //printf("Starting part ghash = %d\n",ghash);
    if(len_genome<5){
        for(i=0;i<len_genome;i++) {
            ghash += prime_powers[i]*genome[i];
            //printf("part ghash = %d\n",ghash);
        }
    }else{
        for(i=0;i<5;i++) {
            ghash += prime_powers[i]*genome[i];
            //printf("part ghash = %d\n",ghash);
        }
        for(;i<len_genome;i++) {
            ghash += genome[i];
            //printf("part ghash = %d\n",ghash);
        }
    }
    //printf("HashGenome(): ghash before modulo = %d\n ",ghash);
    return ghash%HASTHTABLE_SIZE_G; // randomly chosen
}

// from https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
int HashAb(int x){
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x%HASTHTABLE_SIZE_A;
}

void InitialiseHashTabletoNULL(HASHTABLE *ht){
  int i,j;
  for(i=0;i<HASTHTABLE_SIZE_G;i++){
    for(j=0;j<HASTHTABLE_SIZE_A;j++){
      ht->items[i][j]=NULL;
    }
  }
}

int SearchAndInsert_HashTable(HASHTABLE *ht,int gh, int ah, const int * genome, const int len_genome, const int ab){
    int minhd;
    if( ht->items[gh][ah]!=NULL) {
        //printf("SearchAndInsert_HashTable. != NULL, going into SearchItem\n");
        minhd = SearchItem(ht->items[gh][ah], genome, len_genome, ab);
        //printf("SearchAndInsert_HashTable. != NULL, Out of SearchItem, got minhd = %d\n\n",minhd);
        }
    else{
        //printf("SearchAndInsert_HashTable. == NULL, Going to create Item\n");
        minhd = MinHD_hash(genome,len_genome,ab);
        //printf("SearchAndInsert_HashTable. minhd = %d\n",minhd);
        ht->items[gh][ah] = CreateItem(genome,len_genome,ab,minhd);
        //printf("SearchAndInsert_HashTable. Out of CreateItem. ht->items[gh][ah]->genome[0] = %d , ht->items[gh][ah]->ab_minhd->ab = %d\n",ht->items[gh][ah]->genome[0],ht->items[gh][ah]->ab_minhd->ab);
    }
    return minhd;
}

//traverses the binary tree, then looks into AB. Returns minhd.
int SearchItem(HT_ITEM *item, const int * genome, int len_genome,int ab){
    //first test for existence, 
    //if existence  - look for equality
    //      if equality, go down the AB list
    //      if not equality, call recursively checking for left right ... ?     
    // if not existence we are at the end of the tree, and we should place it
    int genome_comp_val;
    int minhd;
    
    //if(item->genome != NULL){
        //0 if equal, <0 if target_genome <= item -> genome, >0 otherwise
        genome_comp_val = CompareGenomes(genome, len_genome, item->genome, item->len_genome); 
        //printf("Comparing genomes, got genome_comp_val = %d\n",genome_comp_val);
        if( genome_comp_val  == 0 ){
            //we got the genome, now check for AB collision
            //printf("We got genome comp val == 0, so we search in SearchAB_linkedlist\n");
            minhd = SearchAB_linkedlist(item->ab_minhd, genome, len_genome, ab);
        }else if(genome_comp_val < 0){
            //printf("We got genome comp val < 0, so we go left\n");
            if(item->left != NULL) {
                //printf("item->left already exists, so we go into it\n");
                minhd = SearchItem(item->left, genome, len_genome, ab);
            }else{
                //printf("item->left does NOT exist, so we make it\n");
                item->left = CreateItem(genome,len_genome,ab,MinHD_hash(genome,len_genome,ab));
                minhd = (item->left)->ab_minhd->minhd;
            }
        }else{
            //printf("We got genome comp val > 0, so we go right\n");
            if(item->right != NULL) {
                //printf("item->right already exists, so we go into it\n");
                minhd = SearchItem(item->right, genome, len_genome, ab);
            }else{
                //printf("item->right does NOT exist, so we make it\n");
                item->right = CreateItem(genome,len_genome,ab,MinHD_hash(genome,len_genome,ab));
                minhd = (item->right)->ab_minhd->minhd;
            }
        }
    // }else{
    //     //make it and set stufff
        
    // }
    return minhd;
}

int SearchAB_linkedlist(AB_BIRTHRATE * ab_minhd, const int * genome, int len_genome, int ab){
    int minhd;
    if (ab_minhd->ab == ab){
        //printf("SearchAB_linkedlist(): found ab_minhd->ab == ab: %d == %d\n", ab_minhd->ab, ab);
        minhd = ab_minhd->minhd;
        //printf("SearchAB_linkedlist(): returning minhd = %d\n", minhd);
    }else{ 
        //printf("SearchAB_linkedlist(): not found here (ab_minhd->ab != ab)), we go to next. Does it exist?\n");
        if( ab_minhd->next == NULL ){
            //end of list, set it
            //printf("SearchAB_linkedlist(): ab_minhd->next == NULL. So we create it\n");
            minhd = MinHD_hash(genome,len_genome,ab);
            ab_minhd->next = (AB_BIRTHRATE*) malloc ( sizeof(AB_BIRTHRATE) );
            (ab_minhd->next)-> ab = ab;
            (ab_minhd->next)-> minhd = minhd;
            (ab_minhd->next)-> next = NULL;
            //printf("SearchAB_linkedlist(): ab_minhd->next created, we can return with minhd = %d\n", (ab_minhd->next)-> minhd);
        }else{
            //printf("SearchAB_linkedlist(): ab_minhd->next != NULL. So we go inside it\n");
            minhd = SearchAB_linkedlist(ab_minhd->next, genome, len_genome, ab);
        }
    }
    
    return minhd;
}

int CompareGenomes(const int * src_genome, int src_len_genome, const int * trg_genome, int trg_len_genome){
    int glength, memc_val;
    if(src_len_genome == trg_len_genome){
        return memcmp(src_genome, trg_genome, src_len_genome*sizeof(src_genome[0]));
    }
    glength = (src_len_genome<trg_len_genome)?src_len_genome:trg_len_genome;
    memc_val = memcmp(src_genome, trg_genome, glength*sizeof(src_genome[0]));
    if(memc_val != 0) return memc_val;
    else{
      //the two genomes are equal in the part that they have in common
      // so if src is shorter we return negative, else positive
      if(src_len_genome<trg_len_genome) return -1;
      else return 1;
    }
}

//CALCULATE AB BIRTHRATE HERE -- JUST MOVE THE FUNCTION HERE
int HammingDistance_hash(int a, int b){
  int h=0;
  int xor = a ^ b;
  while (xor>0){
    h += xor & 1;
    xor >>= 1;
  }
  return h;
}

//MinHD_hash() checks, for every antib at the pos of an individual, 
//the antib with the smallest Hamm Dist
// returns min Hamm Dist
int MinHD_hash(const int * genome, int len_genome, int ab){
    int i;
    int hammdist=999;
    int h=0;
      
    for(i=0;i<len_genome;i++){
        h = HammingDistance_hash( genome[i] , ab );  //finds HD
        if(h<hammdist) hammdist=h; //check if h is the smallest
    }
    return hammdist; 
}