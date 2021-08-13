#ifndef HASHT_H
#define HASHT_H

#define HASTHTABLE_SIZE_G 1499
#define HASTHTABLE_SIZE_A 1499

// [ab,birthrate] list

//we have to declare them, otherwise inner pointers to same type don't work
typedef struct __ab_birthrate AB_BIRTHRATE;
typedef struct __Ht_item HT_ITEM;

struct __ab_birthrate{
        int ab; 
        int minhd;
        AB_BIRTHRATE *next;
    }; 

struct __Ht_item{
    int* genome; //key
    int len_genome;
    HT_ITEM *left,*right; //pointers
    AB_BIRTHRATE *ab_minhd; //ab,birthrate pair
};

// Define the Hash Table here
typedef struct __HashTable {
    // Contains an array of pointers
    // to items
    HT_ITEM* items[HASTHTABLE_SIZE_G][HASTHTABLE_SIZE_A];
    int size;
    int count;
} HASHTABLE;



HT_ITEM* CreateItem(const int* genome, int len_genome, int ab, double birthrate);
int HashGenome(const int* genome, int len_genome);
int HashAb(int x);
void InitialiseHashTabletoNULL(HASHTABLE *ht);
int SearchAndInsert_HashTable(HASHTABLE *ht,int gh,int ah, const int * genome, const int len_genome, const int ab);
int SearchItem(HT_ITEM *item, const int * genome, int len_genome,int ab);
int CompareGenomes(const int * src_genome, int src_len_genome, const int * trg_genome, int trg_len_genome);
int SearchAB_linkedlist(AB_BIRTHRATE * ab_br, const int * genome, int len_genome, int ab);
int HammingDistance_hash(int a, int b);
int MinHD_hash(const int * genome,int len_genome, int ab);




#endif