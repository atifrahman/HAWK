#define HASH_TABLE_LENGTH 50000003
#define SIG_LEVEL 0.05
#define CUTOFF 10000000
#define KMER_LENGTH 31
#define R -1
#define I -2
#define O -3
#define A 0
#define C 1
#define G 2
#define T 3
const char bases[256] = {
  O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, 
  O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O, 
  O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O, 
  O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, 
  O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O
};

