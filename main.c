#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


const char* howToUse = "program \"input_file\" \"output_file\"";


void readInFile(const char*);
void writeOutFile(const char*);

void findMinDominatingSet();

typedef struct Set Set;
typedef struct MultiSet MultiSet;

struct Set
{
  int n;
  int i;
  int* element;
};


struct MultiSet
{
  int n;
  Set** set;
  int* frequency;
};


MultiSet T;
int N = 0;
Set covering;
Set coveringCopy;

void setGraphSize(int N);
void addNeighbourhood(int a, int b);
void cleanT(MultiSet* T);
int largestSet(MultiSet* T);
MultiSet* copy(MultiSet* T);

// Erase all elements in s from each set T
void removeVertexAll(Set* s, MultiSet* T);

// Remove set i from T
void removeSet(int s_i, MultiSet* T);

// Init T ~ Sets of neighbours

int minSet = 0;

int findSubSet(MultiSet* T);


// Minimum Covering Set
int findMCS(MultiSet* T, int U, Set* covering)
{
  //printf("Entering findMCS!\n");

  if (minSet < N && covering->n >= minSet) // TODO ? N
    return 0;

  int setsInT = T->n;

  if (setsInT == 0)
  {
    //printf("returning with U = %d\n", U);
    if (U == 0)
    {
      minSet = covering->n;
      for (int i = 0; i < N; ++i)
      {
        coveringCopy.element[i] = covering->element[i];
      }
      //printf("Covering size %d \n", minSet);


      return 0;
    }
    else
      return 1001; // No (result is not a dominant set)
  }

  // Optimisation #1. (1-frequency)
  for (int i = 0; i < N; ++i)
  {
    if (T->frequency[i] == 1)
    {
      // Find the set that contains v[i]
      for (int j = 0; j < N; ++j)
      {
        if (T->set[j])
        {
          if (T->set[j]->element[i])
          {
            // S_j contains v_i, so take S_j
            MultiSet* Tx = copy(T);
            Set* s = T->set[j];
            int s_n = s->n;
            removeVertexAll(s, Tx);
            removeSet(j, Tx);
            s = 0;
            //addToMCS(s);
            Set covering2;
            covering2.element = malloc(N * sizeof(int));
            memcpy(covering2.element, covering->element, N * sizeof(int));
            covering2.element[j] = 1;
            covering2.n = covering->n + 1;
            int Cx = s_n + findMCS(Tx, U - s_n, &covering2);
            cleanT(Tx);
            return Cx;
          }
        }
      }
    }
  }

  // Optimisation #2. #3.
  int y = findSubSet(T);
  if (y != -1)
  {
    MultiSet* Ty = copy(T);
    removeSet(y, Ty);
    int Cy = findMCS(Ty, U, covering);
    cleanT(Ty);
    return Cy;
  }


  // printf("Sets in T = %d\n", setsInT);

  int s_i = largestSet(T);
  // printf("Set found i = %d\n", s_i);
  Set* s = T->set[s_i];

  // printf("Largest cardinaity = %d\n", s->n);

  // Find C1
  MultiSet* T1 = copy(T);
  removeVertexAll(s, T1);
  removeSet(s_i, T1);
  //addToMCS(s);
  Set covering2;
  covering2.element = malloc(N * sizeof(int));
  for (int i = 0; i < N; ++i)
  {
    covering2.element[i] = covering->element[i];
  }
  covering2.element[s_i] = 1;
  covering2.n = covering->n + 1;
  int C1 = s->n + findMCS(T1, U - s->n, &covering2);
  cleanT(T1);

  // Find C2
  MultiSet* T2 = copy(T);
  removeSet(s_i, T2);
  int C2 = findMCS(T2, U, covering);
  cleanT(T2);

  if (C1 < C2)
  {
    // index_out add s.i
    // printf("add s_i = %d\n", s_i);
    return C1;
  }
  else
  {
    return C2;
  }
}


clock_t begin_time;
clock_t end_time;


int main(int argc, char** argv)
{
  if (argc != 3)
  {
    printf("%s\n", howToUse);
    return 1;
  }

  readInFile(argv[1]);
  begin_time = clock();
  findMinDominatingSet();
  end_time = clock();
  writeOutFile(argv[2]);

  return 0;
}


void readInFile(const char* filename)
{
  printf("Opening %s\n", filename);

  FILE* in_file = fopen(filename, "r");
  char in_num[256];
  //memset(in_num, 0)
  int c = 0;

  // Read first line
  while (N == 0)
  {
    int next = fgetc(in_file);
    char result;
    if (next == EOF)
    {
      printf("End of file encountered, exiting!\n");
      return;
    }
    else
    {
      result = (char)next;
    }

    if (result == '\n')
    {
      N = atoi(in_num);
      printf("%i nodes to read!\n", N);
    }
    else
    {
      in_num[c] = result;
      in_num[c + 1] = 0;
      ++c;
    }
  }

  // Declare graph size
  setGraphSize(N);
  printf("successfully allocated graph size\n");

  // Read subsequent lines
  unsigned long i = 0;
  unsigned long vertex_edges = 0;
  unsigned long total_edges = 0;
  c = 0;
  while (i < N)
  {
    int next = fgetc(in_file);
    char result;
    if (next == EOF)
    {
      printf("End of file encountered, exiting!\n");
      return;
    }
    else
    {
      result = (char)next;
    }

    // Save edge
    if (next == ' ')
    {
      if (c > 0)
      {
        int b = atoi(in_num);
        //addGraphEdge(i, b);
        addNeighbourhood(i, b);
        ++total_edges;
        ++vertex_edges;
      }
      c = 0;
      in_num[c] = 0;
    }
    // Save edge and increment node number
    if (next == '\n')
    {
      if (c > 0)
      {
        int b = atoi(in_num);
        //addGraphEdge(i, b);
        addNeighbourhood(i, b);
        ++total_edges;
        if (vertex_edges == 0)
        {
          // Isolated node
        }
        vertex_edges = 0;
      }
      addNeighbourhood(i, i);
      c = 0;
      in_num[c] = 0;
      ++i;
      if (i % 1000 == 0) printf("...%lu\n", i);
    }
    // Read char to in_num
    else
    {
      in_num[c] = result;
      in_num[c + 1] = 0;
      ++c;
    }
  }

  fclose(in_file);
  printf("%i nodes finished with %lu edges!\n", N, total_edges);

}


void findMinDominatingSet()
{
  //dom_size =
  findMCS(&T, N, &covering);
}


void writeOutFile(const char* filename)
{
  printf("Minimum covering size %d \n", minSet);
  printf("Writing minumum dominant set to %s \n", filename);
  printf("Time spent: %f seconds\n", (float)(end_time - begin_time) / (float)CLOCKS_PER_SEC);

  FILE* out_file = fopen(filename, "w");

  for (int i = 0; i < N; ++i)
  {
    if (coveringCopy.element[i])
    {
      fprintf(out_file, "%d ", i);
    }
  }
  //printf("\n");

  fclose(out_file);
}




void setGraphSize(int N)
{
  T.n = N;
  T.set = malloc(N * sizeof(void*));

  for (int i = 0; i < N; ++i)
  {
    T.set[i] = malloc(sizeof(Set));
    T.set[i]->element = malloc(N * sizeof(int));
    T.set[i]->n = 0;
    T.set[i]->i = i;

    for (int j = 0; j < N; ++j)
    {
      T.set[i]->element[j] = 0;
    }
  }

  T.frequency = malloc(N * sizeof(int));
  memset(T.frequency, 0, N * sizeof(int));

  covering.element = malloc(N * sizeof(int));
  coveringCopy.element = malloc(N * sizeof(int));
  for (int j = 0; j < N; ++j)
  {
    covering.element[j] = 0;
    coveringCopy.element[j] = 0;
  }


  minSet = N;
}


void cleanT(MultiSet* T)
{
  for (int i = 0; i < N; ++i)
  {
    if (T->set[i])
      free(T->set[i]->element);
  }

  free(T->set);

  free(T->frequency);

  free(T);
}


void addNeighbourhood(int a, int b)
{
  if (T.set[a]->element[b] == 0)
  {
    ++T.set[a]->n;
    T.set[a]->element[b] = 1;
  }

  ++T.frequency[b];
}


// Find set with largest size in T
int largestSet(MultiSet* T)
{
  int max = 0;
  int ret = 0;
  for (int i = 0; i < N; ++i)
  {
    if (T->set[i] && T->set[i]->n > max)
    {
        max = T->set[i]->n;
        ret = i;
    }
  }
  return ret;
}


MultiSet* copy(MultiSet* T0)
{
  MultiSet* T1 = malloc(sizeof(MultiSet));
  T1->n = T0->n;
  T1->set = malloc(N * sizeof(void*));

  for (int i = 0; i < N; ++i)
  {
    if (T0->set[i])
    {
      T1->set[i] = malloc(sizeof(Set));
      T1->set[i]->element = malloc(N * sizeof(int));
      T1->set[i]->n = T0->set[i]->n;

      for (int j = 0; j < N; ++j)
      {
        T1->set[i]->element[j] = T0->set[i]->element[j];
      }
    }
    else
    {
      T1->set[i] = 0;
    }
  }

  T1->frequency = malloc(N * sizeof(int));
  memcpy(T1->frequency, T0->frequency, N * sizeof(int));

  return T1;
}


// Remove all elements of s from all sets in T
void removeVertexAll(Set* s, MultiSet* T)
{
  for (int i = 0; i < N; ++i)
  {
    //printf("removeVertexAll %d\n", i);
    if (s->element[i])
    {
      for (int j = 0; j < N; ++j)
      {
        //printf("removeVertexAll %d %d\n", i, j);

        if (T->set[j] && T->set[j]->element[i])
        {
          --T->set[j]->n;
          T->set[j]->element[i] = 0;

          // Has this set ended?
          if (T->set[j]->n == 0)
          {
            free(T->set[j]->element);
            free(T->set[j]);
            T->set[j] = 0;
            --T->n;
          }

          // Frequency
          --T->frequency[i];
        }
      }
    }
  }
}


// Remove set from T
void removeSet(int s_i, MultiSet* T)
{
  if (T->set[s_i])
  {
    --T->n;

    for (int j = 0; j < N; ++j)
    {
      // Frequency
      if (T->set[s_i]->element[j])
        --T->frequency[j];
    }

    T->set[s_i] = 0;
  }
}


// Returns the index of a susbset S in T
int findSubSet(MultiSet* T)
{
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      if (i == j) continue;
      if (T->set[i] == 0) continue;
      if (T->set[j] == 0) continue;

      // is S_i a subset of S_j?
      if (T->set[i]->n > T->set[j]->n) continue;
      int sub = 1;
      for (int k = 0; k < N; ++k)
      {
        if (T->set[i]->element[k] && !T->set[j]->element[k])
        {
          sub = 0;
          break;
        }
      }
      if (sub) return i;
    }
  }

  return -1;
}
