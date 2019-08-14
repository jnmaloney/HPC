#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


//
// Definitions
//
#define OPT_LEVEL_1 1
#define OPT_LEVEL_2 1
//#define OPT_LEVEL_3 1


//
// Struct Types
//
typedef struct Set Set;
typedef struct MultiSet MultiSet;
typedef struct Problem Problem;


struct Set
{
  int n; // Number of elements
  int i; // Identifier
  int* element; // Array of elements
};


struct MultiSet
{
  int n; // Number of sets
  Set** set; // Array of Set pointers
  int* frequency; // Array of ints for counting total element frequency
};


struct Problem
{
  MultiSet* T; // The current state of the problem
  Set* cover; // Elements chosen to include in the set cover
  int U; // Remaining size of the Universe
  double k; // Measure value
  int valid; // The problem has been solved to a complete cover
};


//
//  Function declarations
//


Set* empty_set();
Set* copy_set(Set* a);
void free_set(Set*);

MultiSet* create_multiset();
MultiSet* copy_multiset(MultiSet* T);
void free_multiset(MultiSet* T);

Problem* create_first_problem();
void split_subproblem(Problem* p, int x, Problem** a, Problem** b);
void free_problem(Problem*);

void add_element(Set** set, int element);

int largest_set(MultiSet* T);
int find_subset(MultiSet* T);
void remove_vertex_all(MultiSet* T, int i);
void remove_set(MultiSet* T, int i);

double alpha(int si);
double beta(int fi);
double problem_k_value(Problem* p);
void solve_problem(Problem* p);
void measure_and_conquer(Problem* p);
void solve_and_print();

void read_in_file(const char*);
void write_out_file(const char*);

void setGraphSize(int N);
void addNeighbourhood(int a, int b);

//
// Global vars
//


const char* howToUse = "program \"input_file\" \"output_file\"";

clock_t begin_time;
clock_t end_time;
MultiSet T;
int N = 0;
int g_leastGamma = 0;
Set covering;
Set coveringCopy;

int nProblemSteps = 0;
int nSetsAllocated = 0;
int nSetsDeallocated = 0;
int nMultisetsAllocated = 0;
int nMultisetsDeallocated = 0;
int nProblemsAllocated = 0;
int nProblemsDeallocated = 0;


int findMCS(MultiSet* T, int U, Set* covering);


//
// Entry point
//


int main(int argc, char** argv)
{
  if (argc != 3)
  {
    printf("%s\n", howToUse);
    return 1;
  }

  read_in_file(argv[1]);
  begin_time = clock();
  solve_and_print();
  end_time = clock();
  write_out_file(argv[2]);

  return 0;
}


//
// Function definitions
//


Set* empty_set()
{
  ++nSetsAllocated;
  Set* s = malloc(sizeof(Set));
  memset(s, 0, sizeof(Set));
  s->element = malloc(N * sizeof(int));
  memset(s->element, 0, N * sizeof(int));
  return s;
}


Set* copy_set(Set* a)
{
  if (a == 0) return 0;
  ++nSetsAllocated;
  Set* s = malloc(sizeof(Set));
  s->n = a->n;
  s->element = malloc(N * sizeof(int));
  memcpy(s->element, a->element, N * sizeof(int));
  return s;
}


void free_set(Set* set)
{
  ++nSetsDeallocated;

  free(set->element);
  free(set);
}


MultiSet* create_multiset()
{
  ++nMultisetsAllocated;
  MultiSet* T = malloc(sizeof(MultiSet));
  return T;
}


MultiSet* copy_multiset(MultiSet* T0)
{
  ++nMultisetsAllocated;
  MultiSet* T1 = malloc(sizeof(MultiSet));
  T1->n = T0->n;
  T1->set = malloc(N * sizeof(void*));

  for (int i = 0; i < N; ++i)
  {
    T1->set[i] = copy_set(T0->set[i]);
  }

  T1->frequency = malloc(N * sizeof(int));
  memcpy(T1->frequency, T0->frequency, N * sizeof(int));

  return T1;
}


void free_multiset(MultiSet* T)
{
  ++nMultisetsDeallocated;

  for (int i = 0; i < N; ++i)
  {
    if (T->set[i])
      free_set(T->set[i]);
  }

  free(T->frequency);

  free(T);
}


void free_problem(Problem* P)
{
  ++nProblemsDeallocated;

  free_multiset(P->T);
  free_set(P->cover);
  free(P);
}


void remove_vertex_all(MultiSet* T, int s_i)
{
  Set* s = T->set[s_i];

  for (int i = 0; i < N; ++i)
  {
    if (i == s_i) continue;

    if (s->element[i])
    {
      for (int j = 0; j < N; ++j)
      {
        if (T->set[j] && T->set[j]->element[i])
        {
          --T->set[j]->n;
          T->set[j]->element[i] = 0;

          // Has this set ended?
          if (T->set[j]->n == 0)
          {
            free_set(T->set[j]);
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


void remove_set(MultiSet* T, int i)
{
  if (T->set[i])
  {
    --T->n;

    for (int j = 0; j < N; ++j)
    {
      // Frequency
      if (T->set[i]->element[j])
        --T->frequency[j];
    }

    free_set(T->set[i]);
    T->set[i] = 0;
  }
}


void add_element(Set** set, int element)
{
  if ((*set)->element[element] == 0)
  {
    ++(*set)->n;
  }

  (*set)->element[element] = 1;
}


Problem* create_first_problem()
{
  ++nProblemsAllocated;
  Problem* p = malloc(sizeof(Problem));
  p->T = copy_multiset(&T); // global var
  p->cover = empty_set();
  p->U = N; // global var
  p->valid = 0;
  problem_k_value(p);
}


void split_subproblem(Problem* p, int x, Problem** a, Problem** b)
{
  nProblemsAllocated += 2;

  *a = malloc(sizeof(Problem));
  *b = malloc(sizeof(Problem));

  // a includes S_x
  (*a)->T = copy_multiset(p->T);
  //(*a)->cover = empty_set();
  int u0 = p->T->set[x]->n;
  remove_vertex_all((*a)->T, x);
  remove_set((*a)->T, x);
  // //(*a)->cover = p->cover;
  (*a)->cover = copy_set(p->cover);
  add_element(&((*a)->cover), x);
  (*a)->U = p->U - u0;
  (*a)->valid = 0;
  problem_k_value(*a);
  //printf("Creating A problem\n");
  //printf("a = %i, b = %i\n", a, b);

  // b excludes S_x
  (*b)->T = copy_multiset(p->T);
  //(*b)->cover = empty_set();
  // printf("Copy T\n");
  remove_set((*b)->T, x);
  // printf("rs\n");
  // //(*b)->cover = p->cover;
  (*b)->cover = copy_set(p->cover);
  (*b)->U = p->U;
  (*b)->valid = 0;
  problem_k_value(*b);
  //printf("Creating B problem\n");
}


double alpha(int si)
{
  if (si > 6) si = 6;
  double a[] = {0., 0., 0.377443, 0.754886, 0.909444, 0.976388, 1.};
  return a[si];
}


double beta(int fi)
{
  if (fi > 6) fi = 6;
  double b[] = {0., 0., 0.399418, 0.767579, 0.929850, 0.985614, 0.98232};
  return b[fi];
}


double problem_k_value(Problem* p)
{
  double k = 0.;

  // Frequency elements
  for (int i = 0; i < N; ++i)
  {
    int fi = p->T->frequency[i];
    k += beta(fi);
  }

  // Size elements
  for (int i = 0; i < N; ++i)
  {
    int si = 0.;
    if (p->T->set[i])
    {
      si = p->T->set[i]->n;
    }
    k += alpha(si);
  }

  p->k = k;
}


double multiset_k_value(MultiSet* T)
{
  double k = 0.;

  // Frequency elements
  for (int i = 0; i < N; ++i)
  {
    int fi = T->frequency[i];
    k += beta(fi);
  }

  // Size elements
  for (int i = 0; i < N; ++i)
  {
    int si = 0.;
    if (T->set[i])
    {
      si = T->set[i]->n;
    }
    k += alpha(si);
  }

  return k;
}


void measure_and_conquer(Problem* p)
{
  // //  Find x
  // //x = largest_set(p);
  //
  // // Generate sub problems
  // Problem a;

}


void print_problem(Problem* p)
{
  if (p->valid)
    printf("  Problem, |U| = %i, [VALID]\n", p->U);
  else
    printf("  Problem, |U| = %i\n", p->U);
  printf("MultiSet\n");
  for (int i = 0; i < N; ++i)
  {
    if (p->T->set[i])
    {
      for (int j = 0; j < N; ++j)
      {
        printf("%i ", p->T->set[i]->element[j]);
      }
      printf("\n");
    }
    else
    {
      printf("----\n");
    }
  }
  printf("Frequency\n");
  for (int i = 0; i < N; ++i)
  {
    printf("%i ", p->T->frequency[i]);
  }
  printf("\n");
  printf("Covering\n");
  for (int i = 0; i < N; ++i)
  {
    printf("%i ", p->cover->element[i]);
  }
  printf("\n");

  printf("k = %f\n", p->k);
}


void solve_and_print()
{
  printf("Creating first problem\n");
  Problem* p = create_first_problem();
  //solve_problem(p);
  findMCS(p->T, N, &covering);
  free_problem(p);

  printf("Memory statistics:\n");
  printf("Problems  \t%i allocated \t%i deallocated\n", nProblemsAllocated, nProblemsDeallocated);
  printf("Multisets \t%i allocated \t%i deallocated\n", nMultisetsAllocated, nMultisetsDeallocated);
  printf("Sets      \t%i allocated \t%i deallocated\n", nSetsAllocated, nSetsDeallocated);
  printf("Solution found:\n");
  printf("%i problem branches searched\n", nProblemSteps);
}


void read_in_file(const char* filename)
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


void write_out_file(const char* filename)
{
  printf("Minimum covering size %d \n", g_leastGamma);
  printf("Writing minumum dominant set to %s \n", filename);
  printf("Time spent: %f seconds\n", (float)(end_time - begin_time) / (float)CLOCKS_PER_SEC);

  FILE* out_file = fopen(filename, "w");

  for (int i = 0; i < N; ++i)
  {
    if (coveringCopy.element[i])
    {
      fprintf(out_file, "%d ", i);
      printf("%i ", i);
    }
  }
  printf("\n");

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
int largest_set(MultiSet* T)
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


// Returns the index of a susbset S in T
int find_subset(MultiSet* T)
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


void set_minset(Set* s)
{
  if (g_leastGamma == 0 || s->n < g_leastGamma)
  {
    // Assign new least value
    g_leastGamma = s->n;

    // copy set to least set
    memcpy(coveringCopy.element, s->element, N * sizeof(int));
    coveringCopy.n = s->n;
  }
}


void solve_problem(Problem* p)
{
  // Check search depth and set size
  // ...
  ++nProblemSteps;

  // ...
  if (p->valid)
  {
    return;
  }

  // Check if this path has exceeded the global minimum
  if (g_leastGamma && p->cover->n >= g_leastGamma)
  {
    return;
  }

  // Number of sets in T reaches 0
  // if (p->T->n == 0)
  // {
  //   if (p->U == 0)
  //   {
  //     set_minset(p->cover);
  //     p->valid = 1; // Set Cover (dominant set) found
  //     return;
  //   }
  //   else
  //   {
  //     printf("cccba\n");
  //     return; // No (result is not a dominant set)
  //   }
  // }

  // Optimisation #1. (single-frequency elements)
  for (int i = 0; i < N; ++i)
  {
    if (p->T->frequency[i] == 1)
    {
      // Find the set that contains v[i]
      for (int j = 0; j < N; ++j)
      {
        if (p->T->set[j])
        {
          if (p->T->set[j]->element[i])
          {
            // Assign in-place
            // Set S_j is added to covering
            if (p->cover->element[j] == 0)
            {
              p->cover->element[j] = 1;
              ++p->cover->n;
            }
            int s_n = p->T->set[j]->n;
            p->U -= s_n;
            remove_vertex_all(p->T, j);
            remove_set(p->T, j);
          }
        }
      }
    }
  }

  // Optimisation #2. #3.
  int y = find_subset(p->T);
  while (y != -1)
  {
    // Assign in-place
    // Set S_y is a subset, so remove it
    p->cover->element[y] = 0; // redundant
    remove_set(p->T, y);

    // Repeat
    y = find_subset(p->T);
  }

  // Number of sets in T reaches 0 [repeat]
  if (p->T->n == 0)
  {
    if (p->U == 0)
    {
      set_minset(p->cover);
      p->valid = 1; // Set Cover (dominant set) found
      return;
    }
    else
    {
      return; // No (result is not a dominant set)
    }
  }

  int s_i = largest_set(p->T);

  Problem* p1;
  Problem* p2;
  split_subproblem(p, s_i, &p1, &p2);

  // Choose ordering
  int ab = 0;
  // p1 < p2
  if (p1->k < p2->k) ab = 1;

  if (ab)
  {
    solve_problem(p1);
    free_problem(p1);
    solve_problem(p2);
    free_problem(p2);
  }
  else
  {
    solve_problem(p2);
    free_problem(p2);
    solve_problem(p1);
    free_problem(p1);
  }
}


//
// Recycled code
//

//int findMCS(MultiSet* T, int U, Set* covering);
// int findSubSet(MultiSet* T);
// void removeSet(int s_i, MultiSet* T);
// void removeVertexAll(Set* s, MultiSet* T);
// int largestSet(MultiSet* T);


int findMCS(MultiSet* T, int U, Set* covering)
{
  if (g_leastGamma > 0 && covering->n >= g_leastGamma) // TODO ? N
    return 0;

  int setsInT = T->n;

  if (setsInT == 0)
  {
    if (U == 0)
    {
      g_leastGamma = covering->n;
      for (int i = 0; i < N; ++i)
      {
        coveringCopy.element[i] = covering->element[i];
      }
      return 0;
    }
    else
      return 1001; // No (result is not a dominant set)
  }

#ifdef OPT_LEVEL_1
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
            MultiSet* Tx = copy_multiset(T);
            Set* s = T->set[j];
            int s_n = s->n;
            remove_vertex_all(Tx, j);
            remove_set(Tx, j);
            s = 0;

            Set* covering2 = empty_set();
            memcpy(covering2->element, covering->element, N * sizeof(int));
            covering2->element[j] = 1;
            covering2->n = covering->n + 1;
            int Cx = s_n + findMCS(Tx, U - s_n, covering2);
            free_set(covering2);
            free_multiset(Tx);
            return Cx;
          }
        }
      }
    }
  }
#endif//OPT_LEVEL_1

#ifdef OPT_LEVEL_2
  // Optimisation #2. #3.
  int y = find_subset(T);
  if (y != -1)
  {
    MultiSet* Ty = copy_multiset(T);
    remove_set(Ty, y);
    int Cy = findMCS(Ty, U, covering);
    free_multiset(Ty);
    return Cy;
  }
#endif//OPT_LEVEL_2

  int s_i = largest_set(T);
  Set* s = T->set[s_i];

  int ab = 1;

  // MultiSet* T1 = copy_multiset(T);
  // remove_vertex_all(T1, s_i);
  // remove_set(T1, s_i);

  // MultiSet* T2 = copy_multiset(T);
  // remove_set(T2, s_i);

#ifdef OPT_LEVEL_3
  double k1 = multiset_k_value(T1);
  double k2 = multiset_k_value(T2);
  if (k1 > k2) ab = 0;
#endif//OPT_LEVEL_3

  if (ab) // C1 first, then C2
  {
    // Find C1
    MultiSet* T1 = copy_multiset(T);
    remove_vertex_all(T1, s_i);
    remove_set(T1, s_i);
    Set* covering2 = empty_set();
    covering2->element = malloc(N * sizeof(int));
    for (int i = 0; i < N; ++i)
    {
      covering2->element[i] = covering->element[i];
    }
    covering2->element[s_i] = 1;
    covering2->n = covering->n + 1;
    int C1 = s->n + findMCS(T1, U - s->n, covering2);
    free_set(covering2);
    free_multiset(T1);

    // Find C2
    MultiSet* T2 = copy_multiset(T);
    remove_set(T2, s_i);
    int C2 = findMCS(T2, U, covering);
    free_multiset(T2);
  }
  // else // C2 first, then C1
  // {
  //   // Find C2
  //   int C2 = findMCS(T2, U, covering);
  //   free_multiset(T2);
  //
  //   // Find C1
  //   Set* covering2 = empty_set();
  //   covering2->element = malloc(N * sizeof(int));
  //   for (int i = 0; i < N; ++i)
  //   {
  //     covering2->element[i] = covering->element[i];
  //   }
  //   covering2->element[s_i] = 1;
  //   covering2->n = covering->n + 1;
  //   int C1 = s->n + findMCS(T1, U - s->n, covering2);
  //   free_set(covering2);
  //   free_multiset(T1);
  // }

  return 0;
}
