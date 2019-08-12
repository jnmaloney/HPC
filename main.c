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
typedef struct Problem Problem;

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


struct Problem
{
  MultiSet* T; // Problem sets
  Set* cover; // elements included in our set cover
  int U; // size of Universe
  double k;
  int valid;
};


// Erase all elements in s from each set T
void removeVertexAll(Set* s, MultiSet* T);

// Remove set i from T
void removeSet(int s_i, MultiSet* T);

Problem* create_first_problem();
void split_subproblem(Problem* p, int x, Problem** a, Problem** b);
double problem_k_value(Problem* p);
double alpha(int si);
double beta(int fi);
Set* empty_set();
void remove_vertex_all(MultiSet* T, int i);
void remove_set(MultiSet* T, int i);
void add_element(Set** set, int element);
int largest_set(MultiSet* T);
void measure_and_conquer(Problem* p);
int largestSet(MultiSet* T);
MultiSet* copy(MultiSet* T);

void solve_problem(Problem* p);


MultiSet T;
int N = 0;
int g_leastN;
Set covering;
Set coveringCopy;


Set* empty_set()
{
  Set* s = malloc(sizeof(Set));
  memset(s, 0, sizeof(Set));
  s->element = malloc(N * sizeof(int));
  memset(s->element, 0, N * sizeof(int));
  return s;
}


void remove_vertex_all(MultiSet* T, int i)
{
  removeVertexAll(T->set[i], T);
}


void remove_set(MultiSet* T, int i)
{
  removeSet(i, T);
}


void add_element(Set** set, int element)
{
  //printf("add_element\n");

  if ((*set)->element[element] == 0)
  {
    ++(*set)->n;
  }

  (*set)->element[element] = 1;
}


Problem* create_first_problem()
{
  Problem* p = malloc(sizeof(Problem));
  p->T = copy(&T); // global var
  p->cover = empty_set();
  p->U = N; // global var
  p->valid = 0;
  problem_k_value(p);
}


Set* copy_set(Set* a)
{
  Set* b = malloc(sizeof(Set));
  b->element = malloc(N * sizeof(int));
  memcpy(b->element, a->element, N * sizeof(int));
  b->n = a->n;
  return b;
}


void split_subproblem(Problem* p, int x, Problem** a, Problem** b)
{
  *a = malloc(sizeof(Problem));
  *b = malloc(sizeof(Problem));

  // a includes S_x
  (*a)->T = copy(p->T);
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
  (*b)->T = copy(p->T);
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
  //print_problem(p);

  int n_spaces_in_queue = 12;
  Problem* queue[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  queue[0] = p;
  --n_spaces_in_queue;

  printf("Allocating queue\n");

  //n_spaces_in_queue = 0; //HACK

  // Do until queue is full
  while (n_spaces_in_queue)
  {
    // Find largest k in the pool
    double largest_k = -1.;
    int target_problem = 0;
    for (int i = 0; i < 12; ++i)
    {
      if (queue[i])
      {
        if (largest_k < 0 || largest_k < queue[i]->k)
        {
          largest_k = queue[i]->k;
          target_problem = i;
        }
      }
    }

    // Target problem is found - now find splitting point x
    int x = 0;

    // x is largest set
    int nmax = 0;
    for (int i = 0; i < N; ++i)
    {
      if (queue[target_problem]->T->set[i])
      {
        if (nmax == 0)
        {
          nmax = queue[target_problem]->T->set[i]->n;
          x = i;
        }
        else if (queue[target_problem]->T->set[i]->n > nmax)
        {
          nmax = queue[target_problem]->T->set[i]->n;
          x = i;
        }
      }
    }
    if (nmax == 0)
    {
      //printf("No splitting point is defined!\n");

      // exit the queue
      n_spaces_in_queue = 0;
      break;
    }


    // Split this into 2 subproblems
    //printf("Splitting problem %i on element %i\n", target_problem, x);
    Problem* a = 0;
    Problem* b = 0;
    split_subproblem(queue[target_problem], x, &a, &b);

    // printf("Split problem %i (k=%f) into problem %i (k=%f) and problem %i (k=%f)\n",
    //   target_problem, queue[target_problem]->k, target_problem, a->k, 12 - n_spaces_in_queue, b->k);
    // print_problem(a);
    // print_problem(b);

    queue[target_problem] = a;
    queue[12 - n_spaces_in_queue] = b;
    --n_spaces_in_queue;
  }

  // Print problems
  printf("Problems:\n");
  for (int i = 0; i < 12; ++i)
  {
    if (queue[i])
      printf("%i: \tk = %f\n", i, queue[i]->k);
    else
      printf("%i: ----\n", i);
  }

  // Solve queue
  for (int i = 0; i < 12; ++i)
  {
    if (queue[i])
      solve_problem(queue[i]);
  }

  // Print problems / solutions
  printf("Problems:\n");
  for (int i = 0; i < 12; ++i)
  {
    if (queue[i])
      printf("%i: \tgamma = %i\tvalid = %i\n", i, queue[i]->cover->n, queue[i]->valid);
    else
      printf("%i: ----\n", i);
  }
}


void setGraphSize(int N);
void addNeighbourhood(int a, int b);
void cleanT(MultiSet* T);



// Init T ~ Sets of neighbours

int minSet = 0;

int findSubSet(MultiSet* T);


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
  solve_and_print();

  //dom_size =
  //findMCS(&T, N, &covering);
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
  g_leastN = N;
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
  //printf("copy 0\n");
  MultiSet* T1 = malloc(sizeof(MultiSet));
  T1->n = T0->n;
  T1->set = malloc(N * sizeof(void*));

  //printf("copy 1\n");

  for (int i = 0; i < N; ++i)
  {
    if (T0->set[i])
    {
      // printf("copy s %i / %i\n", i, N);
      //
      // printf("sizeof set %li\n", (long)T1);
      // printf("sizeof set %li\n", (long)T1->set);
      // printf("sizeof set %li\n", sizeof(Set));

      T1->set[i] = malloc(sizeof(Set));
      //printf("x\n");
      T1->set[i]->element = malloc(N * sizeof(int));
      //printf("y\n");
      T1->set[i]->n = T0->set[i]->n;
      //printf("z\n");

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
  //printf("copy 2\n");

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


void set_minset(Set* s)
{
  int leastN = g_leastN;
  if (s->n < leastN)
  {
    // copy set to least set
  }
}


void solve_problem(Problem* p)
{
  // Check search depth and set size
  // ...

  // ...
  if (p->valid) return;

  // Number of sets in T reaches 0
  if (p->T->n == 0)
  {
    if (p->U == 0)
    {
      //set_minset(p->cover);
      p->valid = 1; // useless?
      return;
    }
    else
      return; // No (result is not a dominant set)
  }

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
            removeVertexAll(p->T->set[j], p->T);
            removeSet(j, p->T);
          }
        }
      }
    }
  }

  // // Optimisation #2. #3.
  // int y = findSubSet(p->T);
  // while (y != -1)
  // {
  //   // Assign in-place
  //   // Set S_y is a subset, so remove it
  //   p->cover->element[y] = 0; // redundant
  //   removeSet(y, p->T);
  //
  //   // Repeat
  //   y = findSubSet(p->T);
  // }

  // Number of sets in T reaches 0 [repeat]
  if (p->T->n == 0)
  {
    if (p->U == 0)
    {
      p->valid = 1; // useless?
      return;
    }
    else
      return; // No (result is not a dominant set)
  }

  int s_i = largestSet(p->T);

  Problem* p1;
  Problem* p2;
  split_subproblem(p, s_i, &p1, &p2);

  // Find C1
  solve_problem(p1);
  //cleanT(T1);

  // Find C2
  solve_problem(p2);
  //cleanT(T2);

  // Compare problems
  int n1 = p1->valid ? p1->cover->n : N + p1->cover->n;
  int n2 = p2->valid ? p2->cover->n : N + p2->cover->n;
  // int n1 = p1->cover->n;
  // int n2 = p2->cover->n;
  if (n1 < n2)
  {
    p->cover->n = p1->cover->n;
    p->valid = p1->valid;
  }
  else
  {
    p->cover->n = p2->cover->n;
    p->valid = p2->valid;
  }
}
