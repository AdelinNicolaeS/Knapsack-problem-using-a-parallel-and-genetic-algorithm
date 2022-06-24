#include "genetic_algorithm.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct _thread_argument {
  individual *current_generation;
  individual *next_generation;
  sack_object *objects;
  int sack_capacity;
  int object_count;
  pthread_barrier_t *barrier;
  int count_threads;
  int count_generations;
  int index;
} thread_argument;

int minim(int a, int b) {
  if (a < b) {
    return a;
  }
  return b;
}

int read_input(sack_object **objects, int *object_count, int *sack_capacity,
               int *generations_count, int *threads_number, int argc,
               char *argv[]) {
  FILE *fp;

  if (argc < 4) {
    fprintf(stderr, "Usage:\n\t./tema1 in_file generations_count threads_number\n");
    return 0;
  }

  fp = fopen(argv[1], "r");
  if (fp == NULL) {
    return 0;
  }

  if (fscanf(fp, "%d %d", object_count, sack_capacity) < 2) {
    fclose(fp);
    return 0;
  }

  if (*object_count % 10) {
    fclose(fp);
    return 0;
  }

  sack_object *tmp_objects =
      (sack_object *)calloc(*object_count, sizeof(sack_object));

  for (int i = 0; i < *object_count; ++i) {
    if (fscanf(fp, "%d %d", &tmp_objects[i].profit, &tmp_objects[i].weight) <
        2) {
      free(objects);
      fclose(fp);
      return 0;
    }
  }

  fclose(fp);

  *threads_number = (int)strtol(argv[3], NULL, 10);
  *generations_count = (int)strtol(argv[2], NULL, 10);

  if (*generations_count == 0) {
    free(tmp_objects);

    return 0;
  }

  *objects = tmp_objects;

  return 1;
}

void print_objects(const sack_object *objects, int object_count) {
  for (int i = 0; i < object_count; ++i) {
    printf("%d %d\n", objects[i].weight, objects[i].profit);
  }
}

void print_generation(const individual *generation, int limit) {
  for (int i = 0; i < limit; ++i) {
    for (int j = 0; j < generation[i].chromosome_length; ++j) {
      printf("%d ", generation[i].chromosomes[j]);
    }

    printf("\n%d - %d\n", i, generation[i].fitness);
  }
}

void print_best_fitness(const individual *generation) {
  printf("%d\n", generation[0].fitness);
}

void compute_fitness_function(const sack_object *objects, individual *genome,
                              int object_count, int sack_capacity) {
  int weight = 0;
  int profit = 0;
  int objects_taken = 0;

  for (int j = 0; j < genome->chromosome_length; ++j) {
    if (genome->chromosomes[j]) {
      weight += objects[j].weight;
      profit += objects[j].profit;
      objects_taken++;
    }
  }
  genome->objects_taken = objects_taken;
  genome->fitness = (weight <= sack_capacity) ? profit : 0;
}

int cmpfunc(const void *a, const void *b) {
  individual *first = (individual *)a;
  individual *second = (individual *)b;

  int res = second->fitness - first->fitness; // decreasing by fitness
  if (res == 0) {
    res = first->objects_taken -
          second->objects_taken; // increasing by number of objects in the sack
    if (res == 0) {
      return second->index - first->index;
    }
  }
  return res;
}

void mutate_bit_string_1(const individual *ind, int generation_index) {
  int i, mutation_size;
  int step = 1 + generation_index % (ind->chromosome_length - 2);

  if (ind->index % 2 == 0) {
    // for even-indexed individuals, mutate the first 40% chromosomes by a given
    // step
    mutation_size = ind->chromosome_length * 4 / 10;
    for (i = 0; i < mutation_size; i += step) {
      ind->chromosomes[i] = 1 - ind->chromosomes[i];
    }
  } else {
    // for even-indexed individuals, mutate the last 80% chromosomes by a given
    // step
    mutation_size = ind->chromosome_length * 8 / 10;
    for (i = ind->chromosome_length - mutation_size; i < ind->chromosome_length;
         i += step) {
      ind->chromosomes[i] = 1 - ind->chromosomes[i];
    }
  }
}

void mutate_bit_string_2(const individual *ind, int generation_index) {
  int step = 1 + generation_index % (ind->chromosome_length - 2);

  // mutate all chromosomes by a given step
  for (int i = 0; i < ind->chromosome_length; i += step) {
    ind->chromosomes[i] = 1 - ind->chromosomes[i];
  }
}

void crossover(individual *parent1, individual *child1, int generation_index) {
  individual *parent2 = parent1 + 1;
  individual *child2 = child1 + 1;
  int count = 1 + generation_index % parent1->chromosome_length;

  memcpy(child1->chromosomes, parent1->chromosomes, count * sizeof(int));
  memcpy(child1->chromosomes + count, parent2->chromosomes + count,
         (parent1->chromosome_length - count) * sizeof(int));

  memcpy(child2->chromosomes, parent2->chromosomes, count * sizeof(int));
  memcpy(child2->chromosomes + count, parent1->chromosomes + count,
         (parent1->chromosome_length - count) * sizeof(int));
}

void copy_individual(const individual *from, const individual *to) {
  memcpy(to->chromosomes, from->chromosomes,
         from->chromosome_length * sizeof(int));
}

void free_generation(individual *generation) {
  int i;

  for (i = 0; i < generation->chromosome_length; ++i) {
    free(generation[i].chromosomes);
    generation[i].chromosomes = NULL;
    generation[i].fitness = 0;
  }
}

void do_distribution(int *start, int *end, int range_start, int range_end,
                     int index, int threads_number) {
  int size = range_end - range_start;
  *start = range_start + index * size / threads_number;
  *end = minim(range_start + (index + 1) * size / threads_number, range_end);
}

void *thread_function(void *arg) {
  thread_argument *argument = (thread_argument *)arg;
  pthread_barrier_t *barrier = argument->barrier;
  int count_generations = argument->count_generations;
  int count_threads = argument->count_threads;
  individual *current_generation = argument->current_generation;
  individual *next_generation = argument->next_generation;
  individual *tmp = NULL;
  int index = argument->index;
  int object_count = argument->object_count;
  int sack_capacity = argument->sack_capacity;
  sack_object *objects = argument->objects;

  int start, end;

  int cursor;

  do_distribution(&start, &end, 0, object_count, index, count_threads);

  for (int i = start; i < end; i++) {
    current_generation[i].chromosomes[i] = 1;
  }

  for (int k = 0; k < count_generations; k++) {
    // calculate fitness for each individual
    cursor = 0;
    do_distribution(&start, &end, 0, object_count, index, count_threads);
    for (int i = start; i < end; i++) {
      compute_fitness_function(objects, &current_generation[i], object_count,
                               sack_capacity);
    }
    // wait to finish all calculations
    pthread_barrier_wait(barrier);

    if (index == 0) {
      qsort(current_generation, object_count, sizeof(*current_generation),
            cmpfunc);
    }

    pthread_barrier_wait(barrier);
    // firstly- keep 30% after sort
    do_distribution(&start, &end, 0, object_count * 3 / 10, index,
                    count_threads);
    for (int i = start; i < end; i++) {
      copy_individual(&current_generation[i], &next_generation[cursor + i]);
    }
    cursor += object_count * 3 / 10;

    // first mutation
    do_distribution(&start, &end, 0, object_count * 2 / 10, index,
                    count_threads);
    for (int i = start; i < end; i++) {
      copy_individual(&current_generation[i], &next_generation[cursor + i]);
      mutate_bit_string_1(&next_generation[cursor + i], k);
    }
    cursor += object_count * 2 / 10;

    // second mutation
    do_distribution(&start, &end, object_count * 2 / 10, object_count * 4 / 10,
                    index, count_threads);
    for (int i = start; i < end; i++) {
      copy_individual(current_generation + i,
                      next_generation + i + object_count * 3 / 10);
      mutate_bit_string_2(next_generation + i + object_count * 3 / 10, k);
    }
    cursor += object_count * 2 / 10;

    // crossover
    do_distribution(&start, &end, 0, object_count * 3 / 10, index,
                    count_threads);
    // start always even to not interfer with other individuals' thread
    start -= start % 2;
    if (end % 2 == 1 && index == count_threads - 1) {
      copy_individual(&current_generation[object_count - 1],
                      &next_generation[object_count * 7 / 10 + end - 1]);
    }
    end -= end % 2;

    for (int i = start; i < end; i += 2) {
      crossover(current_generation + i, next_generation + cursor + i, k);
    }

    tmp = current_generation;
    current_generation = next_generation;
    next_generation = tmp;

    pthread_barrier_wait(barrier);
    do_distribution(&start, &end, 0, object_count, index, count_threads);

    for (int i = start; i < end; ++i) {
      current_generation[i].index = i;
    }

    if (k % 5 == 0 && index == 0) {

      print_best_fitness(current_generation);
    }

    pthread_barrier_wait(barrier);
  }
  do_distribution(&start, &end, 0, object_count, index, count_threads);
  for (int i = start; i < end; i++) {
    compute_fitness_function(objects, &current_generation[i], object_count,
                             sack_capacity);
  }
  pthread_barrier_wait(barrier);
  if (index == 0) {
    int get_highest = current_generation[0].fitness;
    for (int i = 1; i < object_count; i++) {
      if (get_highest < current_generation[i].fitness) {
        get_highest = current_generation[i].fitness;
      }
    }
    printf("%d\n", get_highest);
  }
  pthread_exit(NULL);
}
void run_par_genetic_algorithm(sack_object *objects, int object_count,
                               int generations_count, int sack_capacity,
                               int threads_number) {
  pthread_barrier_t barrier;

  thread_argument arguments[threads_number];
  pthread_t threads[threads_number];

  pthread_barrier_init(&barrier, NULL, threads_number);

  individual *current_generation =
      (individual *)calloc(object_count, sizeof(individual));
  individual *next_generation =
      (individual *)calloc(object_count, sizeof(individual));

  // set initial generation (composed of object_count individuals with a single
  // item in the sack)
  for (int i = 0; i < object_count; ++i) {
    current_generation[i].fitness = 0;
    current_generation[i].chromosomes =
        (int *)calloc(object_count, sizeof(int));
    current_generation[i].index = i;
    current_generation[i].chromosome_length = object_count;

    next_generation[i].fitness = 0;
    next_generation[i].chromosomes = (int *)calloc(object_count, sizeof(int));
    next_generation[i].index = i;
    next_generation[i].chromosome_length = object_count;
  }

  for (int i = 0; i < threads_number; i++) {
    arguments[i].barrier = &barrier;
    arguments[i].count_generations = generations_count;
    arguments[i].current_generation = current_generation;
    arguments[i].next_generation = next_generation;
    arguments[i].count_threads = threads_number;
    arguments[i].index = i;
    arguments[i].objects = objects;
    arguments[i].sack_capacity = sack_capacity;
    arguments[i].object_count = object_count;

    pthread_create(&threads[i], NULL, thread_function, &arguments[i]);
  }
  for (int i = 0; i < threads_number; i++) {
    int r = pthread_join(threads[i], NULL);

    if (r) {
      printf("Eroare la asteptarea thread-ului %d\n", i);
      exit(-1);
    }
  }
  free_generation(current_generation);
  free_generation(next_generation);

  free(current_generation);
  free(next_generation);

  pthread_barrier_destroy(&barrier);
}
