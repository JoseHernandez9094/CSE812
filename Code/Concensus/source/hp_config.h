#ifndef HP_CONFIG_H
#define HP_CONFIG_H

#include "config/config.h"

constexpr double mut = .01;

EMP_BUILD_CONFIG( HPConfig,
  GROUP(HARDWARE_GROUP, "Hardware settings"),
  VALUE(NUM_ITER, size_t,     200, "Number of iterations per trial."),
  GROUP(MUTATION_GROUP, "Mutation settings"),
  VALUE(MIN_SIZE_G, size_t,   10, "Minimum number of instructions each function will have."),
  VALUE(MAX_SIZE_G, size_t,  30, "Maximum number of instructions each function will have."),
  VALUE(INST_MUT_RATE, double, mut, "Instruction mustation rate!"),
  VALUE(ARG_MUT_RATE, double, mut, "Argument mutation rate!"),
  VALUE(DEL_MUT_RATE, double, mut, "Deletion rate!"),
  VALUE(INS_MUT_RATE, double, mut, "Insertion rate!"),
  GROUP(EXPERIMENT_GROUP, "Experiment settings"),
  VALUE(POP_SIZE,  size_t, 1000, "Population size."),
  VALUE(NUM_GENS,  size_t, 20000, "Number of generations per experiments."),
  VALUE(RNG_SEED,  int,    1, "Random number seed."),
  VALUE(TOURN_SIZE, size_t,    100, "Number or organims competing in tournament selection."),
  VALUE(SNAP_SHOT,  size_t,   100, "Time that we will take a snapshot of population"),
  VALUE(DIM,  size_t,   3, "Dimension of toroidal grid!"),
  VALUE(DIR, std::string, "./", "Directory where the data is going")
)

#endif
