## Introduction
This project contains the MATLAB source code for the Reinforcement Learning-assisted Two-stage Constrained Multi-objective Evolutionary Algorithm (**RLTS-CMOEA**). Specifically, the algorithm is designed to solve complex constrained multi-objective optimization problems by virtue of a two-stage search strategy and a reinforcement learning mechanism integrated throughout the entire process (including a novel SDQ-learning estimator)


## File Descriptions

### 1. **`RLTS-CMOEA.m`**
- **Purpose**: Implements the main logic of the RLTS-CMOEA algorithm.
- **Key Features**:
  - Initializes three sub-populations (`PopulationA`, `PopulationB` and `PopulationC`).
  - Generates weight vectors and calculates angular relationships between individuals.
  - Executes a two-stage evolutionary process:
    - **Adaptive angle threshold**: Selects angle threshold based on evolutionary progress.
    - **Forward search**: Explores the unconstrained Pareto front (UPF).
    - **Reverse search**: Optimizes the constrained Pareto front (CPF).
  - Executes reinforcement learning-based resource allocation：
    - **SDQ-learing estimator**: Mitigates estimation bias and accelerates convergence speed.
    - **Prior to convergence**: Allocates computational resources to high-value neighborhoods.
    - **Subsequent to convergence**: Allocates computational resources to underdeveloped regions.
  - Calculate constraint boundaries (`VAR`) to handle sparse feasible solutions.
  - Calls `ArchiveUpdate` to maintain the archive of non-dominated solutions.

### 2. **`ArchiveUpdate.m`**
- **Purpose**: Updates the archive by retaining non-dominated solutions.
- **Key Features**:
  - Filters feasible, non-dominated solutions from the population.
  - Removes duplicates and normalizes solutions.
  - Maintains population diversity using crowding distance.
  - Deletes inferior solutions if the archive size exceeds the limit.

### 3. **`CalFitness.m`**
- **Purpose**: Calculates the fitness of each solution.
- **Key Features**:
  - Differentiates between feasible and infeasible solutions using constraint violation (CV).
  - Detects dominance relationships among solutions to calculate dominance counts.
  - Uses Euclidean distance to maintain diversity and ensure even distribution.
  - Returns a fitness value that balances solution quality and diversity.

### 4. **`EnvironmentalSelection.m`**
- **Purpose**: Performs environmental selection to form the next generation population.
- **Key Features**:
  - Computes fitness using `CalFitness`.
  - Sorts and selects the top N solutions based on fitness.
  - Uses a truncation mechanism to reduce the population size if necessary.

### 5. **`EnvironmentalSelection_VAR.m`**
- **Purpose**: An improved environmental selection function with dynamic constraint boundary adjustment.
- **Key Features**:
  - Filters solutions based on a constraint boundary (`VAR`).
  - Complements feasible solutions with infeasible ones if necessary.

## Instructions for Use

1. Ensure the MATLAB environment is installed and properly configured with the required paths.
2. Place all `.m` files in the same working directory.
3. Use the PlatEMO platform to execute the algorithm. Define the optimization problem using the `Problem` class in PlatEMO.
4. Monitor the optimization results, including the evolution of the population and the final Pareto front.

## Notes
- Input data must conform to the multi-objective optimization problem format, including objective values (`objs`) and constraints (`cons`).
- Algorithm performance depends on parameter settings. Key parameters like the number of weight vectors (`WVs`), number of neighbors (`T`), dynamic angle threshold(`beta`), learning rate of RL (`alp`), and discount factor (`gam`) can be fine-tuned in the code.


## Contact
For any questions or suggestions, please contact the author.
