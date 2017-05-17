package vmalloc.algorithm;

import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.Properties;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Population;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.util.TypedProperties;

import vmalloc.Utils;
import vmalloc.evolutionary.VMCwMProblem;

public class BBOAlloc extends EvolutionaryAllocAlgorithm {
    
    private class BBO extends AbstractEvolutionaryAlgorithm {
        
        private Population[] populations;
        private double immigration_rate;
        private double cross_immig_prob;
        private Variation mutation;
        
        public BBO(Problem problem,
                   Population[] populations,
                   double immigration_rate,
                   double cross_immig_prob,
                   Variation mutation,
                   Initialization initialization) {
            super(problem, null, null, initialization);
            assert(immigration_rate >= 0.0 && immigration_rate <= 1.0);
            assert(cross_immig_prob >= 0.0 && cross_immig_prob <= 1.0);
            assert(mutation.getArity() == 1);
            this.populations = populations;
            this.immigration_rate = immigration_rate;
            this.cross_immig_prob = cross_immig_prob;
            this.mutation = mutation;
        }
        
        private int computeConstraintViolation(Solution sol) {
            int violation = 0;
            for (int i = 0; i < sol.getNumberOfConstraints(); ++i) {
                if (sol.getConstraint(i) > 0.0) {
                    violation++;
                }
            }
            return violation;
        }
        
        private Map<Solution, Integer> getConstraintViolations(Population population) {
            Map<Solution, Integer> sol2violation = new IdentityHashMap<Solution, Integer>();
            for (int i = 0; i < population.size(); ++i) {
                Solution sol = population.get(i);
                sol2violation.put(sol, new Integer(computeConstraintViolation(sol)));
            }
            return sol2violation;
        }
        
        private void sortByRank(Population population, final int obj_idx) {
            final Map<Solution, Integer> sol2violation = getConstraintViolations(population);
            population.sort(new Comparator<Solution>() {
                @Override
                public int compare(Solution sol1, Solution sol2) {
                    assert(sol2violation.containsKey(sol1) && sol2violation.containsKey(sol2));
                    int violation1 = sol2violation.get(sol1).intValue();
                    int violation2 = sol2violation.get(sol2).intValue();
                    int violation_diff = violation1 - violation2;
                    if (violation_diff == 0) {
                        double obj_diff = sol1.getObjective(obj_idx) - sol2.getObjective(obj_idx);
                        if (obj_diff > 0.0) {
                            return 1;
                        }
                        else if (obj_diff < 0.0) {
                            return -1;
                        }
                        return 0;
                    }
                    return violation_diff;
                }
            });
        }
        
        private int[] rankSubsystemPopulation(Population population, final int obj_idx) {
            sortByRank(population, obj_idx);
            int[] ranks = new int[population.size()];
            for (int i = 0; i < population.size(); ++i) {
                ranks[i] = i+1;
            }
            return ranks;
        }
        
        private double immigrationRate(int rank, int pop_size) {
            return ((double)rank*(rank+1)) / ((double)pop_size*(pop_size+1));
        }
        
        private double[] computeImmigrationRates(int[] ranks, int pop_size) {
            double[] rates = new double[ranks.length];
            for (int i = 0; i < ranks.length; ++i) {
                rates[i] = immigrationRate(ranks[i], pop_size);
            }
            return rates;
        }
        
        private double immigrationToEmmigrationRate(double rate) {
            return 1.0 - rate;
        }
        
        private double[] computeEmmigrationRates(int[] ranks, int pop_size) {
            double[] rates = computeImmigrationRates(ranks, pop_size);
            for (int i = 0; i < rates.length; ++i) {
                rates[i] = immigrationToEmmigrationRate(rates[i]);
            }
            return rates;
        }
        
        private void normalize(double[] array) {
            double sum = 0.0;
            for (int i = 0; i < array.length; ++i) {
                sum += array[i];
            }
            for (int i = 0; i < array.length; ++i) {
                array[i] /= sum;
            }
        }
        
        private void scale(double[] array) {
            double min = Utils.min(array), max = Utils.max(array), max_min_diff = max-min;
            for (int i = 0; i < array.length; ++i) {
                array[i] = (array[i] - min) / max_min_diff;
            }
        }
        
        private Solution rouletteWheelSelection(Population population, double[] probabilites) {
            Solution result = null;
            double prob_sum = 0.0, roulette_val = PRNG.nextDouble();
            for (int i = 0; i < probabilites.length; ++i) {
                prob_sum += probabilites[i];
                if (roulette_val <= prob_sum) {
                    result = population.get(i);
                    break;
                }
            }
            assert(result != null);
            return result;
        }
        
        private Solution migrate(Solution immig, Solution emmig) {
            Solution new_sol = immig.copy();
            for (int i = 0; i < immig.getNumberOfVariables(); ++i) {
                if (PRNG.nextDouble() <= immigration_rate) {
                    Variable emmig_var = emmig.getVariable(PRNG.nextInt(emmig.getNumberOfVariables()));
                    new_sol.setVariable(i, emmig_var);
                }
                else {
                    new_sol.setVariable(i, immig.getVariable(i));
                }
            }
            return new_sol;
        }
        
        // FIXME: only supports integer encoding
        private int getValue(Variable var) {
            return EncodingUtils.getInt(var);
        }
        
        private double euclideanDistance(Solution sol1, Solution sol2) {
            assert(sol1.getNumberOfVariables() == sol2.getNumberOfVariables());
            double distance = 0.0, npms = instance.getPhysicalMachines().size();
            for (int i = 0; i < sol1.getNumberOfVariables(); ++i) {
                double value1 = getValue(sol1.getVariable(i)), value2 = getValue(sol2.getVariable(i));
                double var_diff = (value1-value2) / npms;
                distance += var_diff * var_diff;
            }
            return Math.sqrt(distance);
        }
        
        private double[] computeDistances(Solution immig, Population emmig_pop) {
            double[] distances = new double[emmig_pop.size()];
            for (int i = 0; i < emmig_pop.size(); ++i) {
                distances[i] = euclideanDistance(immig, emmig_pop.get(i));
            }
            return distances;
        }
        
        private boolean allZeros(double[] array) {
            for (int i = 0; i < array.length; ++i) {
                if (array[i] == 0.0) {
                    return false;
                }
            }
            return true;
        }
        
        private Population crossMigration(Population immig_pop, Population emmig_pop) {
            assert(immig_pop.size() == emmig_pop.size());
            Population offspring = new Population();
            for (int i = 0; i < immig_pop.size(); ++i) {
                Solution immig = immig_pop.get(i);
                if (PRNG.nextDouble() <= cross_immig_prob) {
                    double[] distances = computeDistances(immig, emmig_pop);
                    if (allZeros(distances)) {
                        offspring.add(immig);
                    }
                    else {
                        normalize(distances);
                        Solution emmig = rouletteWheelSelection(emmig_pop, distances);
                        offspring.add(migrate(emmig, emmig));
                    }
                }
                else {
                    offspring.add(immig);
                }
            }
            return offspring;
        }
        
        private boolean solutionEquals(Solution sol1, Solution sol2) {
            assert(sol1.getNumberOfVariables() == sol2.getNumberOfVariables());
            for (int j = 0; j < sol1.getNumberOfVariables(); ++j) {
                if (getValue(sol1.getVariable(j)) != getValue(sol2.getVariable(j))) {
                    return false;
                }
            }
            return true;
        }
        
        private boolean solutionInPopulation(Solution sol, Population population) {
            boolean result = false;
            for (int i = 0; i < population.size() && !result; ++i) {
                result = solutionEquals(sol, population.get(i));
            }
            return result;
        }
        
        @Override
        protected void initialize() {
            initialized = true;
            for (int i = 0; i < populations.length; ++i) {
                Solution[] initial_solutions = initialization.initialize();
                evaluateAll(initial_solutions);
                populations[i].addAll(initial_solutions);
            }
        }

        @Override
        protected void iterate() {
            int n_objectives = problem.getNumberOfObjectives(), n_sub_systems = populations.length;
            Population[] offsprings = new Population[n_sub_systems];
            for (int i = 0; i < n_sub_systems; ++i) {
                offsprings[i] = new Population();
            }
            Solution[] elites = new Solution[n_sub_systems];
            // Perform within-subsystem migration
            for (int i = 0; i < n_sub_systems; ++i) {
                Population population = populations[i];
                Population offspring = offsprings[i];
                int[] ranks = rankSubsystemPopulation(population, i % n_objectives);
                int best_rank = ranks[0];
                elites[i] = population.get(0);
                double[] immig_rates = computeImmigrationRates(ranks, population.size());
                scale(immig_rates);
                for (int j = 0; j < population.size(); ++j) {
                    Solution sol = population.get(j);
                    double immig_rate = immig_rates[j];
                    if (PRNG.nextDouble() <= immig_rate) {
                        double[] emmig_rates = computeEmmigrationRates(ranks, population.size());
                        emmig_rates[j] = 0.0;
                        normalize(emmig_rates);
                        Solution emmig_sol = rouletteWheelSelection(population, emmig_rates);
                        offspring.add(migrate(sol, emmig_sol));
                    }
                    else {
                        offspring.add(sol);
                    }
                    if (ranks[j] < best_rank) {
                        best_rank = ranks[j];
                        elites[i] = sol;
                    }
                }
            }
            // Perform cross-subsystem migration
            for (int i = 0; i < n_sub_systems; ++i) {
                for (int j = i+1; j < n_sub_systems; ++j) {
                    if (    i % n_objectives != j % n_objectives && // only do cross migration for different objectives
                            PRNG.nextDouble() <= 1 / (n_sub_systems - n_sub_systems%n_objectives)) {
                        Population new_i_offspring = crossMigration(offsprings[i], offsprings[j]);
                        Population new_j_offspring = crossMigration(offsprings[j], offsprings[i]);
                        offsprings[i] = new_i_offspring;
                        offsprings[j] = new_j_offspring;
                    }
                }
            }
            // Perform mutation
            for (int i = 0; i < n_sub_systems; ++i) {
                Population offspring = offsprings[i];
                for (int j = 0; j < offspring.size(); ++j) {
                    Solution[] parent = { offspring.get(j) };
                    offspring.replace(j, mutation.evolve(parent)[0]);
                }
            }
            // Perform elitism
            for (int i = 0; i < n_sub_systems; ++i) {
                Population offspring = offsprings[i];
                Solution elite = elites[i];
                if (!solutionInPopulation(elite, offspring)) {
                    evaluateAll(offspring);
                    offspring.add(elite);
                    sortByRank(offspring, i % n_objectives);
                    offspring.remove(offspring.size()-1);
                }
            }
            populations = offsprings;
        }
        
        @Override
        public Population getPopulation() {
            Population population = new Population();
            for (int i = 0; i < populations.length; ++i) {
                Population system_pop = populations[i];
                for (int j = 0; j < system_pop.size(); ++j) {
                    population.add(system_pop.get(j));
                }
            }
            return population;
        }
        
        @Override
        public NondominatedPopulation getResult() {
            return new NondominatedPopulation(getPopulation());
        }
        
    }
    
    private class BBOAllocAlgorithmProvider extends AlgorithmProvider {
        
        @Override
        public Algorithm getAlgorithm(String name, Properties properties, Problem problem) {
            if (name.equals("BBO")) {
                TypedProperties typed_props = new TypedProperties(properties);
                int pop_size = typed_props.getInt("populationSize", 3);
                int n_sub_systems = typed_props.getInt("subSystems", 2*problem.getNumberOfObjectives());
                Population[] populations = new Population[n_sub_systems];
                for (int i = 0; i < n_sub_systems; ++i) {
                    populations[i] = new Population();
                }
                Initialization initialization = makeInitializer(problem, pop_size);
                double immig_rate = typed_props.getDouble("immigrationRate", 0.5);
                double cross_mig_rate = typed_props.getDouble("crossMigrationRate", 0.5);
                Variation mutation =
                        OperatorFactory.getInstance().getVariation("svum", properties, problem);
                return new BBO(problem,
                               populations,
                               immig_rate,
                               cross_mig_rate,
                               mutation,
                               initialization);
            }
            return null;
        }
        
    }

    public BBOAlloc(VMCwMProblem instance) {
        super(instance, "BBO", VMCwMProblem.Encoding.INTEGER);
        AlgorithmFactory.getInstance().addProvider(new BBOAllocAlgorithmProvider());
    }
    
    public void setImmigrationRate(double rate) { exec = exec.withProperty("immigrationRate", rate); }
    public void setMutationRate(double rate) { exec = exec.withProperty("svum.rate", rate); }
    public void setSubSystems(int n) { exec = exec.withProperty("subSystems", n); }
    public void setCrossSystemMigrationRate(double rate) {
        exec = exec.withProperty("crossMigrationRate", rate);
    }

}
