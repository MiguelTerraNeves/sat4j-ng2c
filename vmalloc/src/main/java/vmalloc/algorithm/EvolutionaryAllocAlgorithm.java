package vmalloc.algorithm;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.Set;

import org.moeaframework.Executor;
import org.moeaframework.Instrumenter;
import org.moeaframework.analysis.collector.Accumulator;
import org.moeaframework.analysis.collector.AttachPoint;
import org.moeaframework.analysis.collector.Collector;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.PopulationIO;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.Variation;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.real.UM;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.core.spi.OperatorProvider;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.util.ReferenceSetMerger;
import org.moeaframework.util.TypedProperties;
import org.sat4j.core.Vec;
import org.sat4j.specs.IVec;

import vmalloc.Utils;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.evolutionary.VMCwMProblem;

public abstract class EvolutionaryAllocAlgorithm extends AllocAlgorithm {
    
    public static enum InitializationType {
        RANDOM, RANDOM_PACKING, SHUFFLED_FIRST_FIT, SHUFFLED_VMCWM_HEURISTIC, MIXED;
    }
    
    private String algorithm = null;
    private InitializationType init_type = InitializationType.RANDOM;
    
    protected Executor exec = null;
    
    protected double best_energy;
    protected double best_wastage;
    protected double best_migration;
    
    protected void injectBinPackingResultInSolution(Solution solution,
                                                    MappingVec allocation,
                                                    VirtualMachineVec leftover_vms) {
        IVec<Set<Integer>> anti_coloc_jobs_per_pm =
                new Vec<Set<Integer>>(instance.getPhysicalMachines().size());
        for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
            anti_coloc_jobs_per_pm.push(new HashSet<Integer>());
        }
        // Set solution variables for VMs in partial allocation
        assert(allocation != null);
        for (int i = 0; i < allocation.size(); ++i) {
            VirtualMachine vm = allocation.get(i).getVirtualMachine();
            PhysicalMachine pm = allocation.get(i).getPhysicalMachine();
            int vm_idx = this.instance.getVirtualMachineIndex(vm);
            int pm_idx = this.instance.getPhysicalMachineIndex(pm);
            this.instance.setVariableValue(solution.getVariable(vm_idx), pm_idx);
            if (vm.isAntiColocatable()) {
                anti_coloc_jobs_per_pm.get(pm_idx).add(new Integer(vm.getJobID()));
            }
        }
        // Randomly place leftover VMs
        assert(leftover_vms != null);
        for (int i = 0; i < leftover_vms.size(); ++i) {
            VirtualMachine vm = leftover_vms.get(i);
            IVec<Integer> selectable_idxs = new Vec<Integer>();
            for (int j = 0; j < instance.getPhysicalMachines().size(); ++j) {
                if (    !(vm.isAntiColocatable() &&
                          anti_coloc_jobs_per_pm.get(j).contains(new Integer(vm.getJobID()))) &&
                        vm.canRunInPhysicalMachine(instance.getPhysicalMachines().get(j))) {
                    selectable_idxs.push(new Integer(j));
                }
            }
            int pm_idx = selectable_idxs.get(PRNG.nextInt(selectable_idxs.size())).intValue();
            this.instance.setVariableValue(solution.getVariable(i), pm_idx);
            if (vm.isAntiColocatable()) {
                anti_coloc_jobs_per_pm.get(pm_idx).add(new Integer(vm.getJobID()));
            }
        }
    }
    
    protected class EfficientAggregateConstraintComparator extends AggregateConstraintComparator {

        private static final long serialVersionUID = 6563564748784503050L;

        @Override
        protected double getConstraints(Solution solution) {
            return instance.getConstraintViolation(solution);
        }
    }
    
    protected class EfficientParetoDominanceComparator extends ChainedComparator {

        private static final long serialVersionUID = 5136293013643667614L;

        public EfficientParetoDominanceComparator() {
            super(new EfficientAggregateConstraintComparator(), new ParetoObjectiveComparator());
        }
    }
    
    private abstract class FeasibleInitialization extends RandomInitialization {
        
        public FeasibleInitialization(Problem problem, int populationSize) {
            super(problem, populationSize);
        }
        
        @Override
        public Solution[] initialize() {
            System.out.println("c Attempting feasible initialization");
            Solution[] initial_pop = new Solution[populationSize];
            for (int i = 0; i < populationSize; ++i) {
                initial_pop[i] = problem.newSolution();
                attemptFeasibleInitialization(initial_pop[i]);
            }
            return initial_pop;
        }
        
        protected abstract void attemptFeasibleInitialization(Solution solution);
        
    }
    
    private class RandomBinPackingInitialization extends FeasibleInitialization {

        public RandomBinPackingInitialization(Problem problem, int populationSize) {
            super(problem, populationSize);
        }
        
        @Override
        protected void attemptFeasibleInitialization(Solution solution) {
            // Initialize
            BigInteger[] remaining_cpu_caps = new BigInteger[instance.getPhysicalMachines().size()];
            BigInteger[] remaining_mem_caps = new BigInteger[instance.getPhysicalMachines().size()];
            IVec<Set<Integer>> anti_coloc_jobs_per_pm =
                    new Vec<Set<Integer>>(instance.getPhysicalMachines().size());
            for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
                remaining_cpu_caps[i] = instance.getPhysicalMachines().get(i).getCPU();
                remaining_mem_caps[i] = instance.getPhysicalMachines().get(i).getMemory();
                anti_coloc_jobs_per_pm.push(new HashSet<Integer>());
            }
            // Place VMs randomly while trying to satisfy constraints
            VirtualMachine[] shuffled_vms = new VirtualMachine[instance.getVirtualMachines().size()];
            instance.getVirtualMachines().copyTo(shuffled_vms);
            PRNG.shuffle(shuffled_vms);
            for (int i = 0; i < shuffled_vms.length; ++i) {
                VirtualMachine vm = shuffled_vms[i];
                // Retrieve PM indexes that maintain feasibility
                IVec<Integer> selectable_idxs = new Vec<Integer>();
                for (int j = 0; j < instance.getPhysicalMachines().size(); ++j) {
                    BigInteger rem_cpu = remaining_cpu_caps[j];
                    BigInteger rem_mem = remaining_mem_caps[j];
                    if (    vm.getCPU().compareTo(rem_cpu) <= 0 &&
                            vm.getMemory().compareTo(rem_mem) <= 0 &&
                            !(vm.isAntiColocatable() &&
                              anti_coloc_jobs_per_pm.get(j).contains(new Integer(vm.getJobID()))) &&
                            vm.canRunInPhysicalMachine(instance.getPhysicalMachines().get(j))) {
                        selectable_idxs.push(new Integer(j));
                    }
                }
                // If not possible to maintain feasibility, retrieve allowed PM indexes without an
                // anti-colocatable VM
                if (selectable_idxs.size() == 0) {
                    for (int j = 0; j < instance.getPhysicalMachines().size(); ++j) {
                        if (    !(vm.isAntiColocatable() &&
                                  anti_coloc_jobs_per_pm.get(j).contains(new Integer(vm.getJobID()))) &&
                                vm.canRunInPhysicalMachine(instance.getPhysicalMachines().get(j))) {
                            selectable_idxs.push(new Integer(j));
                        }
                    }
                }
                // Randomly select one of the selected PMs
                int pm_idx = selectable_idxs.get(PRNG.nextInt(selectable_idxs.size())).intValue();
                instance.setVariableValue(solution.getVariable(i), pm_idx);
                remaining_cpu_caps[pm_idx] = remaining_cpu_caps[pm_idx].subtract(vm.getCPU());
                remaining_mem_caps[pm_idx] = remaining_mem_caps[pm_idx].subtract(vm.getMemory());
                if (vm.isAntiColocatable()) {
                    anti_coloc_jobs_per_pm.get(pm_idx).add(new Integer(vm.getJobID()));
                }
            }
        }
        
    }
    
    private class BinPackingInitialization extends FeasibleInitialization {

        private boolean min_migrations;
        private boolean mix;
        
        // Auxiliary parameters for mixed initialization
        private int initialized = 0;
        
        public BinPackingInitialization(Problem problem,
                                        int populationSize,
                                        boolean min_migrations,
                                        boolean mix) {
            super(problem, populationSize);
            this.min_migrations = min_migrations;
            this.mix = mix;
        }
        
        @Override
        protected void attemptFeasibleInitialization(Solution solution) {
            // Attempt to place using shuffled Bin Packing
            BinPackingAllocAlgorithm alloc = null;
            if (this.mix && this.initialized == 0) {
                alloc = new BestFitDecreasingAlloc(instance);
            }
            else if (this.mix && this.initialized == 1) {
                alloc = new BestFitDecreasingAlloc(instance);
                alloc.setMappings(new MappingVec());
            }
            else {
                alloc = new FirstFitDecreasingAlloc(instance);
                if (!this.min_migrations) {
                    alloc.setMappings(new MappingVec());
                }
            }
            alloc.enableVirtualMachineShuffling();
            Utils.stdoutDisable();
            alloc.allocate();
            Utils.stdoutEnable();
            // Generate solution from Bin Packing result
            MappingVec allocation = null;
            VirtualMachineVec leftover_vms = null;
            if (alloc.foundSolution()) {
                allocation = alloc.getSolution();
                leftover_vms = new VirtualMachineVec();
            }
            else {
                allocation = alloc.getPartialAllocation();
                leftover_vms = alloc.getLeftoverVirtualMachines();
            }
            injectBinPackingResultInSolution(solution, allocation, leftover_vms);
            ++this.initialized;
        }
        
    }
    
    private class BestCollector implements Collector {
        
        private final Algorithm algorithm;

        public BestCollector(Algorithm algorithm) {
            super();
            this.algorithm = algorithm;
        }

        @Override
        public Collector attach(Object obj) { return new BestCollector((Algorithm)obj); }

        @Override
        public void collect(Accumulator acc) {
            NondominatedPopulation result = algorithm.getResult();
            double new_best_energy = best_energy, new_best_wastage = best_wastage;
            double new_best_mig = best_migration;
            for (int i = 0; i < result.size(); ++i) {
                if (!result.get(i).violatesConstraints()) {
                    double energy_val = result.get(i).getObjective(VMCwMProblem.ENERGY_OBJ_INDEX);
                    double wastage_val = result.get(i).getObjective(VMCwMProblem.WASTAGE_OBJ_INDEX);
                    if (energy_val < new_best_energy) {
                        new_best_energy = energy_val;
                    }
                    if (wastage_val < new_best_wastage) {
                        new_best_wastage = wastage_val;
                    }
                    if (instance.getMappings().size() > 0) {
                        double migration_val =
                                result.get(i).getObjective(VMCwMProblem.MIGRATION_OBJ_INDEX);
                        if (migration_val < new_best_mig) {
                            new_best_mig = migration_val;
                        }
                    }
                }
            }
            boolean print = false;
            if (new_best_energy < best_energy) {
                best_energy = new_best_energy;
                print = true;
            }
            if (new_best_wastage < best_wastage) {
                best_wastage = new_best_wastage;
                print = true;
            }
            if (instance.getMappings().size() > 0 && new_best_mig < best_migration) {
                best_migration = new_best_mig;
                print = true;
            }
            if (print) {
                if (instance.getMappings().size() > 0) {
                    System.out.printf("e %.5f \tw %.5f \tm %d\n",
                                      new_best_energy, new_best_wastage, (long)new_best_mig);
                }
                else {
                    System.out.printf("e %.5f \tw %.5f\n", new_best_energy, new_best_wastage);
                }
                printElapsedTime();
            }
        }

        @Override
        public AttachPoint getAttachPoint() {
            return AttachPoint.isSubclass(Algorithm.class).and(
                    AttachPoint.not(AttachPoint.isNestedIn(Algorithm.class)));
        }
        
    }
    
    private class SVUM extends UM {
        
        public SVUM(double probability) {
            super(probability);
        }
        
        @Override
        public Solution[] evolve(Solution[] parents) {
            if (PRNG.nextDouble() <= getProbability()) {
                Solution result = parents[0].copy();
                Variable variable = result.getVariable(PRNG.nextInt(result.getNumberOfVariables()));
                if (variable instanceof RealVariable) {
                    evolve((RealVariable)variable);
                }
                return new Solution[] { result };
            }
            return parents;
        }
    
    }
    
    private class SVUMProvider extends OperatorProvider {

        @Override
        public String getMutationHint(Problem problem) { return null; }
        @Override
        public String getVariationHint(Problem problem) { return null; }

        @Override
        public Variation getVariation(String name, Properties properties, Problem problem) {
            if (name.equals("svum")) {
                TypedProperties typed_props = new TypedProperties(properties);
                return new SVUM(typed_props.getDouble("svum.rate", 0.05));
            }
            return null;
        }

    }

    public EvolutionaryAllocAlgorithm(VMCwMProblem instance,
                                      String algorithm,
                                      VMCwMProblem.Encoding encoding) {
        super(instance, encoding);
        this.best_energy = Double.MAX_VALUE;
        this.best_wastage = Double.MAX_VALUE;
        this.best_migration = Double.MAX_VALUE;
        this.exec = new Executor();
        this.algorithm = algorithm;
        OperatorFactory.getInstance().addProvider(new SVUMProvider());
    }
    
    private Executor instrumentExecutor(Executor exec) {
        Instrumenter instrumenter = new Instrumenter()
                .withProblem(VMCWM_PROBLEM)
                .withFrequency(1)
                .attach(new BestCollector(null));
        return exec.withInstrumenter(instrumenter);
    }
    
    private Executor setUpExecutor(Executor exec, String alg) {
        exec = exec.withAlgorithm(alg).withProblem(VMCWM_PROBLEM);
        return instrumentExecutor(exec);
    }
    
    protected Initialization makeInitializer(Problem problem, int pop_size) {
        if (init_type == InitializationType.RANDOM_PACKING) {
            return new RandomBinPackingInitialization(problem, pop_size);
        }
        else if (init_type == InitializationType.SHUFFLED_FIRST_FIT) {
            return new BinPackingInitialization(problem, pop_size, false, false);
        }
        else if (init_type == InitializationType.SHUFFLED_VMCWM_HEURISTIC) {
            return new BinPackingInitialization(problem, pop_size, true, false);
        }
        else if (init_type == InitializationType.MIXED) {
            return new BinPackingInitialization(problem, pop_size, false, true);
        }
        return new RandomInitialization(problem, pop_size);
    }
    
    protected boolean solutionEquals(Solution sol1, Solution sol2) {
        int[] x1 = this.instance.getVariableAssignment(sol1);
        int[] x2 = this.instance.getVariableAssignment(sol2);
        for (int i = 0; i < sol1.getNumberOfVariables(); ++i) {
            if (x1[i] != x2[i]) {
                return false;
            }
        }
        return true;
    }
    
    public void setPopulationSize(int size) { exec = exec.withProperty("populationSize", size); }
    public void setInitializationType(InitializationType type) { this.init_type = type; }
    
    public void dumpReferenceSet(String path) throws IOException {
        ReferenceSetMerger ref_merger = new ReferenceSetMerger();
        for (int i = 0; i < this.results.size(); ++i) {
            ref_merger.add("Seed" + i, this.results.get(i));
        }
        PopulationIO.writeObjectives(new File(path), ref_merger.getCombinedPopulation());
    }
    
    protected void dumpSolution(Solution sol) {
        int[] x = this.instance.getVariableAssignment(sol);
        for (int i = 0; i < sol.getNumberOfVariables(); ++i) {
            System.out.print(x[i] + " ");
        }
    }
    
    protected void printSolution(String label, Solution sol) {
        System.out.print("c Solution " + label + ": ");
        dumpSolution(sol);
        System.out.print("\n");
    }
    
    @Override
    // FIXME: if no time limit is given, algorithm returns prematurely
    protected List<NondominatedPopulation> runMultipleSeeds(int nseeds) {
        System.out.println("c Applying evolutionary optimization");
        System.out.println("c Running with " + nseeds + " different seeds");
        this.exec = setUpExecutor(this.exec, algorithm);
        return this.exec.withMaxTime(1000L*getRemainingTime()).runSeeds(nseeds);
    }
    
    @Override
    // FIXME: if no time limit is given, algorithm returns prematurely
    public void allocate() {
        System.out.println("c Applying evolutionary optimization");
        this.exec = setUpExecutor(this.exec, algorithm);
        this.solutions = this.exec.withMaxTime(1000L*getRemainingTime()).run();
    }

}
