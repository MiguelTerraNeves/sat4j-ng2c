package vmalloc.algorithm;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.moeaframework.Analyzer;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.spi.ProblemFactory;
import org.moeaframework.core.spi.ProblemProvider;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Clock;
import vmalloc.Utils;
import vmalloc.domain.JobVec;
import vmalloc.domain.Mapping;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.evolutionary.VMCwMProblem;

public abstract class AllocAlgorithm {
    
    // =============================================================================================
    // ======================================== PRIVATE API ========================================
    // =============================================================================================

    private long timeout = NO_TIMEOUT;
    
    // Solution related parameters
    private double best_energy = Double.MAX_VALUE;
    private double best_wastage = Double.MAX_VALUE;
    private long best_mig_cost = Long.MAX_VALUE;

    private void copyVecToSet(PhysicalMachineVec vec, Set<PhysicalMachine> set) {
        for (int i = 0; i < vec.size(); ++i) {
            set.add(vec.get(i));
        }
    }
    
    private double[] normalizeRequirements(BigInteger[] reqs, BigInteger cap) {
        double[] norm_reqs = new double[reqs.length];
        for (int i = 0; i < reqs.length; ++i) {
            norm_reqs[i] = Utils.divideBigIntegers(reqs[i], cap, RoundingMode.HALF_EVEN).doubleValue();
        }
        return norm_reqs;
    }
    
    private double[] getObjectiveCosts(int obj_idx) {
        double[] costs = new double[this.solutions.size()];
        for (int i = 0; i < solutions.size(); ++i) {
            costs[i] = solutions.get(i).getObjective(obj_idx);
        }
        return costs;
    }
    
    // =============================================================================================
    // ======================================= PROTECTED API =======================================
    // =============================================================================================
    
    protected static final String VMCWM_PROBLEM = "VMCwMProblem";
    
    protected final VMCwMProblem instance;
    protected NondominatedPopulation solutions = null;
    protected List<NondominatedPopulation> results = null;
    
    protected class VWCwMProblemProvider extends ProblemProvider {
        
        private PhysicalMachineVec pms;
        private JobVec jobs;
        private MappingVec mappings;
        private double max_mig_percentile;
        private VMCwMProblem.Encoding encoding;

        public VWCwMProblemProvider(PhysicalMachineVec pms,
                                    JobVec jobs,
                                    MappingVec mappings,
                                    double max_mig_percentile,
                                    VMCwMProblem.Encoding encoding) {
            this.pms = pms;
            this.jobs = jobs;
            this.mappings = mappings;
            this.max_mig_percentile = max_mig_percentile;
            this.encoding = encoding;
        }
        
        @Override
        public Problem getProblem(String name) {
            if (name.equals(VMCWM_PROBLEM)) {
                return new VMCwMProblem(this.pms,
                                        this.jobs,
                                        this.mappings,
                                        this.max_mig_percentile,
                                        this.encoding);
            }
            return null;
        }
        
        private Solution makeOrigin(int nobj) {
            Solution solution = null;
            if (nobj == 3) {
                solution = new Solution(new double[] {0.0, 0.0, 0.0});
            }
            else {
                solution = new Solution(new double[] {0.0, 0.0});
            }
            return solution;
        }

        @Override
        public NondominatedPopulation getReferenceSet(String name) {
            if (name.equals(VMCWM_PROBLEM)) {
                NondominatedPopulation ref = new NondominatedPopulation();
                Solution ref_sol = makeOrigin(instance.getNumberOfObjectives());
                ref_sol.setObjective(0, instance.getMaxEnergyConsumption());
                ref.add(ref_sol);
                ref_sol = makeOrigin(instance.getNumberOfObjectives());
                ref_sol.setObjective(1, instance.getMaxResourceWastage());
                ref.add(ref_sol);
                if (instance.getNumberOfObjectives() == 3) {
                    ref_sol = makeOrigin(instance.getNumberOfObjectives());
                    ref_sol.setObjective(2, instance.getMaxMigrationCost());
                    ref.add(ref_sol);
                }
                return ref;
            }
            return null;
        }
        
    }
    
    protected double[] getNormalizedCPURequirements(VirtualMachineVec vms, PhysicalMachine pm) {
        return normalizeRequirements(vms.getCPUs(), pm.getCPU());
    }
    protected double[] getNormalizedMemoryRequirements(VirtualMachineVec vms, PhysicalMachine pm) {
        return normalizeRequirements(vms.getMemories(), pm.getMemory());
    }

    protected IVecInt intersectAllowedPhysicalMachineIndexes(PhysicalMachineVec pms,
                                                             VirtualMachine vm1,
                                                             VirtualMachine vm2) {
        IVecInt indexes = new VecInt();
        Set<PhysicalMachine> vm1_unallowed = new HashSet<PhysicalMachine>();
        Set<PhysicalMachine> vm2_unallowed = new HashSet<PhysicalMachine>();
        copyVecToSet(vm1.getUnallowedPhysicalMachines(), vm1_unallowed);
        copyVecToSet(vm2.getUnallowedPhysicalMachines(), vm2_unallowed);
        for (int i = 0; i < pms.size(); ++i) {
            if (!vm1_unallowed.contains(pms.get(i)) && !vm2_unallowed.contains(pms.get(i))) {
                indexes.push(i);
            }
        }
        return indexes;
    }
    
    protected boolean mappedToSamePhysicalMachine(VirtualMachine vm1,
                                                  VirtualMachine vm2,
                                                  Map<VirtualMachine, PhysicalMachine> mapping_map) {
        PhysicalMachine pm1 = null;
        PhysicalMachine pm2 = null;
        if (mapping_map.containsKey(vm1)) {
            pm1 = mapping_map.get(vm1);
        }
        if (mapping_map.containsKey(vm2)) {
            pm2 = mapping_map.get(vm2);
        }
        return (pm1 == null && pm2 == null) ||
               (pm1 != null && pm2 != null && pm1.getID() == pm2.getID());
    }
    
    protected PhysicalMachine[] sortedPhysicalMachines(PhysicalMachineVec pms,
                                                       Comparator<PhysicalMachine> comp) {
        PhysicalMachine[] pm_array = new PhysicalMachine[pms.size()];
        pms.copyTo(pm_array);
        Arrays.sort(pm_array, comp);
        return pm_array;
    }
    
    protected MappingVec solutionToAllocation(Solution solution) {
        MappingVec allocation = new MappingVec();
        int[] x = this.instance.getVariableAssignment(solution);
        for (int i = 0; i < x.length; ++i) {
            allocation.push(new Mapping(this.instance.getVirtualMachines().get(i),
                                        this.instance.getPhysicalMachines().get(x[i])));
        }
        return allocation;
    }
    protected Solution allocationToSolution(MappingVec allocation) {
        Solution solution = this.instance.newSolution();
        for (int i = 0; i < allocation.size(); ++i) {
            PhysicalMachine pm = allocation.get(i).getPhysicalMachine();
            VirtualMachine vm = allocation.get(i).getVirtualMachine();
            int pm_idx = this.instance.getPhysicalMachineIndex(pm);
            int vm_idx = this.instance.getVirtualMachineIndex(vm);
            this.instance.setVariableValue(solution.getVariable(vm_idx), pm_idx);
        }
        return solution;
    }
    
    protected Map<VirtualMachine, PhysicalMachine> allocationToMap(MappingVec allocation) {
        Map<VirtualMachine, PhysicalMachine> map = new HashMap<VirtualMachine, PhysicalMachine>();
        for (int i = 0; i < allocation.size(); ++i) {
            Mapping mapping = allocation.get(i);
            map.put(mapping.getVirtualMachine(), mapping.getPhysicalMachine());
        }
        return map;
    }

    protected long getRemainingTime() {
        if (this.timeout == NO_TIMEOUT) {
            return Long.MAX_VALUE; // FIXME: implicit limit, should change eventually
        }
        assert(this.timeout >= 0);
        return Math.max(this.timeout - (long)Clock.getInstance().getElapsed(), 0);
    }
    
    protected long getTimeout() { return this.timeout; }
    
    protected void saveSolution(Solution solution, boolean print_objectives) {
        this.instance.evaluate(solution);
        if (!solution.violatesConstraints()) {
            this.solutions.add(solution);
            boolean do_print = false;
            if (solution.getObjective(VMCwMProblem.ENERGY_OBJ_INDEX) < this.best_energy) {
                this.best_energy = solution.getObjective(VMCwMProblem.ENERGY_OBJ_INDEX);
                do_print = true;
            }
            if (solution.getObjective(VMCwMProblem.WASTAGE_OBJ_INDEX) < this.best_wastage) {
                this.best_wastage = solution.getObjective(VMCwMProblem.WASTAGE_OBJ_INDEX);
                do_print = true;
            }
            if (    solution.getNumberOfObjectives() == 3 &&
                    solution.getObjective(VMCwMProblem.MIGRATION_OBJ_INDEX) < this.best_mig_cost) {
                this.best_mig_cost = (long)solution.getObjective(VMCwMProblem.MIGRATION_OBJ_INDEX);
                do_print = true;
            }
            if (print_objectives && do_print) {
                if (solution.getNumberOfObjectives() == 3) {
                    System.out.printf("e %.5f \tw %.5f \tm %d\n",
                                      this.best_energy, this.best_wastage, (long)this.best_mig_cost);
                }
                else {
                    System.out.printf("e %.5f \tw %.5f\n", this.best_energy, this.best_wastage);
                }
                printElapsedTime();
            }
        }
    }
    protected void saveSolution(Solution solution) {
        saveSolution(solution, false);
    }
    protected void saveSolution(MappingVec allocation, boolean print_objectives) {
        saveSolution(allocationToSolution(allocation), print_objectives);
    }
    protected void saveSolution(MappingVec allocation) {
        saveSolution(allocation, false);
    }
    
    protected void clearSolutions() { this.solutions.clear(); }
    
    protected List<NondominatedPopulation> runMultipleSeeds(int nseeds) {
        System.out.println("c Running with " + nseeds + " different seeds");
        List<NondominatedPopulation> results = new LinkedList<NondominatedPopulation>();
        long timeout = this.timeout;
        for (int i = 0; i < nseeds; ++i) {
            if (timeout != NO_TIMEOUT) {
                setTimeout(timeout * (i+1));
            }
            allocate();
            results.add(this.solutions);
            this.solutions = new NondominatedPopulation();
        }
        setTimeout(timeout);
        return results;
    }

    protected void printElapsedTime() {
        System.out.println("c Elapsed time: " + Clock.getInstance().getElapsed() + " seconds");
    }

    // =============================================================================================
    // ======================================== PUBLIC API =========================================
    // =============================================================================================

    public static final long NO_TIMEOUT = -1;
    
    public AllocAlgorithm(VMCwMProblem instance) {
        this(instance, VMCwMProblem.Encoding.INTEGER);
    }
    protected AllocAlgorithm(VMCwMProblem instance, VMCwMProblem.Encoding encoding) {
        this.instance = instance;
        this.instance.setEncoding(encoding);
        this.solutions = new @Gen NondominatedPopulation();
        this.results = new @Gen LinkedList<NondominatedPopulation>();
        ProblemFactory.getInstance().addProvider(
                new @Gen VWCwMProblemProvider(this.instance.getPhysicalMachines(),
                                              this.instance.getJobs(),
                                              this.instance.getMappings(),
                                              this.instance.getMaxMigrationPercentile(),
                                              encoding));
    }

    public void setTimeout(long timeout) { this.timeout = timeout; }
    
    public boolean foundSolution() { return this.solutions.size() > 0; }
    
    public IVec<MappingVec> getSolutions() {
        IVec<MappingVec> allocations = new Vec<MappingVec>();
        for (int i = 0; i < this.solutions.size(); ++i) {
            allocations.push(solutionToAllocation(this.solutions.get(i)));
        }
        return allocations;
    }
    
    public void dumpPopulation(String path) throws IOException {
        System.out.println("c Dumping non-dominated population to " + path);
        Analyzer analyzer = new Analyzer().withProblem(VMCWM_PROBLEM);
        if (this.results == null || this.results.size() == 0 ) {
            assert(this.solutions != null);
            analyzer.add("VMCwM", this.solutions);
        }
        else {
            analyzer.addAll("VMCwM", this.results);
        }
        analyzer.saveAs("VMCwM", new File(path));
    }

    public void analyzePopulations(Map<String, String> dataset) throws IOException {
        Analyzer analyzer = new Analyzer()
                .withProblem(VMCWM_PROBLEM)
                .includeInvertedGenerationalDistance()
                .includeHypervolume()
                .showIndividualValues()
                .showStatisticalSignificance();
        for (Iterator<String> it = dataset.keySet().iterator(); it.hasNext();) {
            String label = it.next(), file_path = dataset.get(label);
            analyzer.loadAs(label, new File(file_path));
        }
        analyzer.printAnalysis();
    }
    
    public double[] getSolutionEnergyCosts() {
        return getObjectiveCosts(VMCwMProblem.ENERGY_OBJ_INDEX);
    }
    public double[] getSolutionWastageCosts() {
        return getObjectiveCosts(VMCwMProblem.WASTAGE_OBJ_INDEX);
    }
    public double[] getSolutionMigrationCosts() {
        return getObjectiveCosts(VMCwMProblem.MIGRATION_OBJ_INDEX);
    }

    public abstract void allocate();
    
    public void allocateMultipleSeeds(int nseeds) {
        this.results = runMultipleSeeds(nseeds);
        int seed_idx = 0;
        for (Iterator<NondominatedPopulation> it = this.results.iterator(); it.hasNext();) {
            System.out.println("c Population obtained with seed " + ++seed_idx);
            NondominatedPopulation population = it.next();
            for (int i = 0; i < population.size(); ++i) {
                Solution solution = population.get(i);
                if (!solution.violatesConstraints()) {
                    this.solutions.add(solution);
                    double energy = solution.getObjective(VMCwMProblem.ENERGY_OBJ_INDEX);
                    double wastage = solution.getObjective(VMCwMProblem.WASTAGE_OBJ_INDEX);
                    if (this.instance.getMappings().size() > 0) {
                        double migration = solution.getObjective(VMCwMProblem.MIGRATION_OBJ_INDEX);
                        System.out.printf("e %.5f \tw %.5f \tm %d\n", energy, wastage, (long)migration);
                    }
                    else {
                        System.out.printf("e %.5f \tw %.5f\n", energy, wastage);
                    }
                }
            }
        }
        System.out.println("c Done");
        printElapsedTime();
    }

}
