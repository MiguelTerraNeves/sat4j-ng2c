package vmalloc.algorithm;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Properties;
import java.util.Set;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Population;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.util.TypedProperties;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Utils;
import vmalloc.domain.CPUAggregator;
import vmalloc.domain.Mapping;
import vmalloc.domain.MappingVec;
import vmalloc.domain.MemoryAggregator;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.evolutionary.VMCwMProblem;

public class GGAAlloc extends EvolutionaryAllocAlgorithm {
    
    private class GGA extends AbstractEvolutionaryAlgorithm {
        
        private static final int NO_PM = -1;
        private static final double BETA = 0.5;
        
        private static final String FUZZY_EVAL_ATR = "e";
        private static final String GROUP_ATR = "g";
        
        private double mutation_rate;
        private double crossover_rate;
        
        // Auxiliary parameters
        private final double energy_lb;
        private final double energy_ub;
        private final double wastage_lb;
        private final double wastage_ub;
        private final BigInteger migration_lb;
        private final BigInteger migration_ub;
        
        // Auxiliary structures and statistics
        private IVecInt mapped_pm_idxs;
        private IVec<BigInteger> mapped_mems;
        private BigInteger total_mapped_mem;
        private IVec<IVecInt> pm_groups; // FIXME: Map would be more efficient, but should work fine for few distinct PM types
        private BinPackingAllocAlgorithm bp_alloc;
        
        public GGA(Problem problem,
                   Population population,
                   double mutation_rate,
                   double crossover_rate,
                   Initialization initialization) {
            super(problem, population, null, initialization);
            assert(mutation_rate >= 0.0 && mutation_rate <= 1.0);
            this.mutation_rate = mutation_rate;
            this.crossover_rate = crossover_rate;
            this.bp_alloc = new FirstFitDecreasingAlloc(instance);
            this.bp_alloc.setMaxMemoryMigrationPercentile(0.0);
            this.bp_alloc.enableVirtualMachineShuffling();
            // Compute lower bounds and upper bounds for fuzzy membership computation
            this.energy_lb = energyLowerBound();
            this.energy_ub = energyUpperBound();
            this.wastage_lb = wastageLowerBound();
            this.wastage_ub = wastageUpperBound();
            assert(this.energy_lb < this.energy_ub);
            assert(this.wastage_lb < this.wastage_ub);
            if (problem.getNumberOfObjectives() > 2) {
                this.migration_lb = migrationLowerBound();
                this.migration_ub = migrationUpperBound();
                assert(this.migration_lb.compareTo(this.migration_ub) < 0);
            }
            else {
                this.migration_lb = this.migration_ub = null;
            }
            // Initialize auxiliary structure with sets of VMs mapped to each PM if required
            if (problem.getNumberOfObjectives() > 2) {
                this.mapped_pm_idxs = new VecInt();
                for (int i = 0; i < instance.getVirtualMachines().size(); ++i) {
                    this.mapped_pm_idxs.push(NO_PM);
                }
                this.mapped_mems = new Vec<BigInteger>();
                this.total_mapped_mem = BigInteger.ZERO;
                for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
                    this.mapped_mems.push(BigInteger.ZERO);
                }
                for (int i = 0; i < instance.getMappings().size(); ++i) {
                    PhysicalMachine pm = instance.getMappings().get(i).getPhysicalMachine();
                    VirtualMachine vm = instance.getMappings().get(i).getVirtualMachine();
                    int vm_idx = instance.getVirtualMachineIndex(vm);
                    int pm_idx = instance.getPhysicalMachineIndex(pm);
                    assert(this.mapped_pm_idxs.get(vm_idx) == NO_PM);
                    this.mapped_pm_idxs.set(vm_idx, pm_idx);
                    this.mapped_mems.set(pm_idx, this.mapped_mems.get(pm_idx).add(vm.getMemory()));
                    this.total_mapped_mem = this.total_mapped_mem.add(vm.getMemory());
                }
            }
            // Split PM indexes into groups with equal resource capacities
            this.pm_groups = new Vec<IVecInt>();
            for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
                PhysicalMachine pm = instance.getPhysicalMachines().get(i);
                int j;
                for (j = 0; j < this.pm_groups.size(); ++j) {
                    PhysicalMachine rep_pm =
                            instance.getPhysicalMachines().get(this.pm_groups.get(j).get(0));
                    if (pm.getCPU().equals(rep_pm.getCPU()) && pm.getMemory().equals(rep_pm.getMemory())) {
                        break;
                    }
                }
                if (j >= this.pm_groups.size()) {
                    this.pm_groups.push(new VecInt());
                }
                assert(j < this.pm_groups.size());
                this.pm_groups.get(j).push(i);
            }
        }
        
        private BigDecimal computeEfficiency(BigInteger cap, int max_consumpt, int idle_consumpt) {
            return Utils.divideBigIntegers(cap,
                                           BigInteger.valueOf(max_consumpt - idle_consumpt),
                                           RoundingMode.HALF_DOWN);
        }
        private BigDecimal cpuEfficiencyOf(PhysicalMachine pm) {
            return computeEfficiency(pm.getCPU(), pm.getMaxConsumption(), pm.getIdleConsumption());
        }
        private BigDecimal memEfficiencyOf(PhysicalMachine pm) {
            return computeEfficiency(pm.getMemory(), pm.getMaxConsumption(), pm.getIdleConsumption());
        }
        
        private double energyLowerBound() {
            // Determine most CPU and memory efficient PMs
            PhysicalMachine cpu_efficient_pm = instance.getPhysicalMachines().get(0);
            PhysicalMachine mem_efficient_pm = instance.getPhysicalMachines().get(0);
            BigDecimal best_cpu_efficiency = cpuEfficiencyOf(cpu_efficient_pm);
            BigDecimal best_mem_efficiency = memEfficiencyOf(mem_efficient_pm);
            for (int i = 1; i < instance.getPhysicalMachines().size(); ++i) {
                PhysicalMachine pm = instance.getPhysicalMachines().get(i);
                BigDecimal cpu_efficiency = cpuEfficiencyOf(pm), mem_efficiency = memEfficiencyOf(pm);
                if (cpu_efficiency.compareTo(best_cpu_efficiency) > 0) {
                    best_cpu_efficiency = cpu_efficiency;
                    cpu_efficient_pm = pm;
                }
                if (mem_efficiency.compareTo(best_mem_efficiency) > 0) {
                    best_mem_efficiency = mem_efficiency;
                    mem_efficient_pm = pm;
                }
            }
            // Sort PMs by increasing order of idle energy consumption
            PhysicalMachine[] pm_array = new PhysicalMachine[instance.getPhysicalMachines().size()];
            instance.getPhysicalMachines().copyTo(pm_array);
            Arrays.sort(pm_array, new Comparator<PhysicalMachine>() {
                @Override
                public int compare(PhysicalMachine pm1, PhysicalMachine pm2) {
                    return pm1.getIdleConsumption() - pm2.getIdleConsumption();
                }
            });
            // Compute CPU based LB
            long min_cpu_pms = Utils.floor(Utils.divideBigIntegers(instance.getTotalCPURequirements(),
                                                                   cpu_efficient_pm.getCPU(),
                                                                   RoundingMode.FLOOR)).longValueExact();
            int cpu_consumpt_diff = cpu_efficient_pm.getMaxConsumption() -
                                    cpu_efficient_pm.getIdleConsumption();
            double cpu_lb = min_cpu_pms * (double)cpu_consumpt_diff;
            for (int i = 0; i < min_cpu_pms; ++i) {
                cpu_lb += pm_array[i].getIdleConsumption();
            }
            // Compute memory based LB
            long min_mem_pms = Utils.floor(Utils.divideBigIntegers(instance.getTotalMemoryRequirements(),
                                                                   mem_efficient_pm.getMemory(),
                                                                   RoundingMode.FLOOR)).longValueExact();
            int mem_consumpt_diff = mem_efficient_pm.getMaxConsumption() -
                                    mem_efficient_pm.getIdleConsumption();
            double mem_lb = min_mem_pms * (double)mem_consumpt_diff;
            for (int i = 0; i < min_mem_pms; ++i) {
                mem_lb += pm_array[i].getIdleConsumption();
            }
            // Return highest LB
            return Math.max(cpu_lb, mem_lb);
        }
        
        // FIXME: similar to energyLowerBound
        private double energyUpperBound() {
            // Determine most CPU and memory inefficient PMs
            PhysicalMachine cpu_ineff_pm = instance.getPhysicalMachines().get(0);
            PhysicalMachine mem_ineff_pm = instance.getPhysicalMachines().get(0);
            BigDecimal worst_cpu_eff = cpuEfficiencyOf(cpu_ineff_pm);
            BigDecimal worst_mem_eff = memEfficiencyOf(mem_ineff_pm);
            for (int i = 1; i < instance.getPhysicalMachines().size(); ++i) {
                PhysicalMachine pm = instance.getPhysicalMachines().get(i);
                BigDecimal cpu_efficiency = cpuEfficiencyOf(pm), mem_efficiency = memEfficiencyOf(pm);
                if (cpu_efficiency.compareTo(worst_cpu_eff) < 0) {
                    worst_cpu_eff = cpu_efficiency;
                    cpu_ineff_pm = pm;
                }
                if (mem_efficiency.compareTo(worst_mem_eff) < 0) {
                    worst_mem_eff = mem_efficiency;
                    mem_ineff_pm = pm;
                }
            }
            // Sort PMs by decreasing order of idle energy consumption
            PhysicalMachine[] pm_array = new PhysicalMachine[instance.getPhysicalMachines().size()];
            instance.getPhysicalMachines().copyTo(pm_array);
            Arrays.sort(pm_array, new Comparator<PhysicalMachine>() {
                @Override
                public int compare(PhysicalMachine pm1, PhysicalMachine pm2) {
                    return pm2.getIdleConsumption() - pm1.getIdleConsumption();
                }
            });
            // Compute CPU based UB
            long max_cpu_pms = Utils.ceil(
                    Utils.divideBigIntegers(instance.getTotalCPURequirements(),
                                            cpu_ineff_pm.getCPU(),
                                            RoundingMode.CEILING)).longValueExact();
            max_cpu_pms = Math.min(max_cpu_pms, instance.getPhysicalMachines().size());
            int cpu_consumpt_diff = cpu_ineff_pm.getMaxConsumption() - cpu_ineff_pm.getIdleConsumption(); 
            double cpu_ub = max_cpu_pms * (double)cpu_consumpt_diff;
            for (int i = 0; i < max_cpu_pms; ++i) {
                cpu_ub += pm_array[i].getIdleConsumption();
            }
            // Compute memory based UB
            long max_mem_pms =Utils.ceil(
                    Utils.divideBigIntegers(instance.getTotalMemoryRequirements(),
                                            mem_ineff_pm.getMemory(),
                                            RoundingMode.CEILING)).longValueExact();
            max_mem_pms = Math.min(max_mem_pms, instance.getPhysicalMachines().size());
            int mem_consumpt_diff = mem_ineff_pm.getMaxConsumption() - mem_ineff_pm.getIdleConsumption();
            double mem_ub = max_mem_pms * (double)mem_consumpt_diff;
            for (int i = 0; i < max_mem_pms; ++i) {
                mem_ub += pm_array[i].getIdleConsumption();
            }
            return Math.min(Math.min(cpu_ub, mem_ub),
                            Utils.sum(instance.getPhysicalMachines().getMaxConsumptions()));
        }
        
        private double[] normalizeCPURequirements(PhysicalMachine pm) {
            double[] cpu_reqs = new double[instance.getVirtualMachines().size()];
            for (int i = 0; i < instance.getVirtualMachines().size(); ++i) {
                VirtualMachine vm = instance.getVirtualMachines().get(i);
                cpu_reqs[i] = Utils.divideBigIntegers(vm.getCPU(),
                                                      pm.getCPU(),
                                                      RoundingMode.HALF_EVEN).doubleValue();
            }
            return cpu_reqs;
        }
        
        // FIXME: very similar to normalizeCPURequirements
        private double[] normalizeMemoryRequirements(PhysicalMachine pm) {
            double[] mem_reqs = new double[instance.getVirtualMachines().size()];
            for (int i = 0; i < instance.getVirtualMachines().size(); ++i) {
                VirtualMachine vm = instance.getVirtualMachines().get(i);
                mem_reqs[i] = Utils.divideBigIntegers(vm.getMemory(),
                                                      pm.getMemory(),
                                                      RoundingMode.HALF_EVEN).doubleValue();
            }
            return mem_reqs;
        }
        
        private double wastageLowerBound() {
            // Select PM that can accommodate the VMs in the most balanced way
            BigInteger cpu_req_sum = instance.getTotalCPURequirements();
            BigInteger mem_req_sum = instance.getTotalMemoryRequirements();
            PhysicalMachine best_pm = instance.getPhysicalMachines().get(0);
            BigDecimal best_alloc_diff = Utils.divideBigIntegers(cpu_req_sum,
                                                                 best_pm.getCPU(),
                                                                 RoundingMode.HALF_EVEN).subtract(
                                         Utils.divideBigIntegers(mem_req_sum,
                                                                 best_pm.getMemory(),
                                                                 RoundingMode.HALF_EVEN)).abs();
            for (int i = 1; i < instance.getPhysicalMachines().size(); ++i) {
                PhysicalMachine pm = instance.getPhysicalMachines().get(i);
                BigDecimal new_alloc_diff = Utils.divideBigIntegers(cpu_req_sum,
                                                                    pm.getCPU(),
                                                                    RoundingMode.HALF_EVEN).subtract(
                                            Utils.divideBigIntegers(mem_req_sum,
                                                                    pm.getMemory(),
                                                                    RoundingMode.HALF_EVEN)).abs();
                if (new_alloc_diff.compareTo(best_alloc_diff) < 0) {
                    best_pm = pm;
                    best_alloc_diff = new_alloc_diff;
                }
            }
            // Normalize requirements and compute LB
            double[] cpu_reqs = normalizeCPURequirements(best_pm);
            double[] mem_reqs = normalizeMemoryRequirements(best_pm);
            return Math.abs(Utils.sum(cpu_reqs) - Utils.sum(mem_reqs));
        }
        
        // FIXME: very similar to wastageLowerBound
        private double wastageUpperBound() {
            // Select PM that can accommodate the VMs in the least balanced way
            BigInteger cpu_req_sum = instance.getTotalCPURequirements();
            BigInteger mem_req_sum = instance.getTotalMemoryRequirements();
            PhysicalMachine worst_pm = instance.getPhysicalMachines().get(0);
            BigDecimal worst_alloc_diff = Utils.divideBigIntegers(cpu_req_sum,
                                                                  worst_pm.getCPU(),
                                                                  RoundingMode.HALF_EVEN).subtract(
                                          Utils.divideBigIntegers(mem_req_sum,
                                                                  worst_pm.getMemory(),
                                                                  RoundingMode.HALF_EVEN)).abs();
            for (int i = 1; i < instance.getPhysicalMachines().size(); ++i) {
                PhysicalMachine pm = instance.getPhysicalMachines().get(i);
                BigDecimal new_alloc_diff = Utils.divideBigIntegers(cpu_req_sum,
                                                                    pm.getCPU(),
                                                                    RoundingMode.HALF_EVEN).subtract(
                                            Utils.divideBigIntegers(mem_req_sum,
                                                                    pm.getMemory(),
                                                                    RoundingMode.HALF_EVEN)).abs();
                if (new_alloc_diff.compareTo(worst_alloc_diff) > 0) {
                    worst_pm = pm;
                    worst_alloc_diff = new_alloc_diff;
                }
            }
            // Normalize requirements and compute LB
            double[] cpu_reqs = normalizeCPURequirements(worst_pm);
            double[] mem_reqs = normalizeMemoryRequirements(worst_pm);
            double ub = 0.0;
            for (int i = 0; i < instance.getVirtualMachines().size(); ++i) {
                ub += Math.abs(cpu_reqs[i] - mem_reqs[i]);
            }
            return ub;
        }
        
        private BigInteger migrationLowerBound() { return BigInteger.ZERO; }
        
        private BigInteger migrationUpperBound() {
            VirtualMachineVec mapped_vms = new VirtualMachineVec();
            for (int i = 0; i < instance.getMappings().size(); ++i) {
                mapped_vms.push(instance.getMappings().get(i).getVirtualMachine());
            }
            MemoryAggregator mem_agr = new MemoryAggregator();
            mem_agr.processVirtualMachines(mapped_vms);
            return mem_agr.memorySum();
        }
        
        // FIXME: only supports integer encoding
        private int getValue(Variable var) {
            return EncodingUtils.getInt(var);
        }
        
        private class Group {
            
            private IVecInt var_idxs;
            private int group_idx;
            private double group_eval;
            
            Group(int group_idx, IVecInt var_idxs) {
                this.var_idxs = var_idxs;
                this.group_idx = group_idx;
                // Compute energy efficiency
                PhysicalMachine group_pm = instance.getPhysicalMachines().get(group_idx);
                VirtualMachineVec group_vms = new VirtualMachineVec();
                for (int i = 0; i < var_idxs.size(); ++i) {
                    group_vms.push(instance.getVirtualMachines().get(i));
                }
                CPUAggregator cpu_agr = new CPUAggregator();
                group_vms.accept(cpu_agr);
                double cpu_usage = Utils.divideBigIntegers(cpu_agr.cpuSum(),
                                                           group_pm.getCPU(),
                                                           RoundingMode.HALF_EVEN).doubleValue();
                double energy_efficiency = cpu_usage*group_pm.getMaxConsumption();
                energy_efficiency /= cpu_usage * (group_pm.getMaxConsumption() -
                                                  group_pm.getIdleConsumption()) +
                                     group_pm.getIdleConsumption();
                // Compute resource usage efficiency
                MemoryAggregator mem_agr = new MemoryAggregator();
                group_vms.accept(mem_agr);
                double mem_usage = Utils.divideBigIntegers(mem_agr.memorySum(),
                                                           group_pm.getMemory(),
                                                           RoundingMode.HALF_EVEN).doubleValue();
                double usage_efficiency = cpu_usage * mem_usage;
                this.group_eval = energy_efficiency * usage_efficiency;
                // Compute migration efficiency if required
                if (problem.getNumberOfObjectives() > 2) {
                    BigInteger inside_migged_mem = BigInteger.ZERO;
                    for (int i = 0; i < var_idxs.size(); ++i) {
                        int var_idx = var_idxs.get(i), mapped_pm_idx = mapped_pm_idxs.get(var_idx);
                        if (mapped_pm_idx != NO_PM && mapped_pm_idx != group_idx) {
                            inside_migged_mem = inside_migged_mem.add(group_vms.get(i).getMemory());
                        }
                    }
                    BigInteger mapped_outside_mem =
                            total_mapped_mem.subtract(mapped_mems.get(group_idx));
                    double mig_efficiency =
                            Utils.divideBigIntegers(mapped_outside_mem.subtract(inside_migged_mem),
                                                    mapped_outside_mem,
                                                    RoundingMode.HALF_EVEN).doubleValue();
                    assert(mig_efficiency >= 0.0 && mig_efficiency <= 1.0);
                    this.group_eval *= mig_efficiency;
                }
            }
            
            public int getVariableIndex(int i) { return this.var_idxs.get(i); }
            public int getGroupSize() { return this.var_idxs.size(); }
            public int getGroupIndex() { return this.group_idx; }
            public double getEvaluation() { return this.group_eval; }
            
        }
        
        private Group[] computeGroups(Solution parent) {
            IVec<IVecInt> var_idxs_per_group = new Vec<IVecInt>();
            for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
                var_idxs_per_group.push(new VecInt());
            }
            for (int i = 0; i < parent.getNumberOfVariables(); ++i) {
                int pm_idx = getValue(parent.getVariable(i));
                var_idxs_per_group.get(pm_idx).push(i);
            }
            IVec<Group> groups = new Vec<Group>();
            for (int i = 0; i < var_idxs_per_group.size(); ++i) {
                IVecInt var_idxs = var_idxs_per_group.get(i);
                if (var_idxs.size() > 0) {
                    groups.push(new Group(i, var_idxs));
                }
            }
            return groupsToArray(groups);
        }
        
        private Solution binPackingResultToSolution(BinPackingAllocAlgorithm bp) {
            MappingVec allocation = null;
            VirtualMachineVec leftover_vms = null;
            if (bp.foundSolution()) {
                allocation = bp.getSolution();
                leftover_vms = new VirtualMachineVec();
            }
            else {
                allocation = bp.getPartialAllocation();
                leftover_vms = bp.getLeftoverVirtualMachines();
            }
            assert(allocation.size() + leftover_vms.size() == instance.getVirtualMachines().size());
            Solution solution = problem.newSolution();
            injectBinPackingResultInSolution(solution, allocation, leftover_vms);
            return solution;
        }
        
        private Solution mutate(Solution parent) {
            Group[] groups = getGroups(parent);
            int del_idx = PRNG.nextInt(groups.length);
            MappingVec placement = new MappingVec();
            for (int i = 0; i < groups.length; ++i) {
                if (i != del_idx) {
                    Group group = groups[i];
                    for (int j = 0; j < group.getGroupSize(); ++j) {
                        placement.push(
                                new Mapping(instance.getVirtualMachines().get(group.getVariableIndex(j)),
                                            instance.getPhysicalMachines().get(group.getGroupIndex())));
                    }
                }
            }
            /*FirstFitDecreasingAlloc first_fit = new FirstFitDecreasingAlloc(pms, jobs, placement);
            first_fit.enableVirtualMachineShuffling();
            first_fit.setMaxMemoryMigrationPercentile(0.0);*/
            this.bp_alloc.setMappings(placement);
            Utils.stdoutDisable();
            this.bp_alloc.allocate();
            Utils.stdoutEnable();
            return binPackingResultToSolution(this.bp_alloc);
        }
        
        private Group[] groupsToArray(IVec<Group> groups) {
            Group[] group_array = new Group[groups.size()];
            for (int i = 0; i < groups.size(); ++i) {
                group_array[i] = groups.get(i);
            }
            return group_array;
        }
        
        private Solution rankingCrossover(Solution parent1, Solution parent2) {
            // Sort group vectors in decreasing order of overall efficiency
            Group[] groups1 = getGroups(parent1), groups2 = getGroups(parent2);
            Comparator<Group> comp = new Comparator<Group>() {
                @Override
                public int compare(Group g1, Group g2) {
                    if (g1.getEvaluation() > g2.getEvaluation()) {
                        return 1;
                    }
                    else if (g1.getEvaluation() < g2.getEvaluation()) {
                        return -1;
                    }
                    return 0;
                }
            };
            Arrays.sort(groups1, comp);
            Arrays.sort(groups2, comp);
            // Perform group crossover in decreasing order of efficiency
            Set<Integer> placed_var_idxs = new HashSet<Integer>();
            Set<Integer> used_group_idxs = new HashSet<Integer>();
            IVec<Group> offspring_groups = new Vec<Group>();
            IVecInt leftover_var_idxs = new VecInt();
            for (int i = 0, j = 0; i < groups1.length && j < groups2.length;) {
                // Select next group to inject into offspring
                Group group_to_inject = null;
                if (    i < groups1.length &&
                        (j >= groups2.length ||
                         groups1[i].getEvaluation() > groups2[j].getEvaluation())) {
                    group_to_inject = groups1[i++];
                }
                else {
                    group_to_inject = groups2[j++];
                }
                // Inject group
                int new_group_idx = group_to_inject.getGroupIndex();
                if (used_group_idxs.contains(new Integer(new_group_idx))) {
                    // Retrieve PM group with equal capacities
                    IVecInt pm_group = null;
                    PhysicalMachine pm = instance.getPhysicalMachines().get(new_group_idx);
                    for (int k = 0; pm_group == null; ++k) {
                        PhysicalMachine rep_pm =
                                instance.getPhysicalMachines().get(this.pm_groups.get(i).get(0));
                        if (    pm.getCPU().equals(rep_pm.getCPU()) &&
                                pm.getMemory().equals(rep_pm.getMemory())) {
                            pm_group = this.pm_groups.get(k);
                        }
                    }
                    // Select an equal empty PM
                    int k;
                    for (k = 0; k < pm_group.size(); ++k) {
                        int candidate_idx = pm_group.get(k);
                        if (!used_group_idxs.contains(new Integer(candidate_idx))) {
                            new_group_idx = candidate_idx;
                            break;
                        }
                    }
                    // If no equal empty PM was found, store VMs to place afterwards with First-Fit
                    if (k >= pm_group.size()) {
                        for (k = 0; k < group_to_inject.getGroupSize(); ++k) {
                            leftover_var_idxs.push(group_to_inject.getVariableIndex(k));
                        }
                    }
                }
                assert(!used_group_idxs.contains(new Integer(new_group_idx)));
                IVecInt new_var_idxs = new VecInt();
                for (int k = 0; k < group_to_inject.getGroupSize(); ++k) {
                    int var_idx = group_to_inject.getVariableIndex(k);
                    if (!placed_var_idxs.contains(var_idx)) {
                        new_var_idxs.push(var_idx);
                    }
                }
                if (new_var_idxs.size() > 0) {
                    offspring_groups.push(new Group(new_group_idx, new_var_idxs));
                    used_group_idxs.add(new Integer(new_group_idx));
                }
            }
            // Insert leftover VMs using first fit
            MappingVec placement = new MappingVec();
            for (int i = 0; i < offspring_groups.size(); ++i) {
                Group group = offspring_groups.get(i);
                for (int j = 0; j < group.getGroupSize(); ++j) {
                    placement.push(
                            new Mapping(instance.getVirtualMachines().get(group.getVariableIndex(j)),
                                        instance.getPhysicalMachines().get(group.getGroupIndex())));
                }
            }
            /*FirstFitDecreasingAlloc first_fit = new FirstFitDecreasingAlloc(pms, jobs, placement);
            first_fit.setMaxMemoryMigrationPercentile(0.0);
            first_fit.enableVirtualMachineShuffling();*/
            this.bp_alloc.setMappings(placement);
            Utils.stdoutDisable();
            this.bp_alloc.allocate();
            Utils.stdoutEnable();
            return binPackingResultToSolution(this.bp_alloc);
        }
        
        private double lowEnergyMembership(Solution sol) {
            double energy_cost = Math.min(sol.getObjective(VMCwMProblem.ENERGY_OBJ_INDEX), this.energy_ub);
            assert(energy_cost >= this.energy_lb); // energy_cost may be higher than UB if constraints are violated
            return 1.0 - ((energy_cost - this.energy_lb) / (this.energy_ub - this.energy_lb));
        }
        
        private double lowWastageMembership(Solution sol) {
            BigInteger[] used_cpu_caps = new BigInteger[instance.getPhysicalMachines().size()];
            BigInteger[] used_mem_caps = new BigInteger[instance.getPhysicalMachines().size()];
            for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
                used_cpu_caps[i] = BigInteger.ZERO;
                used_mem_caps[i] = BigInteger.ZERO;
            }
            int[] x = instance.getVariableAssignment(sol);
            for (int i = 0; i < instance.getVirtualMachines().size(); ++i) {
                int pm_idx = x[i];
                used_cpu_caps[pm_idx] =
                        used_cpu_caps[pm_idx].add(instance.getVirtualMachines().get(i).getCPU());
                used_mem_caps[pm_idx] =
                        used_mem_caps[pm_idx].add(instance.getVirtualMachines().get(i).getMemory());
            }
            double total_wastage = 0.0;
            for (int i = 0; i < instance.getPhysicalMachines().size(); ++i) {
                if (    !used_cpu_caps[i].equals(BigInteger.ZERO) && 
                        !used_mem_caps[i].equals(BigInteger.ZERO)) {
                    double used_cpu = used_cpu_caps[i].doubleValue()/
                                      instance.getPhysicalMachines().get(i).getCPU().doubleValue();
                    double used_mem = used_mem_caps[i].doubleValue()/
                                      instance.getPhysicalMachines().get(i).getMemory().doubleValue();
                    double cpu_left = 1.0 - used_cpu;
                    double mem_left = 1.0 - used_mem;
                    total_wastage += Math.abs(cpu_left-mem_left);
                }
            }
            total_wastage = Math.min(total_wastage, this.wastage_ub);
            assert(total_wastage >= this.wastage_lb); // wastage may be higher than UB if constraints are violated
            return 1.0 - ((total_wastage - this.wastage_lb) / (this.wastage_ub - this.wastage_lb));
        }
        
        private double fewMigrationsMembership(Solution sol) {
            BigDecimal mig_cost =
                    new BigDecimal(sol.getObjective(VMCwMProblem.MIGRATION_OBJ_INDEX));
            assert(mig_cost.compareTo(new BigDecimal(this.migration_lb)) >= 0 &&
                   mig_cost.compareTo(new BigDecimal(this.migration_ub)) <= 0);
            BigDecimal mig_lb_diff = mig_cost.subtract(new BigDecimal(this.migration_lb));
            BigDecimal mig_range = new BigDecimal(this.migration_ub.subtract(this.migration_lb));
            return 1.0 - Utils.divideBigDecimals(mig_lb_diff,
                                                 mig_range,
                                                 RoundingMode.HALF_EVEN).doubleValue();
        }
        
        private double fuzzyEvaluation(Solution sol) {
            double[] membership_values = (problem.getNumberOfObjectives() > 2) ? new double[3] :
                                                                                 new double[2];
            membership_values[0] = lowEnergyMembership(sol);
            membership_values[1] = lowWastageMembership(sol);
            if (problem.getNumberOfObjectives() > 2) {
                membership_values[2] = fewMigrationsMembership(sol);
            }
            return BETA * Utils.min(membership_values) + (1.0 - BETA) * Utils.avg(membership_values);
        }
        
        private void ensureFuzzyEvaluation(Solution sol) {
            if (!sol.hasAttribute(FUZZY_EVAL_ATR)) {
                sol.setAttribute(FUZZY_EVAL_ATR, new Double(fuzzyEvaluation(sol)));
            }
        }
        private double getFuzzyEvaluation(Solution sol) {
            ensureFuzzyEvaluation(sol);
            return ((Double)sol.getAttribute(FUZZY_EVAL_ATR)).doubleValue();
        }
        
        private void updatePopulation(Population offspring) {
            evaluateAll(offspring);
            offspring.sort(new Comparator<Solution>() {
                @Override
                public int compare(Solution sol1, Solution sol2) {
                    double violation1 = instance.getConstraintViolation(sol1);
                    double violation2 = instance.getConstraintViolation(sol2);
                    if (violation1 != violation2) {
                        return (violation1 > violation2) ? 1 : -1;
                    }
                    double eval1 = getFuzzyEvaluation(sol1), eval2 = getFuzzyEvaluation(sol2);
                    return Double.compare(eval2, eval1);
                }
            });
            Population next_gen = new Population();
            for (int i = 0; i < offspring.size(); ++i) {
                Solution offspring_sol = offspring.get(i);
                if (    next_gen.size() == 0 ||
                        getFuzzyEvaluation(next_gen.get(next_gen.size()-1)) !=
                                getFuzzyEvaluation(offspring_sol) ||
                        !solutionEquals(next_gen.get(next_gen.size()-1), offspring.get(i))) {
                    next_gen.add(offspring.get(i));
                }
            }
            assert(next_gen.size() >= this.population.size());
            next_gen.truncate(this.population.size(), new Comparator<Solution>() {
                @Override
                public int compare(Solution o1, Solution o2) {
                    return 0; // dummy method
                }
            });
            assert(next_gen.size() == this.population.size());
            this.population.clear();
            this.population.addAll(next_gen);
        }
        
        private void ensureGroups(Solution... solutions) {
            for (int i = 0; i < solutions.length; ++i) {
                Solution sol = solutions[i];
                if (!sol.hasAttribute(GROUP_ATR)) {
                    sol.setAttribute(GROUP_ATR, computeGroups(sol));
                }
            }
        }
        
        private Group[] getGroups(Solution sol) {
            ensureGroups(sol);
            return (Group[])sol.getAttribute(GROUP_ATR);
        }
        
        private void dumpGroupIndexes(IVec<Group> groups) {
            for (int i = 0; i < groups.size(); ++i) {
                System.out.print(groups.get(i).getGroupIndex() + " ");
            }
        }
        
        private void printGroupIndexes(String label, IVec<Group> groups) {
            System.out.print("c Group indexes for " + label + ": ");
            dumpGroupIndexes(groups);
            System.out.print("\n");
        }
        
        @Override
        protected void initialize() {
            initialized = true;
            Solution[] initial_solutions = initialization.initialize();
            evaluateAll(initial_solutions);
            population.addAll(initial_solutions);
        }

        @Override
        protected void iterate() {
            Population offspring = new Population();
            offspring.addAll(this.population);
            // Perform ranking crossover
            for (int i = 0; i < crossover_rate*this.population.size(); ++i) {
                int first_idx = PRNG.nextInt(this.population.size()), second_idx;
                for (second_idx = first_idx;
                     second_idx == first_idx;
                     second_idx = PRNG.nextInt(this.population.size()));
                ensureGroups(this.population.get(first_idx), this.population.get(second_idx));
                offspring.add(rankingCrossover(population.get(first_idx), population.get(second_idx)));
            }
            // Perform mutation
            for (int i = 0; i < mutation_rate*this.population.size(); ++i) {
                int idx = PRNG.nextInt(this.population.size());
                ensureGroups(this.population.get(idx));
                offspring.add(mutate(this.population.get(idx)));
            }
            updatePopulation(offspring);
        }
        
        @Override
        public Population getPopulation() {
            Population ret_pop = new Population();
            ret_pop.addAll(population);
            return ret_pop;
        }
        
        @Override
        public NondominatedPopulation getResult() {
            return new NondominatedPopulation(getPopulation());
        }
        
    }
    
    private class GGAAllocAlgorithmProvider extends AlgorithmProvider {
        
        @Override
        public Algorithm getAlgorithm(String name, Properties properties, Problem problem) {
            if (name.equals("GGA")) {
                TypedProperties typed_props = new TypedProperties(properties);
                int pop_size = typed_props.getInt("populationSize", 100);
                Population population = new Population();
                Initialization initialization = makeInitializer(problem, pop_size);
                double mutation_rate = typed_props.getDouble("gga.mutationRate", 0.0);
                double crossover_rate = typed_props.getDouble("gga.crossoverRate", 0.8);
                return new GGA(problem, population, mutation_rate, crossover_rate, initialization);
            }
            return null;
        }
        
    }

    public GGAAlloc(VMCwMProblem instance) {
        super(instance, "GGA", VMCwMProblem.Encoding.INTEGER);
        AlgorithmFactory.getInstance().addProvider(new GGAAllocAlgorithmProvider());
    }
    
    // FIXME: note that crossover and mutation rate for GGA has different meaning than for DE and GA
    public void setCrossoverRate(double rate) { exec = exec.withProperty("gga.crossoverRate", rate); }
    public void setMutationRate(double rate) { exec = exec.withProperty("gga.mutationRate", rate); }

}
