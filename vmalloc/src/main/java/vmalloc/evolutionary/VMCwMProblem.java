package vmalloc.evolutionary;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.AbstractProblem;
import org.sat4j.core.Vec;
import org.sat4j.specs.IVec;

import vmalloc.Utils;
import vmalloc.domain.CPUAggregator;
import vmalloc.domain.Job;
import vmalloc.domain.JobVec;
import vmalloc.domain.MappingVec;
import vmalloc.domain.MemoryAggregator;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.exception.InvalidEncodingException;

public class VMCwMProblem extends AbstractProblem {
    
    private abstract class VariableEncoding {
        abstract Variable makeVariable();
        abstract void setVariableValue(Variable x, int val);
        abstract int[] getVariableAssignment(Solution solution);
    }
    
    private class IntegerEncoding extends VariableEncoding {
        @Override
        Variable makeVariable() { return EncodingUtils.newInt(0, pms.size()-1); }
        @Override
        void setVariableValue(Variable x, int val) { EncodingUtils.setInt(x, val); }
        @Override
        int[] getVariableAssignment(Solution solution) { return EncodingUtils.getInt(solution); }
    }
    
    private class BinaryIntegerEncoding extends VariableEncoding {
        @Override
        Variable makeVariable() { return EncodingUtils.newBinaryInt(0, pms.size()-1); }
        @Override
        void setVariableValue(Variable x, int val) { EncodingUtils.encode(val, (BinaryVariable)x); }
        @Override
        int[] getVariableAssignment(Solution solution) {
            int[] x = new int[vms.size()];
            for (int i = 0; i < vms.size(); ++i) {
                x[i] = (int)EncodingUtils.decode((BinaryVariable)solution.getVariable(i));
            }
            return x;
        }
    }
    
    private static final double WASTAGE_EPSILON = 0.0001;
    
    public static final int ENERGY_OBJ_INDEX = 0;
    public static final int WASTAGE_OBJ_INDEX = 1;
    public static final int MIGRATION_OBJ_INDEX = 2;
    
    public static final String VIOLATION_ATR = "c";
    
    private final PhysicalMachineVec pms;
    private final JobVec jobs;
    private final VirtualMachineVec vms;
    private final MappingVec mappings;
    private double max_mig_percentile = 1.0;

    // Auxiliary stats and structures
    private final BigInteger total_cpu_req;
    private final BigInteger total_mem_req;
    private final BigInteger total_cpu_cap;
    private final BigInteger total_mem_cap;
    private final IVec<VirtualMachineVec> anti_coloc_vms;
    private final VirtualMachineVec plat_constrained_vms;
    private final Map<String, Integer> vm_id_to_idx;
    private final Map<Integer, Integer> pm_id_to_idx;
    
    private VariableEncoding encoding;
    
    private VariableEncoding makeEncoding(Encoding encoding) {
        if (encoding == Encoding.INTEGER) {
            return new @Gen IntegerEncoding();
        }
        else if (encoding == Encoding.BINARY_INTEGER) {
            return new @Gen BinaryIntegerEncoding();
        }
        throw new InvalidEncodingException();
    }
    
    private Variable makeVariable() { return encoding.makeVariable(); }
    
    private static int calculateNVariables(VirtualMachineVec vms) { return vms.size(); }
    private static int calculateNObjectives(boolean mapping_exists) {
        return mapping_exists ? 3 : 2;
    }
    private static int calculateNConstraints(PhysicalMachineVec pms,
                                             JobVec jobs,
                                             boolean mapping_exists) {
        return 2*pms.size() +
               retrieveAntiColocatableVirtualMachines(jobs).size() +
               retrievePlatformConstrainedVirtualMachines(jobs.flattenJobs()).size() +
               (mapping_exists ? 1 : 0);
    }

    private static BigInteger cpuCapacitySum(PhysicalMachineVec pms) {
        CPUAggregator cpu_agr = new CPUAggregator();
        cpu_agr.processPhysicalMachines(pms);
        return cpu_agr.cpuSum();
    }
    
    private static BigInteger memoryCapacitySum(PhysicalMachineVec pms) {
        MemoryAggregator mem_agr = new MemoryAggregator();
        mem_agr.processPhysicalMachines(pms);
        return mem_agr.memorySum();
    }
    
    private static BigInteger cpuRequirementSum(VirtualMachineVec vms) {
        CPUAggregator cpu_agr = new CPUAggregator();
        cpu_agr.processVirtualMachines(vms);
        return cpu_agr.cpuSum();
    }
    
    private static BigInteger memoryRequirementSum(VirtualMachineVec vms) {
        MemoryAggregator mem_agr = new MemoryAggregator();
        mem_agr.processVirtualMachines(vms);
        return mem_agr.memorySum();
    }
    
    private static IVec<VirtualMachineVec> retrieveAntiColocatableVirtualMachines(JobVec jobs) {
        IVec<VirtualMachineVec> anti_coloc_vms = new @Gen Vec<VirtualMachineVec>();
        for (int i = 0; i < jobs.size(); ++i) {
            Job job = jobs.get(i);
            VirtualMachineVec job_anti_coloc_vms = job.getAntiColocatableVirtualMachines();
            if (job_anti_coloc_vms.size() > 1) {
                anti_coloc_vms.push(job_anti_coloc_vms);
            }
        }
        return anti_coloc_vms;
    }
    
    private static VirtualMachineVec retrievePlatformConstrainedVirtualMachines(VirtualMachineVec vms) {
        VirtualMachineVec plat_constrained_vms = new @Gen VirtualMachineVec();
        for (int i = 0; i < vms.size(); ++i) {
            if (vms.get(i).getUnallowedPhysicalMachines().size() > 0) {
                plat_constrained_vms.push(vms.get(i));
            }
        }
        return plat_constrained_vms;
    }
    
    private int cpuCapacityConstraintIndex(int pm_idx) { return 2*pm_idx; }
    private int memoryCapacityConstraintIndex(int pm_idx) { return 2*pm_idx + 1; }
    private int antiColocationConstraintIndex(int ac_idx) { return 2*this.pms.size() + ac_idx; }
    private int platformConstraintIndex(int pc_idx) {
        return 2*this.pms.size() + this.anti_coloc_vms.size() + pc_idx;
    }
    private int migrationConstraintIndex() {
        return 2*this.pms.size() + this.anti_coloc_vms.size() + this.plat_constrained_vms.size();
    }
    
    private void ensureConstraintViolation(Solution sol) {
        if (!sol.hasAttribute(VIOLATION_ATR)) {
            sol.setAttribute(VIOLATION_ATR, new @Gen Double(Utils.sum(sol.getConstraints())));
        }
    }
    
    public static enum Encoding {
        NONE, INTEGER, BINARY_INTEGER;
    }
    
    public VMCwMProblem(PhysicalMachineVec pms, JobVec jobs, MappingVec mappings) {
        this(pms, jobs, mappings, 1.0);
    }
    public VMCwMProblem(PhysicalMachineVec pms, JobVec jobs, MappingVec mappings, Encoding encoding) {
        this(pms, jobs, mappings, 1.0, encoding);
    }
    public VMCwMProblem(PhysicalMachineVec pms,
                        JobVec jobs,
                        MappingVec mappings,
                        double max_mig_percentile) {
        this(pms, jobs, mappings, max_mig_percentile, Encoding.INTEGER);
    }
    public VMCwMProblem(PhysicalMachineVec pms,
                        JobVec jobs,
                        MappingVec mappings,
                        double max_mig_percentile,
                        Encoding encoding) {
        super(calculateNVariables(jobs.flattenJobs()),
              calculateNObjectives(mappings.size() > 0),
              calculateNConstraints(pms, jobs, mappings.size() > 0));
        this.pms = pms;
        this.jobs = jobs;
        this.vms = this.jobs.flattenJobs();
        this.mappings = mappings;
        assert(max_mig_percentile >= 0.0 && max_mig_percentile <= 1.0);
        this.max_mig_percentile = max_mig_percentile;
        this.anti_coloc_vms = retrieveAntiColocatableVirtualMachines(this.jobs);
        this.plat_constrained_vms = retrievePlatformConstrainedVirtualMachines(this.vms);
        this.pm_id_to_idx = Utils.makePhysicalMachineIDtoIndexMap_Gen(this.pms);
        this.vm_id_to_idx = Utils.makeVirtualMachineIDtoIndexMap_Gen(this.vms);
        this.total_cpu_cap = cpuCapacitySum(this.pms);
        this.total_mem_cap = memoryCapacitySum(this.pms);
        this.total_cpu_req = cpuRequirementSum(this.vms);
        this.total_mem_req = memoryRequirementSum(this.vms);
        this.encoding = makeEncoding(encoding);
    }
    
    public PhysicalMachineVec getPhysicalMachines() { return this.pms; }
    public JobVec getJobs() { return this.jobs; }
    public VirtualMachineVec getVirtualMachines() { return this.vms; }
    public MappingVec getMappings() { return this.mappings; }
    public double getMaxMigrationPercentile() { return this.max_mig_percentile; }
    
    public BigInteger getTotalCPUCapacity() { return this.total_cpu_cap; }
    public BigInteger getTotalMemoryCapacity() { return this.total_mem_cap; }
    public BigInteger getTotalCPURequirements() { return this.total_cpu_req; }
    public BigInteger getTotalMemoryRequirements() { return this.total_mem_req; }
    
    public IVec<VirtualMachineVec> getAntiColocatableVirtualMachines() { return this.anti_coloc_vms; }
    public VirtualMachineVec getPlatformConstrainedVirtualMachines() {
        return this.plat_constrained_vms;
    }
    
    public double getMaxEnergyConsumption() {
        double val = 0.0;
        for (int i = 0; i < pms.size(); ++i) {
            val += (double)pms.get(i).getMaxConsumption();
        }
        return val;
    }
    public double getMaxResourceWastage() {
        //return (double)pms.size()*(1.0+WASTAGE_EPSILON)/2.0;
        return (double)pms.size();
    }
    public double getMaxMigrationCost() {
        double val = 0.0;
        for (int i = 0; i < mappings.size(); ++i) {
            val += mappings.get(i).getVirtualMachine().getMemory().doubleValue();
        }
        return val;
    }

    public int[] getVariableAssignment(Solution solution) {
        return this.encoding.getVariableAssignment(solution);
    }
    
    public double getConstraintViolation(Solution sol) {
        ensureConstraintViolation(sol);
        return ((Double)sol.getAttribute(VMCwMProblem.VIOLATION_ATR)).doubleValue();
    }
    
    public void setVariableValue(Variable var, int val) {
        this.encoding.setVariableValue(var, val);
    }
    
    public int getVirtualMachineIndex(VirtualMachine vm) {
        return this.vm_id_to_idx.get(vm.getID()).intValue();
    }
    public int getPhysicalMachineIndex(PhysicalMachine pm) {
        return this.pm_id_to_idx.get(new Integer(pm.getID())).intValue();
    }
    
    public void setEncoding(VMCwMProblem.Encoding enc) {
        this.encoding = makeEncoding(enc);
    }
    
    public void setMaxMemoryMigrationPercentile(double percentile) {
        assert(percentile >= 0.0 && percentile <= 1.0);
        this.max_mig_percentile = percentile;
    }
    
    @Override
    public void evaluate(Solution solution) {
        // Initialize
        int[] x = getVariableAssignment(solution);
        BigInteger[] used_cpu_caps = new BigInteger[this.pms.size()];
        BigInteger[] used_mem_caps = new BigInteger[this.pms.size()];
        for (int i = 0; i < this.pms.size(); ++i) {
            used_cpu_caps[i] = BigInteger.ZERO;
            used_mem_caps[i] = BigInteger.ZERO;
        }
        for (int i = 0; i < this.vms.size(); ++i) {
            int pm_idx = x[i];
            used_cpu_caps[pm_idx] = used_cpu_caps[pm_idx].add(this.vms.get(i).getCPU());
            used_mem_caps[pm_idx] = used_mem_caps[pm_idx].add(this.vms.get(i).getMemory());
        }
        // Set energy consumption objective function value
        double total_energy = 0.0;
        for (int i = 0; i < this.pms.size(); ++i) {
            if (    !used_cpu_caps[i].equals(BigInteger.ZERO) && 
                    !used_mem_caps[i].equals(BigInteger.ZERO)) {
                double idle_energy = (double)this.pms.get(i).getIdleConsumption();
                double max_energy = (double)this.pms.get(i).getMaxConsumption();
                double used_cpu = used_cpu_caps[i].doubleValue()/
                                  this.pms.get(i).getCPU().doubleValue();
                total_energy += idle_energy;
                total_energy += (max_energy-idle_energy) * used_cpu;
            }
        }
        solution.setObjective(ENERGY_OBJ_INDEX, total_energy);
        // Set resource wastage objective function value
        double total_wastage = 0.0;
        for (int i = 0; i < pms.size(); ++i) {
            if (    !used_cpu_caps[i].equals(BigInteger.ZERO) && 
                    !used_mem_caps[i].equals(BigInteger.ZERO)) {
                double used_cpu = used_cpu_caps[i].doubleValue()/pms.get(i).getCPU().doubleValue();
                double used_mem = used_mem_caps[i].doubleValue()/
                                  pms.get(i).getMemory().doubleValue();
                double cpu_left = 1.0 - used_cpu;
                double mem_left = 1.0 - used_mem;
                /*double wastage = (Math.abs(cpu_left-mem_left) + WASTAGE_EPSILON) /
                                 (2 * (used_cpu+used_mem));*/
                double wastage = Math.abs(cpu_left-mem_left);
                total_wastage += wastage;
            }
        }
        solution.setObjective(WASTAGE_OBJ_INDEX, total_wastage);
        // Set migration objective function value
        BigInteger migged_mem = BigInteger.ZERO;
        if (solution.getNumberOfObjectives() == 3) {
            for (int i = 0; i < mappings.size(); ++i) {
                VirtualMachine vm = mappings.get(i).getVirtualMachine();
                PhysicalMachine pm = mappings.get(i).getPhysicalMachine();
                if (x[getVirtualMachineIndex(vm)] != getPhysicalMachineIndex(pm)) {
                    migged_mem = migged_mem.add(vm.getMemory());
                }
            }
            solution.setObjective(MIGRATION_OBJ_INDEX, migged_mem.doubleValue());
        }
        // Set capacity constraint violations
        for (int i = 0; i < pms.size(); ++i) {
            double diff = used_cpu_caps[i].subtract(pms.get(i).getCPU()).doubleValue();
            solution.setConstraint(cpuCapacityConstraintIndex(i), diff > 0.0 ? diff : 0.0);
            diff = used_mem_caps[i].subtract(pms.get(i).getMemory()).doubleValue();
            solution.setConstraint(memoryCapacityConstraintIndex(i), diff > 0.0 ? diff : 0.0);
        }
        // Set anti-colocation constraint violations
        for (int i = 0; i < anti_coloc_vms.size(); ++i) {
            double violation = 0.0;
            Set<PhysicalMachine> used_pms = new HashSet<PhysicalMachine>();
            for (int j = 0; j < anti_coloc_vms.get(i).size(); ++j) {
                VirtualMachine vm = anti_coloc_vms.get(i).get(j);
                PhysicalMachine pm = pms.get(x[getVirtualMachineIndex(vm)]);
                if (used_pms.contains(pm)) {
                    violation++;
                }
                used_pms.add(pm);
            }
            solution.setConstraint(antiColocationConstraintIndex(i), violation);
        }
        // Set platform constraint violations
        for (int i = 0; i < plat_constrained_vms.size(); ++i) {
            VirtualMachine vm = plat_constrained_vms.get(i);
            if (vm.canRunInPhysicalMachine(pms.get(x[getVirtualMachineIndex(vm)]))) {
                solution.setConstraint(platformConstraintIndex(i), 0.0);
            }
            else {
                solution.setConstraint(platformConstraintIndex(i), 1.0);
            }
        }
        // Set migration constraint violation
        if (solution.getNumberOfObjectives() == 3) {
            BigInteger max_mig_mem =
                    new BigDecimal(total_mem_cap).multiply(
                            new BigDecimal(max_mig_percentile)).toBigInteger();
            double diff = migged_mem.subtract(max_mig_mem).doubleValue();
            solution.setConstraint(migrationConstraintIndex(), diff > 0.0 ? diff : 0.0);
        }
        ensureConstraintViolation(solution);
    }

    @Override
    public Solution newSolution() {
        Solution solution = new @Gen Solution(getNumberOfVariables(),
                                              getNumberOfObjectives(),
                                              getNumberOfConstraints());
        for (int i = 0; i < getNumberOfVariables(); ++i) {
            solution.setVariable(i, makeVariable());
        }
        return solution;
    }

}
