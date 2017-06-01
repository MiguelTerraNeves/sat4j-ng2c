package vmalloc.algorithm;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Utils;
import vmalloc.constraints.ConstraintSolver;
import vmalloc.constraints.PseudoBooleanSolver;
import vmalloc.domain.JobVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.evolutionary.VMCwMProblem;

public abstract class MultiObjectiveConstraintBasedAllocAlgorithm extends ConstraintBasedAllocAlgorithm {

    protected IVecInt pm_vars = null;
    protected IVec<IVec<IVecInt>> job_vars = null;
    protected IVecInt aux_pm_pvars = null;
    protected IVec<IVec<IVecInt>> aux_job_plus_vars = null;
    protected IVec<IVec<IVecInt>> aux_job_minus_vars = null;
    
    // Auxiliary structures for encoding objective functions
    protected IVecInt energy_lits;
    protected IVec<BigDecimal> energy_coeffs;
    protected BigDecimal energy_coeff_sum;
    protected IVecInt wastage_lits;
    protected IVec<BigDecimal> wastage_coeffs;
    protected BigDecimal wastage_coeff_sum;
    protected IVecInt migration_lits;
    protected IVec<BigInteger> migration_coeffs;
    protected BigInteger migration_coeff_sum;
    private Map<Pair<BigInteger,BigInteger>, Pair<IVec<BigDecimal>,IVec<BigDecimal> > > wastage_coeffs_cache;
    
    public MultiObjectiveConstraintBasedAllocAlgorithm(VMCwMProblem instance) {
        this(instance, false);
    }
    public MultiObjectiveConstraintBasedAllocAlgorithm(VMCwMProblem instance, boolean break_symms) {
        super(instance, break_symms);
        this.wastage_coeffs_cache =
                new @Gen HashMap<Pair<BigInteger,BigInteger>, Pair<IVec<BigDecimal>,IVec<BigDecimal> > >();
    }
    

    private Pair<IVec<BigDecimal>,IVec<BigDecimal> > getWastageCoefficients(PhysicalMachine pm,
                                                                            VirtualMachineVec vms) {
        Pair<BigInteger,BigInteger> resource_pair = Pair.of(pm.getCPU(), pm.getMemory());
        if (!this.wastage_coeffs_cache.containsKey(resource_pair)) {
            double[] norm_cpus = getNormalizedCPURequirements(vms, pm);
            double[] norm_mems = getNormalizedMemoryRequirements(vms, pm);
            IVec<BigDecimal> norm_mem_minus_cpu = indexWiseSubtraction(norm_mems, norm_cpus);
            IVec<BigDecimal> norm_cpu_minus_mem = indexWiseSubtraction(norm_cpus, norm_mems);
            this.wastage_coeffs_cache.put(resource_pair,
                                          Pair.of(norm_mem_minus_cpu, norm_cpu_minus_mem));
            assert(norm_mem_minus_cpu.size() == vms.size());
            assert(norm_cpu_minus_mem.size() == vms.size());
        }
        assert(this.wastage_coeffs_cache.containsKey(resource_pair));
        return this.wastage_coeffs_cache.get(resource_pair);
    }
    
    private IVec<BigDecimal> indexWiseSubtraction(double[] a1, double[] a2) {
        assert(a1.length == a2.length);
        IVec<BigDecimal> sub_array = new @Gen Vec<BigDecimal>(a1.length);
        for (int i = 0; i < a1.length; ++i) {
            sub_array.push(new @Gen BigDecimal(a1[i] - a2[i]));
        }
        return sub_array;
    }
    
    protected void addWastageAuxiliaryConstraints(ConstraintSolver solver,
                                                  PhysicalMachineVec pms,
                                                  IVecInt aux_pm_vars,
                                                  JobVec jobs,
                                                  IVec<IVec<IVecInt>> job_vars,
                                                  IVec<IVec<IVecInt>> aux_job_plus_vars,
                                                  IVec<IVec<IVecInt>> aux_job_minus_vars)
                                                          throws ContradictionException {
        assert(pms.size() == aux_pm_vars.size());
        assert(jobs.size() == job_vars.size());
        assert(jobs.size() == aux_job_plus_vars.size());
        assert(jobs.size() == aux_job_minus_vars.size());
        VirtualMachineVec vms = jobs.flattenJobs();
        IVec<IVecInt> vm_vars = flattenJobVars(job_vars);
        assert(vms.size() == vm_vars.size());
        for (int i = 0; i < pms.size(); ++i) {
            IVecInt lits = new VecInt();
            for (int j = 0; j < vms.size(); ++j) {
                lits.push(vm_vars.get(j).get(i));
            }
            PhysicalMachine pm = pms.get(i);
            Pair<IVec<BigDecimal>,IVec<BigDecimal> > coeff_vec_pair = getWastageCoefficients(pm, vms);
            IVec<BigDecimal> norm_mem_minus_cpu = coeff_vec_pair.getLeft();
            IVec<BigDecimal> norm_cpu_minus_mem = coeff_vec_pair.getRight();
            assert(norm_mem_minus_cpu.size() == lits.size());
            assert(norm_cpu_minus_mem.size() == lits.size());
            lits.push(aux_pm_vars.get(i));
            norm_mem_minus_cpu.push(BigDecimal.ONE);
            solver.addGreaterOrEqual(lits, norm_mem_minus_cpu, BigDecimal.ZERO);
            norm_mem_minus_cpu.pop();
            lits.pop();
            lits.push(-aux_pm_vars.get(i));
            norm_cpu_minus_mem.push(BigDecimal.ONE);
            solver.addGreaterOrEqual(lits, norm_cpu_minus_mem, BigDecimal.ZERO);
            norm_cpu_minus_mem.pop();
        }
        IVec<IVecInt> aux_vm_plus_vars = flattenJobVars(aux_job_plus_vars);
        IVec<IVecInt> aux_vm_minus_vars = flattenJobVars(aux_job_minus_vars);
        assert(vms.size() == aux_vm_plus_vars.size());
        assert(vms.size() == aux_vm_minus_vars.size());
        for (int i = 0; i < pms.size(); ++i) {
            IVecInt lits = new VecInt();
            lits.push(-aux_pm_vars.get(i));
            for (int j = 0; j < vms.size(); ++j) {
                lits.push(-aux_vm_plus_vars.get(j).get(i));
                solver.addClause(lits);
                lits.pop();
            }
            lits.pop();
            lits.push(aux_pm_vars.get(i));
            for (int j = 0; j < vms.size(); ++j) {
                lits.push(-aux_vm_minus_vars.get(j).get(i));
                solver.addClause(lits);
                lits.pop();
            }
        }
        for (int i = 0; i < vms.size(); ++i) {
            for (int j = 0; j < pms.size(); ++j) {
                // Encode 'vm_vars - plus_var - minus_var = 0', explicit equality is extremely inefficient
                int vm_var = vm_vars.get(i).get(j);
                int plus_var = aux_vm_plus_vars.get(i).get(j);
                int minus_var = aux_vm_minus_vars.get(i).get(j);
                solver.addClause(new VecInt(new int[] { -vm_var, -plus_var, -minus_var }));
                solver.addClause(new VecInt(new int[] { -vm_var, plus_var, minus_var }));
                solver.addClause(new VecInt(new int[] { vm_var, -plus_var }));
                solver.addClause(new VecInt(new int[] { vm_var, -minus_var }));
                solver.addClause(new VecInt(new int[] { -plus_var, -minus_var }));
            }
        }
    }
    
    protected ConstraintSolver buildSolver() throws ContradictionException {
        System.out.println("c Building formula");
        ConstraintSolver solver = new @Gen PseudoBooleanSolver();
        this.pm_vars = newVarsForPMs(solver, this.instance.getPhysicalMachines());
        this.job_vars = newVarsForJobs(solver,
                                       this.instance.getPhysicalMachines(),
                                       this.instance.getJobs());
        addLowerBoundConstraints(solver, this.instance.getPhysicalMachines(), this.pm_vars);
        addExactlyOnePMConstraints(solver, this.job_vars);
        addPMCapacityConstraints(solver,
                                 this.instance.getPhysicalMachines(),
                                 this.instance.getJobs(),
                                 this.job_vars);
        addVarLinkConstraints(solver, this.pm_vars, this.job_vars);
        addAntiColocationConstraints(solver,
                                     this.instance.getPhysicalMachines(),
                                     this.instance.getJobs(),
                                     this.job_vars);
        addPlatformConstraints(solver,
                               this.instance.getPhysicalMachines(),
                               this.instance.getJobs(),
                               this.job_vars);
        addMigrationConstraint(solver,
                               this.instance.getPhysicalMachines(),
                               this.instance.getJobs(),
                               this.job_vars,
                               this.instance.getMappings(),
                               this.instance.getMaxMigrationPercentile());
        if (this.break_symms) {
            addSymmetryBreakingConstraints(solver,
                                           this.instance.getPhysicalMachines(),
                                           this.pm_vars,
                                           this.instance.getJobs(),
                                           this.job_vars,
                                           this.instance.getMappings());
        }
        this.aux_pm_pvars = newVarsForPMs(solver, this.instance.getPhysicalMachines());
        this.aux_job_plus_vars = newVarsForJobs(solver,
                                                this.instance.getPhysicalMachines(),
                                                this.instance.getJobs());
        this.aux_job_minus_vars = newVarsForJobs(solver,
                                                 this.instance.getPhysicalMachines(),
                                                 this.instance.getJobs());
        addWastageAuxiliaryConstraints(solver,
                                       this.instance.getPhysicalMachines(),
                                       this.aux_pm_pvars,
                                       this.instance.getJobs(),
                                       this.job_vars,
                                       this.aux_job_plus_vars,
                                       this.aux_job_minus_vars);
        printElapsedTime();
        return solver;
    }
    
    private void initializeEnergyObjectiveFunction() {
        VirtualMachineVec vms = this.instance.getJobs().flattenJobs();
        IVec<IVecInt> vm_vars = flattenJobVars(this.job_vars);
        this.energy_lits = new @Gen VecInt();
        this.energy_coeffs = new @Gen Vec<BigDecimal>();
        for (int i = 0; i < this.instance.getPhysicalMachines().size(); ++i) {
            PhysicalMachine pm = this.instance.getPhysicalMachines().get(i);
            this.energy_coeffs.push(new @Gen BigDecimal(pm.getIdleConsumption()));
            this.energy_lits.push(this.pm_vars.get(i));
            int energy_range = pm.getMaxConsumption() - pm.getIdleConsumption();
            double[] norm_cpus = getNormalizedCPURequirements(vms, pm);
            for (int j = 0; j < vms.size(); ++j) {
                this.energy_coeffs.push(new @Gen BigDecimal(energy_range * norm_cpus[j]));
                this.energy_lits.push(vm_vars.get(j).get(i));
            }
        }
        assert(this.energy_coeffs.size() == this.energy_lits.size());
        this.energy_coeff_sum = Utils.bigDecimalVecSum(this.energy_coeffs);
    }
    
    private void initializeWastageObjectiveFunction() {
        VirtualMachineVec vms = this.instance.getJobs().flattenJobs();
        IVec<IVecInt> aux_vm_plus_vars = flattenJobVars(this.aux_job_plus_vars);
        IVec<IVecInt> aux_vm_minus_vars = flattenJobVars(this.aux_job_minus_vars);
        this.wastage_lits = new @Gen VecInt();
        this.wastage_coeffs = new @Gen Vec<BigDecimal>();
        for (int i = 0; i < this.instance.getPhysicalMachines().size(); ++i) {
            PhysicalMachine pm = this.instance.getPhysicalMachines().get(i);
            Pair<IVec<BigDecimal>,IVec<BigDecimal> > coeff_vec_pair = getWastageCoefficients(pm, vms);
            IVec<BigDecimal> norm_mem_minus_cpu = coeff_vec_pair.getLeft();
            IVec<BigDecimal> norm_cpu_minus_mem = coeff_vec_pair.getRight();
            assert(norm_mem_minus_cpu.size() == vms.size());
            assert(norm_cpu_minus_mem.size() == vms.size());
            for (int j = 0; j < vms.size(); ++j) {
                this.wastage_coeffs.push(norm_mem_minus_cpu.get(j));
                this.wastage_lits.push(aux_vm_plus_vars.get(j).get(i));
                this.wastage_coeffs.push(norm_cpu_minus_mem.get(j));
                this.wastage_lits.push(aux_vm_minus_vars.get(j).get(i));
            }
        }
        assert(this.wastage_coeffs.size() == this.wastage_lits.size());
        this.wastage_coeff_sum = Utils.bigDecimalVecSum(this.wastage_coeffs);
    }
    
    private void initializeMigrationObjectiveFunction() {
        VirtualMachineVec vms = this.instance.getJobs().flattenJobs();
        IVec<IVecInt> vm_vars = flattenJobVars(this.job_vars);
        Map<Integer, Integer> pm_id_to_idx =
                Utils.makePhysicalMachineIDtoIndexMap(this.instance.getPhysicalMachines());
        Map<String, Integer> vm_id_to_idx = Utils.makeVirtualMachineIDtoIndexMap(vms);
        this.migration_lits = new @Gen VecInt();
        this.migration_coeffs = new @Gen Vec<BigInteger>();
        for (int i = 0; i < this.instance.getMappings().size(); ++i) {
            VirtualMachine vm = this.instance.getMappings().get(i).getVirtualMachine();
            PhysicalMachine pm = this.instance.getMappings().get(i).getPhysicalMachine();
            int vm_idx = vm_id_to_idx.get(vm.getID()).intValue();
            int pm_idx = pm_id_to_idx.get(new Integer(pm.getID())).intValue();
            this.migration_coeffs.push(vm.getMemory());
            this.migration_lits.push(-vm_vars.get(vm_idx).get(pm_idx));
        }
        assert(this.migration_lits.size() == this.migration_coeffs.size());
        this.migration_coeff_sum = Utils.bigIntegerVecSum(this.migration_coeffs);
    }
    
    protected void initializeObjectiveFunctions() {
        initializeEnergyObjectiveFunction();
        initializeWastageObjectiveFunction();
        if (this.instance.getMappings().size() > 0) {
            initializeMigrationObjectiveFunction();
        }
    }
    
    private BigDecimal computeObjectiveValue(ConstraintSolver solver,
                                             IVecInt lits,
                                             IVec<BigDecimal> coeffs) {
        BigDecimal obj_val = BigDecimal.ZERO;
        for (int i = 0; i < lits.size(); ++i) {
            if (solver.modelValue(lits.get(i))) {
                obj_val = obj_val.add(coeffs.get(i));
            }
        }
        return obj_val;
    }
    protected BigDecimal computeEnergyConsumption(ConstraintSolver solver) {
        return computeObjectiveValue(solver, this.energy_lits, this.energy_coeffs);
    }
    protected BigDecimal computeResourceWastage(ConstraintSolver solver) {
        return computeObjectiveValue(solver, this.wastage_lits, this.wastage_coeffs);
    }
    protected BigInteger computeMigrationCost(ConstraintSolver solver) {
        BigInteger obj_val = BigInteger.ZERO;
        for (int i = 0; i < this.migration_lits.size(); ++i) {
            if (!solver.modelValue(-this.migration_lits.get(i))) {
                obj_val = obj_val.add(this.migration_coeffs.get(i));
            }
        }
        return obj_val;
    }
    
    // used by Pareto-MCS algorithms
    protected IVecInt buildUndefFormulas() {
        IVecInt undef_fmls = new @Gen VecInt();
        for (int i = 0; i < this.energy_lits.size(); ++i) {
            assert(this.energy_coeffs.get(i).compareTo(BigDecimal.ZERO) > 0);
            undef_fmls.push(-this.energy_lits.get(i));
        }
        for (int i = 0; i < this.wastage_lits.size(); ++i) {
            if (this.wastage_coeffs.get(i).compareTo(BigDecimal.ZERO) > 0) {
                undef_fmls.push(-this.wastage_lits.get(i));
            }
            else if (this.wastage_coeffs.get(i).compareTo(BigDecimal.ZERO) < 0) {
                undef_fmls.push(this.wastage_lits.get(i));
            }
        }
        if (this.instance.getMappings().size() > 0) {
            for (int i = 0; i < this.migration_lits.size(); ++i) {
                undef_fmls.push(-this.migration_lits.get(i));
            }
        }
        return undef_fmls;
    }
    
}
