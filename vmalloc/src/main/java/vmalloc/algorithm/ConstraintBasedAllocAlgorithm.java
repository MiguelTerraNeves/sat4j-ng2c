package vmalloc.algorithm;

import static org.sat4j.GlobalDefs.USE_NG2C;
import static org.sat4j.GlobalDefs.ANNOTATE_CONSTRAINTS_EXTERNAL;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.moeaframework.core.PRNG;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Utils;
import vmalloc.constraints.ConstraintID;
import vmalloc.constraints.ConstraintSolver;
import vmalloc.domain.JobVec;
import vmalloc.domain.Mapping;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.evolutionary.VMCwMProblem;

public abstract class ConstraintBasedAllocAlgorithm extends AllocAlgorithm {

    // =============================================================================================
    // ======================================= PRIVATE API =======================================
    // =============================================================================================

    // TODO: turn configurable or adaptive
    private static final double EPSILON = 0.8;
    private static final int KEY_SIZE = 1;
    
    private HashFunctionType hash_type = HashFunctionType.GLOBAL;

    private int mergeSubXOR(ConstraintSolver solver,
                            int left_out,
                            int right_out,
                            IVec<ConstraintID> ids)
            throws ContradictionException {
        int out = newVar(solver);
        ids.push(solver.addRemovableClause(new VecInt(new int[] { -left_out, right_out, out })));
        ids.push(solver.addRemovableClause(new VecInt(new int[] { left_out, -right_out, out })));
        ids.push(solver.addRemovableClause(new VecInt(new int[] { left_out, right_out, -out })));
        ids.push(solver.addRemovableClause(new VecInt(new int[] { -left_out, -right_out, -out })));
        return out;
    }
    
    private int xorToCNF(ConstraintSolver solver, IVecInt lits, IVec<ConstraintID> ids)
            throws ContradictionException {
        assert(lits.size() > 1);
        int split = lits.size() / 2;
        IVecInt left = new VecInt();
        for (int i = 0; i < split; ++i) {
            left.push(lits.get(i));
        }
        IVecInt right = new VecInt();
        for (int i = split; i < lits.size(); ++i) {
            right.push(lits.get(i));
        }
        int left_out, right_out;
        if (left.size() > 1) {
            left_out = xorToCNF(solver, left, ids);
        }
        else {
            left_out = left.get(0);
        }
        if (right.size() > 1) {
            right_out = xorToCNF(solver, right, ids);
        }
        else {
            right_out = right.get(0);
        }
        return mergeSubXOR(solver, left_out, right_out, ids);
    }
    
    private int computeEnumerationThreshold(double epsilon) {
        return (int)(1 + 9.84 * (1 + epsilon/(1+epsilon)) * (1 + 1/epsilon) * (1 + 1/epsilon));
    }
    
    // =============================================================================================
    // ======================================= PROTECTED API =======================================
    // =============================================================================================
    
    protected boolean break_symms = false;
    
    protected int newVar(ConstraintSolver solver) {
        solver.newVars(1);
        return solver.nVars();
    }

    protected IVecInt newVarsForPMs(ConstraintSolver solver, PhysicalMachineVec pms) {
        IVecInt pm_vars = ANNOTATE_CONSTRAINTS_EXTERNAL ? new @Gen VecInt() : new VecInt();
        int nvars = solver.nVars();
        solver.newVars(pms.size());
        for (int i = 1; i <= pms.size(); ++i) {
            pm_vars.push(nvars+i);
        }
        return pm_vars;
    }
    
    protected IVec<IVecInt> newVarsForVMs(ConstraintSolver solver,
                                          PhysicalMachineVec pms,
                                          VirtualMachineVec vms) {
        IVec<IVecInt> vm_vars = ANNOTATE_CONSTRAINTS_EXTERNAL ? new @Gen Vec<IVecInt>() : new Vec<IVecInt>();
        int new_vars = 0, nvars = solver.nVars();
        for (int i = 0; i < vms.size(); ++i) {
            vm_vars.push(ANNOTATE_CONSTRAINTS_EXTERNAL ? new @Gen VecInt() : new VecInt());
            for (int j = 0; j < pms.size(); ++j) {
                vm_vars.get(i).push(++new_vars + nvars);
            }
        }
        solver.newVars(new_vars);
        return vm_vars;
    }
    
    protected IVec<IVec<IVecInt>> newVarsForJobs(ConstraintSolver solver,
                                                 PhysicalMachineVec pms,
                                                 JobVec jobs) {
        IVec<IVec<IVecInt>> job_vars = ANNOTATE_CONSTRAINTS_EXTERNAL ? new @Gen Vec<IVec<IVecInt>>() : new Vec<IVec<IVecInt>>();
        for (int i = 0; i < jobs.size(); ++i) {
            job_vars.push(newVarsForVMs(solver, pms, jobs.get(i).virtualMachinesAsVec()));
            assert(jobs.get(i).nVirtualMachines() == job_vars.get(i).size());
        }
        return job_vars;
    }
    
    protected IVec<IVecInt> flattenJobVars(IVec<IVec<IVecInt>> job_vars) {
        IVec<IVecInt> vm_vars = new Vec<IVecInt>();
        for (int i = 0; i < job_vars.size(); ++i) {
            for (int j = 0; j < job_vars.get(i).size(); ++j) {
                vm_vars.push(job_vars.get(i).get(j));
            }
        }
        return vm_vars;
    }
    
    protected IVecInt flattenLitVectors(IVec<IVecInt> lits) {
        IVecInt flat_lits = new VecInt();
        for (int i = 0; i < lits.size(); ++i) {
            for (int j = 0; j < lits.get(i).size(); ++j) {
                flat_lits.push(lits.get(i).get(j));
            }
        }
        return flat_lits;
    }
    
    protected void addLowerBoundConstraints(ConstraintSolver solver,
                                            PhysicalMachineVec pms,
                                            IVecInt pm_vars) throws ContradictionException {
        solver.addGreaterOrEqual(pm_vars,
                                 new Vec<BigInteger>(pms.getCPUs()),
                                 this.instance.getTotalCPURequirements());
        solver.addGreaterOrEqual(pm_vars,
                                 new Vec<BigInteger>(pms.getMemories()),
                                 this.instance.getTotalMemoryRequirements());
    }
    
    protected void addExactlyOnePMConstraintsForVMs(ConstraintSolver solver, IVec<IVecInt> vm_vars)
            throws ContradictionException {
        for (int i = 0; i < vm_vars.size(); ++i) {
            solver.addExactly(vm_vars.get(i), 1);
        }
    }
    
    protected void addExactlyOnePMConstraints(
            ConstraintSolver solver, IVec<IVec<IVecInt>> job_vars) throws ContradictionException {
        for (int i = 0; i < job_vars.size(); ++i) {
            addExactlyOnePMConstraintsForVMs(solver, job_vars.get(i));
        }
    }
    
    protected void addPMCapacityConstraintsForVMs(ConstraintSolver solver,
                                                  PhysicalMachineVec pms,
                                                  VirtualMachineVec vms,
                                                  IVec<IVecInt> vm_vars)
                                                          throws ContradictionException {
        IVecInt lits = new VecInt();
        for (int i = 0; i < pms.size(); ++i) {
            lits.clear();
            for (int j = 0; j < vm_vars.size(); ++j) {
                lits.push(vm_vars.get(j).get(i));
            }
            solver.addLessOrEqual(lits, new Vec<BigInteger>(vms.getCPUs()), pms.get(i).getCPU());
            solver.addLessOrEqual(lits, new Vec<BigInteger>(vms.getMemories()), pms.get(i).getMemory());
        }
    }
    
    protected void addPMCapacityConstraints(ConstraintSolver solver,
                                            PhysicalMachineVec pms,
                                            JobVec jobs,
                                            IVec<IVec<IVecInt>> job_vars)
                                                    throws ContradictionException {
        addPMCapacityConstraintsForVMs(solver, pms, jobs.flattenJobs(), flattenJobVars(job_vars));
    }
    
    protected void addVarLinkConstraintsForVMs(ConstraintSolver solver,
                                               IVecInt pm_vars,
                                               IVec<IVecInt> vm_vars) throws ContradictionException {
        for (int i = 0; i < pm_vars.size(); ++i) {
            for (int j = 0; j < vm_vars.size(); ++j) {
                IVecInt clause = new VecInt();
                clause.push(-vm_vars.get(j).get(i));
                clause.push(pm_vars.get(i));
                solver.addClause(clause);
            }
        }
    }
    
    protected void addVarLinkConstraints(ConstraintSolver solver,
                                         IVecInt pm_vars,
                                         IVec<IVec<IVecInt>> job_vars)
                                                 throws ContradictionException {
        addVarLinkConstraintsForVMs(solver, pm_vars, flattenJobVars(job_vars));
    }
    
    protected void addAntiColocationConstraints(ConstraintSolver solver,
                                                PhysicalMachineVec pms,
                                                JobVec jobs,
                                                IVec<IVec<IVecInt>> job_vars)
                                                        throws ContradictionException {
        for (int i = 0; i < jobs.size(); ++i) {
            IVec<IVecInt> anti_coloc_vm_vars = new Vec<IVecInt>();
            for (int j = 0; j < jobs.get(i).nVirtualMachines(); ++j) {
                if (jobs.get(i).getVirtualMachine(j).isAntiColocatable()) {
                    anti_coloc_vm_vars.push(job_vars.get(i).get(j));
                }
            }
            for (int j = 0; j < pms.size(); ++j) {
                IVecInt cons_vars = new VecInt();
                for (int k = 0; k < anti_coloc_vm_vars.size(); ++k) {
                    cons_vars.push(anti_coloc_vm_vars.get(k).get(j));
                }
                if (cons_vars.size() > 1) {
                    solver.addAtMost(cons_vars, 1);
                }
            }
        }
    }
    
    protected void addPlatformConstraints(ConstraintSolver solver,
                                          PhysicalMachineVec pms,
                                          JobVec jobs,
                                          IVec<IVec<IVecInt>> job_vars)
                                                  throws ContradictionException {
        Map<Integer, Integer> pm_id_to_idx = Utils.makePhysicalMachineIDtoIndexMap(pms);
        for (int i = 0; i < jobs.size(); ++i) {
            for (int j = 0; j < jobs.get(i).nVirtualMachines(); ++j) {
                VirtualMachine vm = jobs.get(i).getVirtualMachine(j);
                for (int k = 0; k < vm.getUnallowedPhysicalMachines().size(); ++k) {
                    Integer unallowed_id =
                            new Integer(vm.getUnallowedPhysicalMachines().get(k).getID());
                    if (pm_id_to_idx.containsKey(unallowed_id)) {
                        IVecInt unit_cl = new VecInt();
                        unit_cl.push(-job_vars.get(i).get(j).get(pm_id_to_idx.get(unallowed_id)));
                        solver.addClause(unit_cl);
                    }
                }
            }
        }
    }
    
    protected void addMigrationConstraintForVMs(ConstraintSolver solver,
                                                PhysicalMachineVec pms,
                                                VirtualMachineVec vms,
                                                IVec<IVecInt> vm_vars,
                                                MappingVec mappings,
                                                double max_mig_percentile)
                                                        throws ContradictionException {
        if (mappings.size() > 0) {
            assert(max_mig_percentile >= 0.0 && max_mig_percentile <= 1.0);
            Map<Integer, Integer> pm_id_to_idx = Utils.makePhysicalMachineIDtoIndexMap(pms);
            Map<String, Integer> vm_id_to_idx = Utils.makeVirtualMachineIDtoIndexMap(vms);
            IVecInt lits = new VecInt();
            IVec<BigInteger> coeffs = new Vec<BigInteger>();
            BigDecimal total_mem_cap = new BigDecimal(this.instance.getTotalMemoryCapacity());
            BigDecimal mig_percentile = new BigDecimal(max_mig_percentile);
            BigInteger rhs = total_mem_cap.multiply(mig_percentile).toBigInteger();
            for (int i = 0; i < mappings.size(); ++i) {
                VirtualMachine vm = mappings.get(i).getVirtualMachine();
                PhysicalMachine pm = mappings.get(i).getPhysicalMachine();
                int vm_idx = vm_id_to_idx.get(vm.getID()).intValue();
                int pm_idx = pm_id_to_idx.get(new Integer(pm.getID())).intValue();
                coeffs.push(vm.getMemory());
                lits.push(-vm_vars.get(vm_idx).get(pm_idx));
            }
            solver.addLessOrEqual(lits, coeffs, rhs);
        }
    }
    
    protected void addMigrationConstraint(ConstraintSolver solver,
                                          PhysicalMachineVec pms,
                                          JobVec jobs,
                                          IVec<IVec<IVecInt>> job_vars,
                                          MappingVec mappings,
                                          double max_mig_percentile)
                                                  throws ContradictionException {
        addMigrationConstraintForVMs(solver,
                                     pms,
                                     jobs.flattenJobs(),
                                     flattenJobVars(job_vars),
                                     mappings,
                                     max_mig_percentile);
    }
    
    protected void addSymmetryBreakingConstraints(ConstraintSolver solver,
                                                  PhysicalMachineVec pms,
                                                  IVecInt pm_vars,
                                                  JobVec jobs,
                                                  IVec<IVec<IVecInt>> job_vars,
                                                  MappingVec mappings)
                                                          throws ContradictionException {
        Map<VirtualMachine, PhysicalMachine> mapping_map =
                new HashMap<VirtualMachine, PhysicalMachine>();
        for (int i = 0; i < mappings.size(); ++i) {
            assert (!mapping_map.containsKey(mappings.get(i).getVirtualMachine()));
            mapping_map.put(mappings.get(i).getVirtualMachine(), mappings.get(i).getPhysicalMachine());
        }
        // Partition VMs into sets of symmetric VMs
        VirtualMachineVec vms = jobs.flattenJobs();
        IVec<IVecInt> vm_vars = flattenJobVars(job_vars);
        IVec<VirtualMachineVec> sym_groups = new Vec<VirtualMachineVec>();
        IVec<IVec<IVecInt>> sym_group_vars = new Vec<IVec<IVecInt>>();
        for (int i = 0; i < vms.size(); ++i) {
            VirtualMachine vm = vms.get(i);
            int j;
            for (j = 0; j < sym_groups.size(); ++j) {
                VirtualMachine rep_vm = sym_groups.get(j).get(0);
                if (((rep_vm.getJobID() == vm.getJobID() && rep_vm.isAntiColocatable() == vm.isAntiColocatable())
                        || (!rep_vm.isAntiColocatable() && !vm.isAntiColocatable()))
                        && rep_vm.getCPU().equals(vm.getCPU()) && rep_vm.getMemory().equals(vm.getMemory())
                        && mappedToSamePhysicalMachine(rep_vm, vm, mapping_map)) {
                    sym_groups.get(j).push(vm);
                    sym_group_vars.get(j).push(vm_vars.get(i));
                    break;
                }
            }
            if (j == sym_groups.size()) {
                VirtualMachineVec new_group = new VirtualMachineVec();
                IVec<IVecInt> new_group_vars = new Vec<IVecInt>();
                new_group.push(vm);
                new_group_vars.push(vm_vars.get(i));
                sym_groups.push(new_group);
                sym_group_vars.push(new_group_vars);
            }
        }
        // Add symmetry breaking constraints for VMs
        for (int i = 0; i < sym_groups.size(); ++i) {
            VirtualMachineVec sym_group = sym_groups.get(i);
            IVec<IVecInt> group_vars = sym_group_vars.get(i);
            boolean simplify = sym_group.get(0).isAntiColocatable();
            for (int j = 0; j < sym_group.size()-1; ++j) {
                VirtualMachine vm1 = sym_group.get(j);
                VirtualMachine vm2 = sym_group.get(j+1);
                IVecInt vm1_vars = group_vars.get(j);
                IVecInt vm2_vars = group_vars.get(j+1);
                IVecInt sym_pm_indexes = intersectAllowedPhysicalMachineIndexes(pms, vm1, vm2);
                simplify = simplify && sym_pm_indexes.size() == pms.size(); // FIXME: possible to simplify in more cases, but not trivial
                for (int k1 = 0; k1 < sym_pm_indexes.size(); ++k1) {
                    for (int k2 = 0; k2 <= k1; ++k2) {
                        if (k2 < k1 || vm1.isAntiColocatable()) {
                            IVecInt clause = new VecInt();
                            clause.push(-vm1_vars.get(sym_pm_indexes.get(k1)));
                            clause.push(-vm2_vars.get(sym_pm_indexes.get(k2)));
                            solver.addClause(clause);
                        }
                    }
                }
            }
            if (simplify) {
                for (int j = 0; j < sym_group.size(); ++j) {
                    IVecInt vmj_vars = group_vars.get(j);
                    for (int k = 0; k < j; ++k) {
                        IVecInt unit_cl = new VecInt();
                        unit_cl.push(-vmj_vars.get(k));
                        solver.addClause(unit_cl);
                    }
                    for (int k = 0; k < sym_group.size()-j-1; ++k) {
                        IVecInt unit_cl = new VecInt();
                        unit_cl.push(-vmj_vars.get(pms.size()-k-1));
                        solver.addClause(unit_cl);
                    }
                }
            }
        }
    }
    
    protected void addConjunction(ConstraintSolver solver, IVecInt lits) throws ContradictionException {
        for (int i = 0; i < lits.size(); ++i) {
            IVecInt unit_cl = new VecInt();
            unit_cl.push(lits.get(i));
            solver.addClause(unit_cl);
        }
    }
    
    // FIXME: could be encapsulated in the ConstraintSolver class
    protected void addRemovableConjunction(ConstraintSolver solver,
                                           IVecInt conjunction,
                                           IVec<ConstraintID> to_remove) throws ContradictionException {
        for (int i = 0; i < conjunction.size(); ++i) {
            IVecInt unit_cl = new VecInt();
            unit_cl.push(conjunction.get(i));
            to_remove.push(solver.addRemovableClause(unit_cl));
        }
    }

    // FIXME: could be encapsulated in the ConstraintSolver class
    protected int addRemovableXORConstraint(ConstraintSolver solver,
                                            IVecInt lits,
                                            boolean negated,
                                            IVec<ConstraintID> ids)
            throws ContradictionException {
        ids.clear();
        int activator;
        if (lits.size() > 1) {
            int out = xorToCNF(solver, lits, ids);
            activator = (negated) ? -out : out;
            //ids.push(solver.addRemovableClause(new VecInt(new int[] { (negated) ? -out : out })));
        }
        else if (lits.size() == 1) {
            activator = (negated) ? -lits.get(0) : lits.get(0);
            //ids.push(solver.addRemovableClause(
            //        new VecInt(new int[] { (negated) ? -lits.get(0) : lits.get(0) })));
        }
        else {
            activator = newVar(solver);
            ids.push(solver.addRemovableClause(new VecInt(new int[] { -activator })));
        }
        return activator;
    }

    protected IVec<ConstraintID> setHashFunction(ConstraintSolver solver,
                                                 IVecInt lits,
                                                 IVecInt asms) {
        IVec<ConstraintID> ids = new Vec<ConstraintID>();
        boolean done = false;
        while (!done) {
            asms.clear();
            try {
                for (int i = 0; i < KEY_SIZE; ++i) {
                    boolean key_bit = PRNG.nextBoolean(), alpha_0 = PRNG.nextBoolean();
                    IVecInt bit_lits = new VecInt();
                    for (int j = 0; j < lits.size(); ++j) {
                        if (PRNG.nextBoolean()) {
                            bit_lits.push(lits.get(j));
                        }
                    }
                    asms.push(addRemovableXORConstraint(solver, bit_lits, alpha_0 == key_bit, ids));
                }
                done = true;
            }
            catch (ContradictionException ce) {
                System.out.println("c Hash function led to empty cell, generating another one");
                solver.removeConstraints(ids);
                ids.clear();
            }
        }
        return ids;
    }
    
    protected HashFunctionType getHashType() { return this.hash_type; }

    protected IVec<ConstraintID> setHashFunction(ConstraintSolver solver, IVecInt lits) {
        IVec<ConstraintID> ids = new Vec<ConstraintID>();
        IVecInt activators = new VecInt();
        boolean done = false;
        while (!done) {
            try {
                ids = setHashFunction(solver, lits, activators);
                for (int i = 0; i < activators.size(); ++i) {
                    solver.addClause(new VecInt(new int[] { activators.get(i) }));
                }
                done = true;
            }
            catch (ContradictionException ce) {
                System.out.println("c Hash function led to empty cell, generating another one");
                solver.removeConstraints(ids);
                ids.clear();
                activators.clear();
            }
        }
        return ids;
    }
    
    protected int getEnumerationThreshold() { return computeEnumerationThreshold(EPSILON); }
    
    // used by MCS algorithms
    protected IVecInt extractSatisfied(ConstraintSolver solver, IVecInt undef_fmls) {
        IVecInt sat = new VecInt();
        int i = 0;
        while (i < undef_fmls.size()) {
            if (    (undef_fmls.get(i) < 0 && !solver.modelValue(-undef_fmls.get(i)) ||
                    (undef_fmls.get(i) > 0 && solver.modelValue(undef_fmls.get(i))))) {
                sat.push(undef_fmls.get(i));
                undef_fmls.set(i, undef_fmls.get(undef_fmls.size()-1));
                undef_fmls.pop();
            }
            else {
                i++;
            }
        }
        return sat;
    }
    
    protected Set<Integer> getUsedPhysicalMachineIndexesForVMs(ConstraintSolver solver,
                                                               PhysicalMachineVec pms,
                                                               IVec<IVecInt> vm_vars) {
        Set<Integer> indexes = new HashSet<Integer>();
        for (int i = 0; i < pms.size(); ++i) {
            for (int j = 0; j < vm_vars.size(); ++j) {
                if (solver.modelValue(vm_vars.get(j).get(i))) {
                    indexes.add(new Integer(i));
                    break;
                }
            }
        }
        return indexes;
    }
    
    protected Set<Integer> getUsedPhysicalMachineIndexes(ConstraintSolver solver,
                                                         PhysicalMachineVec pms,
                                                         IVec<IVec<IVecInt>> job_vars) {
        return getUsedPhysicalMachineIndexesForVMs(solver, pms, flattenJobVars(job_vars));
    }
    
    protected MappingVec modelToAllocationForVMs(ConstraintSolver solver,
                                                 PhysicalMachineVec pms,
                                                 VirtualMachineVec vms,
                                                 IVec<IVecInt> vm_vars) {
        MappingVec allocation = new MappingVec();
        for (int i = 0; i < vm_vars.size(); ++i) {
            for (int j = 0; j < vm_vars.get(i).size(); ++j) {
                if (solver.modelValue(vm_vars.get(i).get(j))) {
                    allocation.push(new Mapping(vms.get(i), pms.get(j)));
                    break;
                }
            }
        }
        return allocation;
    }
    
    protected MappingVec modelToAllocation(ConstraintSolver solver,
                                           PhysicalMachineVec pms,
                                           JobVec jobs,
                                           IVec<IVec<IVecInt>> job_vars) {
        return modelToAllocationForVMs(solver, pms, jobs.flattenJobs(), flattenJobVars(job_vars));
    }

    protected int getUsedPMsCount(ConstraintSolver solver,
                                  PhysicalMachineVec pms,
                                  IVec<IVec<IVecInt>> job_vars) {
        return getUsedPhysicalMachineIndexes(solver, pms, job_vars).size();
    }

    protected void checkSAT(ConstraintSolver solver) { checkSAT(solver, new VecInt()); }
    
    protected void checkSAT(ConstraintSolver solver, IVecInt asms) {
        if (getTimeout() != NO_TIMEOUT) {
            assert(getTimeout() > 0);
            solver.setTimeout((int)getRemainingTime());
        }
        solver.solve(asms);
    }
    
    protected void printUnsatisfiable() {
        System.out.println("c Instance is unsatisfiable");
        printElapsedTime();
    }
    
    protected void printOptimum() {
        System.out.println("c Proved optimality");
        printElapsedTime();
    }
    
    protected void printTimeoutMessage() {
        System.out.println("c TIMEOUT!");
        printElapsedTime();
    }
    
    protected void printUsedPMsCount(int used) {
        System.out.println("c Solution using " + used + " PMs found");
        printElapsedTime();
    }
    
    // =============================================================================================
    // ======================================== PUBLIC API =========================================
    // =============================================================================================
    
    public static enum HashFunctionType {
        NONE, GLOBAL, SEPARATED;
    }
    
    public ConstraintBasedAllocAlgorithm(VMCwMProblem instance) { super(instance); }
    public ConstraintBasedAllocAlgorithm(VMCwMProblem instance, boolean break_symms) {
        this(instance);
        this.break_symms = break_symms;
    }
    
    public void setHashFunctionType(HashFunctionType type) { this.hash_type = type; }
    
}
