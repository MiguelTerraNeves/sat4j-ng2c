package vmalloc.algorithm;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.constraints.ConstraintID;
import vmalloc.constraints.ConstraintSolver;
import vmalloc.constraints.PseudoBooleanSolver;
import vmalloc.evolutionary.VMCwMProblem;
import vmalloc.exception.InvalidHashFunctionTypeException;

// TODO: new abstract class for algorithms that use hash functions with more generic methods
public class HashEnumAlloc extends ConstraintBasedAllocAlgorithm {

    protected IVecInt pm_vars = null;
    protected IVec<IVec<IVecInt>> job_vars = null;
    
    public HashEnumAlloc(VMCwMProblem instance) { super(instance); }
    public HashEnumAlloc(VMCwMProblem instance, boolean break_symms) {
        super(instance, break_symms);
    }
    
    protected ConstraintSolver buildSolver() throws ContradictionException {
        System.out.println("c Building formula");
        ConstraintSolver solver = new PseudoBooleanSolver();
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
                                     this.instance.getJobs(), this.job_vars);
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
        printElapsedTime();
        return solver;
    }
    
    private void blockSolution(ConstraintSolver solver, IVec<IVec<IVecInt>> job_vars)
            throws ContradictionException {
        IVecInt or_lits = new VecInt();
        for (int i = 0; i < job_vars.size(); ++i) {
            for (int j = 0; j < job_vars.get(i).size(); ++j) {
                IVecInt vm_vars = job_vars.get(i).get(j);
                for (int k = 0; k < vm_vars.size(); ++k) {
                    if (solver.modelValue(vm_vars.get(k))) {
                        or_lits.push(-vm_vars.get(k));
                        break;
                    }
                }
            }
        }
        solver.addClause(or_lits);
    }
    
    // FIXME: should be in ConstraintBasedAllocAlgorithm
    protected IVec<ConstraintID> setHashFunction(ConstraintSolver solver) {
        if (getHashType() == HashFunctionType.GLOBAL) {
            return setGlobalHashFunction(solver);
        }
        else if (getHashType() == HashFunctionType.SEPARATED) {
            return setSeparatedHashFunction(solver);
        }
        else if (getHashType() == HashFunctionType.NONE) {
            return new Vec<ConstraintID>();
        }
        else {
            throw new InvalidHashFunctionTypeException();
        }
    }
    
    // FIXME: should be in ConstraintBasedAllocAlgorithm
    private IVec<ConstraintID> setGlobalHashFunction(ConstraintSolver solver) {
        System.out.println("c Setting global hash function");
        IVecInt vars = flattenLitVectors(flattenJobVars(this.job_vars));
        IVec<ConstraintID> ids = setHashFunction(solver, vars);
        System.out.println("c Global hash function set");
        printElapsedTime();
        return ids;
    }
    
    // FIXME: should be in ConstraintBasedAllocAlgorithm
    private IVec<ConstraintID> setSeparatedHashFunction(ConstraintSolver solver) {
        System.out.println("c Setting separated hash functions");
        IVec<IVecInt> vm_vars = flattenJobVars(this.job_vars);
        IVec<ConstraintID> ids = new Vec<ConstraintID>();
        for (int i = 0; i < this.instance.getPhysicalMachines().size(); ++i) {
            IVecInt vars = new VecInt(vm_vars.size());
            for (int j = 0; j < vm_vars.size(); ++j) {
                vars.unsafePush(vm_vars.get(j).get(i));
            }
            setHashFunction(solver, vars).copyTo(ids);;
        }
        System.out.println("c Separated hash functions set");
        printElapsedTime();
        return ids;
    }
    
    @Override
    public void allocate() {
        System.out.println("c Initializing");
        ConstraintSolver solver = null;
        try {
            solver = buildSolver();
        }
        catch (ContradictionException e) {
            printUnsatisfiable();
            return;
        }
        int enum_threshold = getEnumerationThreshold();
        IVec<ConstraintID> to_remove = setHashFunction(solver);
        int nsols_in_cell = 0;
        while (true) {
            checkSAT(solver);
            if (!solver.isSolved()) {
                printTimeoutMessage();
                return;
            }
            else if (solver.isSatisfiable()) {
                ++nsols_in_cell;
                saveSolution(modelToAllocation(solver,
                                               this.instance.getPhysicalMachines(),
                                               this.instance.getJobs(),
                                               this.job_vars),
                             true);
                try {
                    blockSolution(solver, this.job_vars);
                }
                catch (ContradictionException ce) {
                    nsols_in_cell = enum_threshold; // force new hash function
                }
            }
            if (!solver.isSatisfiable() && getHashType() == HashFunctionType.NONE) {
                return;
            }
            if (    !solver.isSatisfiable() ||
                    (getHashType() != HashFunctionType.NONE && nsols_in_cell >= enum_threshold)) {
                assert(to_remove.size() > 0);
                solver.removeConstraints(to_remove);
                to_remove = setHashFunction(solver);
                nsols_in_cell = 0;
            }
        }
    }

}
