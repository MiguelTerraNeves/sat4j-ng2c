package vmalloc.algorithm;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Utils;
import vmalloc.constraints.ConstraintID;
import vmalloc.constraints.ConstraintSolver;
import vmalloc.constraints.PseudoBooleanSolver;
import vmalloc.evolutionary.VMCwMProblem;

public class MCSAlloc extends HashEnumAlloc {
    
    public MCSAlloc(VMCwMProblem instance) { super(instance); }
    public MCSAlloc(VMCwMProblem instance, boolean break_symms) {
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
        printElapsedTime();
        return solver;
    }
    
    private IVecInt buildUndefFormulas(IVecInt pm_vars) {
        IVecInt undef_fmls = new VecInt();
        for (int i = 0; i < pm_vars.size(); ++i) {
            undef_fmls.push(-pm_vars.get(i));
        }
        return undef_fmls;
    }
    
    // FIXME: only checks if hash function type is NONE
    // FIXME: it'll run forever if no time limit is given when hash functions are enabled
    @Override
    public void allocate() {
        System.out.println("c WARNING: MCS minimizes number of servers, not energy consumption");
        System.out.println("c Initializing");
        ConstraintSolver solver = null;
        try {
            solver = buildSolver();
        }
        catch (ContradictionException e) {
            printUnsatisfiable();
            return;
        }
        int ub = this.instance.getPhysicalMachines().size()+1;
        IVecInt asms = new VecInt();
        IVecInt undef_fmls = buildUndefFormulas(this.pm_vars);
        boolean mcs_exists = false, next_mcs = false;
        IVec<ConstraintID> to_remove = new Vec<ConstraintID>();
        IVec<ConstraintID> hash_ids = new Vec<ConstraintID>(); // dummy vec
        IVecInt hash_asms = new VecInt();
        if (getHashType() != HashFunctionType.NONE) {
            System.out.println("c Setting initial hash function");
            //hash_ids = setHashFunction(solver);
            hash_ids = setHashFunction(solver,
                                       flattenLitVectors(flattenJobVars(this.job_vars)),
                                       hash_asms);
            assert(hash_asms.size() > 0);
        }
        while (true) {
            System.out.println("c Computing mapping");
            checkSAT(solver, hash_asms);
            if (!solver.isSolved()) {
                printTimeoutMessage();
                return;
            }
            else if (solver.isSatisfiable()) {
                mcs_exists = true;
                int new_ub = getUsedPMsCount(solver, this.instance.getPhysicalMachines(), this.job_vars);
                if (new_ub < ub) {
                    ub = new_ub;
                    printUsedPMsCount(ub);
                    saveSolution(modelToAllocation(solver,
                                                   this.instance.getPhysicalMachines(),
                                                   this.instance.getJobs(),
                                                   this.job_vars));
                }
                //asms.clear(); FIXME: sat4j bugs if I don't re-add the asms after removing hash function
                IVecInt new_asms = extractSatisfied(solver, undef_fmls);
                try {
                    //addRemovableConjunction(solver, asms, to_remove);
                    addRemovableConjunction(solver, new_asms, to_remove);
                    to_remove.push(solver.addRemovableClause(undef_fmls));
                    for (int i = 0; i < new_asms.size(); ++i) {
                        asms.push(new_asms.get(i));
                    }
                }
                catch (ContradictionException e) {
                    next_mcs = true;
                }
            }
            else {
                if (!mcs_exists) {
                    printOptimum();
                    return;
                }
                if (    hash_asms.size() > 0 &&
                        Utils.valuesIntersect(hash_asms, solver.unsatExplanation())) {
                    System.out.println("c Removing hash function");
                    IVecInt hash_neg_clause = new VecInt();
                    for (int i = 0; i < hash_asms.size(); ++i) {
                        hash_neg_clause.push(-hash_asms.get(i));
                    }
                    hash_asms.clear();
                    try {
                        hash_ids.push(solver.addRemovableClause(hash_neg_clause));
                    }
                    catch (ContradictionException e) {
                        next_mcs = true;
                    }
                }
                /*if (hash_ids.size() > 0) {
                    System.out.println("c Removing hash function");
                    solver.removeConstraints(hash_ids);
                    hash_ids.clear();
                    try { // FIXME: sat4j bugs if I don't re-add the asms after removing hash function
                        addRemovableConjunction(solver, asms, to_remove);
                    }
                    catch (ContradictionException ce) { // these asms were added previously; exception never thrown
                        throw new RuntimeException(ce);
                    }
                }*/
                else {
                    hash_asms.clear();
                    next_mcs = true;
                }
            }
            if (mcs_exists && next_mcs) {
                System.out.println("c MCS computed, generating another one");
                mcs_exists = next_mcs = false;
                solver.removeConstraints(to_remove);
                to_remove.clear();
                asms.clear();
                try {
                    solver.addClause(undef_fmls); // block MCS
                    undef_fmls = buildUndefFormulas(this.pm_vars);
                }
                catch (ContradictionException e) {
                    printOptimum();
                    return;
                }
                if (getHashType() != HashFunctionType.NONE) {
                    System.out.println("c Generating new hash function");
                    //hash_ids = setHashFunction(solver);
                    solver.removeConstraints(hash_ids);
                    hash_ids = setHashFunction(solver,
                                               flattenLitVectors(flattenJobVars(this.job_vars)),
                                               hash_asms);
                }
            }
        }
    }

}
