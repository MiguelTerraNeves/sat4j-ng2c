package vmalloc.algorithm;

import static org.sat4j.GlobalDefs.USE_NG2C;
import static org.sat4j.GlobalDefs.ANNOTATE_SOLVER_STRUCTS;

import org.sat4j.core.Vec;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.constraints.ConstraintID;
import vmalloc.constraints.ConstraintSolver;
import vmalloc.evolutionary.VMCwMProblem;

public class ParetoCLD extends MultiObjectiveConstraintBasedAllocAlgorithm {

    public ParetoCLD(VMCwMProblem instance) { this(instance, false); }
    public ParetoCLD(VMCwMProblem instance, boolean break_symms) {
        super(instance, break_symms);
    }

    @Override
    public void allocate() {
        System.out.println("c Initializing");
        if (USE_NG2C) {
            System.newAllocGen();
        }
        ConstraintSolver solver = null;
        try {
            solver = buildSolver();
        }
        catch (ContradictionException e) {
            printUnsatisfiable();
            return;
        }
        System.out.println("c Initializing objective functions");
        initializeObjectiveFunctions();
        printElapsedTime();
        IVecInt undef_fmls = buildUndefFormulas();
        boolean mcs_exists = false, next_mcs = false;
        IVec<ConstraintID> to_remove = ANNOTATE_SOLVER_STRUCTS ? new @Gen Vec<ConstraintID>() : new Vec<ConstraintID>();
        int perm_gen;
        if (USE_NG2C) {
            perm_gen = System.getAllocGen();
            System.newAllocGen();
        }
        while (true) {
            System.out.println("c Computing mapping");
            checkSAT(solver);
            System.out.println("c Solver returned");
            printElapsedTime();
            if (!solver.isSolved()) {
                printTimeoutMessage();
                return;
            }
            else if (solver.isSatisfiable()) {
                mcs_exists = true;
                int current_gen;
                if (USE_NG2C) {
                    current_gen = System.getAllocGen();
                    System.setAllocGen(perm_gen);
                }
                saveSolution(modelToAllocation(solver,
                                               this.instance.getPhysicalMachines(),
                                               this.instance.getJobs(),
                                               this.job_vars),
                             true);
                if (USE_NG2C) {
                    System.setAllocGen(current_gen);
                }
                IVecInt asms = extractSatisfied(solver, undef_fmls);
                try {
                    addRemovableConjunction(solver, asms, to_remove);
                    to_remove.push(solver.addRemovableClause(undef_fmls));
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
                else {
                    next_mcs = true;
                }
            }
            if (mcs_exists && next_mcs) {
                System.out.println("c MCS computed, generating another one");
                mcs_exists = next_mcs = false;
                solver.removeConstraints(to_remove);
                to_remove.clear();
                try {
                    int current_gen;
                    if (USE_NG2C) {
                        current_gen = System.getAllocGen();
                        System.setAllocGen(perm_gen);
                    }
                    solver.addClause(undef_fmls); // block MCS
                    undef_fmls = buildUndefFormulas();
                    if (USE_NG2C) {
                        System.setAllocGen(current_gen);
                    }
                }
                catch (ContradictionException e) {
                    printOptimum();
                    return;
                }
                if (USE_NG2C) {
                    System.newAllocGen();
                }
            }
        }
    }

}
