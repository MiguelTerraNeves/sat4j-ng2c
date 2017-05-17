package vmalloc.algorithm;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.constraints.ConstraintID;
import vmalloc.constraints.ConstraintSolver;
import vmalloc.evolutionary.VMCwMProblem;

public class ParetoLBX extends MultiObjectiveConstraintBasedAllocAlgorithm {

    public ParetoLBX(VMCwMProblem instance) { this(instance, false); }
    public ParetoLBX(VMCwMProblem instance, boolean break_symms) {
        super(instance, break_symms);
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
        System.out.println("c Initializing objective functions");
        initializeObjectiveFunctions();
        printElapsedTime();
        IVecInt undef_fmls = buildUndefFormulas(), asms = new VecInt(), mcs = new VecInt();
        IVec<ConstraintID> to_remove = new Vec<ConstraintID>();
        while (true) {
            System.out.println("c Computing mapping");
            checkSAT(solver, asms);
            if (!solver.isSolved()) {
                printTimeoutMessage();
                return;
            }
            else if (solver.isSatisfiable()) {
                saveSolution(modelToAllocation(solver,
                                               this.instance.getPhysicalMachines(),
                                               this.instance.getJobs(),
                                               this.job_vars),
                             true);
                IVecInt satisfied = extractSatisfied(solver, undef_fmls);
                try {
                    addRemovableConjunction(solver, satisfied, to_remove); // FIXME: could use assumptions instead (consider giving it a try)
                }
                catch (ContradictionException e) {
                    assert(false); // should not happen
                }
            }
            else {
                if (asms.size() == 0) {
                    assert(mcs.size() == 0);
                    printOptimum();
                    return;
                }
                else {
                    try {
                        to_remove.push(solver.addRemovableClause(new VecInt(new int[] { -asms.get(0) })));
                        mcs.push(asms.get(0));
                    }
                    catch (ContradictionException e) {
                        assert(false); // should not happen
                    }
                }
            }
            if (undef_fmls.size() > 0) {
                asms = new VecInt(new int[] { undef_fmls.last() });
                undef_fmls.pop();
            }
            else {
                assert(mcs.size() > 0);
                System.out.println("c MCS computed, generating another one");
                solver.removeConstraints(to_remove);
                to_remove.clear();
                asms.clear();
                try {
                    solver.addClause(mcs);
                    undef_fmls = buildUndefFormulas();
                    mcs.clear();
                }
                catch (ContradictionException e) {
                    printOptimum();
                    return;
                }
            }
        }
    }

}
