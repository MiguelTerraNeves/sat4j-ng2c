package vmalloc.algorithm;

import java.math.BigInteger;

import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.constraints.INewBestHandler;
import vmalloc.constraints.PBOSolver;
import vmalloc.evolutionary.VMCwMProblem;

public class PBOAlloc extends ConstraintBasedAllocAlgorithm {

    protected IVecInt pm_vars = null;
    protected IVec<IVec<IVecInt>> job_vars = null;
    
    public PBOAlloc(VMCwMProblem instance) { this(instance, false); }
    public PBOAlloc(VMCwMProblem instance, boolean break_symms) {
        super(instance, break_symms);
    }
    
    protected PBOSolver buildSolver() throws ContradictionException {
        System.out.println("c Building formula");
        PBOSolver solver = new PBOSolver();
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
        solver.setObjectiveFunction(this.pm_vars);
        printElapsedTime();
        return solver;
    }
    
    @Override
    public void allocate() {
        System.out.println("c WARNING: PBO minimizes number of servers, not energy consumption");
        PBOSolver solver;
        try {
            solver = buildSolver();
        }
        catch (ContradictionException e) {
            printUnsatisfiable();
            return;
        }
        solver.setTimeout(getRemainingTime());
        solver.solve(new INewBestHandler() {
            @Override
            public void handleNewBest(BigInteger best) {
                printUsedPMsCount(best.intValue());
            }
        });
        if (solver.isSolved()) {
            if (solver.isSatisfiable()) {
                saveSolution(modelToAllocation(solver,
                                               this.instance.getPhysicalMachines(),
                                               this.instance.getJobs(),
                                               job_vars));
            }
            if (solver.foundOptimum()) {
                printOptimum();
            }
        }
    }

}
