package vmalloc.algorithm;

import java.math.BigInteger;

import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.constraints.ExternalWBOSolver;
import vmalloc.constraints.INewBestHandler;
import vmalloc.evolutionary.VMCwMProblem;

public class WBOAlloc extends ConstraintBasedAllocAlgorithm {
    
    protected IVecInt pm_vars = null;
    protected IVec<IVec<IVecInt>> job_vars = null;

    public WBOAlloc(VMCwMProblem instance) { this(instance, false); }
    public WBOAlloc(VMCwMProblem instance, boolean break_symms) {
        super(instance, break_symms);
    }
    
    private void addSoftPMOffConstraints(ExternalWBOSolver solver, IVecInt pm_vars) {
        for (int i = 0; i < pm_vars.size(); ++i) {
            IVecInt vars = new VecInt();
            vars.push(pm_vars.get(i));
            solver.addSoftAtMost(vars, 0, BigInteger.ONE);
        }
    }
    
    protected ExternalWBOSolver buildSolver() throws ContradictionException {
        System.out.println("c Building formula");
        ExternalWBOSolver solver = new ExternalWBOSolver();
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
        addSoftPMOffConstraints(solver, this.pm_vars);
        printElapsedTime();
        return solver;
    }

    @Override
    public void allocate() {
        System.out.println("c WARNING: WBO minimizes number of servers, not energy consumption");
        ExternalWBOSolver solver;
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
