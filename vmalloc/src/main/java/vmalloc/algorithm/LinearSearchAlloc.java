package vmalloc.algorithm;

import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.constraints.ConstraintSolver;
import vmalloc.constraints.PseudoBooleanSolver;
import vmalloc.evolutionary.VMCwMProblem;

public class LinearSearchAlloc extends ConstraintBasedAllocAlgorithm {
    
    protected IVecInt pm_vars = null;
    protected IVec<IVec<IVecInt>> job_vars = null;
    
    public LinearSearchAlloc(VMCwMProblem instance) { super(instance); }
    public LinearSearchAlloc(VMCwMProblem instance, boolean break_symms) {
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
    
    @Override
    public void allocate() {
        System.out.println("c WARNING: LS minimizes number of servers, not energy consumption");
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
        while (true) {
            System.out.println("c Computing mapping");
            checkSAT(solver);
            if (!solver.isSolved()) {
                printTimeoutMessage();
                return;
            }
            else if (solver.isSatisfiable()) {
                ub = getUsedPMsCount(solver, this.instance.getPhysicalMachines(), this.job_vars);
                printUsedPMsCount(ub);
                saveSolution(modelToAllocation(solver,
                                               this.instance.getPhysicalMachines(),
                                               this.instance.getJobs(),
                                               this.job_vars));
                try {
                    solver.addAtMost(this.pm_vars, ub-1);
                }
                catch (ContradictionException e) {
                    printOptimum();
                    return;
                }
            }
            else {
                printOptimum();
                return;
            }
        }
    }

}
