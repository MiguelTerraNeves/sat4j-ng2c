package vmalloc.constraints;

import java.math.BigInteger;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.pb.IPBSolver;
import org.sat4j.pb.SolverFactory;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IConstr;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;
import org.sat4j.specs.TimeoutException;

/**
 * Constraint solver wrapper around SAT4J for solving instances of the Pseudo-Boolean Satisfaction
 * (PBS) problem. Works incrementally.
 * @author Miguel Terra-Neves
 */
public class PseudoBooleanSolver extends ConstraintSolver {
    
    /**
     * Abstract superclass for constraint objects. This class and respective subclasses are
     * auxiliary and are required because SAT4J does not add constraints satisfied by unit
     * propagation (this is not documented). If constraints are removed using
     * {@link PseudoBooleanSolver#removeConstraint(ConstraintID)} or
     * {@link PseudoBooleanSolver#removeConstraints(IVec)}, then we need to add those that were
     * satisfied by unit propagation, or else the solver might produce incorrect results in
     * future calls to {@link PseudoBooleanSolver#solve()} or
     * {@link PseudoBooleanSolver#solve(IVecInt)}.
     */
    private abstract class Constraint {
        
        /**
         * Literals in the constraint.
         */
        private IVecInt lits = null;
        
        /**
         * Creates an instance of a constraint.
         * @param lits The literals in the constraint.
         */
        Constraint(IVecInt lits) {
            this.lits = new @Gen VecInt(lits.size());
            lits.copyTo(this.lits);
        }
        
        /**
         * Retrieves the literals in the constraint.
         * @return The constraint's literals.
         */
        protected IVecInt getLits() { return this.lits; }
        
        /**
         * Adds the constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        abstract IConstr addRemovable(IPBSolver solver) throws ContradictionException;
        
    }
    
    /**
     * Abstract superclass for constraint objects with a right-hand side.
     */
    private abstract class ConstraintWithRHS extends Constraint {
        
        /**
         * The right-hand side of the constraint.
         */
        private BigInteger rhs = null;
        
        /**
         * Creates an instance of a constraint with a right-hand side.
         * @param lits The constraint's literals.
         * @param rhs The constraint's right-hand side.
         */
        ConstraintWithRHS(IVecInt lits, BigInteger rhs) {
            super(lits);
            this.rhs = rhs;
        }
        
        /**
         * Retrieves the right-hand side of the constraint.
         * @return The constraint's right-hand side.
         */
        protected BigInteger getRHS() { return this.rhs; }
        
    }
    
    /**
     * Abstract superclass for constraint objects with coefficients (like Pseudo-Boolean
     * constraints) and a right-hand side.
     */
    private abstract class ConstraintWithCoeffsAndRHS extends ConstraintWithRHS {
        
        /**
         * The coefficients in the constraint.
         */
        private IVec<BigInteger> coeffs = null;
        
        /**
         * Creates an instance of a constraint with coefficients and a right-hand side.
         * @param lits The literals in the constraint.
         * @param coeffs The coefficients in the constraint.
         * @param rhs The constraint's right-hand side.
         */
        ConstraintWithCoeffsAndRHS(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) {
            super(lits, rhs);
            this.coeffs = new @Gen Vec<BigInteger>(coeffs.size());
            coeffs.copyTo(this.coeffs);
        }
        
        /**
         * Retrieves the coefficients of the constraint.
         * @return The constraint's coefficients.
         */
        protected IVec<BigInteger> getCoeffs() { return this.coeffs; }
        
    }
    
    /**
     * Constraint object that represents a disjunction of literals.
     */
    private class Clause extends Constraint {
        
        /**
         * Creates an instance of a clause.
         * @param lits The clause's literals.
         */
        Clause(IVecInt lits) { super(lits); }
        
        /**
         * Adds the clause to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the clause.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addClause(getLits());
        }
        
    }
    
    /**
     * Constraint object that represents the constraint that the sum of a given set of literals
     * must be less or equal to a given value.
     */
    private class AtMost extends ConstraintWithRHS {
        
        /**
         * Creates an instance of an at-most constraint.
         * @param lits The literals in the at-most constraint.
         * @param rhs The constraint's right-hand side.
         */
        AtMost(IVecInt lits, int rhs) { super(lits, BigInteger.valueOf(rhs)); }
        
        /**
         * Adds the at-most constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the at-most constraint.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addAtMost(getLits(), getRHS().intValueExact());
        }
        
    }
    
    /**
     * Constraint object that represents the constraint that the sum of a given set of literals
     * must be equal or greater than a given value.
     */
    private class AtLeast extends ConstraintWithRHS {
        
        /**
         * Creates an instance of an at-least constraint.
         * @param lits The literals in the at-least constraint.
         * @param rhs The constraint's right-hand side.
         */
        AtLeast(IVecInt lits, int rhs) { super(lits, BigInteger.valueOf(rhs)); }
        
        /**
         * Adds the at-least constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the at-least constraint.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addAtLeast(getLits(), getRHS().intValueExact());
        }
        
    }
    
    /**
     * Constraint object that represents the constraint that the sum of a given set of literals
     * must be equal to a given value.
     */
    private class Exactly extends ConstraintWithRHS {
        
        /**
         * Creates an instance of an exactly constraint.
         * @param lits The literals in the exactly constraint.
         * @param rhs The constraint's right-hand side.
         */
        Exactly(IVecInt lits, int rhs) { super(lits, BigInteger.valueOf(rhs)); }
        
        /**
         * Adds the exactly constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the exactly constraint.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addExactly(getLits(), getRHS().intValueExact());
        }
        
    }
    
    /**
     * Constraint object that represents the constraint that the sum of a given set of literals
     * times the corresponding coefficients must be greater or equal to a given value.
     */
    private class GreaterOrEqual extends ConstraintWithCoeffsAndRHS {
        
        /**
         * Creates an instance of a greater or equal constraint.
         * @param lits The constraint's literals.
         * @param coeffs The constraint's coefficients.
         * @param rhs The constraint's right-hand side.
         */
        public GreaterOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) {
            super(lits, coeffs, rhs);
        }
        
        /**
         * Adds the greater or equal constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the greater or equal constraint.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addPseudoBoolean(getLits(), getCoeffs(), true, getRHS());
        }
        
    }
    
    /**
     * Constraint object that represents the constraint that the sum of a given set of literals
     * times the corresponding coefficients must be less or equal to a given value.
     */
    private class LessOrEqual extends ConstraintWithCoeffsAndRHS {
        
        /**
         * Creates an instance of a less or equal constraint.
         * @param lits The constraint's literals.
         * @param coeffs The constraint's coefficients.
         * @param rhs The constraint's right-hand side.
         */
        public LessOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) {
            super(lits, coeffs, rhs);
        }
        
        /**
         * Adds the less or equal constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the less or equal constraint.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addPseudoBoolean(getLits(), getCoeffs(), false, getRHS());
        }
        
    }
    
    /**
     * Constraint object that represents the constraint that the sum of a given set of literals
     * times the corresponding coefficients must be equal to a given value.
     */
    private class Equal extends ConstraintWithCoeffsAndRHS {
        
        /**
         * Creates an instance of an equal constraint.
         * @param lits The constraint's literals.
         * @param coeffs The constraint's coefficients.
         * @param rhs The constraint's right-hand side.
         */
        public Equal(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) {
            super(lits, coeffs, rhs);
        }
        
        /**
         * Adds the equal constraint to a SAT4J PBS solver.
         * @param solver The SAT4J PBS solver.
         * @return A SAT4J constraint object for the equal constraint.
         * @throws ContradictionException If the solver detects that the addition of this constraint
         * would result in a contradiction.
         */
        @Override
        IConstr addRemovable(IPBSolver solver) throws ContradictionException {
            return solver.addExactly(getLits(), getCoeffs(), getRHS());
        }
        
    }
    
    /**
     * An instance of the underlying SAT4J PBS solver.
     */
    protected IPBSolver solver = null;
    
    /**
     * Boolean used to store the satisfiability of the PBS instance of the last call to
     * {@link #solve()} or {@link #solve(IVecInt)}. True if the instance is satisfiable, false
     * otherwise.
     */
    protected boolean is_sat = false;
    
    /**
     * Boolean used to store if the PBS instance was solved successfully on the last call to
     * {@link #solve()} or {@link #solve(IVecInt)}. True if so, false otherwise.
     */
    protected boolean is_solved = false;
    
    /**
     * Mapping of constraint ids to SAT4J's {@link IConstr} objects. In SAT4J, constraints are
     * removed by providing the {@link IConstr} object, returned by the add methods, to the
     * {@link IPBSolver#removeConstr(IConstr)} method.
     */
    private Map<ConstraintID, IConstr> rem_map = new @Gen HashMap<ConstraintID, IConstr>();
    
    /**
     * Map used to store constraints that SAT4J satisfied by unit propagation.
     * @see #removeConstraint(ConstraintID)
     */
    private Map<ConstraintID, Constraint> unit_sat_constraints = new @Gen HashMap<ConstraintID, Constraint>();
    
    /**
     * Creates an instance of a Pseudo-Boolean Satisfaction solver.
     */
    public PseudoBooleanSolver() { solver = SolverFactory.newDefault(); }

    @Override
    public void newVar() { newVars(1); }

    @Override
    public void newVars(int nvars) { solver.newVar(solver.nVars() + nvars); }

    @Override
    public int nVars() { return solver.nVars(); }
    
    /**
     * Creates a fresh constraint id for a given removable constraint object and stores a mapping
     * between that id and the corresponding SAT4J {@link IConstr} object. If SAT4J satisfied the
     * constraint by unit propagation, then the constraint is stored as well. If some of the
     * constraints are removed in the future, then the constraints satisfied by unit propagation
     * must be added again to the SAT4J solver for correctness.
     * @param constr The SAT4J {@link IConstr} object. If this is null, then the constraint was
     * satisfied by unit propagation.
     * @param constraint The constraint object to be stored if SAT4J satisfied the constraint by
     * unit propagation.
     * @return The id assigned to the constraint.
     * @see {@link #removeConstraint(ConstraintID)} {@link #removeConstraints(IVec)}
     */
    private ConstraintID storeRemovable(IConstr constr, Constraint constraint) {
        ConstraintID id = ConstraintID.makeFresh();
        rem_map.put(id, constr);
        if (constr == null) {
            unit_sat_constraints.put(id, constraint);
        }
        return id;
    }

    @Override
    public ConstraintID addRemovableExactly(IVecInt lits, int rhs) throws ContradictionException {
        return storeRemovable(solver.addExactly(lits, rhs), new @Gen Exactly(lits, rhs));
    }

    @Override
    public ConstraintID addRemovableAtMost(IVecInt lits, int rhs) throws ContradictionException {
        return storeRemovable(solver.addAtMost(lits, rhs), new @Gen AtMost(lits, rhs));
    }

    @Override
    public ConstraintID addRemovableAtLeast(IVecInt lits, int rhs) throws ContradictionException {
        return storeRemovable(solver.addAtLeast(lits, rhs), new @Gen AtLeast(lits, rhs));
    }

    @Override
    public void addExactly(IVecInt lits, int rhs) throws ContradictionException {
        solver.addExactly(lits, rhs);
    }

    @Override
    public void addAtMost(IVecInt lits, int rhs) throws ContradictionException {
        solver.addAtMost(lits, rhs);
    }

    @Override
    public void addAtLeast(IVecInt lits, int rhs) throws ContradictionException {
        solver.addAtLeast(lits, rhs);
    }

    @Override
    public ConstraintID addRemovableEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        return storeRemovable(solver.addExactly(lits, coeffs, rhs), new @Gen Equal(lits, coeffs, rhs));
    }

    @Override
    public ConstraintID addRemovableGreaterOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        return storeRemovable(solver.addPseudoBoolean(lits, coeffs, true, rhs),
                              new @Gen GreaterOrEqual(lits, coeffs, rhs));
    }

    @Override
    public ConstraintID addRemovableLessOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        return storeRemovable(solver.addPseudoBoolean(lits, coeffs, false, rhs),
                              new @Gen LessOrEqual(lits, coeffs, rhs));
    }

    @Override
    public void addEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        solver.addExactly(lits, coeffs, rhs);
    }

    @Override
    public void addGreaterOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        solver.addPseudoBoolean(lits, coeffs, true, rhs);
    }

    @Override
    public void addLessOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        solver.addPseudoBoolean(lits, coeffs, false, rhs);
    }

    @Override
    public ConstraintID addRemovableClause(IVecInt lits) throws ContradictionException {
        return storeRemovable(solver.addClause(lits), new @Gen Clause(lits));
    }

    @Override
    public void addClause(IVecInt lits) throws ContradictionException {
        solver.addClause(lits);
    }

    // FIXME: removeConstraints' implementation is inefficient: the set of constraints satisfied by
    // unit propagation is satisfied every time a single constraint is removed; not problematic so
    // far, but should be optimized if it becomes a bottleneck
    @Override
    public void removeConstraint(ConstraintID id) {
        assert(rem_map.containsKey(id));
        IConstr constr_obj = rem_map.get(id);
        if (constr_obj != null) { // not documented, but sat4j may return null if the constraint was satisfied by unit propagation
            solver.removeConstr(constr_obj);
            // Re-add constraints previously satisfied by unit propagation
            Iterator<ConstraintID> it = unit_sat_constraints.keySet().iterator();
            while (it.hasNext()) {
                ConstraintID sat_id = it.next();
                Constraint constraint = unit_sat_constraints.get(sat_id);
                try {
                    IConstr constr = constraint.addRemovable(solver);
                    if (constr != null) {
                        rem_map.put(sat_id, constr);
                        it.remove();
                    }
                }
                catch (ContradictionException ce) { // this will never happen
                    throw new RuntimeException(ce); // if formula was satisfied with the constraint, then it also is without
                }
            }
        }
        else {
            unit_sat_constraints.remove(id);
        }
        rem_map.remove(id);
    }
    
    @Override
    protected long getRemainingTime() {
        long timeout = super.getRemainingTime();
        return (timeout > Integer.MAX_VALUE) ? Integer.MAX_VALUE : timeout; // FIXME: implicit limit, should change eventually
    }

    @Override
    public void solve(IVecInt asms) {
        solver.setTimeout((int)getRemainingTime());
        try {
            is_sat = solver.isSatisfiable(asms);
            is_solved = true;
        }
        catch (TimeoutException e) {
            is_solved = false;
        }
    }

    @Override
    public boolean isSolved() { return is_solved; }

    @Override
    public boolean isSatisfiable() { return is_sat; }

    @Override
    public boolean isUnsatisfiable() { return !isSatisfiable(); }

    @Override
    public boolean modelValue(int var) { return solver.model(var); }
    
    @Override
    public IVecInt unsatExplanation() {
        try {
            IVecInt explanation = solver.unsatExplanation();
            return explanation == null ? new VecInt() : explanation;
        }
        catch (NullPointerException npe) { /* left empty on purpose, bug in sat4j may cause a NullPointerException */ }
        return new VecInt();
    }

}
