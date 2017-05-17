package vmalloc.constraints;

import java.math.BigDecimal;
import java.math.BigInteger;

import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Clock;

/**
 * Abstract superclass for constraint solvers.
 * @author Miguel Terra-Neves
 */
public abstract class ConstraintSolver {
    
    /**
     * Integer constant used to represent undefined timeouts.
     */
    protected static final int NO_TIMEOUT = -1;
    
    /**
     * Time instant from which calls to {@link #solve()} or {@link #solve(IVecInt)} must terminate.
     * {@code timeout} is specified in seconds starting from the moment at which the application was
     * launched. If {@code timeout} is set to {@link #NO_TIMEOUT}, then no time limit is imposed.
     */
    private long timeout = NO_TIMEOUT;
    
    /**
     * Creates a new Boolean variable in the solver.
     * @see #nVars()
     */
    public abstract void newVar();
    
    /**
     * Creates multiple new Boolean variables in the solver.
     * @param nvars The number of variables to create.
     * @see #nVars()
     */
    public abstract void newVars(int nvars);
    
    /**
     * Retrieves the number of variables created in the solver. Constraints added to the solver
     * cannot have variables with indexes larger than the value returned by this method.
     * @return The number of variables created in the solver.
     */
    public abstract int nVars();
    
    /**
     * Adds a removable constraint stating that the sum of the literals in {@code lits} must be
     * equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract ConstraintID addRemovableExactly(IVecInt lits, int rhs)
            throws ContradictionException;
    
    /**
     * Adds a removable constraint stating that the sum of the literals in {@code lits} must be less
     * or equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract ConstraintID addRemovableAtMost(IVecInt lits, int rhs)
            throws ContradictionException;
    
    /**
     * Adds a removable constraint stating that the sum of the literals in {@code lits} must be
     * equal or larger than {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract ConstraintID addRemovableAtLeast(IVecInt lits, int rhs)
            throws ContradictionException;
    
    /**
     * Adds a constraint stating that the sum of the literals in {@code lits} must be equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract void addExactly(IVecInt lits, int rhs) throws ContradictionException;
    
    /**
     * Adds a constraint stating that the sum of the literals in {@code lits} must be less or equal
     * to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract void addAtMost(IVecInt lits, int rhs) throws ContradictionException;
    
    /**
     * Adds a constraint stating that the sum of the literals in {@code lits} must be equal or
     * larger than {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract void addAtLeast(IVecInt lits, int rhs) throws ContradictionException;
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract ConstraintID addRemovableEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException;
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract ConstraintID addRemovableGreaterOrEqual(
            IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) throws ContradictionException;
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract ConstraintID addRemovableLessOrEqual(
            IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) throws ContradictionException;
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract void addEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException;
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract void addGreaterOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException;
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public abstract void addLessOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException;
    
    /**
     * Auxiliary class used to store the result of applying a scaling factor to a set of
     * non-integral coefficients and rhs of a Pseudo-Boolean constraint.
     * @see ConstraintSolver#scaleToInteger(IVec, BigDecimal)
     */
    private class ScaledResult {
        
        /**
         * The coefficients after applying the scaling factor.
         */
        private IVec<BigInteger> coeffs;
        
        /**
         * The right-hand side after applying the scaling factor.
         */
        private BigInteger rhs;
        
        /**
         * Creates an instance of a {@code ScaledResult} object to hold the results of a scaling
         * operation.
         * @param coeffs The scaled coefficients.
         * @param rhs The scaled right-hand side.
         */
        ScaledResult(IVec<BigInteger> coeffs, BigInteger rhs) {
            this.coeffs = coeffs;
            this.rhs = rhs;
        }
        
        /**
         * Retrieves the scaled coefficients.
         * @return The scaled coefficients.
         */
        IVec<BigInteger> getCoefficients() { return this.coeffs; }
        
        /**
         * Retrieves the scaled right-hand side.
         * @return The scaled right-hand side.
         */
        BigInteger getRightHandSide() { return this.rhs; }
        
    }
    
    /**
     * Computes the scaling factor to be provided as an argument to
     * {@link BigDecimal#scaleByPowerOfTen(int)} in order to turn all coefficients and the
     * right-hand side into integers.
     * @param coeffs The coefficients to scale.
     * @param rhs The right-hand side to scale.
     * @return The scaling factor.
     */
    private int computeScaleFactorExponent(IVec<BigDecimal> coeffs, BigDecimal rhs) {
        assert(coeffs.size() > 0);
        int factor = rhs.scale();
        for (int i = 0; i < coeffs.size(); ++i) {
            int scale = coeffs.get(i).scale();
            factor = (scale > factor) ? scale : factor;
        }
        return factor;
    }
    
    /**
     * Given a set of possible non-integral coefficients and the right-hand side of a
     * Pseudo-Boolean constraint, converts those coefficients and right-hand side to integers in a
     * way that preserves the constraint's model set.
     * @param coeffs The coefficients to scale.
     * @param rhs The right-hand side to scale.
     * @return A {@link ScaledResult} object with the scaled integer coefficients and right-hand
     * side.
     */
    private ScaledResult scaleToInteger(IVec<BigDecimal> coeffs, BigDecimal rhs) {
        int factor = computeScaleFactorExponent(coeffs, rhs);
        IVec<BigInteger> int_coeffs = new Vec<BigInteger>();
        for (int i = 0; i < coeffs.size(); ++i) {
            int_coeffs.push(coeffs.get(i).scaleByPowerOfTen(factor).toBigIntegerExact());
        }
        assert(int_coeffs.size() == coeffs.size());
        return new ScaledResult(int_coeffs, rhs.scaleByPowerOfTen(factor).toBigIntegerExact());
    }
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public ConstraintID addRemovableEqual(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        return addRemovableEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide());
    }
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public ConstraintID addRemovableGreaterOrEqual(
            IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs) throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        return addRemovableGreaterOrEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide());
    }
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public ConstraintID addRemovableLessOrEqual(
            IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs) throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        return addRemovableLessOrEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide());
    }
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public void addEqual(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        addEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide());
    }
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public void addGreaterOrEqual(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        addGreaterOrEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide());
    }
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public void addLessOrEqual(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        addLessOrEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide());
    }
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater than
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public ConstraintID addRemovableGreater(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        return addRemovableGreaterOrEqual(lits,
                                          scaled.getCoefficients(),
                                          scaled.getRightHandSide().add(BigInteger.ONE));
    }
    
    /**
     * Adds a removable Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less than
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @return The constraint's id.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public ConstraintID addRemovableLess(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        return addRemovableLessOrEqual(lits,
                                       scaled.getCoefficients(),
                                       scaled.getRightHandSide().subtract(BigInteger.ONE));
    }
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater than
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public void addGreater(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        addGreaterOrEqual(lits, scaled.getCoefficients(), scaled.getRightHandSide().add(BigInteger.ONE));
    }
    
    /**
     * Adds a Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less than
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @throws ContradictionException If the solver detects that the addition of this constraint
     * would result in a contradiction.
     */
    public void addLess(IVecInt lits, IVec<BigDecimal> coeffs, BigDecimal rhs)
            throws ContradictionException {
        ScaledResult scaled = scaleToInteger(coeffs, rhs);
        addLessOrEqual(lits,
                       scaled.getCoefficients(),
                       scaled.getRightHandSide().subtract(BigInteger.ONE));
    }
    
    /**
     * Adds a removable constraint stating that at least one of the literals in {@code lits} must be
     * satisfied.
     * @param lits The clause's literals.
     * @return The clause's constraint id.
     * @throws ContradictionException If the solver detects that the addition of the clause would
     * result in a contradiction.
     */
    public abstract ConstraintID addRemovableClause(IVecInt lits) throws ContradictionException;
    
    /**
     * Adds a constraint stating that at least one of the literals in {@code lits} must be
     * satisfied.
     * @param lits The clause's literals.
     * @throws ContradictionException If the solver detects that the addition of the clause would
     * result in a contradiction.
     */
    public abstract void addClause(IVecInt lits) throws ContradictionException;
    
    /**
     * Removes a constraint from the solver.
     * @param id The id of the constraint to be removed.
     */
    public abstract void removeConstraint(ConstraintID id);
    
    /**
     * Removes a set of constraints from the solver.
     * @param ids The ids of the constraints to be removed.
     */
    public void removeConstraints(IVec<ConstraintID> ids) {
        for (int i = 0; i < ids.size(); ++i) {
            removeConstraint(ids.get(i));
        }
    }
    
    /**
     * Creates a vector of unit coefficients (of {@link BigInteger#ONE}).
     * @param ncoeffs The number of coefficients to create.
     * @return The vector of unit coefficients.
     */
    protected IVec<BigInteger> makeUnitCoeffs(int ncoeffs) {
        IVec<BigInteger> coeffs = new Vec<BigInteger>(ncoeffs);
        for (int i = 0; i < ncoeffs; ++i) {
            coeffs.push(BigInteger.ONE);
        }
        return coeffs;
    }
    
    /**
     * Sets the time seconds allotted for future calls to {@link #solve()} or
     * {@link #solve(IVecInt)}.
     * @param timeout The allotted time in seconds.
     */
    public void setTimeout(long timeout) {
        this.timeout = timeout + (long)Clock.getInstance().getElapsed();
    }
    
    /**
     * Retrieves the remaining time allotted for calls to {@link #solve()} or
     * {@link #solve(IVecInt)}. Default value is {@link Long#MAX_VALUE}.
     * @return The remaining time in seconds.
     */
    protected long getRemainingTime() {
        if (timeout == NO_TIMEOUT) {
            return Long.MAX_VALUE; // FIXME: implicit limit, should change eventually
        }
        assert(timeout >= 0);
        return timeout - (long)Clock.getInstance().getElapsed();
    }
    
    /**
     * Solves the constraint satisfaction problem under a given set of assumptions.
     * @param asms A set of literals that must be satisfied by the solution.
     */
    public abstract void solve(IVecInt asms);
    
    /**
     * Solves the constraint satisfaction problem.
     */
    public void solve() { solve(new VecInt()); }
    
    /**
     * Checks if the solver was able to solve the constraint satisfaction problem. Must be called
     * after {@link #solve()} or {@link #solve(IVecInt)}.
     * @return True if the problem was solved, false otherwise (in this case, usually due to
     * timeout).
     */
    public abstract boolean isSolved();
    
    /**
     * Checks if the constraint set is satisfiable. Must be called after {@link #solve()} or
     * {@link #solve(IVecInt)}.
     * @return True if it is satisfiable, false otherwise.
     */
    public abstract boolean isSatisfiable();
    
    /**
     * Checks if the constraint set is unsatisfiable. Must be called after {@link #solve()} or
     * {@link #solve(IVecInt)}.
     * @return True if it is unsatisfiable, false otherwise.
     */
    public abstract boolean isUnsatisfiable();
    
    /**
     * If the constraint set is satisfiable, returns the value of a given variable in the solution
     * found by the solver. Must be called after {@link #solve()} or {@link #solve(IVecInt)}.
     * @param var The variable.
     * @return True if {@code var} has value 1, false otherwise.
     * @see #isSatisfiable()
     */
    public abstract boolean modelValue(int var);
    
    /**
     * If the constraint set is unsatisfiable, returns a subset of the assumptions that are
     * responsible for unsatisfiability. Must be called after {@link #solve(IVecInt)}.
     * @return A vector with the assumptions literals responsible for unsatisfiability.
     * @see #isUnsatisfiable()
     */
    public abstract IVecInt unsatExplanation();
    
}
