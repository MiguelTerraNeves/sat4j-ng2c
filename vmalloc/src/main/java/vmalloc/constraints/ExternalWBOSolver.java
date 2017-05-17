package vmalloc.constraints;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringWriter;
import java.io.Writer;
import java.math.BigInteger;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.sat4j.core.Vec;
import org.sat4j.core.VecInt;
import org.sat4j.specs.ContradictionException;
import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.Configuration;
import vmalloc.Utils;
import vmalloc.exception.NotSupportedException;

/**
 * Constrained optimization solver wrapper around an external CLASP solver for solving Weighted
 * Boolean Optimization (WBO) problems. Does not support constraint removal and explanations of
 * unsatisfiability.
 * @author Miguel Terra-Neves
 */
public class ExternalWBOSolver extends ConstraintSolver {

    /**
     * Prefix for temporary files containing WBO instances to be solved by CLASP.
     */
    private static final String CLASP_INPUT_PREFIX = "vmalloc";
    
    /**
     * CLASP's WBO representation of greater or equal.
     */
    private static final String GREATER_OR_EQUAL = ">=";
    
    /**
     * CLASP's WBO representation of equal.
     */
    private static final String EQUAL = "=";
    
    /**
     * CLASP's representation that the optimum of an instance was found.
     */
    private static final String OPTIMUM = "s OPTIMUM";
    
    /**
     * CLASP's representation that a given instance is unsatisfiable.
     */
    private static final String UNSATISFIABLE = "s UNSATISFIABLE";
    
    /**
     * The WBO format includes a header with the minimum, maximum and total cost of the soft
     * constraints. {@code NO_COST} is a {@link BigInteger} constant that represents that these
     * costs have not been initialized yet.
     */
    private static final BigInteger NO_COST = BigInteger.ONE.negate();
    
    /**
     * {@link BigInteger} constant used to represent the case where a solution has not been found
     * yet.
     */
    private static final BigInteger NO_BEST = BigInteger.ONE.negate();
    
    /**
     * Counts the number of variables in the problem instance.
     */
    private int nvars = 0;
    
    /**
     * Counts the number of hard constraints in the problem instance.
     */
    private int nhard = 0;
    
    /**
     * Counts the number of soft constraints in the problem instance.
     */
    private int nsoft = 0;
    
    /**
     * The minimum cost among soft constraints.
     */
    private BigInteger mincost = NO_COST;
    
    /**
     * The maximum cost among soft constraints.
     */
    private BigInteger maxcost = NO_COST;
    
    /**
     * The total cost of the soft constraints.
     */
    private BigInteger sumcost = BigInteger.ZERO;
    
    /**
     * {@link StringWriter} object that accumulates the hard constraints, in the WBO file format, to
     * be included in the WBO instance file.
     */
    private StringWriter hard_writer = new StringWriter();
    
    /**
     * {@link StringWriter} object that accumulates the soft constraints, in the WBO file format, to
     * be included in the WBO instance file.
     */
    private StringWriter soft_writer = new StringWriter();
    
    /**
     * Boolean used to store the satisfiability of the WBO instance of the last call to
     * {@link #solve()} or {@link #solve(INewBestHandler)}. True if the instance is satisfiable,
     * false otherwise.
     */
    private boolean is_sat = false;
    
    /**
     * Boolean used to store if at least one solution was found on the last call to
     * {@link #solve()} or {@link #solve(INewBestHandler)}. True if so, false otherwise.
     */
    private boolean found_solution = false;
    
    /**
     * Boolean used to store if the solver found an optimal solution on the last call to
     * {@link #solve()} or {@link #solve(INewBestHandler)}. True if so, false otherwise. 
     */
    private boolean found_optimum = false;
    
    /**
     * Stores the cost of the best solution found on the last call to {@link #solve()} or
     * {@link #solve(INewBestHandler)}. If no solution was found, the default value is
     * {@link #NO_BEST}.
     */
    private BigInteger best = NO_BEST;
    
    /**
     * Stores the model of the best solution found on the last call to {@link #solve()} or
     * {@link #solve(INewBestHandler)}.
     */
    private Set<Integer> model = new HashSet<Integer>();
    
    @Override
    public void newVar() { ++nvars; }

    @Override
    public void newVars(int nvars) { this.nvars += nvars; }

    @Override
    public int nVars() { return nvars; }
    
    /**
     * Converts a variable to its string representation.
     * @param var The variable.
     * @return The string representation of variable {@code var}.
     */
    private String varToString(int var) { return "x" + var; }
    
    /**
     * Writes a Pseudo-Boolean constraint to a given writer object.
     * @param writer The writer object where the constraint is to be written.
     * @param lits The constraint's literals.
     * @param coeffs The constraint's coefficients.
     * @param rhs The constraint's right-hand side.
     * @param op The relation operator of the constraint (like "=" or ">=").
     */
    private void writePseudoBoolean(Writer writer,
                                    IVecInt lits,
                                    IVec<BigInteger> coeffs,
                                    BigInteger rhs,
                                    String op) {
        assert(lits.size() == coeffs.size());
        assert(op.equals(GREATER_OR_EQUAL) || op.equals(EQUAL));
        IVecInt vars = new VecInt(lits.size());
        IVec<BigInteger> norm_coeffs = new Vec<BigInteger>(coeffs.size());
        BigInteger norm_rhs = rhs;
        for (int i = 0; i < lits.size(); ++i) {
            int lit = lits.get(i);
            if (lit < 0) {
                vars.push(-lit);
                norm_coeffs.push(coeffs.get(i).negate());
                norm_rhs = norm_rhs.subtract(coeffs.get(i));
            }
            else {
                vars.push(lit);
                norm_coeffs.push(coeffs.get(i));
            }
        }
        assert(lits.size() == vars.size());
        try {
            for (int i = 0; i < vars.size(); ++i) {
                writer.write(norm_coeffs.get(i) + " " + varToString(vars.get(i)) + " ");
            }
            writer.write(op + " " + norm_rhs + " ;\n");
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to write PseudoBoolean constraint", e);
        }
    }
    
    /**
     * Writes an equal constraint to a given writer object.
     * @param writer The writer object where the constraint is to be written.
     * @param lits The constraint's literals.
     * @param coeffs The constraint's coefficients.
     * @param rhs The constraint's right-hand side.
     */
    private void writeEqual(Writer writer, IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs) {
        writePseudoBoolean(writer, lits, coeffs, rhs, EQUAL);
    }
    
    /**
     * Writes a greater or equal constraint to a given writer object.
     * @param writer The writer object where the constraint is to be written.
     * @param lits The constraint's literals.
     * @param coeffs The constraint's coefficients.
     * @param rhs The constraint's right-hand side.
     */
    private void writeGreaterOrEqual(Writer writer,
                                     IVecInt lits,
                                     IVec<BigInteger> coeffs,
                                     BigInteger rhs) {
        writePseudoBoolean(writer, lits, coeffs, rhs, GREATER_OR_EQUAL);
    }
    
    /**
     * Writes a less or equal constraint to a given writer object.
     * @param writer The writer object where the constraint is to be written.
     * @param lits The constraint's literals.
     * @param coeffs The constraint's coefficients.
     * @param rhs The constraint's right-hand side.
     */
    private void writeLessOrEqual(Writer writer,
                                  IVecInt lits,
                                  IVec<BigInteger> coeffs,
                                  BigInteger rhs) {
        IVec<BigInteger> neg_coeffs = new Vec<BigInteger>(coeffs.size());
        for (int i = 0; i < coeffs.size(); ++i) {
            neg_coeffs.push(coeffs.get(i).negate());
        }
        writeGreaterOrEqual(writer, lits, neg_coeffs, rhs.negate());
    }
    
    @Override
    public void addExactly(IVecInt lits, int rhs) throws ContradictionException {
        addEqual(lits, makeUnitCoeffs(lits.size()), BigInteger.valueOf(rhs));
    }

    @Override
    public void addAtMost(IVecInt lits, int rhs) throws ContradictionException {
        addLessOrEqual(lits, makeUnitCoeffs(lits.size()), BigInteger.valueOf(rhs));
    }

    @Override
    public void addAtLeast(IVecInt lits, int rhs) throws ContradictionException {
        addGreaterOrEqual(lits, makeUnitCoeffs(lits.size()), BigInteger.valueOf(rhs));
    }

    @Override
    public void addEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        writeEqual(hard_writer, lits, coeffs, rhs);
        ++nhard;
    }

    @Override
    public void addGreaterOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        writeGreaterOrEqual(hard_writer, lits, coeffs, rhs);
        ++nhard;
    }

    @Override
    public void addLessOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        writeLessOrEqual(hard_writer, lits, coeffs, rhs);
        ++nhard;
    }

    @Override
    public void addClause(IVecInt lits) throws ContradictionException {
        addGreaterOrEqual(lits, makeUnitCoeffs(lits.size()), BigInteger.ONE);
    }
    
    /**
     * Prepares the CLASP solver for execution and retrieves the path to the CLASP executable.
     * @return The path to the CLASP executable.
     */
    private String prepareSolver() {
        // Create temporary executable
        String clasp_path = Configuration.getInstance().getCLASPExecutablePath();
        File tmp_exe = Utils.makeTemporaryFile(clasp_path, "", true);
        tmp_exe.setExecutable(true);
        // Copy executable
        try {
            FileUtils.copyInputStreamToFile(ClassLoader.getSystemResourceAsStream(clasp_path), tmp_exe);
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to create temporary clasp executable", e);
        }
        return tmp_exe.getAbsolutePath();
    }
    
    /**
     * Dumps the WBO instance to a temporary file in the WBO format.
     * @return The path to the WBO instance file.
     */
    private String writeWBOFile() {
        File tmp_wbo = Utils.makeTemporaryFile(CLASP_INPUT_PREFIX, ".wbo", true);
        try {
            BufferedWriter tmp_wbo_writer = new BufferedWriter(new FileWriter(tmp_wbo));
            tmp_wbo_writer.write("* #variable= " + nvars + " #constraint= " + (nhard+nsoft));
            tmp_wbo_writer.write(" #soft= " + nsoft + " mincost= " + mincost + " maxcost= " + maxcost);
            tmp_wbo_writer.write(" sumcost= " + sumcost + "\n");
            tmp_wbo_writer.write("soft: " + sumcost.add(BigInteger.ONE) + " ;\n");
            tmp_wbo_writer.write(soft_writer.toString());
            tmp_wbo_writer.write(hard_writer.toString());
            tmp_wbo_writer.close();
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to write temporary .wbo file", e);
        }
        return tmp_wbo.getAbsolutePath();
    }
    
    /**
     * Builds the command to execute in order to call the CLASP solver on the WBO instance.
     * @param exe_path The path to the CLASP executable.
     * @param wbo_path The path to the WBO instance file.
     * @return A {@link String} array with the command.
     */
    private String[] buildCLASPCommand(String exe_path, String wbo_path) {
        String[] cmd = new String[4];
        cmd[0] = exe_path;
        cmd[1] = "--time-limit=" + getRemainingTime();
        cmd[2] = "--quiet=1,0,2";
        cmd[3] = wbo_path;
        return cmd;
    }
    
    /**
     * Parses a CLASP output line with a portion of the model of a solution and stores the variables
     * that are assigned truth value 1.
     * @param line The output line.
     * @param true_vars A set in which to store the variables that are assigned truth value 1.
     */
    private void parseModels(String line, Set<Integer> true_vars) {
        assert(line.charAt(0) == 'v');
        String[] tokens = line.split(" ");
        for (int i = 1; i < tokens.length; ++i) {
            assert(tokens[i].charAt(0) == 'x' ||
                   tokens[i].charAt(0) == '-' ||
                   tokens[i].charAt(0) == '0');
            if (tokens[i].charAt(0) == 'x') {
                Integer lit = new Integer(tokens[i].substring(1));
                true_vars.add(lit);
            }
        }
    }
    
    /**
     * Processes the cost of a new, possibly better, solution. If the new solution has cost less
     * than {@link #best}, then the new cost is stored and passed to the handler.
     * @param new_best The cost of the new solution.
     * @param handler The new best cost handler.
     * @see #solve(INewBestHandler)
     */
    private void processNewBest(BigInteger new_best, INewBestHandler handler) {
        if (best.equals(NO_BEST) || best.compareTo(new_best) > 0) {
            best = new_best;
            handler.handleNewBest(new_best);
        }
    }
    
    /**
     * Solves the Weighted Boolean Optimization instance.
     */
    @Override
    public void solve() {
        solve(new INewBestHandler() { // dummy new best handler
            @Override
            public void handleNewBest(BigInteger best) {}
        });
    }
    
    /**
     * Solves the Weighted Boolean Optimization instance. {@link ExternalWBOSolver} does not support
     * assumptions. Calling this method with assumptions will result in a
     * {@link NotSupportedException}. Therefore, {@link #solve()} or {@link #solve(INewBestHandler)}
     * should be used instead.
     * @param asms A set of literals that must be satisfied by the solution.
     */
    @Deprecated
    @Override
    public void solve(IVecInt asms) {
        if (asms.size() > 0) {
            throw new NotSupportedException("PBOSolver does not support assumptions.");
        }
        solve();
    }
    
    /**
     * Solves the Weighted Boolean Optimization instance, calling a given handler whenever a better
     * solution is found.
     * @param handler The callback object to be invoked whenever a better solution is found.
     */
    public void solve(INewBestHandler handler) {
        System.out.println("c Making temporary CLASP executable");
        String exe_path = prepareSolver();
        System.out.println("c Encoding instance into temporary .wbo file");
        String wbo_path = writeWBOFile();
        ProcessBuilder proc_builder = new ProcessBuilder(buildCLASPCommand(exe_path, wbo_path));
        proc_builder.redirectErrorStream(true);
        best = sumcost.add(BigInteger.ONE);
        try {
            Process clasp_proc = proc_builder.start();
            BufferedReader clasp_out =
                    new BufferedReader(new InputStreamReader(clasp_proc.getInputStream()));
            String line = null;
            while ((line = clasp_out.readLine()) != null) {
                if (line.length() == 0) {
                    continue;
                }
                if (line.charAt(0) == 'o') {
                    is_sat = found_solution = true;
                    processNewBest(new BigInteger(line.substring(2)), handler);
                }
                else if (line.charAt(0) == 'v') {
                    parseModels(line, model);
                }
                else if (line.charAt(0) == 's') {
                    if (line.startsWith(OPTIMUM)) {
                        found_optimum = true;
                    }
                    else if (line.startsWith(UNSATISFIABLE)) {
                        found_solution = found_optimum = true;
                    }
                    break;
                }
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to launch CLASP", e);
        }
    }

    @Override
    public boolean isSolved() { return found_solution; }

    @Override
    public boolean isSatisfiable() { return is_sat; }

    @Override
    public boolean isUnsatisfiable() { return !is_sat; }

    @Override
    public boolean modelValue(int var) { return model.contains(new Integer(var)); }
    
    /**
     * Checks if an optimal solution was found. Must be called after {@link #solve()} or
     * {@link #solve(INewBestHandler)}.
     * @return True if the optimal solution was found, false otherwise.
     */
    public boolean foundOptimum() { return found_optimum; }
    
    /**
     * If a solution was found on the last call to {@link #solve()} or
     * {@link #solve(INewBestHandler)}, retrieves the cost of the best solution found by the solver.
     * @return The cost of the best solution.
     * @see ConstraintSolver#isSatisfiable()
     */
    public BigInteger getBest() { return best; }
    
    /**
     * Writes the weight of a soft constraint to a given writer object.
     * @param writer The writer object where the constraint is to be written.
     * @param weight The weight of the soft constraint.
     */
    private void writeWeight(Writer writer, BigInteger weight) {
        try {
            writer.write("[" + weight + "] ");
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to write soft constraint weight", e);
        }
    }
    
    /**
     * Updates the minimum, maximum and total costs of soft constraints after adding a new soft
     * constraint.
     * @param weight The weight of the new soft constraint.
     * @see {@link #addSoftAtLeast(IVecInt, int, BigInteger)}
     * {@link #addSoftAtMost(IVecInt, int, BigInteger)}
     * {@link #addSoftGreaterOrEqual(IVecInt, IVec, BigInteger, BigInteger)}
     * {@link #addSoftLessOrEqual(IVecInt, IVec, BigInteger, BigInteger)}
     */
    private void updateCosts(BigInteger weight) {
        if (mincost.equals(NO_COST)) {
            assert(maxcost.equals(NO_COST));
            mincost = weight;
            maxcost = weight;
        }
        else if (mincost.compareTo(weight) > 0) {
            mincost = weight;
        }
        else if (maxcost.compareTo(weight) < 0) {
            maxcost = weight;
        }
        sumcost = sumcost.add(weight);
    }
    
    /**
     * Adds a soft constraint stating that the sum of the literals in {@code lits} must be less or
     * equal to {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @param weight The cost of unsatisfying the soft constraint.
     */
    public void addSoftAtMost(IVecInt lits, int rhs, BigInteger weight) {
        writeWeight(soft_writer, weight);
        writeLessOrEqual(soft_writer, lits, makeUnitCoeffs(lits.size()), BigInteger.valueOf(rhs));
        updateCosts(weight);
        ++nsoft;
    }
    
    /**
     * Adds a soft constraint stating that the sum of the literals in {@code lits} must be equal or
     * larger than {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @param weight The cost of unsatisfying the soft constraint.
     */
    public void addSoftAtLeast(IVecInt lits, int rhs, BigInteger weight) {
        writeWeight(soft_writer, weight);
        writeGreaterOrEqual(soft_writer, lits, makeUnitCoeffs(lits.size()), BigInteger.valueOf(rhs));
        updateCosts(weight);
        ++nsoft;
    }
    
    /**
     * Adds a soft Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be greater or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @param weight The cost of unsatisfying the soft constraint.
     */
    public void addSoftGreaterOrEqual(IVecInt lits,
                                      IVec<BigInteger> coeffs,
                                      BigInteger rhs,
                                      BigInteger weight) {
        writeWeight(soft_writer, weight);
        writeGreaterOrEqual(soft_writer, lits, coeffs, rhs);
        updateCosts(weight);
        ++nsoft;
    }
    
    /**
     * Adds a soft Pseudo-Boolean constraint stating that the sum of a set of coefficients in
     * {@code coeffs} times the corresponding literals in {@code lits} must be less or equal to
     * {@code rhs}.
     * @param lits The literals in the left-hand side of the constraint.
     * @param coeffs The coefficients in the left-hand side of the constraint.
     * @param rhs The value in the right-hand side of the constraint.
     * @param weight The cost of unsatisfying the soft constraint.
     */
    public void addSoftLessOrEqual(IVecInt lits,
                                   IVec<BigInteger> coeffs,
                                   BigInteger rhs,
                                   BigInteger weight) {
        writeWeight(soft_writer, weight);
        writeLessOrEqual(soft_writer, lits, coeffs, rhs);
        updateCosts(weight);
        ++nsoft;
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableExactly(IVecInt lits, int rhs) throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableAtMost(IVecInt lits, int rhs) throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableAtLeast(IVecInt lits, int rhs) throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableGreaterOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableLessOrEqual(IVecInt lits, IVec<BigInteger> coeffs, BigInteger rhs)
            throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public ConstraintID addRemovableClause(IVecInt lits) throws ContradictionException {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }

    // TODO: implement
    @Override
    public void removeConstraint(ConstraintID id) {
        throw new NotSupportedException("ExternalWBOSolver does not support removable constraints.");
    }
    
    @Override
    public IVecInt unsatExplanation() {
        throw new NotSupportedException("ExternalWBOSolver does not support unsat explanations.");
    }

}
