package vmalloc;

import static org.sat4j.GlobalDefs.USE_NG2C;
import static org.sat4j.GlobalDefs.ANNOTATE_SOLVER_STRUCTS;
import static org.sat4j.GlobalDefs.ANNOTATE_INSTANCE;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.moeaframework.core.FrameworkException;
import org.sat4j.specs.IVec;

import vmalloc.algorithm.AllocAlgorithm;
import vmalloc.algorithm.BBOAlloc;
import vmalloc.algorithm.BestFitDecreasingAlloc;
import vmalloc.algorithm.ConstraintBasedAllocAlgorithm;
import vmalloc.algorithm.ConstraintBasedAllocAlgorithm.HashFunctionType;
import vmalloc.algorithm.DBEAAlloc;
import vmalloc.algorithm.DEAlloc;
import vmalloc.algorithm.EvolutionaryAllocAlgorithm;
import vmalloc.algorithm.EvolutionaryAllocAlgorithm.InitializationType;
import vmalloc.algorithm.FirstFitDecreasingAlloc;
import vmalloc.algorithm.GAAlloc;
import vmalloc.algorithm.GGAAlloc;
import vmalloc.algorithm.GIAAlloc;
import vmalloc.algorithm.HashEnumAlloc;
import vmalloc.algorithm.LinearSearchAlloc;
import vmalloc.algorithm.MCSAlloc;
import vmalloc.algorithm.PBOAlloc;
import vmalloc.algorithm.PSOAlloc;
import vmalloc.algorithm.ParetoCLD;
import vmalloc.algorithm.ParetoLBX;
import vmalloc.algorithm.WBOAlloc;
import vmalloc.domain.Job;
import vmalloc.domain.JobVec;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.evolutionary.VMCwMProblem;
import vmalloc.exception.HeuristicReductionFailedException;
import vmalloc.preprocess.HeuristicReducer;

public class VMAlloc {

    private static final String LINEAR_SEARCH = "LS";
    private static final String MCS = "MCS";
    private static final String WBO = "WBO";
    private static final String PBO = "PBO";
    private static final String FFD = "FFD";
    private static final String BFD = "BFD";
    private static final String DE = "DE";
    private static final String PSO = "PSO";
    private static final String GA = "GA";
    private static final String DBEA = "DBEA";
    private static final String BBO = "BBO";
    private static final String GGA = "GGA";
    private static final String GIA = "GIA";
    private static final String HASH = "HE";
    private static final String PARETO_CLD = "PCLD";
    private static final String PARETO_LBX = "PLBX";
    
    private static final String RAND_INIT = "RAND";
    private static final String RAND_PACKING_INIT = "RBP";
    private static final String SHUFFLED_FIRST_FIT_INIT = "SFF";
    private static final String SHUFFLED_VMCWM_INIT = "SVMCWM";
    private static final String MIXED_INIT = "MIXED";
    
    private static final String NO_HASH = "NONE";
    private static final String GLOBAL_HASH = "GBL";
    private static final String SEPARATED_HASH = "SEP";
    
    private static final String DEFAULT_POPSIZE = "100";
    private static final String DEFAULT_DE_CROSSOVER_RATE = "0.1";
    private static final String DEFAULT_DE_STEP_SIZE = "0.5";
    private static final String DEFAULT_PSO_ARCHIVE_SIZE = "100";
    private static final String DEFAULT_PSO_MUTATION_RATE = "0.1";
    private static final String DEFAULT_PSO_DISTRIBUTION_INDEX = "20.0";
    private static final String DEFAULT_BBO_POPSIZE = "3";
    private static final String DEFAULT_BBO_IMMIGRATION_RATE = "0.5";
    private static final String DEFAULT_BBO_MUTATION_RATE = "0.05";
    private static final String DEFAULT_BBO_CROSS_MIGRATION_RATE = "0.5";
    private static final String DEFAULT_GA_MUTATION_RATE = "0.05";
    private static final String DEFAULT_GA_CROSSOVER_RATE = "0.5";
    private static final String DEFAULT_GGA_MUTATION_RATE = "0.0";
    private static final String DEFAULT_GGA_CROSSOVER_RATE = "0.8";
    
    private static final String DEFAULT_MIG_PERCENTILE = "1.0";
    
    // TODO: print proper help message
    private static void printHelpMessage(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("vmalloc", options);
    }
    
    private static void discardPlatformRestrictions(JobVec jobs) {
        for (int i = 0; i < jobs.size(); ++i) {
            Job job = jobs.get(i);
            for (int j = 0; j < job.nVirtualMachines(); ++j) {
                job.getVirtualMachine(j).clearUnallowedPhysicalMachines();
            }
        }
    }
    
    private static void discardAntiColocationConstraints(JobVec jobs) {
        for (int i = 0; i < jobs.size(); ++i) {
            Job job = jobs.get(i);
            for (int j = 0; j < job.nVirtualMachines(); ++j) {
                job.getVirtualMachine(j).setFullyColocatable();
            }
        }
    }
    
    private static void printMapping(MappingVec maps) {
        for (int i = 0; i < maps.size(); ++i) {
            PhysicalMachine pm = maps.get(i).getPhysicalMachine();
            VirtualMachine vm = maps.get(i).getVirtualMachine();
            System.out.println("p " + vm.getJobID() + "-" + vm.getIndex() + " -> " + pm.getID());
        }
    }
    
    // FIXME: HACK!! requires an instance of an allocation algorithm
    private static void analyze(AllocAlgorithm alloc, String dataset_desc) throws IOException {
        String[] label_fpath_pairs = dataset_desc.split(";");
        Map<String, String> dataset = new HashMap<String, String>();
        for (int i = 0; i < label_fpath_pairs.length; ++i) {
            String[] label_fpath_pair = label_fpath_pairs[i].split(":");
            String label = label_fpath_pair[0], path = label_fpath_pair[1];
            dataset.put(label, path);
        }
        try {
            alloc.analyzePopulations(dataset);
        }
        catch(FrameworkException fe) {
            throw new RuntimeException("Not possible to analyze populations for " + dataset_desc, fe);
        }
    }
    
    public static void main(String[] args) {
        if (USE_NG2C) {
            System.newAllocGen();
        }
        Clock.getInstance().reset();
        System.out.println("c Parsing");
        Options options = ANNOTATE_SOLVER_STRUCTS ? new @Gen Options() : new Options();
        options.addOption("a", "algorithm", true,
                          "Choose the allocation algorithm. Options are Linear Search (" +
                          LINEAR_SEARCH + "), Minimum Correction Set (" + MCS + "), " +
                          "Weighted Boolean Optimization (" + WBO + "), Pseudo-Boolean " +
                          "Optimization (" + PBO + "), First-Fit Decreasing (" + FFD +
                          "), Best-Fit Decreasing (" + BFD + "), Differential Evolution (" +
                          DE + "), Particle Swarm Optimization (" + PSO + "), Genetic Algorithm (" +
                          GA + "), Biogeography-Based Optimization (" + BBO + "), Grouping Genetic " +
                          "Algorithm (" + GGA + "), Guided Improvement Algorithm (" + GIA + "), " +
                          "Hash-based Enumeration (" + HASH + "), Pareto CLD (" + PARETO_CLD + ") and " +
                          "Pareto LBX (" + PARETO_LBX + "). Default is " + WBO + ".");
        options.addOption("t", "time-limit", true,
                          "Set the time limit in seconds. No limit by default.");
        options.addOption("m", "migration-percentile", true,
                          "Set the fraction of total memory capacity that can be used up in " +
                          "migrations. Default is " + DEFAULT_MIG_PERCENTILE + ".");
        options.addOption("r", "heuristic-reduction", false, "Enable heuristic reduction.");
        options.addOption("ra", "reduction-algorithm", true,
                          "Choose the allocation algorithm used for heuristic reduction. Options " +
                          "are " + FFD + " and " + BFD + ". Default is " + BFD + ".");
        options.addOption("s", "break-symmetries", false, "Enable symmetry breaking.");
        options.addOption("ip", "ignore-platform", false,
                          "Ignore platform specific restrictions.");
        options.addOption("ht", "hash-type", true, "Sets hash function type to use for Hash-based " +
                          "Enumeration. Options are None (" + NO_HASH + "), Global (" + GLOBAL_HASH +
                          ") or Separated (" + SEPARATED_HASH + "). Default is " + GLOBAL_HASH + ".");
        options.addOption("ic", "ignore-colocation", false, "Ignore anti-colocation constraints.");
        options.addOption("pa", "print-allocations", false, "Enable solution printing.");
        options.addOption("it", "initialization-type", true,
                          "Attempts to force initialization of population in evolutionary " +
                          "algorithms with feasible solutions using the given approach. Options " +
                          "are Random (" + RAND_INIT + "), Random Bin-Packing (" + RAND_PACKING_INIT +
                          "), Shuffled First-Fit (" + SHUFFLED_FIRST_FIT_INIT + "), Shuffled " +
                          "VMCwM Heuristic (" + SHUFFLED_VMCWM_INIT + ") and Mixed Initialization " +
                          "(" + MIXED_INIT + "). Default is " + RAND_INIT + ".");
        options.addOption("dp", "dump-population", true,
                          "Dump the evolutionary algorithm's final population to the given file. " +
                          "This option exists mostly for evaluation purposes.");
        options.addOption("ap", "analyze-populations", true,
                          "Enables analysis of sets of populations stored in files generated using " +
                          "the 'dp' option instead of computing an allocation. The population " +
                          "files are specified by providing a list of <label>:<file_path> pairs, " +
                          "separated by ';'. The last entry is assumed to be just of the form " +
                          "<file_path> and the path to the file with the reference set for " +
                          "hypervolume computation. This option exists mostly for evaluation " +
                          "purposes.");
        options.addOption("ms", "multiple-seeds", true,
                          "Runs the selected algorithm multiple times with the given number of " +
                          "different seeds. The timeout is reset between runs. This option " +
                          "exists mostly for evaluation purposes.");
        options.addOption("ps", "population-size", true,
                          "Set the population size for evolutionary algorithms. Default is " +
                          DEFAULT_POPSIZE + ". For " + BBO + ", this sets the population size per " +
                          "subsystem. In this case, default is " + DEFAULT_BBO_POPSIZE + ".");
        options.addOption("cr", "crossover-rate", true,
                          "Set the crossover rate for the " + DE + ", " + GA + " and " + GGA +
                          " algorithms. Default for " + DE + " is " + DEFAULT_DE_CROSSOVER_RATE +
                          ", for " + GA + " is " + DEFAULT_GA_CROSSOVER_RATE + " and for " + GGA +
                          " is " + DEFAULT_GGA_CROSSOVER_RATE + ".");
        options.addOption("ss", "step-size", true,
                          "Set the step size for the " + DE + " algorithm. Default is " +
                          DEFAULT_DE_STEP_SIZE + ".");
        options.addOption("as", "archive-size", true,
                          "Set the archive size for the " + PSO + " algorithm. Default is " +
                          DEFAULT_PSO_ARCHIVE_SIZE + ".");
        options.addOption("mr", "mutation-rate", true,
                          "Set the mutation rate for the " + PSO + ", " + BBO + ", " + GA + " and " +
                          GGA + " algorithms. " + PSO + " uses polinomial mutation and the default " +
                          "value is " + DEFAULT_PSO_MUTATION_RATE + ". " + BBO + " and " + GA +
                          " use uniform single variable mutation and the default values are " +
                          DEFAULT_BBO_MUTATION_RATE + " and " + DEFAULT_GA_MUTATION_RATE +
                          " respectively. " + GGA + " uses group mutation and the default value " +
                          "is " + DEFAULT_GGA_MUTATION_RATE);
        options.addOption("di", "distribution-index", true,
                          "Set the distribution index for polynomial mutation in the " + PSO +
                          " algorithm. Default is " + DEFAULT_PSO_DISTRIBUTION_INDEX + ".");
        options.addOption("ir", "immigration-rate", true,
                          "Set the immigration rate for the " + BBO + " algorithm. Default is " +
                          DEFAULT_BBO_IMMIGRATION_RATE + ".");
        options.addOption("cmr", "cross-migration-rate", true,
                          "Set the cross migration rate for the " + BBO + " algorithm. Default is " +
                          DEFAULT_BBO_CROSS_MIGRATION_RATE + ".");
        CommandLineParser cl_parser = ANNOTATE_SOLVER_STRUCTS ? new @Gen DefaultParser() : new DefaultParser();
        try {
            CommandLine cl = cl_parser.parse(options, args);
            InputParser in_parser = ANNOTATE_SOLVER_STRUCTS ? new @Gen InputParser(cl.getArgs()[0]) : new InputParser(cl.getArgs()[0]);
            in_parser.parse();
            System.out.println("c Parsing time: " + Clock.getInstance().getElapsed() + " seconds");
            PhysicalMachineVec orig_pms = in_parser.getPhysicalMachines();
            JobVec orig_jobs = in_parser.getJobs();
            MappingVec orig_mappings = in_parser.getMappings();
            double orig_max_mig_percentile =
                    Double.parseDouble(cl.getOptionValue("m", DEFAULT_MIG_PERCENTILE));
            PhysicalMachineVec pms = orig_pms;
            JobVec jobs = orig_jobs;
            MappingVec mappings = orig_mappings;
            double max_mig_percentile = orig_max_mig_percentile;
            ProblemStatistics stats = ANNOTATE_INSTANCE ? new @Gen ProblemStatistics(pms, jobs, mappings, max_mig_percentile) : new ProblemStatistics(pms, jobs, mappings, max_mig_percentile);
            stats.printStatistics();
            boolean break_symms = false;
            if (cl.hasOption("s")) {
                System.out.println("c Symmetry breaking enabled");
                break_symms = true;
            }
            if (cl.hasOption("ip")) {
                System.out.println("c Discarding platform constraints");
                discardPlatformRestrictions(jobs);
            }
            if (cl.hasOption("ic")) {
                System.out.println("c Discarding anti-colocation constraints");
                discardAntiColocationConstraints(jobs);
            }
            VMCwMProblem instance = ANNOTATE_INSTANCE ? new @Gen VMCwMProblem(pms, jobs, mappings, max_mig_percentile) : new VMCwMProblem(pms, jobs, mappings, max_mig_percentile);
            if (cl.hasOption("r")) {
                AllocAlgorithm reduction_alg = null;
                if (!cl.hasOption("ra") || cl.getOptionValue("ra").equals(BFD)) {
                    reduction_alg = ANNOTATE_SOLVER_STRUCTS ? new @Gen BestFitDecreasingAlloc(instance) : new BestFitDecreasingAlloc(instance);
                }
                else if (cl.getOptionValue("ra").equals(FFD)) {
                    reduction_alg = ANNOTATE_SOLVER_STRUCTS ? new @Gen FirstFitDecreasingAlloc(instance) : new FirstFitDecreasingAlloc(instance);
                }
                else {
                    printHelpMessage(options);
                    return;
                }
                System.out.println("c Applying heuristic reduction");
                HeuristicReducer reducer =
                        ANNOTATE_SOLVER_STRUCTS ? new @Gen HeuristicReducer(pms, jobs, mappings, max_mig_percentile, reduction_alg) : new HeuristicReducer(pms, jobs, mappings, max_mig_percentile, reduction_alg);
                try {
                    reducer.apply();
                    System.out.println("c Solution using " + reducer.getPhysicalMachines().size() +
                                       " PMs found");
                    if (reducer.getPhysicalMachines().size() < pms.size()) {
                        pms = reducer.getPhysicalMachines();
                        mappings = reducer.getMappings();
                        max_mig_percentile = reducer.getMaximumMigrationPercentile();
                    }
                    System.out.println("c Elapsed time: " + Clock.getInstance().getElapsed() +
                                       " seconds");
                    stats = ANNOTATE_INSTANCE ? new @Gen ProblemStatistics(pms, jobs, mappings, max_mig_percentile) : new ProblemStatistics(pms, jobs, mappings, max_mig_percentile);
                    stats.printStatistics();
                }
                catch (HeuristicReductionFailedException e) {
                    System.out.println("c Heuristic reduction failed");
                }
            }
            AllocAlgorithm alloc = null;
            if (!cl.hasOption("a") || cl.getOptionValue("a").equals(WBO)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen WBOAlloc(instance, break_symms) : new WBOAlloc(instance, break_symms);
            }
            else if (cl.getOptionValue("a").equals(LINEAR_SEARCH)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen LinearSearchAlloc(instance, break_symms) : new LinearSearchAlloc(instance, break_symms);
            }
            else if (cl.getOptionValue("a").equals(PBO)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen PBOAlloc(instance, break_symms) : new PBOAlloc(instance, break_symms);
            }
            else if (cl.getOptionValue("a").equals(PARETO_CLD)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen ParetoCLD(instance, break_symms) : new ParetoCLD(instance, break_symms);
            }
            else if (cl.getOptionValue("a").equals(PARETO_LBX)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen ParetoLBX(instance, break_symms) : new ParetoLBX(instance, break_symms);
            }
            else if (cl.getOptionValue("a").equals(FFD)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen FirstFitDecreasingAlloc(instance) : new FirstFitDecreasingAlloc(instance);
            }
            else if (cl.getOptionValue("a").equals(BFD)) {
                alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen BestFitDecreasingAlloc(instance) : new BestFitDecreasingAlloc(instance);
            }
            else if (cl.getOptionValue("a").equals(HASH) ||
                     cl.getOptionValue("a").equals(GIA) ||
                     cl.getOptionValue("a").equals(MCS)) {
                ConstraintBasedAllocAlgorithm cb_alloc = null;
                if (cl.getOptionValue("a").equals(HASH)) {
                    cb_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen HashEnumAlloc(instance, break_symms) : new HashEnumAlloc(instance, break_symms);
                    System.out.println("c ========= HE Configuration =========");
                }
                else if (cl.getOptionValue("a").equals(GIA)) {
                    cb_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen GIAAlloc(instance, break_symms) : new GIAAlloc(instance, break_symms);
                    System.out.println("c ======== GIA Configuration =========");
                }
                else {
                    cb_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen MCSAlloc(instance, break_symms) : new MCSAlloc(instance, break_symms);
                    System.out.println("c ======== MCS Configuration =========");
                }
                String hash_type = cl.getOptionValue("ht", GLOBAL_HASH);
                if (hash_type.equals(GLOBAL_HASH)) {
                    cb_alloc.setHashFunctionType(HashFunctionType.GLOBAL);
                }
                else if (hash_type.equals(SEPARATED_HASH)) {
                    cb_alloc.setHashFunctionType(HashFunctionType.SEPARATED);
                }
                else if (hash_type.equals(NO_HASH)) {
                    cb_alloc.setHashFunctionType(HashFunctionType.NONE);
                }
                else {
                    printHelpMessage(options);
                    return;
                }
                System.out.println("c  Hash Function Type:    " + hash_type);
                System.out.println("c ====================================");
                alloc = cb_alloc;
            }
            else {
                // Evolutionary approaches
                EvolutionaryAllocAlgorithm ea_alloc = null;
                if (cl.getOptionValue("a").equals(DE)) {
                    DEAlloc de_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen DEAlloc(instance) : new DEAlloc(instance);
                    double cr = Double.parseDouble(cl.getOptionValue("cr", DEFAULT_DE_CROSSOVER_RATE));
                    double ss = Double.parseDouble(cl.getOptionValue("ss", DEFAULT_DE_STEP_SIZE));
                    de_alloc.setCrossoverRate(cr);
                    de_alloc.setStepSize(ss);
                    ea_alloc = de_alloc;
                    System.out.println("c ========= DE Configuration =========");
                    System.out.println("c  Crossover rate:       " + cr);
                    System.out.println("c  Step size:            " + ss);
                }
                else if (cl.getOptionValue("a").equals(PSO)) {
                    PSOAlloc pso_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen PSOAlloc(instance) : new PSOAlloc(instance);
                    System.out.println("c ========= PSO Configuration ========");
                    double mr = Double.parseDouble(cl.getOptionValue("mr", DEFAULT_PSO_MUTATION_RATE));
                    double di = Double.parseDouble(cl.getOptionValue("di",
                                                                     DEFAULT_PSO_DISTRIBUTION_INDEX));
                    int as = Integer.parseInt(cl.getOptionValue("as", DEFAULT_PSO_ARCHIVE_SIZE));
                    pso_alloc.setMutationRate(mr);
                    pso_alloc.setDistributionIndex(di);
                    pso_alloc.setArchiveSize(as);
                    ea_alloc = pso_alloc;
                    System.out.println("c  Mutation rate:        " + mr);
                    System.out.println("c  Distribution index:   " + di);
                    System.out.println("c  Archive size:         " + as);
                }
                else if (cl.getOptionValue("a").equals(BBO)) {
                    BBOAlloc bbo_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen BBOAlloc(instance) : new BBOAlloc(instance);
                    System.out.println("c ========= BBO Configuration ========");
                    double mr = Double.parseDouble(cl.getOptionValue("mr", DEFAULT_BBO_MUTATION_RATE));
                    double ir = Double.parseDouble(
                            cl.getOptionValue("ir", DEFAULT_BBO_IMMIGRATION_RATE));
                    double cmr = Double.parseDouble(
                            cl.getOptionValue("cmr", DEFAULT_BBO_CROSS_MIGRATION_RATE));
                    bbo_alloc.setMutationRate(mr);
                    bbo_alloc.setImmigrationRate(ir);
                    bbo_alloc.setCrossSystemMigrationRate(cmr);
                    ea_alloc = bbo_alloc;
                    System.out.println("c  Immigration rate:     " + ir);
                    System.out.println("c  Mutation rate:        " + mr);
                    System.out.println("c  Cross migration rate: " + cmr);
                    // FIXME: needs to be corrected if BBO implementation changes
                    System.out.println("c  Number of subsystems: " +
                                       ((mappings.size() > 0) ? "6" : "4"));
                }
                else if (cl.getOptionValue("a").equals(GA)) {
                    GAAlloc ga_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen GAAlloc(instance) : new GAAlloc(instance);
                    double mr = Double.parseDouble(cl.getOptionValue("mr", DEFAULT_GA_MUTATION_RATE));
                    double cr = Double.parseDouble(cl.getOptionValue("cr", DEFAULT_GA_CROSSOVER_RATE));
                    ga_alloc.setMutationRate(mr);
                    ga_alloc.setCrossoverRate(cr);
                    ea_alloc = ga_alloc;
                    System.out.println("c ========= GA Configuration =========");
                    System.out.println("c  Crossover rate:       " + cr);
                    System.out.println("c  Mutation rate:        " + mr);
                }
                else if (cl.getOptionValue("a").equals(GGA)) {
                    GGAAlloc gga_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen GGAAlloc(instance) : new GGAAlloc(instance);
                    double mr = Double.parseDouble(cl.getOptionValue("mr", DEFAULT_GGA_MUTATION_RATE));
                    double cr = Double.parseDouble(cl.getOptionValue("cr", DEFAULT_GGA_CROSSOVER_RATE));
                    gga_alloc.setMutationRate(mr);
                    gga_alloc.setCrossoverRate(cr);
                    ea_alloc = gga_alloc;
                    System.out.println("c ========= GGA Configuration ========");
                    System.out.println("c  Crossover rate:       " + cr);
                    System.out.println("c  Mutation rate:        " + mr);
                }
                else if (cl.getOptionValue("a").equals(DBEA)) {
                    ea_alloc = ANNOTATE_SOLVER_STRUCTS ? new @Gen DBEAAlloc(instance) : new DBEAAlloc(instance);
                    System.out.println("c ======== DBEA Configuration ========");
                }
                else {
                    printHelpMessage(options);
                    return;
                }
                String default_pop_size = null;
                if (cl.getOptionValue("a").equals(BBO)) {
                    default_pop_size = DEFAULT_BBO_POPSIZE;
                }
                else {
                    default_pop_size = DEFAULT_POPSIZE;
                }
                int ps = Integer.parseInt(cl.getOptionValue("ps", default_pop_size));
                if (!cl.getOptionValue("a").equals(DBEA)) {
                    ea_alloc.setPopulationSize(ps);
                    System.out.println("c  Population size:      " + ps);
                }
                System.out.print("c  Initialization type:  ");
                InitializationType init_type;
                if (!cl.hasOption("it") || cl.getOptionValue("it").equals(RAND_INIT)) {
                    init_type = InitializationType.RANDOM;
                    System.out.print(RAND_INIT + Utils.NEW_LINE);
                }
                else if (cl.getOptionValue("it").equals(RAND_PACKING_INIT)) {
                    init_type = InitializationType.RANDOM_PACKING;
                    System.out.print(RAND_PACKING_INIT + Utils.NEW_LINE);
                }
                else if (cl.getOptionValue("it").equals(SHUFFLED_FIRST_FIT_INIT)) {
                    init_type = InitializationType.SHUFFLED_FIRST_FIT;
                    System.out.print(SHUFFLED_FIRST_FIT_INIT + Utils.NEW_LINE);
                }
                else if (cl.getOptionValue("it").equals(SHUFFLED_VMCWM_INIT)) {
                    init_type = InitializationType.SHUFFLED_VMCWM_HEURISTIC;
                    System.out.print(SHUFFLED_VMCWM_INIT + Utils.NEW_LINE);
                }
                else if (cl.getOptionValue("it").equals(MIXED_INIT)) {
                    init_type = InitializationType.MIXED;
                    System.out.print(MIXED_INIT + Utils.NEW_LINE);
                }
                else {
                    printHelpMessage(options);
                    return;
                }
                System.out.println("c ====================================");
                ea_alloc.setInitializationType(init_type);
                alloc = ea_alloc;
            }
            if (cl.hasOption("t")) {
                alloc.setTimeout(Integer.parseInt(cl.getOptionValue("t")));
            }
            if (cl.hasOption("ap")) {
                analyze(alloc, cl.getOptionValue("ap"));
                return;
            }
            System.setAllocGen(0);
            if (cl.hasOption("ms")) {
                alloc.allocateMultipleSeeds(Integer.parseInt(cl.getOptionValue("ms")));
            }
            else {
                alloc.allocate();
            }
            if (alloc.foundSolution()) {
                System.out.println("s SUCCESS");
                double[] energy_costs = alloc.getSolutionEnergyCosts();
                double[] wastage_costs = alloc.getSolutionWastageCosts();
                double[] migration_costs = null;
                if (mappings.size() > 0) {
                    migration_costs = alloc.getSolutionMigrationCosts();
                    assert(energy_costs.length == migration_costs.length);
                }
                assert(energy_costs.length == wastage_costs.length);
                for (int i = 0; i < energy_costs.length; ++i) {
                    if (mappings.size() > 0) {
                        System.out.printf("e %.5f \tw %.5f \tm %d\n",
                                          energy_costs[i], wastage_costs[i], (long)migration_costs[i]);
                    }
                    else {
                        System.out.printf("e %.5f \tw %.5f\n",
                                          energy_costs[i], wastage_costs[i]);
                    }
                }
                if (cl.hasOption("pa")) {
                    IVec<MappingVec> solutions = alloc.getSolutions();
                    for (int i = 0; i < solutions.size(); ++i) {
                        assert(Utils.allocationIsValid(orig_pms,
                                orig_jobs,
                                orig_mappings,
                                orig_max_mig_percentile,
                                solutions.get(i)));
                        System.out.println("s SOLUTION " + i);
                        printMapping(solutions.get(i));
                    }
                }
            }
            else {
                System.out.println("s FAILURE");
            }
            // FIXME: maybe should check if solution was found first?
            if (cl.hasOption("dp")) {
                alloc.dumpPopulation(cl.getOptionValue("dp"));
            }
        }
        catch (ParseException e) {
            printHelpMessage(options);
        }
        catch (NumberFormatException e) {
            e.printStackTrace();
            printHelpMessage(options);
        }
        catch (IOException e) {
            System.out.println("c PARSING ERROR!");
        }
    }

}
