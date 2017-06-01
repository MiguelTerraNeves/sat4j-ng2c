package vmalloc.preprocess;

import java.math.BigInteger;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.sat4j.core.Vec;
import org.sat4j.specs.IVec;

import vmalloc.Utils;
import vmalloc.algorithm.AllocAlgorithm;
import vmalloc.domain.JobVec;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.exception.HeuristicReductionFailedException;

// TODO: allow multiple options
public class HeuristicReducer {

    private PhysicalMachineVec pms;
    private JobVec jobs;
    private MappingVec mappings;
    private double max_mig_percentile;
    private AllocAlgorithm reduction_alg;

    public HeuristicReducer(PhysicalMachineVec pms,
                            JobVec jobs,
                            MappingVec mappings,
                            double max_mig_percentile,
                            AllocAlgorithm reduction_alg) {
        assert(max_mig_percentile >= 0.0 && max_mig_percentile <= 1.0);
        this.pms = pms;
        this.jobs = jobs;
        this.mappings = mappings;
        this.max_mig_percentile = max_mig_percentile;
        this.reduction_alg = reduction_alg;
    }
    
    public PhysicalMachineVec getPhysicalMachines() { return this.pms; }
    public JobVec getJobs() { return this.jobs; }
    public MappingVec getMappings() { return this.mappings; }
    public double getMaximumMigrationPercentile() { return this.max_mig_percentile; }
    
    public void apply() throws HeuristicReductionFailedException {
        System.out.println("c WARNING: timeout ignored when performing heuristic reduction");
        Utils.stdoutDisable();
        reduction_alg.allocate();
        Utils.stdoutEnable();
        if (reduction_alg.foundSolution()) {
            assert(Utils.allocationIsValid(this.pms,
                                           this.jobs,
                                           this.mappings,
                                           this.max_mig_percentile,
                                           reduction_alg.getSolutions().get(0)));
            BigInteger total_mem_cap = BigInteger.ZERO;
            Map<Integer, Integer> pm_id_to_idx = new HashMap<Integer, Integer>();
            for (int i = 0; i < this.pms.size(); ++i) {
                pm_id_to_idx.put(new Integer(this.pms.get(i).getID()), new Integer(i));
                total_mem_cap = total_mem_cap.add(this.pms.get(i).getMemory());
            }
            MappingVec allocation = reduction_alg.getSolutions().get(0);
            IVec<Set<VirtualMachine>> placement = new Vec<Set<VirtualMachine>>(this.pms.size());
            for (int i = 0; i < this.pms.size(); ++i) {
                placement.push(new HashSet<VirtualMachine>());
            }
            for (int i = 0; i < allocation.size(); ++i) {
                PhysicalMachine pm = allocation.get(i).getPhysicalMachine();
                VirtualMachine vm = allocation.get(i).getVirtualMachine();
                int pm_idx = pm_id_to_idx.get(new Integer(pm.getID())).intValue();
                placement.get(pm_idx).add(vm);
            }
            double new_mig_budget = this.max_mig_percentile;
            MappingVec new_maps = new MappingVec();
            for (int i = 0; i < this.mappings.size(); ++i) {
                PhysicalMachine pm = this.mappings.get(i).getPhysicalMachine();
                VirtualMachine vm = this.mappings.get(i).getVirtualMachine();
                int pm_idx = pm_id_to_idx.get(new Integer(pm.getID())).intValue();
                if (placement.get(pm_idx).contains(vm)) {
                    new_maps.push(this.mappings.get(i));
                }
                else {
                    new_mig_budget -= Utils.toPercentile(vm.getMemory(), total_mem_cap);
                }
            }
            assert(new_mig_budget >= 0.0);
            BigInteger new_total_mem_cap = BigInteger.ZERO;
            PhysicalMachineVec new_pms = new PhysicalMachineVec();
            for (int i = 0; i < this.pms.size(); ++i) {
                if (placement.get(i).size() > 0) {
                    new_pms.push(this.pms.get(i));
                    new_total_mem_cap = new_total_mem_cap.add(this.pms.get(i).getMemory());
                }
            }
            new_mig_budget = Utils.scalePercentile(total_mem_cap, new_mig_budget, new_total_mem_cap);
            new_mig_budget = Utils.normalizePercentile(new_mig_budget);
            this.pms = new_pms;
            this.mappings = new_maps;
            this.max_mig_percentile = new_mig_budget;
        }
        else {
            throw new HeuristicReductionFailedException();
        }
    }
    
}
