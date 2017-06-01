package vmalloc.algorithm;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import org.moeaframework.core.PRNG;
import org.sat4j.core.Vec;
import org.sat4j.specs.IVec;

import vmalloc.Utils;
import vmalloc.domain.Machine;
import vmalloc.domain.Mapping;
import vmalloc.domain.MappingVec;
import vmalloc.domain.MemoryAggregator;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;
import vmalloc.evolutionary.VMCwMProblem;

public abstract class BinPackingAllocAlgorithm extends AllocAlgorithm {
    
    // =============================================================================================
    // ======================================== PRIVATE API ========================================
    // =============================================================================================
    
    private MappingVec mappings = null;
    private double max_mig_percentile;
    private MappingVec partial_allocation = null;
    private VirtualMachineVec leftover_vms = null;
    private boolean shuffle_vms = false;
    private PriorityQueue<UsageInfo> pm_heap = null;
    private VirtualMachine[] vm_array = null;
    
    private void placeMappings(UsageInfo[] infos, Map<Integer,Integer> id2idx, MappingVec mappings) {
        for (int i = 0; i < mappings.size(); ++i) {
            PhysicalMachine pm = mappings.get(i).getPhysicalMachine();
            int pm_idx = id2idx.get(new Integer(pm.getID())).intValue();
            infos[pm_idx].placeVirtualMachine(mappings.get(i).getVirtualMachine());
        }
    }
    
    private void addAllUsageInfos(PriorityQueue<UsageInfo> heap, UsageInfo[] infos) {
        for (int i = 0; i < infos.length; ++i) {
            heap.add(infos[i]);
        }
    }

    private PriorityQueue<UsageInfo> buildPhysicalMachineHeap(PhysicalMachineVec pms,
                                                              MappingVec mappings) {
        PriorityQueue<UsageInfo> heap =
                new PriorityQueue<UsageInfo>(pms.size(), makeUsageInfoComparator());
        UsageInfo[] use_infos = new UsageInfo[pms.size()];
        Map<Integer,Integer> id2idx = Utils.makePhysicalMachineIDtoIndexMap(pms);
        for (int i = 0; i < pms.size(); ++i) {
            PhysicalMachine pm = pms.get(i);
            int pm_idx = id2idx.get(new Integer(pm.getID())).intValue();
            use_infos[pm_idx] = new UsageInfo(pms.get(i));
        }
        placeMappings(use_infos, id2idx, mappings);
        addAllUsageInfos(heap, use_infos);
        return heap;
    }
    
    private IVec<UsageInfo> popAllUsageInfos(PriorityQueue<UsageInfo> heap) {
        IVec<UsageInfo> infos = new Vec<UsageInfo>();
        while (!pm_heap.isEmpty()) {
            infos.push(heap.poll());
        }
        return infos;
    }
    
    private void clearHeapUsages(PriorityQueue<UsageInfo> heap) {
        IVec<UsageInfo> infos = popAllUsageInfos(heap);
        for (int i = 0; i < infos.size(); ++i) {
            UsageInfo info = infos.get(i);
            info.clear();
            heap.add(info);
        }
    }
    
    private void placeMappings(PriorityQueue<UsageInfo> heap,
                               PhysicalMachineVec pms,
                               MappingVec mappings) {
        IVec<UsageInfo> infos = popAllUsageInfos(heap);
        UsageInfo[] infos_array = new UsageInfo[pms.size()];
        Map<Integer,Integer> id2idx = Utils.makePhysicalMachineIDtoIndexMap(pms);
        for (int i = 0; i < infos.size(); ++i) {
            PhysicalMachine pm = infos.get(i).getPhysicalMachine();
            int pm_idx = id2idx.get(new Integer(pm.getID())).intValue();
            infos_array[pm_idx] = infos.get(i);
        }
        placeMappings(infos_array, id2idx, mappings);
        addAllUsageInfos(heap, infos_array);
    }
    
    private PriorityQueue<UsageInfo> duplicatePhysicalMachineHeap(PriorityQueue<UsageInfo> heap) {
        PriorityQueue<UsageInfo> copy =
                new PriorityQueue<UsageInfo>(heap.size(), makeUsageInfoComparator());
        for (Iterator<UsageInfo> it = heap.iterator(); it.hasNext();) {
            copy.add(new UsageInfo(it.next()));
        }
        return copy;
    }
    
    private Set<VirtualMachine> getMappedVirtualMachines(MappingVec mappings) {
        Set<VirtualMachine> mapped_vms = new HashSet<VirtualMachine>();
        for (int i = 0; i < mappings.size(); ++i) {
            mapped_vms.add(mappings.get(i).getVirtualMachine());
        }
        return mapped_vms;
    }
    
    private PhysicalMachine placeVirtualMachine(PriorityQueue<UsageInfo> pm_heap, VirtualMachine vm) {
        PhysicalMachine chosen_pm = null;
        IVec<UsageInfo> popped_infos = new Vec<UsageInfo>();
        while (chosen_pm == null && !pm_heap.isEmpty()) {
            UsageInfo info = pm_heap.poll();
            if (info.canHostVirtualMachine(vm)) {
                info.placeVirtualMachine(vm);
                chosen_pm = info.getPhysicalMachine();
            }
            popped_infos.push(info);
        }
        for (int i = 0; i < popped_infos.size(); ++i) {
            pm_heap.add(popped_infos.get(i));
        }
        return chosen_pm;
    }
    
    private PhysicalMachineVec getUsedPhysicalMachines(PriorityQueue<UsageInfo> pm_heap) {
        PhysicalMachineVec used_pms = new PhysicalMachineVec();
        for (Iterator<UsageInfo> it = pm_heap.iterator(); it.hasNext();) {
            UsageInfo info = it.next();
            if (!info.isEmpty()) {
                used_pms.push(info.getPhysicalMachine());
            }
        }
        return used_pms;
    }
    
    private UsageInfo[] getSortedInfos(PriorityQueue<UsageInfo> pm_heap) {
        UsageInfo[] infos = new UsageInfo[pm_heap.size()];
        pm_heap.toArray(infos);
        Arrays.sort(infos, makeUsageInfoComparator());
        return infos;
    }
    
    private MappingVec retrieveMappings(PriorityQueue<UsageInfo> pm_heap) {
        MappingVec mappings = new MappingVec();
        for (Iterator<UsageInfo> it = pm_heap.iterator(); it.hasNext();) {
            UsageInfo info = it.next();
            PhysicalMachine pm = info.getPhysicalMachine();
            for (Iterator<VirtualMachine> vm_it = info.iterator(); vm_it.hasNext();) {
                VirtualMachine vm = vm_it.next();
                mappings.push(new Mapping(vm, pm));
            }
        }
        return mappings;
    }
    
    private BigInteger notMigratedMemoryRequirementSum(VirtualMachine[] vms,
                                                         Set<VirtualMachine> to_migrate) {
        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < vms.length; ++i) {
            if (to_migrate.contains(vms[i])) {
                sum = sum.add(vms[i].getMemory());
            }
        }
        return sum;
    }
    
    private Machine[] physicalMachineArrayToMachineArray(PhysicalMachine[] pms) {
        Machine[] ms = new Machine[pms.length];
        for (int i = 0; i < pms.length; ++i) {
            ms[i] = pms[i];
        }
        return ms;
    }
    
    private Machine[] virtualMachineArrayToMachineArray(VirtualMachine[] vms) {
        Machine[] ms = new Machine[vms.length];
        for (int i = 0; i < vms.length; ++i) {
            ms[i] = vms[i];
        }
        return ms;
    }
    
    private Comparator<Machine> makeDecreasingMachineComparator() {
        return new Comparator<Machine>() {
            @Override
            public int compare(Machine m0, Machine m1) {
                BigInteger diff = m1.getCPU().subtract(m0.getCPU());
                if (diff.equals(BigInteger.ZERO)) {
                    diff = m1.getMemory().subtract(m0.getMemory());
                }
                return diff.intValue();
            }
        };
    }

    private void sortPhysicalMachineArray(PhysicalMachine[] pms) {
        Arrays.sort(pms, makeDecreasingPhysicalMachineComparator());
    }
    
    private void sortVirtualMachineArray(VirtualMachine[] vms) {
        Arrays.sort(vms, makeDecreasingVirtualMachineComparator());
    }
    
    private void savePartialAllocation(MappingVec allocation, VirtualMachineVec leftover_vms) {
        assert(allocation != null);
        assert(leftover_vms != null);
        this.partial_allocation = allocation;
        this.leftover_vms = leftover_vms;
    }
    
    private void reset() {
        clearSolutions();
        clearHeapUsages(this.pm_heap);
        placeMappings(this.pm_heap, this.instance.getPhysicalMachines(), this.mappings);
    }
    
    // =============================================================================================
    // ======================================= PROTECTED API =======================================
    // =============================================================================================

    protected class UsageInfo implements Iterable<VirtualMachine>, Comparable<UsageInfo> {
        
        private PhysicalMachine pm = null;
        private Set<VirtualMachine> placed_vms = new HashSet<VirtualMachine>();
        private BigInteger used_cpu = BigInteger.ZERO;
        private BigInteger used_mem = BigInteger.ZERO;
        private Set<Integer> placed_anti_coloc_job_ids = new HashSet<Integer>();
        
        UsageInfo(PhysicalMachine pm) {
            this.pm = pm;
        }
        
        UsageInfo(UsageInfo info) {
            this.pm = info.pm;
            this.placed_vms = new HashSet<VirtualMachine>(info.placed_vms);
            this.used_cpu = info.used_cpu;
            this.used_mem = info.used_mem;
            this.placed_anti_coloc_job_ids = new HashSet<Integer>(info.placed_anti_coloc_job_ids);
        }
        
        PhysicalMachine getPhysicalMachine() {
            return pm;
        }
        
        boolean isEmpty() {
            return placed_vms.isEmpty();
        }
        
        BigInteger getLeftoverCPU() {
            return pm.getCPU().subtract(used_cpu);
        }
        
        BigInteger getLeftoverMemory() {
            return pm.getMemory().subtract(used_mem);
        }
        
        double getLeftoverCPUPercentile() {
            return Utils.toPercentile(getLeftoverCPU(), pm.getCPU());
        }
        
        double getLeftoverMemoryPercentile() {
            return Utils.toPercentile(getLeftoverMemory(), pm.getMemory());
        }
        
        VirtualMachine[] getPlacedVirtualMachines() {
            VirtualMachine[] vms = new VirtualMachine[placed_vms.size()];
            placed_vms.toArray(vms);
            return vms;
        }
        
        boolean canHostVirtualMachine(VirtualMachine vm) {
            return getLeftoverCPU().compareTo(vm.getCPU()) >= 0 &&
                   getLeftoverMemory().compareTo(vm.getMemory()) >= 0 &&
                   vm.canRunInPhysicalMachine(pm) &&
                   (!vm.isAntiColocatable() ||
                    !placed_anti_coloc_job_ids.contains(new Integer(vm.getJobID())));
        }
        
        void placeVirtualMachine(VirtualMachine vm) {
            placed_vms.add(vm);
            used_cpu = used_cpu.add(vm.getCPU());
            used_mem = used_mem.add(vm.getMemory());
            if (vm.isAntiColocatable()) {
                placed_anti_coloc_job_ids.add(vm.getJobID());
            }
        }
        
        void clear() {
            this.placed_vms.clear();
            this.placed_anti_coloc_job_ids.clear();
            this.used_cpu = BigInteger.ZERO;
            this.used_mem = BigInteger.ZERO;
        }

        @Override
        public Iterator<VirtualMachine> iterator() {
            return placed_vms.iterator();
        }
        
        // FIXME: should use an equals, but it is not implemented for PhysicalMachine
        @Override
        public boolean equals(Object obj) {
            if (obj instanceof UsageInfo) {
                UsageInfo info = (UsageInfo)obj;
                return info.getPhysicalMachine().getID() == pm.getID();
            }
            return false;
        }

        @Override
        public int compareTo(UsageInfo other) {
            PhysicalMachine other_pm = other.pm;
            int diff = pm.getCPU().compareTo(other_pm.getCPU());
            if (diff == 0) {
                diff = pm.getMemory().compareTo(other_pm.getMemory());
            }
            return diff;
        }
        
        public String toString() {
            return String.format("ServerID=%d, UsedCPU=%d, LeftoverCPU=%f, UsedMemory=%d, " +
                                 "LeftoverMemory=%f",
                                 pm.getID(),
                                 used_cpu, getLeftoverCPUPercentile(),
                                 used_mem, getLeftoverMemoryPercentile());
        }
        
    }

    protected Comparator<PhysicalMachine> makeDecreasingPhysicalMachineComparator() {
        return new Comparator<PhysicalMachine>() {
            Comparator<Machine> cap_comp = makeDecreasingMachineComparator();
            @Override
            public int compare(PhysicalMachine m0, PhysicalMachine m1) {
                int diff = cap_comp.compare(m0, m1);
                if (diff == 0) {
                    diff = m0.getID() - m1.getID();
                }
                return diff;
            }
        };
    }
    
    protected Comparator<VirtualMachine> makeDecreasingVirtualMachineComparator() {
        return new Comparator<VirtualMachine>() {
            Comparator<Machine> cap_comp = makeDecreasingMachineComparator();
            @Override
            public int compare(VirtualMachine m0, VirtualMachine m1) {
                int diff = cap_comp.compare(m0, m1);
                if (diff == 0) {
                    diff = m0.getJobID() - m1.getJobID();
                    if (diff == 0) {
                        diff = m0.getIndex() - m1.getIndex();
                    }
                }
                return diff;
            }
        };
    }
    
    protected abstract Comparator<UsageInfo> makeUsageInfoComparator();
    
    // =============================================================================================
    // ======================================== PACKAGE API ========================================
    // =============================================================================================
    
    void enableVirtualMachineShuffling() { this.shuffle_vms = true; }
    
    MappingVec getPartialAllocation() {
        assert(this.leftover_vms != null);
        return partial_allocation;
    }
    
    VirtualMachineVec getLeftoverVirtualMachines() {
        assert(this.leftover_vms != null);
        return leftover_vms;
    }
    
    void setMappings(MappingVec mappings) {
        this.mappings = mappings;
        reset();
    }
    
    void setMaxMemoryMigrationPercentile(double max_mig_percentile) {
        this.max_mig_percentile = max_mig_percentile;
        reset();
    }
    
    // =============================================================================================
    // ======================================== PUBLIC API =========================================
    // =============================================================================================
    
    public BinPackingAllocAlgorithm(VMCwMProblem instance) {
        super(instance);
        this.mappings = this.instance.getMappings();
        this.max_mig_percentile = this.instance.getMaxMigrationPercentile();
        this.pm_heap = buildPhysicalMachineHeap(this.instance.getPhysicalMachines(),
                                                this.mappings);
        this.vm_array = new VirtualMachine[this.instance.getVirtualMachines().size()];
        this.instance.getVirtualMachines().copyTo(this.vm_array);
    }
    
    public MappingVec getSolution() { return getSolutions().get(0); }
    
    @Override
    public void allocate() {
        System.out.println("c WARNING: BinPacking minimizes number of servers, not energy consumption");
        System.out.println("c WARNING: timeout ignored when applying BinPacking");
        if (this.shuffle_vms) {
            PRNG.shuffle(this.vm_array);
        }
        else {
            sortVirtualMachineArray(this.vm_array);
        }
        Set<VirtualMachine> mapped_vms = getMappedVirtualMachines(this.mappings);
        // Do First Fit for the not mapped VMs remaining
        for (int i = 0; i < this.vm_array.length; ++i) {
            VirtualMachine vm = this.vm_array[i];
            if (!mapped_vms.contains(vm)) {
                PhysicalMachine chosen_pm = placeVirtualMachine(this.pm_heap, vm);
                if (chosen_pm == null) {
                    VirtualMachineVec leftover_vms = new VirtualMachineVec();
                    for (int j = i; j < this.vm_array.length; ++j) {
                        leftover_vms.push(this.vm_array[j]);
                    }
                    savePartialAllocation(retrieveMappings(this.pm_heap), leftover_vms);
                    System.out.println("c BinPacking failed to find a feasible placement");
                    return;
                }
            }
        }
        // Discard unnecessary PMs
        PhysicalMachineVec reduced_pms = getUsedPhysicalMachines(this.pm_heap);
        MemoryAggregator mem_agr = new MemoryAggregator();
        reduced_pms.accept(mem_agr);
        BigInteger reduced_mem_cap = mem_agr.memorySum();
        double mig_budget = 
                Utils.scalePercentile(this.instance.getTotalMemoryCapacity(),
                                      this.max_mig_percentile,
                                      reduced_mem_cap);
        mig_budget = Utils.normalizePercentile(mig_budget);
        // Migrate VMs from smaller PMs
        if (this.max_mig_percentile > 0.0) { // FIXME: could be better
            this.pm_heap = buildPhysicalMachineHeap(reduced_pms, retrieveMappings(this.pm_heap));
            Set<VirtualMachine> to_migrate = new HashSet<VirtualMachine>(mapped_vms);
            UsageInfo[] sorted_infos = getSortedInfos(this.pm_heap); // FIXME: another heap in reverse order would be more efficient
            int i;
            for (i = sorted_infos.length-1; i >= 0; --i) {
                UsageInfo info = sorted_infos[i];
                VirtualMachine[] placed_vms = info.getPlacedVirtualMachines();
                if (this.shuffle_vms) {
                    PRNG.shuffle(placed_vms);
                }
                else {
                    sortVirtualMachineArray(placed_vms);
                }
                BigInteger used_mem_cap_to_mig = notMigratedMemoryRequirementSum(placed_vms, to_migrate);
                double used_percentile = Utils.toPercentile(used_mem_cap_to_mig, reduced_mem_cap);
                if (used_percentile <= mig_budget) {
                    PriorityQueue<UsageInfo> tmp_pm_heap = duplicatePhysicalMachineHeap(this.pm_heap);
                    tmp_pm_heap.remove(info);
                    int j;
                    for (j = 0; j < placed_vms.length; ++j) {
                        VirtualMachine vm = placed_vms[j];
                        PhysicalMachine chosen_pm = placeVirtualMachine(tmp_pm_heap, vm);
                        if (chosen_pm == null) {
                            break;
                        }
                    }
                    if (j == placed_vms.length) {
                        for (j = 0; j < placed_vms.length; ++j) {
                            to_migrate.remove(placed_vms[j]);
                        }
                        this.pm_heap = tmp_pm_heap;
                        sorted_infos = getSortedInfos(this.pm_heap); // FIXME: infos for which migration was attempted already could be skipped
                        i = sorted_infos.length;
                        mig_budget -= used_percentile;
                    }
                }
            }
        }
        // Save solution
        MappingVec allocation = retrieveMappings(pm_heap);
        saveSolution(allocation);
        if (foundSolution()) {
            System.out.println("c Solution using " + pm_heap.size() + " PMs found");
        }
        else {
            System.out.println("c Solution found violates constraints");
            savePartialAllocation(allocation, new VirtualMachineVec());
        }
        printElapsedTime();
    }
    
}
