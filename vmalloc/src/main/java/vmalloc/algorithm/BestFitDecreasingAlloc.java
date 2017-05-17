package vmalloc.algorithm;

import java.util.Comparator;

import vmalloc.domain.PhysicalMachine;
import vmalloc.evolutionary.VMCwMProblem;

public class BestFitDecreasingAlloc extends BinPackingAllocAlgorithm {

    public BestFitDecreasingAlloc(VMCwMProblem instance) { super(instance); }

    @Override
    protected Comparator<UsageInfo> makeUsageInfoComparator() {
        return new Comparator<UsageInfo>() {
            Comparator<PhysicalMachine> cap_comp = makeDecreasingPhysicalMachineComparator();
            @Override
            public int compare(UsageInfo ui1, UsageInfo ui2) {
                double diff = ui1.getLeftoverCPUPercentile() - ui2.getLeftoverCPUPercentile();
                if (diff == 0.0) {
                    diff = ui1.getLeftoverMemoryPercentile() - ui2.getLeftoverMemoryPercentile();
                    if (diff == 0.0) {
                        return cap_comp.compare(ui1.getPhysicalMachine(), ui2.getPhysicalMachine());
                    }
                }
                assert(diff != 0.0);
                return (diff > 0.0) ? 1 : -1;
            }
        };
    }
    
}
