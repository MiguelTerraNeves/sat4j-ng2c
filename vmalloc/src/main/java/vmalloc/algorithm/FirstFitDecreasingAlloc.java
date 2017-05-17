package vmalloc.algorithm;

import java.util.Comparator;

import vmalloc.domain.PhysicalMachine;
import vmalloc.evolutionary.VMCwMProblem;

public class FirstFitDecreasingAlloc extends BinPackingAllocAlgorithm {
    
    public FirstFitDecreasingAlloc(VMCwMProblem instance) { super(instance); }
    
    @Override
    protected Comparator<UsageInfo> makeUsageInfoComparator() {
        return new Comparator<UsageInfo>() {
            Comparator<PhysicalMachine> cap_comp = makeDecreasingPhysicalMachineComparator();
            @Override
            public int compare(UsageInfo ui1, UsageInfo ui2) {
                return cap_comp.compare(ui1.getPhysicalMachine(), ui2.getPhysicalMachine());
            }
        };
    }

}
