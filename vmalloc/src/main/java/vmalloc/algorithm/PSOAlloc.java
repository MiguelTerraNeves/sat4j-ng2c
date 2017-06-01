package vmalloc.algorithm;

import vmalloc.evolutionary.VMCwMProblem;

public class PSOAlloc extends EvolutionaryAllocAlgorithm {

    // TODO: implement provider for feasible initialization
    
    public PSOAlloc(VMCwMProblem instance) {
        super(instance, "SMPSO", VMCwMProblem.Encoding.INTEGER);
    }
    
    public void setArchiveSize(int size) { exec = exec.withProperty("archiveSize", size); }
    public void setMutationRate(double rate) { exec = exec.withProperty("pm.rate", rate); }
    public void setDistributionIndex(double index) {
        exec = exec.withProperty("pm.distributionIndex", index);
    }
    
}
