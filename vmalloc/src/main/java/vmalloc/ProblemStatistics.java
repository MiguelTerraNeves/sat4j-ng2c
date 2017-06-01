package vmalloc;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.HashSet;
import java.util.Set;

import org.sat4j.specs.IVec;

import vmalloc.domain.Job;
import vmalloc.domain.Mapping;
import vmalloc.domain.PhysicalMachine;

public class ProblemStatistics {

    private IVec<PhysicalMachine> pms;
    private IVec<Job> jobs;
    private IVec<Mapping> mappings;
    private double max_mig_percentile = 1.0;
    
    // PM stats
    private BigInteger total_cpu_cap = BigInteger.ZERO;
    private BigInteger total_mem_cap = BigInteger.ZERO;
    private BigInteger min_cpu_cap = BigInteger.ZERO;
    private BigInteger min_mem_cap = BigInteger.ZERO;
    private BigInteger max_cpu_cap = BigInteger.ZERO;
    private BigInteger max_mem_cap = BigInteger.ZERO;
    
    // Job stats
    private int total_vms = 0;
    private int min_vms = 0;
    private int max_vms = 0;
    
    // VM stats
    private BigInteger total_cpu_req = BigInteger.ZERO;
    private BigInteger total_mem_req = BigInteger.ZERO;
    private BigInteger min_cpu_req = BigInteger.ZERO;
    private BigInteger min_mem_req = BigInteger.ZERO;
    private BigInteger max_cpu_req = BigInteger.ZERO;
    private BigInteger max_mem_req = BigInteger.ZERO;
    
    // Mapping stats
    private BigInteger total_cpu_map_req = BigInteger.ZERO;
    private BigInteger total_mem_map_req = BigInteger.ZERO;
    private BigInteger min_cpu_map_req = BigInteger.ZERO;
    private BigInteger min_mem_map_req = BigInteger.ZERO;
    private BigInteger max_cpu_map_req = BigInteger.ZERO;
    private BigInteger max_mem_map_req = BigInteger.ZERO;
    private int used_pms = 0;
    
    public ProblemStatistics(IVec<PhysicalMachine> pms,
                             IVec<Job> jobs,
                             IVec<Mapping> mappings,
                             double max_mig_percentile) {
        this.pms = pms;
        this.jobs = jobs;
        this.mappings = mappings;
        this.max_mig_percentile = max_mig_percentile;
        initStatistics();
    }
    
    private void initStatistics() {
        for (int i = 0; i < this.pms.size(); ++i) {
            BigInteger cpu_cap = this.pms.get(i).getCPU();
            BigInteger mem_cap = this.pms.get(i).getMemory();
            this.total_cpu_cap = this.total_cpu_cap.add(cpu_cap);
            this.total_mem_cap = this.total_mem_cap.add(mem_cap);
            if (this.min_cpu_cap.compareTo(cpu_cap) > 0 ||
                this.min_cpu_cap.equals(BigInteger.ZERO)) {
                this.min_cpu_cap = cpu_cap;
            }
            if (this.max_cpu_cap.compareTo(cpu_cap) < 0 ||
                this.max_cpu_cap.equals(BigInteger.ZERO)) {
                this.max_cpu_cap = cpu_cap;
            }
            if (this.min_mem_cap.compareTo(mem_cap) > 0 ||
                this.min_mem_cap.equals(BigInteger.ZERO)) {
                this.min_mem_cap = mem_cap;
            }
            if (this.max_mem_cap.compareTo(mem_cap) < 0 ||
                this.max_mem_cap.equals(BigInteger.ZERO)) {
                this.max_mem_cap = mem_cap;
            }
        }
        for (int i = 0; i < this.jobs.size(); ++i) {
            Job job = this.jobs.get(i);
            this.total_vms += job.nVirtualMachines();
            if (this.min_vms > job.nVirtualMachines() || this.min_vms == 0) {
                this.min_vms = job.nVirtualMachines();
            }
            if (this.max_vms < job.nVirtualMachines() || this.max_vms == 0) {
                this.max_vms = job.nVirtualMachines();
            }
            for (int j = 0; j < job.nVirtualMachines(); ++j) {
                BigInteger cpu_req = job.getVirtualMachine(j).getCPU();
                BigInteger mem_req = job.getVirtualMachine(j).getMemory();
                this.total_cpu_req = this.total_cpu_req.add(cpu_req);
                this.total_mem_req = this.total_mem_req.add(mem_req);
                if (this.min_cpu_req.compareTo(cpu_req) > 0 ||
                    this.min_cpu_req.equals(BigInteger.ZERO)) {
                    this.min_cpu_req = cpu_req;
                }
                if (this.max_cpu_req.compareTo(cpu_req) < 0 ||
                    this.max_cpu_req.equals(BigInteger.ZERO)) {
                    this.max_cpu_req = cpu_req;
                }
                if (this.min_mem_req.compareTo(mem_req) > 0 ||
                    this.min_mem_req.equals(BigInteger.ZERO)) {
                    this.min_mem_req = mem_req;
                }
                if (this.max_mem_req.compareTo(mem_req) < 0 ||
                    this.max_mem_req.equals(BigInteger.ZERO)) {
                    this.max_mem_req = mem_req;
                }
            }
        }
        Set<PhysicalMachine> used_pms = new HashSet<PhysicalMachine>();
        for (int i = 0; i < this.mappings.size(); ++i) {
            BigInteger cpu_req = this.mappings.get(i).getVirtualMachine().getCPU();
            BigInteger mem_req = this.mappings.get(i).getVirtualMachine().getMemory();
            this.total_cpu_map_req = this.total_cpu_map_req.add(cpu_req);
            this.total_mem_map_req = this.total_mem_map_req.add(mem_req);
            if (this.min_cpu_map_req.compareTo(cpu_req) > 0 ||
                this.min_cpu_map_req.equals(BigInteger.ZERO)) {
                this.min_cpu_map_req = cpu_req;
            }
            if (this.max_cpu_map_req.compareTo(cpu_req) < 0 ||
                this.max_cpu_map_req.equals(BigInteger.ZERO)) {
                this.max_cpu_map_req = cpu_req;
            }
            if (this.min_mem_map_req.compareTo(mem_req) > 0 ||
                this.min_mem_map_req.equals(BigInteger.ZERO)) {
                this.min_mem_map_req = mem_req;
            }
            if (this.max_mem_map_req.compareTo(mem_req) < 0 ||
                this.max_mem_map_req.equals(BigInteger.ZERO)) {
                this.max_mem_map_req = mem_req;
            }
            used_pms.add(this.mappings.get(i).getPhysicalMachine());
        }
        this.used_pms = used_pms.size();
    }
    
    private BigDecimal divideBigIntegers(BigInteger int1, BigInteger int2) {
        return Utils.divideBigIntegers(int1, int2, RoundingMode.HALF_UP);
    }
    
    public int nPhysicalMachines() { return this.pms.size(); }
    public int nJobs() { return this.jobs.size(); }
    public int nVirtualMachines() { return this.total_vms; }
    public int nMappings() { return this.mappings.size(); }
    public int nUsedPhysicalMachines() { return this.used_pms; }
    
    public double getAllowedMigrationPercentile() { return this.max_mig_percentile; }
    public BigInteger getMigrationCapacityBudget() {
        BigDecimal budget =
                new BigDecimal(this.total_mem_cap).multiply(new BigDecimal(this.max_mig_percentile));
        return budget.toBigInteger();
    }
    
    public BigInteger getTotalCPUCapacity() { return this.total_cpu_cap; }
    public BigInteger getTotalMemoryCapacity() { return this.total_mem_cap; }
    public BigInteger getMinimumCPUCapacity() { return this.min_cpu_cap; }
    public BigInteger getMinimumMemoryCapacity() { return this.min_mem_cap; }
    public BigInteger getMaximumCPUCapacity() { return this.max_cpu_cap; }
    public BigInteger getMaximumMemoryCapacity() { return this.max_mem_cap; }
    public BigDecimal getAverageCPUCapacity() {
        return divideBigIntegers(this.total_cpu_cap, BigInteger.valueOf(nPhysicalMachines()));
    }
    public BigDecimal getAverageMemoryCapacity() {
        return divideBigIntegers(this.total_mem_cap, BigInteger.valueOf(nPhysicalMachines()));
    }
    
    public int getMinimumVirtualMachinesPerJob() { return this.min_vms; }
    public int getMaximumVirtualMachinesPerJob() { return this.max_vms; }
    public double getAverageVirtualMachinesPerJob() { return this.total_vms/this.jobs.size(); }
    
    public BigInteger getTotalCPURequirements() { return this.total_cpu_req; }
    public BigInteger getTotalMemoryRequirements() { return this.total_mem_req; }
    public BigInteger getMinimumCPURequirements() { return this.min_cpu_req; }
    public BigInteger getMinimumMemoryRequirements() { return this.min_mem_req; }
    public BigInteger getMaximumCPURequirements() { return this.max_cpu_req; }
    public BigInteger getMaximumMemoryRequirements() { return this.max_mem_req; }
    public BigDecimal getAverageCPURequirements() {
        return divideBigIntegers(this.total_cpu_req, BigInteger.valueOf(nVirtualMachines()));
    }
    public BigDecimal getAverageMemoryRequirements() {
        return divideBigIntegers(this.total_mem_req, BigInteger.valueOf(nVirtualMachines()));
    }
    
    public BigInteger getMappingTotalCPURequirements() { return this.total_cpu_map_req; }
    public BigInteger getMappingTotalMemoryRequirements() { return this.total_mem_map_req; }
    public BigInteger getMappingMinimumCPURequirements() { return this.min_cpu_map_req; }
    public BigInteger getMappingMinimumMemoryRequirements() { return this.min_mem_map_req; }
    public BigInteger getMappingMaximumCPURequirements() { return this.max_cpu_map_req; }
    public BigInteger getMappingMaximumMemoryRequirements() { return this.max_mem_map_req; }
    public BigDecimal getMappingAverageCPURequirements() {
        return divideBigIntegers(
                this.total_cpu_map_req, BigInteger.valueOf(nUsedPhysicalMachines()));
    }
    public BigDecimal getMappingAverageMemoryRequirements() {
        return divideBigIntegers(
                this.total_mem_map_req, BigInteger.valueOf(nUsedPhysicalMachines()));
    }
    
    public double getCPUUsagePercentile() {
        return divideBigIntegers(this.total_cpu_req, this.total_cpu_cap).doubleValue();
    }
    public double getMemoryUsagePercentile() {
        return divideBigIntegers(this.total_mem_req, this.total_mem_cap).doubleValue();
    }
    public double getMappingCPUUsagePercentile() {
        return divideBigIntegers(this.total_cpu_map_req, this.total_cpu_cap).doubleValue();
    }
    public double getMappingMemoryUsagePercentile() {
        return divideBigIntegers(this.total_mem_map_req, this.total_mem_cap).doubleValue();
    }
    
    public void printStatistics() {
        System.out.println("c =========================== Problem Statistics ===========================");
        System.out.println("c          | PMs        | Jobs       | VMs        | Mappings   | Used PMs");
        System.out.println("c  Number  | " +
                           String.format("%10d", nPhysicalMachines()) + " | " + 
                           String.format("%10d", nJobs()) + " | " +
                           String.format("%10d", nVirtualMachines()) + " | " +
                           String.format("%10d", nMappings()) + " | " +
                           String.format("%10d", nUsedPhysicalMachines()));
        System.out.println("c  -------------------------------------------------------------------------");
        System.out.println("c  -------------------------------------------------------------------------");
        System.out.println("c                               | CPU                | RAM");
        System.out.println("c  -------------------------------------------------------------------------");
        System.out.println("c  Total Capacity               | " +
                           String.format("%18d", getTotalCPUCapacity()) + " | " +
                           String.format("%18d", getTotalMemoryCapacity()));
        System.out.println("c  Minimum Capacity             | " +
                              String.format("%18d", getMinimumCPUCapacity()) + " | " +
                              String.format("%18d", getMinimumMemoryCapacity()));
        System.out.println("c  Maximum Capacity             | " +
                              String.format("%18d", getMaximumCPUCapacity()) + " | " +
                              String.format("%18d", getMaximumMemoryCapacity()));
        System.out.println("c  Average Capacity             | " +
                              String.format("%18.2f", getAverageCPUCapacity()) + " | " +
                              String.format("%18.2f", getAverageMemoryCapacity()));
        System.out.println("c  -------------------------------------------------------------------------");
        System.out.println("c  Total Requirements           | " +
                           String.format("%18d", getTotalCPURequirements()) + " | " +
                           String.format("%18d", getTotalMemoryRequirements()));
        System.out.println("c  Minimum Requirements         | " +
                              String.format("%18d", getMinimumCPURequirements()) + " | " +
                              String.format("%18d", getMinimumMemoryRequirements()));
        System.out.println("c  Maximum Requirements         | " +
                              String.format("%18d", getMaximumCPURequirements()) + " | " +
                              String.format("%18d", getMaximumMemoryRequirements()));
        System.out.println("c  Average Requirements         | " +
                              String.format("%18.2f", getAverageCPURequirements()) + " | " +
                              String.format("%18.2f", getAverageMemoryRequirements()));
        System.out.println("c  Usage Percentile             | " +
                              String.format("%18.16f", getCPUUsagePercentile()) + " | " +
                              String.format("%18.16f", getMemoryUsagePercentile()));
        if (nMappings() > 0) {
            System.out.println("c  -------------------------------------------------------------------------");
            System.out.println("c  Total Mapping Requirements   | " +
                    String.format("%18d", getMappingTotalCPURequirements()) + " | " +
                    String.format("%18d", getMappingTotalMemoryRequirements()));
            System.out.println("c  Minimum Mapping Requirements | " +
                    String.format("%18d", getMappingMinimumCPURequirements()) + " | " +
                    String.format("%18d", getMappingMinimumMemoryRequirements()));
            System.out.println("c  Maximum Mapping Requirements | " +
                    String.format("%18d", getMappingMaximumCPURequirements()) + " | " +
                    String.format("%18d", getMappingMaximumMemoryRequirements()));
            System.out.println("c  Average Mapping Requirements | " +
                    String.format("%18.2f", getMappingAverageCPURequirements()) + " | " +
                    String.format("%18.2f", getMappingAverageMemoryRequirements()));
            System.out.println("c  Mapping Usage Percentile     | " +
                    String.format("%18.16f", getMappingCPUUsagePercentile()) + " | " +
                    String.format("%18.16f", getMappingMemoryUsagePercentile()));
            System.out.println("c  Migration Budget Percentile  |                  - | " +
                    String.format("%18.16f", getAllowedMigrationPercentile()));
            System.out.println("c  Migration Capacity Budget    |                  - | " +
                    String.format("%18d", getMigrationCapacityBudget()));
        }
        System.out.println("c ==========================================================================");
    }
    
}
