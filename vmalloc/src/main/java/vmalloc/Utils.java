package vmalloc;

import static org.sat4j.GlobalDefs.USE_NG2C;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.sat4j.specs.IVec;
import org.sat4j.specs.IVecInt;

import vmalloc.domain.JobVec;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;
import vmalloc.domain.VirtualMachineVec;

public class Utils {
    
    public static final String NEW_LINE = System.getProperty("line.separator");
    
    private static final PrintStream STDOUT = System.out;
    private static final PrintStream DUMMY_OUT = new PrintStream(new OutputStream() {
        @Override
        public void write(int b) throws IOException {}
    });
    
    public static void stdoutDisable() { System.setOut(DUMMY_OUT); }
    public static void stdoutEnable() { System.setOut(STDOUT); }

    public static File makeTemporaryFile(String prefix, String suffix, boolean unique, File tmp_dir) {
        if (unique) {
            prefix += "_" + Long.toString(System.currentTimeMillis());
        }
        File tmp_file = null;
        try {
            if (tmp_dir != null) {
                tmp_file = File.createTempFile(prefix, suffix, tmp_dir);
            }
            else {
                tmp_file = File.createTempFile(prefix, suffix);
            }
        }
        catch (IOException e) {
            throw new RuntimeException("Failed to create temporary file " + prefix + suffix, e);
        }
        tmp_file.deleteOnExit();
        return tmp_file;
    }
    
    public static File makeTemporaryFile(String prefix, String suffix, boolean unique) {
        return makeTemporaryFile(prefix, suffix, unique, null);
    }
    
    public static double normalizePercentile(double percentile) {
        if (percentile > 1.0) {
            percentile = 1.0;
        }
        else if (percentile < 0.0) {
            percentile = 0.0;
        }
        return percentile;
    }
    
    public static BigDecimal divideBigIntegers(BigInteger int1, BigInteger int2, RoundingMode mode) {
        return new BigDecimal(int1).divide(new BigDecimal(int2), 10, mode);
    }
    
    public static BigDecimal divideBigDecimals(BigDecimal d1, BigDecimal d2, RoundingMode mode) {
        return d1.divide(d2, 10, mode);
    }
    
    public static BigDecimal floor(BigDecimal n) {
        return n.setScale(0, RoundingMode.FLOOR);
    }
    
    public static BigDecimal ceil(BigDecimal n) {
        return n.setScale(0, RoundingMode.CEILING);
    }
    
    public static double toPercentile(BigInteger q, BigInteger d) {
        return divideBigIntegers(q, d, RoundingMode.UP).doubleValue();
    }
    public static double toPercentile(BigDecimal q, BigDecimal d) {
        return divideBigDecimals(q, d, RoundingMode.UP).doubleValue();
    }
    
    public static double scalePercentile(BigInteger before_ref, double before_per, BigInteger new_ref) {
        BigDecimal inv_k = new BigDecimal(before_per).multiply(new BigDecimal(before_ref));
        return divideBigDecimals(inv_k, new BigDecimal(new_ref), RoundingMode.DOWN).doubleValue();
    }
    
    public static double min(double[] array) {
        double min = Double.MAX_VALUE;
        for (int i = 0; i < array.length; ++i) {
            min = (array[i] < min) ? array[i] : min;
        }
        return min;
    }
    
    public static double max(double[] array) {
        double min = Double.MIN_VALUE;
        for (int i = 0; i < array.length; ++i) {
            min = (array[i] > min) ? array[i] : min;
        }
        return min;
    }
    
    // FIXME: these two are very similar
    public static double sum(double[] array) {
        double sum = 0.0;
        for (int i = 0; i < array.length; ++i) {
            sum += array[i];
        }
        return sum;
    }
    public static long sum(long[] array) {
        long sum = 0;
        for (int i = 0; i < array.length; ++i) {
            sum += array[i];
        }
        return sum;
    }
    
    public static boolean valueInVec(int val, IVecInt vec) {
        boolean result = false;
        for (int i = 0; !result && i < vec.size(); ++i) {
            if (val == vec.get(i)) {
                result = true;
            }
        }
        return result;
    }
    
    public static boolean valuesIntersect(IVecInt vec1, IVecInt vec2) {
        boolean result = false;
        for (int i = 0; !result && i < vec1.size(); ++i) {
            if (valueInVec(vec1.get(i), vec2)) {
                result = true;
            }
        }
        return result;
    }
    
    // FIXME: these two are very similar
    public static BigInteger bigIntegerVecSum(IVec<BigInteger> vec) {
        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < vec.size(); ++i) {
            sum = sum.add(vec.get(i));
        }
        return sum;
    }
    public static BigDecimal bigDecimalVecSum(IVec<BigDecimal> vec) {
        BigDecimal sum = BigDecimal.ZERO;
        for (int i = 0; i < vec.size(); ++i) {
            sum = sum.add(vec.get(i));
        }
        return sum;
    }
    
    public static double avg(double[] array) {
        return sum(array) / array.length;
    }
    
    public static long[] bigIntegerArrayToLongArray(BigInteger[] vec) {
        long[] long_vec = new long[vec.length];
        for (int i = 0; i < vec.length; ++i) {
            long_vec[i] = vec[i].longValue();
        }
        return long_vec;
    }

    public static Map<Integer, Integer> makePhysicalMachineIDtoIndexMap_core(PhysicalMachineVec pms, Map<Integer, Integer> pm_id_to_idx) {
        for (int i = 0; i < pms.size(); ++i) {
            pm_id_to_idx.put(USE_NG2C ? new @Gen Integer(pms.get(i).getID()) : new Integer(pms.get(i).getID()),
                             USE_NG2C ? new @Gen Integer(i) : new Integer(i));
        }
        return pm_id_to_idx;
    }
    public static Map<Integer, Integer> makePhysicalMachineIDtoIndexMap(PhysicalMachineVec pms) {
        return makePhysicalMachineIDtoIndexMap_core(pms, new HashMap<Integer, Integer>());
    }
    public static Map<Integer, Integer> makePhysicalMachineIDtoIndexMap_Gen(PhysicalMachineVec pms) {
        return makePhysicalMachineIDtoIndexMap_core(pms, USE_NG2C ? new @Gen HashMap<Integer, Integer>() : new HashMap<Integer, Integer>());
    }
    
    public static Map<String, Integer> makeVirtualMachineIDtoIndexMap_core(VirtualMachineVec vms, Map<String, Integer> vm_id_to_idx) {
        for (int i = 0; i < vms.size(); ++i) {
            vm_id_to_idx.put(vms.get(i).getID(), USE_NG2C ? new @Gen Integer(i) : new Integer(i));
        }
        return vm_id_to_idx;
    }
    public static Map<String, Integer> makeVirtualMachineIDtoIndexMap(VirtualMachineVec vms) {
        return makeVirtualMachineIDtoIndexMap_core(vms, new HashMap<String, Integer>());
    }
    public static Map<String, Integer> makeVirtualMachineIDtoIndexMap_Gen(VirtualMachineVec vms) {
        return makeVirtualMachineIDtoIndexMap_core(vms, USE_NG2C ? new @Gen HashMap<String, Integer>() : new HashMap<String, Integer>());
    }
    
    public static boolean allocationIsValid(PhysicalMachineVec pms,
                                            JobVec jobs,
                                            MappingVec mappings,
                                            double max_mig_percentile,
                                            MappingVec allocation) {
        // Check if VMs are mapped to a single PM
        Set<String> mapped_vms = new HashSet<String>();
        for (int i = 0; i < allocation.size(); ++i) {
            String vm_id = allocation.get(i).getVirtualMachine().getID();
            if (mapped_vms.contains(vm_id)) {
                System.out.println("ERROR: VM " + vm_id + " placed in multiple PMs");
                return false;
            }
            mapped_vms.add(vm_id);
        }
        // Check if all VMs are mapped
        int nvms = 0;
        for (int i = 0; i < jobs.size(); ++i) {
            nvms += jobs.get(i).nVirtualMachines();
        }
        if (nvms != mapped_vms.size()) {
            System.out.println("ERROR: only " + mapped_vms.size() + " mapped out of a total of " + nvms);
            return false;
        }
        // Check if capacities are not exceeded
        Map<Integer, VirtualMachineVec> vms_per_pm = new HashMap<Integer, VirtualMachineVec>();
        for (int i = 0; i < pms.size(); ++i) {
            vms_per_pm.put(new Integer(pms.get(i).getID()), new VirtualMachineVec());
        }
        for (int i = 0; i < allocation.size(); ++i) {
            Integer pm_id = new Integer(allocation.get(i).getPhysicalMachine().getID());
            VirtualMachine vm = allocation.get(i).getVirtualMachine();
            vms_per_pm.get(pm_id).push(vm);
        }
        for (int i = 0; i < pms.size(); ++i) {
            PhysicalMachine pm = pms.get(i);
            BigInteger cpu_usage = BigInteger.ZERO;
            BigInteger mem_usage = BigInteger.ZERO;
            IVec<VirtualMachine> vms = vms_per_pm.get(new Integer(pm.getID()));
            for (int j = 0; j < vms.size(); ++j) {
                cpu_usage = cpu_usage.add(vms.get(j).getCPU());
                mem_usage = mem_usage.add(vms.get(j).getMemory());
            }
            if (cpu_usage.compareTo(pm.getCPU()) > 0 || mem_usage.compareTo(pm.getMemory()) > 0) {
                System.out.println("ERROR: capacity exceeded for PM " + pm.getID());
                return false;
            }
        }
        // Check if anti-colocation is not violated
        for (int i = 0; i < pms.size(); ++i) {
            VirtualMachineVec vms = vms_per_pm.get(new Integer(pms.get(i).getID()));
            Set<Integer> anti_coloc_jobs = new HashSet<Integer>();
            for (int j = 0; j < vms.size(); ++j) {
                VirtualMachine vm = vms.get(j);
                Integer job_id = new Integer(vm.getJobID());
                if (vm.isAntiColocatable()) {
                    if (anti_coloc_jobs.contains(job_id)) {
                        System.out.println("ERROR: anti-colocation violated for job " + job_id);
                        return false;
                    }
                    anti_coloc_jobs.add(job_id);
                }
            }
        }
        // Check if migration percentile is not exceeded
        Map<String, Integer> vm_to_pm = new HashMap<String, Integer>();
        for (int i = 0; i < allocation.size(); ++i) {
            vm_to_pm.put(allocation.get(i).getVirtualMachine().getID(),
                         new Integer(allocation.get(i).getPhysicalMachine().getID()));
        }
        BigDecimal total_mem_cap = BigDecimal.ZERO;
        for (int i = 0; i < pms.size(); ++i) {
            total_mem_cap = total_mem_cap.add(new BigDecimal(pms.get(i).getMemory()));
        }
        BigInteger max_mig_mem =
                total_mem_cap.multiply(new BigDecimal(max_mig_percentile)).toBigInteger();
        BigInteger migged_mem = BigInteger.ZERO;
        for (int i = 0; i < mappings.size(); ++i) {
            VirtualMachine vm = mappings.get(i).getVirtualMachine();
            Integer prev_pm_id = new Integer(mappings.get(i).getPhysicalMachine().getID());
            if (!vm_to_pm.get(vm.getID()).equals(prev_pm_id)) {
                migged_mem = migged_mem.add(vm.getMemory());
            }
        }
        if (migged_mem.compareTo(max_mig_mem) > 0) {
            System.out.println("ERROR: migration budget exceeded");
            return false;
        }
        return true;
    }
    
}
