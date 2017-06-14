package vmalloc;

import static org.sat4j.GlobalDefs.USE_NG2C;
import static org.sat4j.GlobalDefs.ANNOTATE_INSTANCE;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;

import org.omg.PortableInterceptor.USER_EXCEPTION;

import vmalloc.domain.Job;
import vmalloc.domain.JobVec;
import vmalloc.domain.Mapping;
import vmalloc.domain.MappingVec;
import vmalloc.domain.PhysicalMachine;
import vmalloc.domain.PhysicalMachineVec;
import vmalloc.domain.VirtualMachine;

public class InputParser {

    private static final String TRUE_STRING = "True";
    private static final String FALSE_STRING = "False";
    
    private String fname;
    private PhysicalMachineVec pms = null;
    private JobVec jobs = null;
    private MappingVec mappings = null;
    
    public InputParser(String fname) {
        this.fname = fname;
    }
    
    private void parsePhysicalMachines(BufferedReader reader) throws IOException {
        int npms = Integer.parseInt(reader.readLine().trim());
        for (int i = 0; i < npms; ++i) {
            String line = reader.readLine();
            String[] tokens = line.split(" ");
            assert(tokens.length == 5);
            int id = Integer.parseInt(tokens[0]);
            BigInteger cpu = ANNOTATE_INSTANCE ? new @Gen BigInteger(tokens[1]) : new BigInteger(tokens[1]);
            BigInteger mem = ANNOTATE_INSTANCE ? new @Gen BigInteger(tokens[2]) : new BigInteger(tokens[2]);
            int idle_consume = Integer.parseInt(tokens[3]);
            int max_consume = Integer.parseInt(tokens[4]);
            this.pms.push(ANNOTATE_INSTANCE ? new @Gen PhysicalMachine(id, cpu, mem, idle_consume, max_consume)
                                            : new PhysicalMachine(id, cpu, mem, idle_consume, max_consume));
        }
    }
    
    private PhysicalMachineVec parseUnallowedPhysicalMachines(String str) {
        String[] tokens = str.split(",");
        PhysicalMachineVec unallowed_pms = ANNOTATE_INSTANCE ? new @Gen PhysicalMachineVec() : new PhysicalMachineVec();
        assert(tokens.length < this.pms.size());
        for (String token_id : tokens) {
            int pm_id = Integer.parseInt(token_id);
            unallowed_pms.push(this.pms.get(pm_id));
        }
        return unallowed_pms;
    }
    
    private void parseVirtualMachines(BufferedReader reader) throws IOException {
        int nvms = Integer.parseInt(reader.readLine().trim());
        for (int i = 0; i < nvms; ++i) {
            String line = reader.readLine();
            String[] tokens = line.split(" ");
            assert(tokens.length == 5 || tokens.length == 6);
            int job_id = Integer.parseInt(tokens[0]);
            while (job_id >= this.jobs.size()) {
                this.jobs.push(ANNOTATE_INSTANCE ? new @Gen Job(job_id) : new Job(job_id));
            }
            int vm_idx = Integer.parseInt(tokens[1]);
            BigInteger cpu = ANNOTATE_INSTANCE ? new @Gen BigInteger(tokens[2]) : new BigInteger(tokens[2]);
            BigInteger mem = ANNOTATE_INSTANCE ? new @Gen BigInteger(tokens[3]) : new BigInteger(tokens[3]);
            boolean anti_coloc = TRUE_STRING.equals(tokens[4]);
            assert(anti_coloc || FALSE_STRING.equals(tokens[4]));
            if (tokens.length > 5) {
                PhysicalMachineVec unallowed_pms = parseUnallowedPhysicalMachines(tokens[5]);
                this.jobs.get(job_id).addVirtualMachine(
                        ANNOTATE_INSTANCE ? new @Gen VirtualMachine(job_id, vm_idx, cpu, mem, anti_coloc, unallowed_pms)
                                          : new VirtualMachine(job_id, vm_idx, cpu, mem, anti_coloc, unallowed_pms));
            }
            else {
                this.jobs.get(job_id).addVirtualMachine(
                        ANNOTATE_INSTANCE ? new @Gen VirtualMachine(job_id, vm_idx, cpu, mem, anti_coloc)
                                          : new VirtualMachine(job_id, vm_idx, cpu, mem, anti_coloc));
            }
        }
    }
    
    private void parseMappings(BufferedReader reader) throws IOException {
        int nmaps = Integer.parseInt(reader.readLine().trim());
        for (int i = 0; i < nmaps; ++i) {
            String line = reader.readLine();
            String[] tokens = line.split(" ");
            assert(tokens.length == 3);
            int job_id = Integer.parseInt(tokens[0]);
            int vm_idx = Integer.parseInt(tokens[1]);
            int pm_id = Integer.parseInt(tokens[2]);
            this.mappings.push(ANNOTATE_INSTANCE ? new @Gen Mapping(this.jobs.get(job_id).getVirtualMachine(vm_idx),
                                                                    this.pms.get(pm_id))
                                                 : new Mapping(this.jobs.get(job_id).getVirtualMachine(vm_idx),
                                                               this.pms.get(pm_id)));
        }
    }
    
    // FIXME: assumes complying input file
    // FIXME: assumes IDs and indexes conforming with the object order in the file
    public void parse() throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(this.fname));
        if (USE_NG2C) {
            System.newAllocGen();
        }
        this.pms = ANNOTATE_INSTANCE ? new @Gen PhysicalMachineVec() : new PhysicalMachineVec();
        this.jobs = ANNOTATE_INSTANCE ? new @Gen JobVec() : new JobVec();
        this.mappings = ANNOTATE_INSTANCE ? new @Gen MappingVec() : new MappingVec();
        parsePhysicalMachines(reader);
        parseVirtualMachines(reader);
        parseMappings(reader);
        if (USE_NG2C) {
            System.setAllocGen(0);
        }
        reader.close();
    }
    
    public PhysicalMachineVec getPhysicalMachines() { return this.pms; }
    public JobVec getJobs() { return this.jobs; }
    public MappingVec getMappings() { return this.mappings; }
    
}
