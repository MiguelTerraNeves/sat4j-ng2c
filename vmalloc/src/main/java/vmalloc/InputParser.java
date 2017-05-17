package vmalloc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;

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
    private PhysicalMachineVec pms;
    private JobVec jobs;
    private MappingVec mappings;
    
    public InputParser(String fname) {
        this.fname = fname;
        this.pms = new PhysicalMachineVec();
        this.jobs = new JobVec();
        this.mappings = new MappingVec();
    }
    
    private void parsePhysicalMachines(BufferedReader reader) throws IOException {
        int npms = Integer.parseInt(reader.readLine().trim());
        for (int i = 0; i < npms; ++i) {
            String line = reader.readLine();
            String[] tokens = line.split(" ");
            assert(tokens.length == 5);
            int id = Integer.parseInt(tokens[0]);
            BigInteger cpu = new BigInteger(tokens[1]);
            BigInteger mem = new BigInteger(tokens[2]);
            int idle_consume = Integer.parseInt(tokens[3]);
            int max_consume = Integer.parseInt(tokens[4]);
            this.pms.push(new PhysicalMachine(id, cpu, mem, idle_consume, max_consume));
        }
    }
    
    private PhysicalMachineVec parseUnallowedPhysicalMachines(String str) {
        String[] tokens = str.split(",");
        PhysicalMachineVec unallowed_pms = new PhysicalMachineVec();
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
                this.jobs.push(new Job(job_id));
            }
            int vm_idx = Integer.parseInt(tokens[1]);
            BigInteger cpu = new BigInteger(tokens[2]);
            BigInteger mem = new BigInteger(tokens[3]);
            boolean anti_coloc = TRUE_STRING.equals(tokens[4]);
            assert(anti_coloc || FALSE_STRING.equals(tokens[4]));
            if (tokens.length > 5) {
                PhysicalMachineVec unallowed_pms = parseUnallowedPhysicalMachines(tokens[5]);
                this.jobs.get(job_id).addVirtualMachine(
                        new VirtualMachine(job_id, vm_idx, cpu, mem, anti_coloc, unallowed_pms));
            }
            else {
                this.jobs.get(job_id).addVirtualMachine(
                        new VirtualMachine(job_id, vm_idx, cpu, mem, anti_coloc));
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
            this.mappings.push(new Mapping(this.jobs.get(job_id).getVirtualMachine(vm_idx),
                                           this.pms.get(pm_id)));
        }
    }
    
    // FIXME: assumes complying input file
    // FIXME: assumes IDs and indexes conforming with the object order in the file
    public void parse() throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(this.fname));
        parsePhysicalMachines(reader);
        parseVirtualMachines(reader);
        parseMappings(reader);
        reader.close();
    }
    
    public PhysicalMachineVec getPhysicalMachines() { return this.pms; }
    public JobVec getJobs() { return this.jobs; }
    public MappingVec getMappings() { return this.mappings; }
    
}
