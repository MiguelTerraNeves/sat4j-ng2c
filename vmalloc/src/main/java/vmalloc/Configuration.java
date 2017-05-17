package vmalloc;

import java.io.IOException;
import java.util.Properties;

public class Configuration {

    private static Configuration instance = null;
    
    private Properties props = null;
    
    private Configuration() throws IOException {
        props = new Properties();
        props.load(ClassLoader.getSystemResourceAsStream("config.properties"));
    }
    
    public static Configuration getInstance() {
        if (instance == null) {
            try {
                instance = new Configuration();
            }
            catch (IOException e) {
                throw new RuntimeException("Failed to open configuration properties file", e);
            }
        }
        return instance;
    }
    
    public String getCLASPExecutablePath() { return props.getProperty("claspexe"); }
    
}
