package vmalloc;

public class Clock {

    private static Clock instance = null;
    
    private long start;
    
    private Clock() { }
    
    public static Clock getInstance() {
        if (instance == null) {
            instance = new Clock();
        }
        return instance;
    }
    
    public void reset() { this.start = System.nanoTime(); }
    public double getElapsed() { return (double)(System.nanoTime() - this.start) / 1000000000.0; }
    
}
