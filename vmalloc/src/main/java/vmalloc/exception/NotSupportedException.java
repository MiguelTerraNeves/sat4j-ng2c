package vmalloc.exception;

public class NotSupportedException extends RuntimeException {
    
    private static final long serialVersionUID = 1L;
    
    public NotSupportedException() { super(); }
    
    public NotSupportedException(String string) { super(string); }
    
}
