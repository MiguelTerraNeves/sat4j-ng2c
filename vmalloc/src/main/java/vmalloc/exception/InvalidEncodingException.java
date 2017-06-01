package vmalloc.exception;

public class InvalidEncodingException extends RuntimeException {

    private static final long serialVersionUID = 1L;

    public InvalidEncodingException() { super(); }
    
    public InvalidEncodingException(String string) { super(string); }
    
}
