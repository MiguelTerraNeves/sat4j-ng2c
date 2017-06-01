package vmalloc.exception;

public class InvalidHashFunctionTypeException extends RuntimeException {

    private static final long serialVersionUID = 1L;
    
    public InvalidHashFunctionTypeException() { super(); }
    
    public InvalidHashFunctionTypeException(String string) { super(string); }

}
