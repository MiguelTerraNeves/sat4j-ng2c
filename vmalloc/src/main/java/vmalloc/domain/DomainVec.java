package vmalloc.domain;

import static org.sat4j.GlobalDefs.USE_NG2C;

import java.util.Comparator;
import java.util.Iterator;

import org.omg.PortableInterceptor.USER_EXCEPTION;
import org.sat4j.core.Vec;
import org.sat4j.specs.IVec;

/**
 * Abstract superclass for vectors that contain domain objects.
 * @author Miguel Terra-Neves
 * @param <Type> The type of object contained in the vector.
 */
public abstract class DomainVec<Type> implements IVec<Type> {

    private static final long serialVersionUID = -1183845827636061971L;
    
    /**
     * The actual vector object used to contain domain objects.
     */
    private IVec<Type> vec = USE_NG2C ? new @Gen Vec<Type>() : new Vec<Type>();
    
    /**
     * Creates an instance of a domain object vector.
     */
    public DomainVec() { this(new Vec<Type>()); }
    
    /**
     * Creates an instance of a domain object vector with the contents of another vector.
     * @param vec A vector implementing the {@link IVec} interface.
     */
    public DomainVec(IVec<Type> vec) { vec.copyTo(this.vec); }
    
    /**
     * Creates an instance of domain object vector with the contents of an array.
     * @param array The array.
     */
    public DomainVec(Type[] array) { this.vec = USE_NG2C ? new @Gen Vec<Type>(array) : new Vec<Type>(array); }
    
    /* 
     * Implementation of all methods in the IVec interface. All calls redirected to the actual Vec
     * object. Vec is a final class, so extending Vec was not possible.
     */
    @Override
    public int size() { return this.vec.size(); }
    @Override
    public void shrink(int nofelems) { this.vec.shrink(nofelems); }
    @Override
    public void shrinkTo(int newsize) { this.vec.shrinkTo(newsize); }
    @Override
    public void pop() { this.vec.pop(); }
    @Override
    public void growTo(int newsize, Type pad) { this.vec.growTo(newsize, pad); }
    @Override
    public void ensure(int nsize) { this.vec.ensure(nsize); }
    @Override
    public IVec<Type> push(Type elem) { return this.vec.push(elem); }
    @Override
    public void unsafePush(Type elem) { this.vec.unsafePush(elem); }
    @Override
    public void insertFirst(Type elem) { this.vec.insertFirst(elem); }
    @Override
    public void insertFirstWithShifting(Type elem) { this.vec.insertFirstWithShifting(elem); }
    @Override
    public void clear() { this.vec.clear(); }
    @Override
    public Type last() { return this.vec.last(); }
    @Override
    public Type get(int i) { return this.vec.get(i); }
    @Override
    public void set(int i, Type o) { this.vec.set(i, o); }
    @Override
    public void remove(Type elem) { this.vec.remove(elem); }
    @Override
    public Type delete(int i) { return this.vec.delete(i); }
    @Override
    public void copyTo(IVec<Type> copy) { this.vec.copyTo(copy); }
    @Override
    public <E> void copyTo(E[] dest) { this.vec.copyTo(dest); }
    @Override
    public Type[] toArray() { return this.vec.toArray(); }
    @Override
    public void moveTo(IVec<Type> dest) { this.vec.moveTo(dest); }
    @Override
    public void moveTo(int dest, int source) { this.vec.moveTo(dest, source); }
    @Override
    public void sort(Comparator<Type> comparator) { this.vec.sort(comparator); }
    @Override
    public void sortUnique(Comparator<Type> comparator) { this.vec.sortUnique(comparator); }
    @Override
    public boolean isEmpty() { return this.vec.isEmpty(); }
    @Override
    public Iterator<Type> iterator() { return this.vec.iterator(); }
    @Override
    public boolean contains(Type element) { return this.vec.contains(element); }
    @Override
    public int indexOf(Type element) { return this.vec.indexOf(element); }
   
}
