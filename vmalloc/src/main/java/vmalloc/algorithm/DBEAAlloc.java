package vmalloc.algorithm;

import java.util.Properties;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.algorithm.DBEA;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Variation;
import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.core.spi.OperatorFactory;

import vmalloc.evolutionary.VMCwMProblem;

public class DBEAAlloc extends EvolutionaryAllocAlgorithm {

    private class DBEAAllocAlgorithmProvider extends AlgorithmProvider {
        
        @Override
        public Algorithm getAlgorithm(String name, Properties properties, Problem problem) {
            if (name.equals("DBEA")) {
                int divisionsInner = 0, divisionsOuter = 99;
                if (problem.getNumberOfObjectives() == 3) {
                    divisionsOuter = 12;
                }
                assert(divisionsOuter == 2);
                int pop_size = (int)(CombinatoricsUtils.binomialCoefficient(
                        problem.getNumberOfObjectives() + divisionsOuter - 1, divisionsOuter));
                Initialization initialization = makeInitializer(problem, pop_size);
                Variation variation = OperatorFactory.getInstance().getVariation(
                        "sbx+svum", properties, problem);
                return new DBEA(problem, initialization, variation, divisionsOuter, divisionsInner);
            }
            return null;
        }
        
    }

    public DBEAAlloc(VMCwMProblem instance) {
        super(instance, "DBEA", VMCwMProblem.Encoding.BINARY_INTEGER);
        AlgorithmFactory.getInstance().addProvider(new DBEAAllocAlgorithmProvider());
    }

    // TODO: property set methods
    
}
