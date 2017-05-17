package vmalloc.algorithm;

import java.util.Properties;

import org.moeaframework.algorithm.GDE3;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedSortingPopulation;
import org.moeaframework.core.Problem;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.operator.real.DifferentialEvolution;
import org.moeaframework.core.operator.real.DifferentialEvolutionSelection;
import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.util.TypedProperties;

import vmalloc.evolutionary.VMCwMProblem;

public class DEAlloc extends EvolutionaryAllocAlgorithm {

    private class DEAllocAlgorithmProvider extends AlgorithmProvider {
        
        @Override
        public Algorithm getAlgorithm(String name, Properties properties, Problem problem) {
            if (name.equals("GDE3")) {
                TypedProperties typed_props = new TypedProperties(properties);
                int pop_size = typed_props.getInt("populationSize", 100);
                DominanceComparator comparator = new EfficientParetoDominanceComparator();
                NondominatedSortingPopulation population = new NondominatedSortingPopulation(comparator);
                Initialization initialization = makeInitializer(problem, pop_size);
                DifferentialEvolutionSelection selection =  new DifferentialEvolutionSelection();
                DifferentialEvolution variation = 
                        (DifferentialEvolution)OperatorFactory.getInstance().getVariation("de",
                                                                                          properties,
                                                                                          problem);
                return new GDE3(problem, population, comparator, selection, variation, initialization);
            }
            return null;
        }
        
    }
    
    public DEAlloc(VMCwMProblem instance) {
        super(instance, "GDE3", VMCwMProblem.Encoding.INTEGER);
        AlgorithmFactory.getInstance().addProvider(new DEAllocAlgorithmProvider());
    }
    
    public void setCrossoverRate(double rate) { exec = exec.withProperty("de.crossoverRate", rate); }
    public void setStepSize(double step) { exec = exec.withProperty("de.stepSize", step); }

}
