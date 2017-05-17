package vmalloc.algorithm;

import java.util.Properties;

import org.moeaframework.algorithm.NSGAII;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedSortingPopulation;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Variation;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.CrowdingComparator;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.comparator.ParetoDominanceComparator;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.spi.AlgorithmFactory;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.util.TypedProperties;

import vmalloc.evolutionary.VMCwMProblem;

public class GAAlloc extends EvolutionaryAllocAlgorithm {
    
    private class GAAllocAlgorithmProvider extends AlgorithmProvider {
        
        @Override
        public Algorithm getAlgorithm(String name, Properties properties, Problem problem) {
            if (name.equals("NSGAII")) {
                TypedProperties typed_props = new TypedProperties(properties);
                int pop_size = typed_props.getInt("populationSize", 100);
                Initialization initialization = makeInitializer(problem, pop_size);
                DominanceComparator comparator = new EfficientParetoDominanceComparator();
                NondominatedSortingPopulation population = new NondominatedSortingPopulation(comparator);
                TournamentSelection selection =
                        new TournamentSelection(2, new ChainedComparator(new ParetoDominanceComparator(),
                                                                         new CrowdingComparator()));
                Variation variation = OperatorFactory.getInstance().getVariation("ux+svum",
                                                                                 properties,
                                                                                 problem);
                return new NSGAII(problem, population, null, selection, variation, initialization);
            }
            return null;
        }
        
    }

    public GAAlloc(VMCwMProblem instance) {
        super(instance, "NSGAII", VMCwMProblem.Encoding.BINARY_INTEGER);
        AlgorithmFactory.getInstance().addProvider(new GAAllocAlgorithmProvider());
    }

    public void setCrossoverRate(double rate) { exec = exec.withProperty("ux.rate", rate); }
    public void setMutationRate(double rate) { exec = exec.withProperty("svum.rate", rate); }

}
