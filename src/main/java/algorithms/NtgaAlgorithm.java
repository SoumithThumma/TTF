package algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import model.Solution;
import model.TravelingThiefProblem;


// an improved NTGA (Non-dominated Tournament Genetic Algorithm) implementation for TTP (traveling thief problem)
public class NtgaAlgorithm implements Algorithm{
    int populationSize;
    /**
     * the hyper-parameters below can by changed logically
     * */
    int epochs = 10;  // how many iteration will the algorithm run
    double initPackingRate = 0.05;
    int tournamentSize = 50;
    double orderCrossoverRate = 0.01;
    double uniformCrossoverRate = 0.01;
    double mutationRate = 0.05;

    // Initiate the number of solutions from the problem
    public NtgaAlgorithm(int numOfSolutions) {
        this.populationSize = numOfSolutions;
    }

    @Override
    public List<Solution> solve(TravelingThiefProblem problem) {
        // init population
        List<Solution> population = initPopulation(problem, populationSize, initPackingRate);

        // generation limitation
        for (int epoch = 0; epoch < epochs; ++epoch) {
            // init a new generation and individual index
            List<Solution> newGeneration = new ArrayList<>();
            int solutionIndex = 0;

            // non-dominated sorting
            nonDominatedSorting(population, false, true);

            while (newGeneration.size() < populationSize) {
                List<Solution> parents = new ArrayList<>();
                // select two individuals
                for (int n = 0; n < 2; ++n) {
                    // tournament selection
                    Solution parent = tournamentSelect(population, tournamentSize, populationSize);
                    parents.add(parent);
                }
                // order crossover (OX)
                List<Solution> offspring = orderCrossover(problem, parents, orderCrossoverRate, uniformCrossoverRate);
                // in-place mutation
                mutate(offspring, mutationRate, false);
                // in-place clone prevent - if a child is cloned from original population then mutate it
                clonePrevent(offspring, population, mutationRate, false);

                // improved here - offspring must not worse than parents to make sure population optimized
                // if you do not want, just comment the while
                while (parentsBetterThanOffspring(problem, offspring, parents)){
                    // generate offspring
                    offspring = orderCrossover(problem, parents, orderCrossoverRate, uniformCrossoverRate);
                    mutate(offspring, mutationRate, false);
                    clonePrevent(offspring, population, mutationRate, false);
                }

                // evaluate offspring and add into new generation
                for (Solution child : offspring) {
                    child = problem.evaluate(child.pi, child.z, true);  // update objectives
                    child.index = solutionIndex++;
                    newGeneration.add(child);
                }
            }
            // reset population
            population = new ArrayList<>(newGeneration);

            // show epoch number
            if (epochs > 10 && epoch % (epochs / 10) == 0) {
                System.out.println("epoch: " + epoch);

                // show objectives for each mini-epoch
//                for (Solution test : population){
//                    System.out.println(test.objectives);
//                }
            }
        }

        // show objective for the latest generation
//        for (Solution last : population){
//            System.out.println(last.objectives);
//        }

        return population;
    }

    /**
     * whether parents better than offspring
     */
    private boolean parentsBetterThanOffspring(TravelingThiefProblem problem, List<Solution> offspring, List<Solution> parents){
        List<Solution> tempPopulation = new ArrayList<>();
        int individualIndex = 0;
        for (Solution parent : parents){
            List<Integer> tempParentPI = new ArrayList<>(parent.pi);
            List<Boolean> tempParentZ = new ArrayList<>(parent.z);
            Solution tempParent = problem.evaluate(tempParentPI, tempParentZ, true);
            tempParent.index = individualIndex++;
            tempPopulation.add(tempParent);
        }
        for (Solution child : offspring){
            List<Integer> tempChildPI = new ArrayList<>(child.pi);
            List<Boolean> tempChildZ = new ArrayList<>(child.z);
            Solution tempChild = problem.evaluate(tempChildPI, tempChildZ, true);
            tempChild.index = individualIndex++;
            tempPopulation.add(tempChild);
        }
        List<Integer> tempRanks = nonDominatedSorting(tempPopulation, false, false);
        int minOffspringRank = Collections.min(tempRanks.subList(2, 4));
        int maxParentsRank = Collections.max(tempRanks.subList(0, 2));
        return maxParentsRank < minOffspringRank;
    }

    /**
     * clone prevention - check whether original population contains individual of new population
     * @param newPopulation the population to be checked
     * @param originalPopulation source population
     */
    private void clonePrevent(List<Solution> newPopulation, List<Solution> originalPopulation, double mutationRate, boolean showNotice){
        for (Solution child : newPopulation){
            while (isCloned(child, originalPopulation)){
                if (showNotice)
                    System.out.println("clone prevent execution");
                // place child into a list
                List<Solution> childInList = new ArrayList<>();
                childInList.add(child);
                mutate(childInList, mutationRate, false);  // in-place mutate the individual
            }
        }
    }

    /**
     * M-gene / Swap individual mutation (in-place)
     * @param IND either an individual or a population
     * @param mutationRate the probability of mutation
     * @param useSwapMutate whether use Swap Mutation for pi
     */
    private void mutate(List<Solution> IND, double mutationRate, boolean useSwapMutate){
        Random rand = new Random();

        int percentMutationRate = (int) (mutationRate * 100);  // mutation rate in hundred percent
        // M-gene Mutation Z
        for (Solution individual : IND) {
            for (int i = 0; i < individual.z.size(); ++i) {
                if (rand.nextInt(100) < percentMutationRate) {
                    if (individual.z.get(i)) {
                        individual.z.set(i, false);
                    } else {
                        individual.z.set(i, true);
                    }
                }
            }
        }

        // decide which type of mutation for pi
        if (useSwapMutate) {
            // Swap Mutation PI
            for (Solution individual : IND) {
                for (int i = 1; i < individual.pi.size(); ++i) {  // the first tour should not be swap mutated
                    if (rand.nextInt(100) < percentMutationRate) {
                        // Generate integers in the interval [1, size)
                        int swapPosition = rand.nextInt(individual.pi.size() - 1) + 1;
                        while (swapPosition == i)  // make sure the position to be swapped is different from current
                            swapPosition = rand.nextInt(individual.pi.size() - 1) + 1;
                        // swap
                        Collections.swap(individual.pi, i, swapPosition);
                    }
                }
            }
        }
        else {
            int reverseStart = rand.nextInt(IND.get(0).pi.size() - 1) + 1;  // random select a point from [1, piSize)
            int piMutateNumber = (int) (mutationRate * IND.get(0).pi.size());  // number of pi gene to be mutated
            int reverseEnd = reverseStart + piMutateNumber;  // reverse end point index
            if (reverseEnd > IND.get(0).pi.size()){  // index limitation
                reverseEnd = IND.get(0).pi.size();
            }
            for (Solution individual : IND){
                Collections.reverse(individual.pi.subList(reverseStart, reverseEnd));  // mutate execution
            }

        }
    }

    /**
     * OX order / uniform crossover
     * @param population a subset of population which only have two individuals
     * @param uniformCrossoverRate uniform crossover rate
     * @param orderCrossoverRate the percentage of parent gene not be reserved
     * @return generated offspring by order crossover operation
     */
    private List<Solution> orderCrossover(TravelingThiefProblem problem, List<Solution> population, double orderCrossoverRate, double uniformCrossoverRate){
        Random rand = new Random();
        // get two parents
        Solution parent1 = population.get(0);
        Solution parent2 = population.get(1);

        // pi size
        int size = parent1.pi.size();
        int sublistLength = (int)((1.00 - orderCrossoverRate) * size);  // sublist length
        // choose two random numbers for the start and end indices of the slice
        int start = rand.nextInt(size - sublistLength);
        int end = start + sublistLength;
//        System.out.println("slice: " + start + ", " + end);

        // add the sublist in between the start and end points
        List<Integer> sublist1 = new ArrayList<>(parent1.pi.subList(start, end));
        List<Integer> sublist2 = new ArrayList<>(parent2.pi.subList(start, end));
        // init children
        List<Integer> firstChildPI = new ArrayList<>();
        List<Integer> secondChildPI = new ArrayList<>();

        // iterate over each city in the parent tours
        for (int i = 0; i < size; ++i){
            // get the city at the current index in each of the two parent tours
            int currentCityInTour1 = parent1.pi.get(i);
            int currentCityInTour2 = parent2.pi.get(i);
            // if sublist1 does not already contain the current city in parent2, add it
            if (!sublist1.contains(currentCityInTour2))
                firstChildPI.add(currentCityInTour2);
            // if sublist2 does not already contain the current city in parent1, add it
            if (!sublist2.contains(currentCityInTour1))
                secondChildPI.add(currentCityInTour1);
        }
        // add sublist into child
        firstChildPI.addAll(start, sublist1);
        secondChildPI.addAll(start, sublist2);

        // crossover rate in hundred percent
        int percentUniformCrossoverRate = (int) (uniformCrossoverRate * 100);

        List<Boolean> firstChildZ = new ArrayList<>();
        List<Boolean> secondChildZ = new ArrayList<>();

        // perform uniform crossover for Z
        for (int sizeIndex = 0; sizeIndex < parent1.z.size(); ++sizeIndex){
            if (rand.nextInt(100) < percentUniformCrossoverRate){
                firstChildZ.add(parent2.z.get(sizeIndex));
                secondChildZ.add(parent1.z.get(sizeIndex));
            }
            else {
                firstChildZ.add(parent1.z.get(sizeIndex));
                secondChildZ.add(parent2.z.get(sizeIndex));
            }
        }

        // init two children
        Solution child1 = new Solution();
//        child1.pi = firstChildPI;
//        child1.z = firstChildZ;
        child1.pi = new ArrayList<>(firstChildPI);
        child1.z = new ArrayList<>(firstChildZ);
        Solution child2 = new Solution();
//        child2.pi = secondChildPI;
//        child2.z = secondChildZ;
        child2.pi = new ArrayList<>(secondChildPI);
        child2.z = new ArrayList<>(secondChildZ);
        // init children population
        List<Solution> children = new ArrayList<>();
        // add into children
        children.add(child1);
        children.add(child2);

        return children;
    }

    /**
     * tournament selection
     * comparison operator: I ≥r J if (Irank < Jrank) - the formula should be only used for NTGA
     * @param tournamentSize the number of individuals will be compare by comparison operator
     * @return the best individual
     * */
    private Solution tournamentSelect(List<Solution> population, int tournamentSize, int populationSize){
        Random rand = new Random();
        Solution best = population.get(rand.nextInt(populationSize));  // random select a individual
        for (int i = 1; i < tournamentSize; ++i){
            Solution individual = population.get(rand.nextInt(populationSize));  // random select a individual
            if (individual.rank < best.rank)  // comparison operator
                best = individual;
        }
        return best;
    }

    /**
     * fast non-dominated sorting with time complexity O(MN^2) which M is objectives and N is individual number
     * @param updateRank whether update the individual rank
     * */
    private List<Integer> nonDominatedSorting(List<Solution> population, boolean showInfo, boolean updateRank){
        // how many other individuals can dominate an individual - 多少个体能支配它
        List<Integer> dominated = new ArrayList<>(population.size());
        // list of individuals that an individual dominates - 支配哪些个体
        List<List<Integer>> dominates = new ArrayList<>(population.size());
        List<List<Integer>> paretoFront = new ArrayList<>();  // Pareto front

        List<Integer> ranks = new ArrayList<>();  // a list is used to store rank

        for (Solution s : population){
            if (updateRank) {
                s.rank = Integer.MAX_VALUE;  // reset rank for each individual
            }
            // init list
            dominated.add(0);
            dominates.add(new ArrayList<>());
        }
        paretoFront.add(new ArrayList<>());

        for (Solution p : population){
            for (Solution q : population){
                if (p == q)
                    continue;
                else if (p.getRelation(q) == 1){  // p dominates q
                    dominates.get(p.index).add(q.index);
                }
                else if (p.getRelation(q) == -1){  // p is dominated by q
                    dominated.set(p.index, dominated.get(p.index) + 1);
                }
            }
            if (dominated.get(p.index) == 0){  // pareto front
                if (updateRank) {
                    p.rank = 0;
                }
                else {
                    ranks.add(0);
                }
                paretoFront.get(0).add(p.index);
            }
        }

        if (showInfo) {
            System.out.println("dominated: " + dominated);
            System.out.println("dominates: " + dominates);
        }

        int i = 0;
        while (paretoFront.get(i).size() > 0){
            List<Integer> temp = new ArrayList<>();
            for (int m = 0; m < paretoFront.get(i).size(); ++m){
                for (int n = 0; n < dominates.get(paretoFront.get(i).get(m)).size(); ++n){
                    int otherIndividual = dominates.get(paretoFront.get(i).get(m)).get(n);
                    dominated.set(otherIndividual, dominated.get(otherIndividual)-1);
                    if (dominated.get(otherIndividual) == 0){
                        if (updateRank) {
                            population.get(otherIndividual).rank = i + 1;  // sorting
                        }
                        else {
                            ranks.add(i + 1);
                        }
                        temp.add(otherIndividual);
                    }
                }
            }
            ++i;
            paretoFront.add(temp);
        }
        if (showInfo)
            System.out.println("paretoFront: " + paretoFront);

        if (updateRank){
            return null;
        }
        else {
            return ranks;
        }
    }

    /**
     * Initialise population
     * @param initPackingRate initialised packing rate (Z)
     */
    private List<Solution> initPopulation(TravelingThiefProblem problem, int populationSize, double initPackingRate){
        Random rand = new Random();

        List<Solution> population = new ArrayList<>();  // init a population
        int individualIndex = 0;  // init individual index
        int packingRate = (int) (initPackingRate * 100);  // convert packing rate into hundred percent (int)
//        System.out.println("packing rate: " + packingRate);

        while(population.size() < populationSize){
            // random init tour (pi) scheme
            List<Integer> pi = getIndex(1, problem.numOfCities);
            Collections.shuffle(pi);
            pi.add(0,0);

            // init packing scheme
            List<Boolean> z = new ArrayList<>(problem.numOfItems);
            for (int i = 0; i < problem.numOfItems; ++i){
                if (rand.nextInt(100) < packingRate){
                    z.add(true);
                }
                else {
                    z.add(false);
                }
            }

            // evaluate the individual
            Solution s = problem.evaluate(pi, z, true);
            if (s != null) {
                s.index = individualIndex++;
                population.add(s);  // add the individual into population
            }
        }
        return population;
    }

    /**
     * the function is used to judge whether the individual is cloned from population
     */
    private boolean isCloned(Solution individual, List<Solution> population){
        for (Solution other : population){
            if (individual.equalsInDesignSpace(other))  // compare the genotype between individual and other
                return true;
        }
        return false;
    }

    private List<Integer> getIndex(int low, int high) {
        List<Integer> l = new ArrayList<>();
        for (int j = low; j < high; j++) {
            l.add(j);
        }
        return l;
    }
}
