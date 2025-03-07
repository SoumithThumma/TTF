package algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.stream.Collectors;

import model.Solution;
import model.TravelingThiefProblem;

// NSGA II Algorithm for team Lambda to solve the traveling thief problem
public class LambdaAlgorithm implements Algorithm {

	// Initial population size
	int populationSize;

	// Tournament size
	int tournamentSize = 2;

	// Mutation Factor
	double mutationFactor = 0.01;

	// Clone Prevention
	boolean clonePrevention = false;

	// Number of generations to run the algorithm for
	int totalGenerations = 10000;

	public List<Solution> entries = new LinkedList<>();

	// Initiate the number of solutions from the problem
	public LambdaAlgorithm(int numOfSolutions) {
		this.populationSize = numOfSolutions;
	}

	@Override
	public List<Solution> solve(TravelingThiefProblem problem) {

		Random rand = new Random();

		// Calculate the total weight of the bags in the problem
		Double totalWeight = 0.0;
		for (Double weight : problem.weight) {
			totalWeight += weight;
		}

		// Set the probability ratio out of 100
		int ratio = (int) ((problem.maxWeight / totalWeight) * 90);

		// List to hold initial population
		List<Solution> population = new ArrayList<>();

		// Init i for indices
		int i = 0;

		// Create initial population
		while (population.size() < populationSize) {
			List<Integer> pi;

			// Create a random permutation
			pi = getIndex(1, problem.numOfCities);
			Collections.shuffle(pi);
			pi.add(0, 0);

			// Create a packing plan
			List<Boolean> z = new ArrayList<>(problem.numOfItems);
			for (int j = 0; j < problem.numOfItems; j++) {

				if (rand.nextInt(100) <= ratio) {
					z.add(true);
				} else {
					z.add(false);
				}
			}

			// Evaluate for this random tour
			Solution s = problem.evaluate(pi, z, true);
			if (s != null) {
				s.index = i;
				i++;
				population.add(s);
			}

		}
		entries.addAll(population);

		calculatePopulationFitness(population);

		// Starting point for generations
		int gen = 0;

		// Comparator for Pareto front and crowding distance
		SolutionComparator solutionComparator = new SolutionComparator();

		// Start time
		long startTime = System.currentTimeMillis();

		// Run the algorithm for the required number of generations
		while (gen <= totalGenerations) {
			List<Solution> childPopulation = new ArrayList<>();

			while (childPopulation.size() < populationSize) {
				// Starting point of the NSGA II Algorithm
				tournamentSelection(population, problem, childPopulation);
			}

			population.addAll(childPopulation);

			// Reset indices
			for (int j = 0; j < population.size(); j++) {
				population.get(j).index = j;
			}
			// Calculate the fitness with the children in the population
			calculatePopulationFitness(population);
			// Sort population according to fitness
			Collections.sort(population, solutionComparator);

			// Get the desired population size after evaluating fitness
			population = population.subList(0, populationSize);

			// Reset indices
			for (int j = 0; j < population.size(); j++) {
				population.get(j).index = j;
			}

			// Reset entries
			entries = new LinkedList<>();
			entries.addAll(population);

			gen++;

			System.out.println(gen);

		}

		// End time
		long endTime = System.currentTimeMillis();

		System.out.println("That took " + (endTime - startTime) / 1000 + " seconds");

		return population;
	}

	private void tournamentSelection(List<Solution> population, TravelingThiefProblem problem,
			List<Solution> childPopulation) {
		// Initiate random class
		Random rand = new Random();
		// List to keep track of solutions that have been selected
		List<Integer> prevNumbers = new ArrayList<>();

		// Get a random number from population and make it parent A
		int randomNumber = rand.nextInt(populationSize);
		prevNumbers.add(randomNumber);

		// Set values for the parent A object
		Solution parentA = new Solution();
		parentA.index = randomNumber;
		parentA.rank = population.get(randomNumber).rank;
		parentA.crowdingDistance = population.get(randomNumber).crowdingDistance;
		// Initiate t = 1 and keep incrementing until t reaches the desired tournament
		// size
		int t = 1;

		// Run tournament through random selection and select the parent with the best
		// fitness
		while (t != tournamentSize) {
			// Pick a number that has not already been picked
			randomNumber = rand.nextInt(populationSize);
			while (prevNumbers.contains(randomNumber)) {
				randomNumber = randomNumber + 1;
				if (randomNumber == populationSize) {
					randomNumber = 0;
				}
			}
			prevNumbers.add(randomNumber);

			// Evaluate fitness of the new random selection to the previously selected
			// parent. Update values of parent A if needed
			if (!isFirstSolutionBetter(parentA, population.get(randomNumber))) {
				parentA.index = randomNumber;
				parentA.rank = population.get(randomNumber).rank;
				parentA.crowdingDistance = population.get(randomNumber).crowdingDistance;
			}
			t++;
		}

		// Reset prev Numbers list to hold only parent A while picking a parent B
		prevNumbers = new ArrayList<>();
		prevNumbers.add(parentA.index);

		// Get parent B
		randomNumber = rand.nextInt(populationSize);
		while (prevNumbers.contains(randomNumber)) {
			randomNumber = randomNumber + 1;
			if (randomNumber == populationSize) {
				randomNumber = 0;
			}
		}
		prevNumbers.add(randomNumber);
		// Set values of parent B
		Solution parentB = new Solution();
		parentB.index = randomNumber;
		parentB.rank = population.get(randomNumber).rank;
		parentB.crowdingDistance = population.get(randomNumber).crowdingDistance;

		t = 1;
		// Run tournament selection for parent B same as for parent A
		while (t != tournamentSize) {
			randomNumber = rand.nextInt(populationSize);

			while (prevNumbers.contains(randomNumber)) {
				randomNumber = randomNumber + 1;
				if (randomNumber == populationSize) {
					randomNumber = 0;
				}
			}

			prevNumbers.add(randomNumber);

			if (!isFirstSolutionBetter(parentB, population.get(randomNumber))) {
				parentB.index = randomNumber;
				parentB.rank = population.get(randomNumber).rank;
				parentB.crowdingDistance = population.get(randomNumber).crowdingDistance;
			}
			t++;
		}
		// After tournament selection is done run a crossover
		crossover(parentA, parentB, population, problem, childPopulation);
	}

	// Order crossover function
	private void odercrossover(Solution parentA, Solution parentB, List<Solution> population,
			TravelingThiefProblem problem, List<Solution> childPopulation) {
		Random rand = new Random();

		Solution parentAfromPop = population.get(parentA.index);

		// Lists to hold paths of crossover

		List<Integer> childCPi = new ArrayList<>();
		List<Integer> childDPi = new ArrayList<>();

		Solution parentBfromPop = population.get(parentB.index);

		// Choosing a random crossover point

		int size = parentAfromPop.pi.size();
		int sublistLength = (int) ((1.00 - 0.01) * size); // sublist length
		// choose two random numbers for the start and end indices of the slice
		int start = rand.nextInt(size - sublistLength);
		int end = start + sublistLength;

		// add the sublist in between the start and end points
		List<Integer> sublist1 = new ArrayList<>(parentAfromPop.pi.subList(start, end));
		List<Integer> sublist2 = new ArrayList<>(parentBfromPop.pi.subList(start, end));

		// iterate over each city in the parent tours
		for (int i = 0; i < size; ++i) {
			// get the city at the current index in each of the two parent tours
			int currentCityInTour1 = parentAfromPop.pi.get(i);
			int currentCityInTour2 = parentBfromPop.pi.get(i);
			// if sublist1 does not already contain the current city in parent2, add it
			if (!sublist1.contains(currentCityInTour2))
				childCPi.add(currentCityInTour2);
			// if sublist2 does not already contain the current city in parent1, add it
			if (!sublist2.contains(currentCityInTour1))
				childDPi.add(currentCityInTour1);
		}
		// add sublist into child
		childCPi.addAll(start, sublist1);
		childDPi.addAll(start, sublist2);

		// Uniform crossover between child C and D for z
		List<Boolean> childCz = new ArrayList<>();
		List<Boolean> childDz = new ArrayList<>();
		for (int i = 0; i < parentAfromPop.z.size(); i++) {
			if (rand.nextInt(2) == 0) {
				childCz.add(parentAfromPop.z.get(i));
				childDz.add(parentBfromPop.z.get(i));
			} else {
				childCz.add(parentBfromPop.z.get(i));
				childDz.add(parentAfromPop.z.get(i));
			}
		}

		// Start the mutation
		this.mutate(childCPi, childDPi, childCz, childDz, population, childPopulation, problem);

	}

	private void crossover(Solution parentA, Solution parentB, List<Solution> population, TravelingThiefProblem problem,
			List<Solution> childPopulation) {
		Random rand = new Random();

		// Choosing a random crossover point
		Solution parentAfromPop = population.get(parentA.index);
		int piSize = parentAfromPop.pi.size();
		int crossOverPoint = rand.nextInt(piSize);

		// Lists to hold paths of crossover

		List<Integer> parentARightPath = new ArrayList<>();

		List<Integer> parentBRightPath = new ArrayList<>();

		List<Integer> childCPi = new ArrayList<>();
		List<Integer> childDPi = new ArrayList<>();

		Solution parentBfromPop = population.get(parentB.index);
		// Add the path to the right of the cross over point
		for (int i = crossOverPoint; i < piSize; i++) {
			parentARightPath.add(parentAfromPop.pi.get(i));
			parentBRightPath.add(parentBfromPop.pi.get(i));
		}

		// Add points to child C from parent B if they are not in the right side of
		// Parent A
		for (int i = 0; i < piSize; i++) {

			Integer parentBpi = parentBfromPop.pi.get(i);
			if (!parentARightPath.contains(parentBpi)) {
				childCPi.add(parentBpi);
			}

		}
		// Add all points in the right path of parent A to child C
		childCPi.addAll(parentARightPath);

		// Repeat the same as above from child D from parent As
		for (int i = 0; i < piSize; i++) {

			Integer parentApi = parentAfromPop.pi.get(i);
			if (!parentBRightPath.contains(parentApi)) {
				childDPi.add(parentApi);
			}

		}
		childDPi.addAll(parentBRightPath);

		// Uniform crossover between child C and D for z
		List<Boolean> childCz = new ArrayList<>();
		List<Boolean> childDz = new ArrayList<>();
		for (int i = 0; i < parentAfromPop.z.size(); i++) {
			if (rand.nextInt(2) == 0) {
				childCz.add(parentAfromPop.z.get(i));
				childDz.add(parentBfromPop.z.get(i));
			} else {
				childCz.add(parentBfromPop.z.get(i));
				childDz.add(parentAfromPop.z.get(i));
			}
		}

		// Start the mutation
		this.mutate(childCPi, childDPi, childCz, childDz, population, childPopulation, problem);

	}

	private void mutate(List<Integer> childCPi, List<Integer> childDPi, List<Boolean> childCz, List<Boolean> childDz,
			List<Solution> population, List<Solution> childPopulation, TravelingThiefProblem problem) {

		Random rand = new Random();

		// Run mutation for path
		int piSize = childCPi.size();
		int reversePoint = rand.nextInt(piSize);
		if (reversePoint == 0) {
			reversePoint = 1;
		}
		int numberOfPathMutations = (int) (mutationFactor * piSize);
		int reverseEndPoint = numberOfPathMutations + reversePoint;
		if (reverseEndPoint > piSize) {
			reverseEndPoint = piSize - 1;
		}
		Collections.reverse(childCPi.subList(reversePoint, reverseEndPoint));

		int reversePointforD = rand.nextInt(piSize);
		if (reversePointforD == 0) {
			reversePointforD = 1;
		}
		int reverseEndPointforD = numberOfPathMutations + reversePointforD;
		if (reverseEndPointforD > piSize) {
			reverseEndPointforD = piSize - 1;
		}
		Collections.reverse(childDPi.subList(reversePointforD, reverseEndPointforD));

		// Lists to hold previous mutating points to avoid repetition
		List<Integer> prevCmutations = new ArrayList<>();
		List<Integer> prevDmutations = new ArrayList<>();

		// Run the desired number of mutations
		int mutation = 0;

		int numberOfMutations = (int) (mutationFactor * childCz.size());
		while (mutation < numberOfMutations) {

			int indexToMutateforC = rand.nextInt(childCz.size());
			// If previous mutations does not contain this point flip the solution
			while (prevCmutations.contains(indexToMutateforC)) {
				indexToMutateforC = rand.nextInt(childCz.size());
			}
			prevCmutations.add(indexToMutateforC);

			if (childCz.get(indexToMutateforC)) {
				childCz.set(indexToMutateforC, false);
			} else {
				childCz.set(indexToMutateforC, true);
			}

			int indexToMutateforD = rand.nextInt(childDz.size());
			// If previous mutations does not contain this point flip the solution
			while (prevDmutations.contains(indexToMutateforD)) {
				indexToMutateforD = rand.nextInt(childDz.size());
			}
			prevDmutations.add(indexToMutateforD);
			if (childDz.get(indexToMutateforD)) {
				childDz.set(indexToMutateforD, false);
			} else {
				childDz.set(indexToMutateforD, true);
			}

			mutation++;
		}

		// Evaluate the new children with mutations
		Solution childC = problem.evaluate(childCPi, childCz, true);
		Solution childD = problem.evaluate(childDPi, childDz, true);
		if (!clonePrevention) {
			childPopulation.add(childC);
			childPopulation.add(childD);
			entries.add(childC);
			entries.add(childD);
		} else {
			if (checkForClones(childC)) {
				childPopulation.add(childC);
				entries.add(childC);
			}

			if (checkForClones(childD)) {
				childPopulation.add(childD);
				entries.add(childD);
			}
		}

	}

	private boolean checkForClones(Solution s) {

		boolean isAdded = true;

		for (Iterator<Solution> it = entries.iterator(); it.hasNext();) {
			Solution other = it.next();

			// If equal in design space
			if (s.equalsInDesignSpace(other)) {
				System.out.println("Prevented");
				isAdded = false;
				break;
			}
		}

		return isAdded;

	}

	// Evaluate which one of the solutions s1 or s2 is better using rank and
	// crowding distance
	private boolean isFirstSolutionBetter(Solution s1, Solution s2) {
		if (s1.rank < s2.rank) {
			return true;
		} else if (s1.rank == s2.rank && s1.crowdingDistance > s2.crowdingDistance) {
			return true;
		}

		return false;
	}

	// Determine the fitness of all the solutions in the population
	private void calculatePopulationFitness(List<Solution> population) {

		// Reset rank and crowding distance
		for (Solution s : population) {
			s.rank = -1;
			s.crowdingDistance = -1;
		}

		// Normalize values
		normalizeFitnessValues(population);
		List<Solution> remainingToBeRanked = new ArrayList<>(population);

		// Set initial rank
		int rank = 1;

		while (remainingToBeRanked.size() > 0) {

			List<Integer> solutionsIndexInRank = new ArrayList<>();

			// If a solution is non dominated set the current rank
			for (int i = 0; i < remainingToBeRanked.size(); i++) {
				if (isNotDominated(remainingToBeRanked, remainingToBeRanked.get(i))) {
					remainingToBeRanked.get(i).rank = rank;
					solutionsIndexInRank.add(remainingToBeRanked.get(i).index);
				}
			}

			Iterator<Solution> i = remainingToBeRanked.iterator();

			// Remove solutions that were ranked
			while (i.hasNext()) {
				Solution sol = i.next();
				if (sol.time == Double.MAX_VALUE) {
					sol.rank = Integer.MAX_VALUE;
					i.remove();
				} else if (solutionsIndexInRank.contains(sol.index)) {
					i.remove();
				}

			}
			rank++;
		}

		// Map to hold solutions grouped by ranks
		Map<Integer, List<Solution>> mapByRanks = population.stream().collect(Collectors.groupingBy(Solution::getRank));

		// Calculate crowding distance of solutions in every front
		for (Entry<Integer, List<Solution>> entry : mapByRanks.entrySet()) {
			if (entry.getKey() != Integer.MAX_VALUE) {
				calculateCrowdingDistance(entry.getValue());
			}

		}

	}

	// Method to calculate crowding distance
	private void calculateCrowdingDistance(List<Solution> solutions) {
		solutions.sort(Comparator.comparing(a -> a.normalizedTime));

		for (int i = 0; i < solutions.size(); i++) {
			if (i == 0 || i == solutions.size() - 1) {
				solutions.get(i).crowdingDistance = Double.MAX_VALUE;
			} else {
				Solution current = solutions.get(i);
				Solution left = solutions.get(i - 1);
				Solution right = solutions.get(i + 1);

				double distanceLeft = euclideanDistance(current, left);
				double distanceRight = euclideanDistance(current, right);
				solutions.get(i).crowdingDistance = distanceLeft + distanceRight;

			}
		}

	}

	// Method to calculate euclidean distance
	public double euclideanDistance(Solution a, Solution b) {
		return Math.sqrt(Math.pow(a.normalizedTime - b.normalizedTime, 2)
				+ Math.pow(a.normalizedProfit - b.normalizedProfit, 2));
	}

	// Check whether a solution is dominated or not in the population to assign rank
	private boolean isNotDominated(List<Solution> remainingToBeRanked, Solution solution) {
		for (Solution s : remainingToBeRanked) {
			if (s == solution) {
				continue;
			}

			if ((s.objectives.get(0) <= solution.objectives.get(0) && s.objectives.get(1) <= solution.objectives.get(1))
					&& (s.objectives.get(0) < solution.objectives.get(0)
							|| s.objectives.get(1) < solution.objectives.get(1))) {
				return false;
			}

		}
		return true;
	}

	// Normalize time and profit when calculating fitness
	private void normalizeFitnessValues(List<Solution> population) {
		Double maxTime = 0.0, maxProfit = 0.0;
		for (Solution s : population) {
			if (s.time != Double.MAX_VALUE && s.time > maxTime) {
				maxTime = s.time;
			}
			if (s.profit != Double.MAX_VALUE && s.profit > maxProfit) {
				maxProfit = s.profit;
			}
		}
		for (Solution s : population) {
			s.normalizedTime = s.time / maxTime;
			s.normalizedProfit = s.profit / maxProfit;
		}

	}

	private List<Integer> getIndex(int low, int high) {
		List<Integer> l = new ArrayList<>();
		for (int j = low; j < high; j++) {
			l.add(j);
		}
		return l;
	}

}
