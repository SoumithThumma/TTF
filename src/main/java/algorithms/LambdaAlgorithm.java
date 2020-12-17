package algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.stream.Collectors;

import model.Solution;
import model.TravelingThiefProblem;

public class LambdaAlgorithm implements Algorithm {

	private int numOfSolutions;

	public LambdaAlgorithm(int numOfSolutions) {
		this.numOfSolutions = numOfSolutions;
	}

	@Override
	public List<Solution> solve(TravelingThiefProblem problem) {

		Random rand = new Random();
		Double totalWeight = 0.0;
		for (Double weight : problem.weight) {
			totalWeight += weight;

		}

		int ratio = (int) ((problem.maxWeight / totalWeight) * 90);

		// List to hold initial population
		List<Solution> population = new ArrayList<>();

		// Initial population size
		int populationSize = numOfSolutions;

		// Tournament size
		int tournamentSize = 4;

		int mutationFactor = 3;

		int totalGenerations = 100000;

		int i = 0;
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

			// evaluate for this random tour
			Solution s = problem.evaluate(pi, z, true);
			if (s != null) {
				s.index = i;
				i++;
				population.add(s);
			}

		}
		calculatePopulationFitness(population);

		int gen = 0;

		SolutionComparator solutionComparator = new SolutionComparator();

		while (gen <= totalGenerations) {
			List<Solution> childPopulation = new ArrayList<>();

			while (childPopulation.size() < populationSize) {
				tournamentSelection(tournamentSize, population, populationSize, problem, mutationFactor,
						childPopulation);
			}
			population.addAll(childPopulation);

			for (int j = 0; j < population.size(); j++) {
				population.get(j).index = j;
			}
			calculatePopulationFitness(population);
			Collections.sort(population, solutionComparator);

			population = population.subList(0, populationSize);

			for (int j = 0; j < population.size(); j++) {
				population.get(j).index = j;
			}
			gen++;
			System.out.println(gen);
		}

		return population;
	}

	private void tournamentSelection(int tournamentSize, List<Solution> population, int populationSize,
			TravelingThiefProblem problem, int mutationFactor, List<Solution> childPopulation) {
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
		// Initiate t =1 and keep incrementing until t reaches the desired tournament
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
		crossover(parentA, parentB, population, problem, mutationFactor, childPopulation);
	}

	private void crossover(Solution parentA, Solution parentB, List<Solution> population, TravelingThiefProblem problem,
			int mutationFactor, List<Solution> childPopulation) {
		Random rand = new Random();

		// Choosing a random crossover point
		Solution parentAfromPop = population.get(parentA.index);
		int piSize = parentAfromPop.pi.size();
		int crossOverPoint = rand.nextInt(piSize);

		List<Integer> parentARightPath = new ArrayList<>();

		List<Integer> parentBRightPath = new ArrayList<>();

		List<Integer> childCPi = new ArrayList<>();
		List<Integer> childDPi = new ArrayList<>();

		Solution parentBfromPop = population.get(parentB.index);
		for (int i = crossOverPoint; i < piSize; i++) {
			parentARightPath.add(parentAfromPop.pi.get(i));
			parentBRightPath.add(parentBfromPop.pi.get(i));
		}

		for (int i = 0; i < piSize; i++) {

			Integer parentBpi = parentBfromPop.pi.get(i);
			if (!parentARightPath.contains(parentBpi)) {
				childCPi.add(parentBpi);
			}

		}
		childCPi.addAll(parentARightPath);

		for (int i = 0; i < piSize; i++) {

			Integer parentApi = parentAfromPop.pi.get(i);
			if (!parentBRightPath.contains(parentApi)) {
				childDPi.add(parentApi);
			}

		}
		childDPi.addAll(parentBRightPath);

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

		this.mutate(childCPi, childDPi, childCz, childDz, population, mutationFactor, childPopulation, problem);

	}

	private void mutate(List<Integer> childCPi, List<Integer> childDPi, List<Boolean> childCz, List<Boolean> childDz,
			List<Solution> population, int mutationFactor, List<Solution> childPopulation,
			TravelingThiefProblem problem) {

		Random rand = new Random();

		// Lists to hold previous mutating points to avoid repetition
		List<Integer> prevCmutations = new ArrayList<>();
		List<Integer> prevDmutations = new ArrayList<>();

		// Run the desired number of mutations
		int mutation = 0;

		while (mutation < mutationFactor) {

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

		Solution childC = problem.evaluate(childCPi, childCz, true);
		Solution childD = problem.evaluate(childDPi, childDz, true);

		childPopulation.add(childC);
		childPopulation.add(childD);

	}

	private boolean isFirstSolutionBetter(Solution s1, Solution s2) {
		if (s1.rank < s2.rank) {
			return true;
		} else if (s1.rank == s2.rank && s1.crowdingDistance > s2.crowdingDistance) {
			return true;
		}

		return false;
	}

	private void calculatePopulationFitness(List<Solution> population) {
		for (Solution s : population) {
			s.rank = -1;
			s.crowdingDistance = -1;
		}

		normalizeFitnessValues(population);
		List<Solution> remainingToBeRanked = new ArrayList<>(population);

		int rank = 1;

		while (remainingToBeRanked.size() > 0) {

			List<Integer> solutionsIndexInRank = new ArrayList<>();

			for (int i = 0; i < remainingToBeRanked.size(); i++) {
				if (isNotDominated(remainingToBeRanked, remainingToBeRanked.get(i))) {
					remainingToBeRanked.get(i).rank = rank;
					solutionsIndexInRank.add(remainingToBeRanked.get(i).index);
				}
			}

			Iterator<Solution> i = remainingToBeRanked.iterator();

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

		Map<Integer, List<Solution>> mapByRanks = population.stream().collect(Collectors.groupingBy(Solution::getRank));

		for (Entry<Integer, List<Solution>> entry : mapByRanks.entrySet()) {
			if (entry.getKey() != Integer.MAX_VALUE) {
				calculateCrowdingDistance(entry.getValue());
			}

		}

	}

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

	public double euclideanDistance(Solution a, Solution b) {
		return Math.sqrt(Math.pow(a.normalizedTime - b.normalizedTime, 2)
				+ Math.pow(a.normalizedProfit - b.normalizedProfit, 2));
	}

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
