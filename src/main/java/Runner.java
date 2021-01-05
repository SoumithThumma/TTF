import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import algorithms.Algorithm;
import algorithms.LambdaAlgorithm;
import model.Solution;
import model.TravelingThiefProblem;

class Runner {

	static final ClassLoader LOADER = Runner.class.getClassLoader();

	public static void main(String[] args) throws IOException {

		// List<String> instanceToRun = Arrays.asList("a280-n279");
		List<String> instanceToRun = Competition.INSTANCES;

		for (String instance : instanceToRun) {

			// readProblem the problem from the file
			String fname = String.format("resources/%s.txt", instance);
			InputStream is = LOADER.getResourceAsStream(fname);

			TravelingThiefProblem problem = Util.readProblem(is);
			problem.name = instance;

			// number of solutions that will be finally necessary for submission - not used
			// here
			int numOfSolutions = Competition.numberOfSolutions(problem);
			System.out.println(problem.name + " ");
			// initialize your algorithm
			Algorithm algorithm = new LambdaAlgorithm(numOfSolutions);
			// Algorithm algorithm = new ExhaustiveSearch();

			// use it to to solve the problem and return the non-dominated set
			List<Solution> nds = algorithm.solve(problem);

			double leastTime = nds.get(0).time;
			double profit = nds.get(0).profit;
			double bestProfit = nds.get(0).profit;
			double time = nds.get(0).time;

			// sort by time and printSolutions it
			// nds.sort(Comparator.comparing(a -> a.time));

			for (Solution s : nds) {
				if (s.profit > bestProfit) {
					bestProfit = s.profit;
					time = s.time;
				}
				if (s.time < leastTime) {
					leastTime = s.time;
					profit = s.profit;
				}
			}

			System.out.println("Least Time: " + leastTime + " , Proft: " + profit);
			System.out.println("Best Profit: " + bestProfit + " Time: " + time);

//			Util.printSolutions(nds, true);

			File dir = new File("results");
			if (!dir.exists())
				dir.mkdirs();
			Util.writeSolutions("results", Competition.TEAM_NAME, problem, nds);

		}

	}

}