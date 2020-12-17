package algorithms;

import java.util.Comparator;

import model.Solution;

public class SolutionComparator implements Comparator<Solution> {

	@Override
	public int compare(Solution o1, Solution o2) {

		int c = Integer.compare(o1.rank, o2.rank);
		if (c == 0) {
			return Double.compare(o1.crowdingDistance, o2.crowdingDistance);
		}

		return c;
	}

}
