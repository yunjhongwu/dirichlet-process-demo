/*
 * Created on Nov 18, 2014
 * JAVA SE 8
 * Auther: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu 
 */

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Collectors;

public class BlockGibbsSampler extends GibbsSampler {
	HashMap<Integer, Integer> amended;

	public BlockGibbsSampler(double alpha, double theta, double beta,
			double xi, int initClusters, int maxNumClusters,
			ArrayList<Point2D> data) {
		super(alpha, theta, beta, xi, initClusters, maxNumClusters, data);
	}

	public void nextIter(int i) {
		int nextCluster = MHKernel();
		if (nextCluster != -1 || clusters.size() <= maxNumClusters)
			if (nextCluster != labels[i]
					&& sampler.nextExponential(1) > MHThreshold(i, labels[i],
							nextCluster)) {
				nextCluster = (nextCluster == -1) ? emptyClusters.peek()
						: nextCluster;

				if (!clusters.containsKey(nextCluster))
					clusters.put(nextCluster, new double[5]);

				if (!amended.containsKey(labels[i]))
					amended.put(labels[i],
							(int) (clusters.get(labels[i])[0] - 1));
				else
					amended.put(labels[i], amended.get(labels[i]) - 1);
				if (!amended.containsKey(nextCluster))
					amended.put(nextCluster,
							(int) (clusters.get(nextCluster)[0] + 1));
				else
					amended.put(nextCluster, amended.get(nextCluster) + 1);

				labels[i] = nextCluster;
			}
	}

	public void next(int k) {
		amended = new HashMap<Integer, Integer>();

	    ord.forEach(i -> nextIter(i));

		for (Integer c : amended.keySet())
			if (amended.get(c) > 0)
				clusters.get(c)[0] = amended.get(c);
			else {
				clusters.remove(c);
				emptyClusters.push(c);
			}

		if (clusters.containsKey(emptyClusters.peek()))
			emptyClusters.pop();
		updateAllMoments(amended.keySet().stream()
				.filter(i -> amended.get(i) > 0).collect(Collectors.toSet()));

	}
}