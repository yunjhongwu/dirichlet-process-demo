import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class VectorGibbsSampler extends GibbsSampler {
	HashMap<Integer, AtomicInteger> modifiedClusters;

	public VectorGibbsSampler(double alpha, double theta, double beta,
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

				if (!modifiedClusters.containsKey(labels[i]))
					modifiedClusters.put(labels[i], new AtomicInteger(
							(int) (clusters.get(labels[i])[0] - 1)));
				else
					modifiedClusters.get(labels[i]).decrementAndGet();
				if (!modifiedClusters.containsKey(nextCluster))
					modifiedClusters.put(nextCluster, new AtomicInteger(
							(int) (clusters.get(nextCluster)[0] + 1)));
				else
					modifiedClusters.get(nextCluster).incrementAndGet();

				labels[i] = nextCluster;
			}
	}

	public void next(int k) {
		modifiedClusters = new HashMap<Integer, AtomicInteger>();

		Collections.shuffle(ord);
		ord.forEach(i -> nextIter(i));

		for (Integer c : modifiedClusters.keySet())
			if (modifiedClusters.get(c).intValue() > 0)
				clusters.get(c)[0] = modifiedClusters.get(c).doubleValue();
			else {
				clusters.remove(c);
				emptyClusters.push(c);
			}

		if (clusters.containsKey(emptyClusters.peek()))
			emptyClusters.pop();
		updateAllMoments(modifiedClusters.keySet().stream()
				.filter(i -> modifiedClusters.get(i).intValue() > 0)
				.collect(Collectors.toSet()));

	}
}