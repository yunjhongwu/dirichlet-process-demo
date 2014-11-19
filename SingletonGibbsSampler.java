import java.util.ArrayList;
import java.util.Collections;

public class SingletonGibbsSampler extends GibbsSampler {
	public SingletonGibbsSampler(double alpha, double theta, double beta,
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
				nextCluster = (nextCluster == -1) ? emptyClusters.pop()
						: nextCluster;

				if (!clusters.containsKey(nextCluster))
					clusters.put(nextCluster, new double[5]);

				clusters.get(nextCluster)[0]++;
				if (clusters.get(labels[i])[0] > 1)
					clusters.get(labels[i])[0]--;
				else {
					clusters.remove(labels[i]);
					emptyClusters.push(labels[i]);
				}
				updateMoments(i, nextCluster, labels[i]);
				labels[i] = nextCluster;
			}
	}

	public void next(int k) {
		Collections.shuffle(ord);
		ord.forEach(i -> nextIter(i));

		if (k % 1000 == 0)
			updateAllMoments();
	}

}
