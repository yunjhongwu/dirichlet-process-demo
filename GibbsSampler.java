/*
 * Created on Nov 18, 2014
 * JAVA SE 8
 * Auther: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu 
 */

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.Stack;

import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.random.RandomDataGenerator;

public abstract class GibbsSampler {
	final int n;
	final int maxNumClusters;
	final double alpha, theta, beta, xi;
	final ArrayList<Point2D> data;
	int[] labels;

	HashMap<Integer, double[]> clusters = new HashMap<Integer, double[]>();
	Stack<Integer> emptyClusters = new Stack<Integer>();
	ArrayList<Integer> ord = new ArrayList<Integer>();
	RandomDataGenerator sampler = new RandomDataGenerator();

	public GibbsSampler(double alpha, double theta, double beta, double xi,
			int initClusters, int maxNumClusters, ArrayList<Point2D> data) {
		this.data = data;
		this.n = data.size();
		this.labels = new int[n];
		this.alpha = alpha;
		this.theta = theta;
		this.beta = beta;
		this.xi = xi;
		this.maxNumClusters = maxNumClusters;
		initClusters(initClusters);
	}

	public void initClusters(int numClusters) {
		for (int i = 0; i < n; i++) {
			ord.add(i);
			if (i < numClusters)
				clusters.put(i, new double[5]);
			else
				emptyClusters.push(n - 1 - i);
		}
		Collections.shuffle(ord);
		ord.stream().forEach(i -> {
			clusters.get(i % numClusters)[0]++;
			labels[i] = i % numClusters;
		});
		clusters.put(-1, new double[5]);
		clusters.get(-1)[0] = -1;
		updateAllMoments();
	}

	private double posteriorVariance(double mu, double m2, double size) {
		return (2 * theta + size)
				/ (2 * beta + (m2 - Math.pow(mu, 2) / (size + xi)));
		// return sampler.nextGamma(theta + size / 2.0,
		// 2 / (2 * beta + (m2 - Math.pow(mu, 2) / (size + xi))));
	}

	private double posteriorMean(double mu, double size) {
		return mu / (size + xi);
		// return sampler.nextGaussian(mu / (size + xi), Math.sqrt((1 +
		// xi/ size) * sigma2 ));
	}

	private double logNormalLikelihood(int i, int j) {
		double s2x = posteriorVariance(clusters.get(j)[1], clusters.get(j)[3],
				clusters.get(j)[0]);
		double s2y = posteriorVariance(clusters.get(j)[2], clusters.get(j)[4],
				clusters.get(j)[0]);
		return -0.5
				* (Math.log(s2x)
						+ Math.log(s2y)
						+ Math.pow(
								data.get(i).x
										- posteriorMean(clusters.get(j)[1],
												clusters.get(j)[0]), 2) / s2x + Math
						.pow(data.get(i).y
								- posteriorMean(clusters.get(j)[2],
										clusters.get(j)[0]), 2)
						/ s2y);
	}

	protected double MHThreshold(int i, int j, int k) {
		if (clusters.get(j)[0] == 1)
			return 0;
		double t = Math.log(clusters.get(j)[0] - 1) + logNormalLikelihood(i, j);
		t -= Math.log((k == -1) ? alpha / clusters.size() : clusters.get(k)[0]);
		t -= logNormalLikelihood(i, k);
		return t;
	}

	protected int MHKernel() {
		int s = 0;
		int r = sampler.nextInt(0, clusters.size() - 1);
		for (Integer c : clusters.keySet())
			if (r == s++)
				return c;
		return -1;
	}

	protected void updateMoments(int i, int j, int k) {
		clusters.get(j)[1] += data.get(i).x;
		clusters.get(j)[2] += data.get(i).y;
		clusters.get(j)[3] += Math.pow(data.get(i).x, 2);
		clusters.get(j)[4] += Math.pow(data.get(i).y, 2);
		if (clusters.containsKey(k)) {
			clusters.get(k)[1] -= data.get(i).x;
			clusters.get(k)[2] -= data.get(i).y;
			clusters.get(k)[3] -= Math.pow(data.get(i).x, 2);
			clusters.get(k)[4] -= Math.pow(data.get(i).y, 2);
		}
	}

	protected void updateAllMoments() {
		updateAllMoments(clusters.keySet());
	}

	protected void updateAllMoments(Set<Integer> amended) {
		for (Integer c : amended) {
			clusters.get(c)[1] = 0;
			clusters.get(c)[2] = 0;
			clusters.get(c)[3] = 0;
			clusters.get(c)[4] = 0;
		}

		for (int i = 0; i < n; i++) {
			clusters.get(labels[i])[1] += data.get(i).x;
			clusters.get(labels[i])[2] += data.get(i).y;
			clusters.get(labels[i])[3] += Math.pow(data.get(i).x, 2);
			clusters.get(labels[i])[4] += Math.pow(data.get(i).y, 2);
		}
	}

	public double getResidual(int n, final double[] p, final double[] mux,
			final double[] muy) {
		double[] res = new double[(clusters.size() - 1) * mux.length];
		double[] q = new double[clusters.size() - 1];
		int i = 0;
		for (Integer c : clusters.keySet())
			if (c > -1) {
				for (int j = 0; j < mux.length; j++)
					res[i * mux.length + j] = Math.pow(
							mux[j]
									- posteriorMean(clusters.get(c)[1],
											clusters.get(c)[0]), 2)
							+ Math.pow(
									muy[j]
											- posteriorMean(clusters.get(c)[2],
													clusters.get(c)[0]), 2);
				q[i++] = clusters.get(c)[0] / (double) n;
			}

		return getOptResidual(res, q, p);
	}

	private double getOptResidual(final double[] res, final double[] cols,
			final double[] rows) {
		LinearObjectiveFunction f = new LinearObjectiveFunction(res, 0);
		Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		for (int i = 0; i < cols.length; i++) {
			double[] b = new double[cols.length * rows.length];
			for (int j = i * rows.length; j < (i + 1) * rows.length; j++)
				b[j]++;
			constraints.add(new LinearConstraint(b, Relationship.EQ, cols[i]));
		}

		for (int i = 0; i < rows.length; i++) {
			double[] b = new double[cols.length * rows.length];
			for (int j = i; j < b.length; j += rows.length)
				b[j]++;
			constraints.add(new LinearConstraint(b, Relationship.EQ, rows[i]));
		}
		SimplexSolver solver = new SimplexSolver();

		return solver.optimize(new MaxIter(10000), f,
				new LinearConstraintSet(constraints), GoalType.MINIMIZE,
				new NonNegativeConstraint(true)).getSecond();
	}

	public abstract void nextIter(int i);

	public abstract void next(int k);

}