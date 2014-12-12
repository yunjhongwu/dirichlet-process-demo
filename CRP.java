/*
 * Created on Nov 18 2014
 * JAVA SE 8
 * Author: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu
 */

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomDataGenerator;

public class CRP {
	// X \sim N(\mu_c, \sigma^2)
	// \mu_c \sim N(0, \gamma\sigma^2)
	// \sigma^2 \sim Gamma(\theta, \beta)

	final double alpha, theta, eta, gamma;
	int n = 0;

	RandomDataGenerator sampler = new RandomDataGenerator();
	ArrayList<Integer> size = new ArrayList<Integer>();
	ArrayList<double[]> moments = new ArrayList<double[]>();

	public CRP(double alpha, double theta, double beta, double gamma, int seed) {
		sampler.reSeed(seed);
		this.alpha = alpha;
		this.theta = theta;
		this.eta = 1 / beta;
		this.gamma = gamma;
		initCentroids(10000);
	}

	public CRP(double alpha, double theta, double beta, double gamma) {
		this(alpha, theta, beta, gamma, (new RandomDataGenerator()).nextInt(0,
				Integer.MAX_VALUE));

	}

	private void initCentroids(int K) {
		for (int k = 0; k < K; k++) {
			double[] m = new double[4];
			m[2] = sampler.nextGamma(theta, eta);
			m[3] = sampler.nextGamma(theta, eta);
			m[0] = sampler.nextGaussian(0, Math.sqrt(gamma * m[2]));
			m[1] = sampler.nextGaussian(0, Math.sqrt(gamma * m[3]));
			moments.add(m);
		}
	}

	private Point2D getPosition(int cluster) {
		return new Point2D(cluster, (float) sampler.nextGaussian(
				moments.get(cluster)[0], moments.get(cluster)[2]),
				(float) sampler.nextGaussian(moments.get(cluster)[1],
						moments.get(cluster)[3]));
	}

	public Point2D next() {
		int cluster = -1;
		double sample = sampler.nextUniform(0, 1) * (n + alpha);
		if (sample < n) {
			while (sample > 0)
				sample -= size.get(++cluster);
			size.set(cluster, size.get(cluster) + 1);
		} else {
			size.add(1);
			cluster++;
		}
		n++;

		return getPosition(cluster);
	}
}