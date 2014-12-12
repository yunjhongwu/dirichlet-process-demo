/*
 * Created on Nov 12 2014
 * JAVA SE 8
 * Author: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.ImageIO;

public class DPSimulator {
	@SuppressWarnings("unused")
	public static void main(String[] args) throws InterruptedException,
			IOException {
		/* Parameters */
		final int n = 50000;
		final int maxIters = Integer.MAX_VALUE;
		final double alpha0 = 1; // Alpha for the generative model
		final double alpha = 1; // Alpha for the Gibbs sampler
		final double theta = 20; // Shape parameter of the gamma distribution
		final double beta = 20; // Scale parameter of the gamma distribution
		final double gamma = 20; // Separability of clusters
		final int maxNumClusters = n;
		int initClusters = 20;

		final boolean singleton = true;
		final int seed = 1;
		final int visual = 1;
		final int eval = 0;
		final boolean fout = false;
		final int saveChart = Integer.MAX_VALUE;

		/* Generating data */
		System.out.print("Generating data...");
		final ArrayList<Point2D> data = new ArrayList<Point2D>();
		int[] labels = new int[n];
		CRP crp = new CRP(alpha0, theta, beta, gamma, seed);
		for (int i = 0; i < n; i++)
			data.add(crp.next());

		Collections.sort(data);
		for (int i = 0; i < n; i++)
			labels[i] = data.get(i).cluster;

		/* Getting parameters */
		double[] proportion = new double[crp.size.size()];
		double[] centroidx = new double[crp.size.size()];
		double[] centroidy = new double[crp.size.size()];
		for (int c = 0; c < proportion.length; c++) {
			proportion[c] = crp.size.get(c) / (double) n;
			centroidx[c] = crp.moments.get(c)[0];
			centroidy[c] = crp.moments.get(c)[1];
		}
		ArrayList<Integer> freq = crp.size;
		Collections.sort(freq, Collections.reverseOrder());
		System.out.println("done.");
		System.out.println(n + " data points and " + crp.size.size()
				+ " clusters with size " + freq.toString() + " generated.");

		crp = null;
		initClusters = (initClusters < 1) ? freq.size() : initClusters;

		/* Simulation */
		GibbsSampler gibbs = (singleton) ? new SingletonGibbsSampler(alpha,
				theta, beta, gamma, initClusters, maxNumClusters, data)
				: new BlockGibbsSampler(alpha, theta, beta, gamma,
						initClusters, maxNumClusters, data);
		gibbs.setTrueCentroids(proportion, centroidx, centroidy);

		ScatterPlot currentPlot = null;
		EvaluationPlot numsPlot = null;
		EvaluationPlot residualPlot = null;
		if (visual > 0) {
			new DistributionPlot(freq, proportion.length);
			numsPlot = new EvaluationPlot(initClusters, "Number of clusters");
			currentPlot = ScatterPlot.initPlots(data, labels, gibbs.labels,
					freq.size(), singleton);

			if (eval > 0)
				residualPlot = new EvaluationPlot(0, "Wasserstein distance");
		}

		final String filename = "num_of_clusters_" + n + "_" + initClusters
				+ "_alpha" + alpha + "_beta" + beta + "_theta" + theta
				+ "_gamma" + gamma + "_"
				+ ((singleton) ? "singleton" : "block");
		PrintWriter file = null;
		PrintWriter fileWD = null;
		if (fout) {
			file = new PrintWriter(new BufferedWriter(new FileWriter(new File(
					filename + ".txt"), true)));
			fileWD = new PrintWriter(new BufferedWriter(new FileWriter(
					new File(filename + "_WD.txt"), true)));
		}

		long startTime = System.nanoTime();
		for (int k = 1; k <= maxIters; k++) {
			gibbs.next(k);
			String output = String
					.format("Iteration %d; %f milliseconds per iteration; %d clusters. ",
							k, (System.nanoTime() - startTime) / 1000000.0 / k,
							gibbs.clusters.size() - 1);
			System.out.println(output);

			if (fout)
				file.println(output);
			if (visual > 0 && k % visual == 0) {
				currentPlot.updateColors(gibbs.labels);
				numsPlot.updateSeries(k, gibbs.clusters.size() - 1);
				if (eval > 0 && k % eval == 0)
					residualPlot.updateSeries(k, gibbs.getResidual());

			}
			if (fout && eval > 0 && k % eval == 0)
				fileWD.println(gibbs.getResidual());

			if (k == saveChart && visual > 0)
				ImageIO.write(numsPlot.chart.createBufferedImage(650, 400),
						"png", new File(filename + ".png"));

		}

		if (fout) {
			file.close();
			fileWD.close();
		}
	}
}
