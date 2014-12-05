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
	public static void main(String[] args) throws InterruptedException,
			IOException {
		final int n = 10000;
		final int maxIters = Integer.MAX_VALUE;
		final double alpha = 1;
		final double theta = 100;
		final double beta = 20;
		final double xi = 20;
		final int maxNumClusters = n;
		final boolean singleton = true;
		int initClusters = 0;

		final int seed = 3;
		final int visual = 1;
		final int eval = 10;
		final boolean fout = true;
		final int saveChart = 3000;

		/* Generating data */
		System.out.print("Generating data...");
		final ArrayList<Point2D> data = new ArrayList<Point2D>();
		int[] labels = new int[n];
		CRP crp = new CRP(alpha, theta, beta, xi, seed);
		for (int i = 0; i < n; i++)
			data.add(crp.next());

		Collections.sort(data);
		for (int i = 0; i < n; i++)
			labels[i] = data.get(i).cluster;

		/* Getting parameters */
		double[] proportion = new double[crp.size.size()];
		double[] centroidx = new double[crp.moments.size()];
		double[] centroidy = new double[crp.moments.size()];

		for (int c = 0; c < proportion.length; c++) {
			proportion[c] = crp.size.get(c) / (double) n;
			centroidx[c] = crp.moments.get(c)[0];
			centroidy[c] = crp.moments.get(c)[1];
		}
		ArrayList<Integer> freq = crp.size;
		Collections.sort(freq, Collections.reverseOrder());
		System.out.println("done.");
		System.out.println(n + " data points and " + crp.moments.size()
				+ " clusters with size " + freq.toString() + " generated.");

		crp = null;
		initClusters = (initClusters < 1) ? freq.size() : initClusters;

		/* Simulation */
		GibbsSampler gibbs = (singleton) ? new SingletonGibbsSampler(alpha,
				theta, beta, xi, initClusters, maxNumClusters, data)
				: new BlockGibbsSampler(alpha, theta, beta, xi, initClusters,
						maxNumClusters, data);

		ScatterPlot currentPlot = null;
		EvaluationPlot numsPlot = null;
		EvaluationPlot residualPlot = null;
		if (visual > 0) {
			new DistributionPlot(freq, proportion.length);
			numsPlot = new EvaluationPlot(initClusters, "Number of clusters");
			residualPlot = new EvaluationPlot(0, "Wasserstein distance");
			if (eval > 0)
				currentPlot = ScatterPlot.initPlots(data, labels, gibbs.labels,
						freq.size(), singleton);
		}

		final String filename = "num_of_clusters_" + n + "_" + initClusters
				+ "_" + ((singleton) ? "singleton" : "block");
		PrintWriter file = null;
		if (fout)
			file = new PrintWriter(new BufferedWriter(new FileWriter(new File(
					filename + ".txt"), true)));

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
					residualPlot
							.updateSeries(k, gibbs.getResidual(proportion,
									centroidx, centroidy));
			}
			if (k == saveChart)
				ImageIO.write(numsPlot.chart.createBufferedImage(650, 400),
						"png", new File(filename + ".png"));
		}

		if (fout)
			file.close();
	}
}
