/*
 * Created on Mon Nov 12 2014
 * JAVA SE 8
 * Author: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu
 */

import java.awt.Color;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Collections;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.ui.RefineryUtilities;

public class DPSimulator {
	final static Color[] colors = { Color.red, Color.yellow, Color.blue,
			Color.green, Color.white, Color.magenta, Color.orange, Color.pink,
			Color.cyan, Color.lightGray };

	protected static void updateColors(XYPlot plot, int[] labels) {
		plot.setRenderer(new XYLineAndShapeRenderer(false, true) {
			private static final long serialVersionUID = 1L;

			@Override
			public Paint getItemPaint(int row, int col) {
				Color baseColor = colors[labels[col] % 10];
				switch (labels[col] / 10 % 2) {
				case 1:
					return baseColor.darker();
				default:
					return baseColor;
				}
			}

			@Override
			public Shape getSeriesShape(int series) {
				return new Ellipse2D.Float(0f, 0f, 1f, 1f);
			}
		});
	}

	public static ScatterPlot initPlots(ScatterPlot truePlot,
			ScatterPlot currentPlot, ArrayList<Point2D> data,
			final int[] labels, final int[] glabels) {
		truePlot = new ScatterPlot(data, labels, colors, "Data");
		truePlot.pack();
		RefineryUtilities.centerFrameOnScreen(truePlot);
		truePlot.setSize(685, 650);
		truePlot.setLocation(0, 20);
		truePlot.setVisible(true);

		currentPlot = new ScatterPlot(data, glabels, colors, "DP Model");
		currentPlot.pack();
		RefineryUtilities.centerFrameOnScreen(currentPlot);
		currentPlot.setSize(685, 650);
		currentPlot.setLocation(685, 20);
		currentPlot.setVisible(true);
		return currentPlot;
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) throws InterruptedException {
		final int n = 10000;
		final int maxIters = Integer.MAX_VALUE;
		final double alpha = 1;
		final double theta = 100;
		final double beta = 20;
		final double xi = 20;
		final int initClusters = 1;
		final int maxNumClusters;
		final int visual = 1;
		final int eval = 0;
		final boolean singleton = true;

		/* Generating data */
		System.out.print("Generating data...");
		final ArrayList<Point2D> data = new ArrayList<Point2D>();
		int[] labels = new int[n];
		CRP crp = new CRP(alpha, theta, beta, xi, 0);
		for (int i = 0; i < n; i++)
			data.add(crp.next());

		Collections.sort(data);
		for (int i = 0; i < n; i++)
			labels[i] = data.get(i).cluster;

		maxNumClusters = n; // crp.moments.size();

		/* Getting parameters */
		double[] proportion = new double[crp.size.size()];
		double[] centroidx = new double[crp.moments.size()];
		double[] centroidy = new double[crp.moments.size()];

		for (int c = 0; c < proportion.length; c++) {
			proportion[c] = crp.size.get(c) / (double) n;
			centroidx[c] = crp.moments.get(c)[0];
			centroidy[c] = crp.moments.get(c)[1];
		}

		System.out.println("done.");
		System.out.println(n + " data points and " + crp.moments.size()
				+ " clusters with size " + crp.size.toString() + " generated.");
		crp = null;

		/* Simulation */
		GibbsSampler gibbs = (singleton) ? new SingletonGibbsSampler(alpha,
				theta, beta, xi, initClusters, maxNumClusters, data)
				: new VectorGibbsSampler(alpha, theta, beta, xi, initClusters,
						maxNumClusters, data);

		ScatterPlot truePlot = null;
		ScatterPlot currentPlot = null;
		if (visual > 0)
			currentPlot = initPlots(truePlot, currentPlot, data, labels,
					gibbs.labels);

		long startTime = System.nanoTime();
		for (int k = 0; k < maxIters; k++) {
			gibbs.next(k);
			System.out
					.format("Iteration %d; %f milliseconds per iteration; %d clusters. ",
							k, (System.nanoTime() - startTime) / 1000000.0
									/ (k + 1), (gibbs.clusters.size() - 1));
			System.out.println((eval > 0 && k % eval == 0) ? "residual = "
					+ gibbs.getResidual(n, proportion, centroidx, centroidy)
					: "");

			if (visual > 0 && k % visual == 0)
				updateColors(currentPlot.plot, gibbs.labels);
		}
	}
}
