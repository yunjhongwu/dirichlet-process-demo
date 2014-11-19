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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

public class DPSimulator {
	final static Color[] colors = { Color.red, Color.yellow, Color.blue,
			Color.green, Color.white, Color.magenta, Color.orange, Color.pink,
			Color.cyan, Color.lightGray };

	public static class Point2D implements Comparable<Point2D> {
		final int cluster;
		final float x;
		final float y;

		public Point2D(int cluster, float x, float y) {
			this.cluster = cluster;
			this.x = x;
			this.y = y;
		}

		public int compareTo(Point2D that) {
			if (this.x != that.x)
				return Float.compare(this.x, that.x);
			else
				return Float.compare(this.y, that.y);
		}

		public String toString() {
			return cluster + " (" + x + ", " + y + ")";
		}
	}

	public static class CRP {
		// X \sim N(\mu_c, \sigma^2)
		// \mu_c \sim N(0, \xi\sigma^2)
		// \sigma^2 \sim Gamma(\theta, \beta)

		final double alpha, theta, eta, xi;
		int n = 0;

		RandomDataGenerator sampler = new RandomDataGenerator();
		ArrayList<Integer> size = new ArrayList<Integer>();
		ArrayList<double[]> moments = new ArrayList<double[]>();

		public CRP(double alpha, double theta, double beta, double xi) {
			sampler.reSeed();
			this.alpha = alpha;
			this.theta = theta;
			this.eta = 1 / beta;
			this.xi = xi;
		}

		public CRP(double alpha, double theta, double beta, double xi, int i) {
			this(alpha, theta, beta, xi);
			sampler.reSeed(i);
		}

		public Point2D getPosition(int cluster) {
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
				double[] m = new double[4];
				m[2] = sampler.nextGamma(theta, eta);
				m[3] = sampler.nextGamma(theta, eta);
				m[0] = sampler.nextGaussian(0, Math.sqrt(xi * m[2]));
				m[1] = sampler.nextGaussian(0, Math.sqrt(xi * m[3]));
				moments.add(m);
				cluster++;
			}
			n++;

			return getPosition(cluster);
		}
	}

	public static abstract class GibbsSampler {
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

		public double posteriorVariance(double mu, double m2, double size) {
			return (2 * theta + size)
					/ (2 * beta + (m2 - Math.pow(mu, 2) / (size + xi)));
			// return sampler.nextGamma(theta + size / 2.0,
			// 2 / (2 * beta + (m2 - Math.pow(mu, 2) / (size + xi))));
		}

		public double posteriorMean(double mu, double size) {
			return mu / (size + xi);
			// return sampler.nextGaussian(mu / (size + xi), Math.sqrt((1 +
			// xi/ size) * sigma2 ));
		}

		public double logNormalLikelihood(int i, int j) {
			double s2x = posteriorVariance(clusters.get(j)[1],
					clusters.get(j)[3], clusters.get(j)[0]);
			double s2y = posteriorVariance(clusters.get(j)[2],
					clusters.get(j)[4], clusters.get(j)[0]);
			return -0.5
					* (Math.log(s2x)
							+ Math.log(s2y)
							+ Math.pow(
									data.get(i).x
											- posteriorMean(clusters.get(j)[1],
													clusters.get(j)[0]), 2)
							/ s2x + Math.pow(
							data.get(i).y
									- posteriorMean(clusters.get(j)[2],
											clusters.get(j)[0]), 2)
							/ s2y);
		}

		public double MHThreshold(int i, int j, int k) {
			if (clusters.get(j)[0] == 1)
				return 0;
			double t = Math.log(clusters.get(j)[0] - 1)
					+ logNormalLikelihood(i, j);
			t -= Math.log((k == -1) ? alpha / clusters.size()
					: clusters.get(k)[0]);
			t -= logNormalLikelihood(i, k);
			return t;
		}

		public int MHKernel() {
			int s = 0;
			int r = sampler.nextInt(0, clusters.size() - 1);
			for (Integer c : clusters.keySet())
				if (r == s++)
					return c;
			return -1;
		}

		public void updateMoments(int i, int j, int k) {
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

		public void updateAllMoments() {
			updateAllMoments(clusters.keySet());
		}

		public void updateAllMoments(Set<Integer> modifiedClusters) {
			for (Integer c : modifiedClusters) {
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
												- posteriorMean(
														clusters.get(c)[2],
														clusters.get(c)[0]), 2);
					q[i++] = clusters.get(c)[0] / (double) n;
				}

			return getOptResidual(res, q, p);
		}

		public double getOptResidual(final double[] res, final double[] cols,
				final double[] rows) {
			LinearObjectiveFunction f = new LinearObjectiveFunction(res, 0);
			Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
			for (int i = 0; i < cols.length; i++) {
				double[] b = new double[cols.length * rows.length];
				for (int j = i * rows.length; j < (i + 1) * rows.length; j++)
					b[j]++;
				constraints.add(new LinearConstraint(b, Relationship.EQ,
						cols[i]));
			}

			for (int i = 0; i < rows.length; i++) {
				double[] b = new double[cols.length * rows.length];
				for (int j = i; j < b.length; j += rows.length)
					b[j]++;
				constraints.add(new LinearConstraint(b, Relationship.EQ,
						rows[i]));
			}
			SimplexSolver solver = new SimplexSolver();

			return solver.optimize(new MaxIter(10000), f,
					new LinearConstraintSet(constraints), GoalType.MINIMIZE,
					new NonNegativeConstraint(true)).getSecond();
		}

		public abstract void nextIter(int i);

		public abstract void next(int k);

	}

	public static class SingletonGibbsSampler extends GibbsSampler {
		public SingletonGibbsSampler(double alpha, double theta, double beta,
				double xi, int initClusters, int maxNumClusters,
				ArrayList<Point2D> data) {
			super(alpha, theta, beta, xi, initClusters, maxNumClusters, data);
		}

		public void nextIter(int i) {
			int nextCluster = MHKernel();
			if (nextCluster != -1 || clusters.size() <= maxNumClusters)
				if (nextCluster != labels[i]
						&& sampler.nextExponential(1) > MHThreshold(i,
								labels[i], nextCluster)) {
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

	public static class VectorGibbsSampler extends GibbsSampler {
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
						&& sampler.nextExponential(1) > MHThreshold(i,
								labels[i], nextCluster)) {
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

	public static class ScatterPlot extends ApplicationFrame {
		private static final long serialVersionUID = 1L;
		public XYPlot plot;

		public ScatterPlot(final ArrayList<Point2D> data, int[] labels,
				Color[] colors, String title) {
			super(title);
			XYSeriesCollection xyseries = getData(data);

			final JFreeChart chart = ChartFactory.createScatterPlot("DP Model",
					"X", "Y", xyseries, PlotOrientation.VERTICAL, false, false,
					false);

			final ChartPanel panel = new ChartPanel(chart);
			plot = (XYPlot) chart.getPlot();
			updateColors(plot, labels);
			plot.setDomainCrosshairVisible(false);
			plot.setRangeCrosshairVisible(false);
			plot.setDomainGridlinesVisible(false);
			plot.setRangeGridlinesVisible(false);
			plot.setBackgroundPaint(Color.black);

			chart.setBackgroundPaint(Color.black);
			chart.setBorderVisible(false);
			panel.setBackground(Color.black);
			setContentPane(panel);
		}

		private XYSeriesCollection getData(final ArrayList<Point2D> data) {
			XYSeriesCollection xyseries = new XYSeriesCollection();
			XYSeries series = new XYSeries("data");
			for (int i = 0; i < data.size(); i++)
				series.add(data.get(i).x, data.get(i).y);
			xyseries.addSeries(series);
			return xyseries;
		}
	}

	private static void updateColors(XYPlot plot, int[] labels) {
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
