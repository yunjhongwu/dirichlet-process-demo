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
import java.util.HashSet;
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
		ArrayList<Double> mux = new ArrayList<Double>();
		ArrayList<Double> muy = new ArrayList<Double>();
		ArrayList<Double> s2x = new ArrayList<Double>();
		ArrayList<Double> s2y = new ArrayList<Double>();

		public CRP(double alpha, double theta, double beta, double xi) {
			this.alpha = alpha;
			this.theta = theta;
			this.eta = 1 / beta;
			this.xi = xi;

		}

		public Point2D getPosition(int cluster) {
			return new Point2D(cluster, (float) sampler.nextGaussian(
					mux.get(cluster), s2x.get(cluster)),
					(float) sampler.nextGaussian(muy.get(cluster),
							s2y.get(cluster)));
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
				double s1 = sampler.nextGamma(theta, eta);
				double s2 = sampler.nextGamma(theta, eta);
				s2x.add(s1);
				s2y.add(s2);
				mux.add(sampler.nextGaussian(0, Math.sqrt(xi * s1)));
				muy.add(sampler.nextGaussian(0, Math.sqrt(xi * s2)));
				cluster++;
			}
			n++;

			return getPosition(cluster);
		}

	}

	public static class GibbsSampler {

		final int n;
		final int maxNumClusters;
		final double alpha, theta, beta, xi;
		final ArrayList<Point2D> data;
		int[] labels;

		// size, mux, muy, s2x, s2y
		HashMap<Integer, double[]> clusters = new HashMap<Integer, double[]>();
		HashSet<Integer> modifiedClusters = new HashSet<Integer>();
		Stack<Integer> emptyClusters = new Stack<Integer>();
		ArrayList<Integer> ord = new ArrayList<Integer>();
		RandomDataGenerator sampler = new RandomDataGenerator();

		public GibbsSampler(double alpha, double theta, double beta, double xi, int initClusters,
				int maxNumClusters, ArrayList<Point2D> data) {
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
			double xsigma2 = (xi + 1) * theta / beta;
			for (int i = 0; i < n; i++) {
				ord.add(i);
				if(i < numClusters)
					clusters.put(i, new double[5]);
				else
					emptyClusters.push(i);
			}
			Collections.shuffle(ord);
			ord.forEach(i -> {
				clusters.get(i % numClusters)[0]++;
				labels[i] = i % numClusters;
			});
			clusters.put(-1, new double[5]);
			clusters.get(-1)[0] = -1;
			clusters.get(-1)[3] = xsigma2;
			clusters.get(-1)[4] = xsigma2;
			updateMoments(clusters.keySet());
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
			return -0.5
					* (Math.log(clusters.get(j)[3])
							+ Math.log(clusters.get(j)[4])
							+ Math.pow(data.get(i).x - clusters.get(j)[1], 2)
							/ clusters.get(j)[3] + Math.pow(data.get(i).y
							- clusters.get(j)[2], 2)
							/ clusters.get(j)[4]);
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

		public void updateMoments(Set<Integer> modifiedClusters) {
			for (Integer c : modifiedClusters) {
				clusters.get(c)[1] = 0;
				clusters.get(c)[2] = 0;
				clusters.get(c)[3] = 0;
				clusters.get(c)[4] = 0;
			}

			for (int i = 0; i < n; i++)
				if (modifiedClusters.contains(labels[i])) {
					clusters.get(labels[i])[1] += data.get(i).x;
					clusters.get(labels[i])[2] += data.get(i).y;
					clusters.get(labels[i])[3] += Math.pow(data.get(i).x, 2);
					clusters.get(labels[i])[4] += Math.pow(data.get(i).y, 2);
				}

			for (Integer c : modifiedClusters) {
				clusters.get(c)[3] = posteriorVariance(clusters.get(c)[1],
						clusters.get(c)[3], clusters.get(c)[0]);
				clusters.get(c)[4] = posteriorVariance(clusters.get(c)[2],
						clusters.get(c)[4], clusters.get(c)[0]);
				clusters.get(c)[1] = posteriorMean(clusters.get(c)[1],
						clusters.get(c)[0]);
				clusters.get(c)[2] = posteriorMean(clusters.get(c)[2],
						clusters.get(c)[0]);
			}

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

					clusters.get(labels[i])[0]--;
					clusters.get(nextCluster)[0]++;
					modifiedClusters.add(labels[i]);
					modifiedClusters.add(nextCluster);
					labels[i] = nextCluster;
				}
		}

		public void next() {
			modifiedClusters = new HashSet<Integer>();
			Collections.shuffle(ord);
			ord.forEach(i -> nextIter(i));

			if (!emptyClusters.isEmpty()
					&& clusters.containsKey(emptyClusters.peek()))
				emptyClusters.pop();

			for (Integer c : new HashSet<Integer>(modifiedClusters))
				if (clusters.get(c)[0] == 0) {
					clusters.remove(c);
					emptyClusters.push(c);
					modifiedClusters.remove(c);
				}
			updateMoments(modifiedClusters);
		}
	}

	public static double getResidual(int n,
			HashMap<Integer, double[]> clusters, final double[] p,
			final double[] mux, final double[] muy) {
		double[] res = new double[(clusters.size() - 1) * mux.length];
		double[] q = new double[clusters.size() - 1];
		int i = 0;
		for (Integer c : clusters.keySet())
			if (c > -1) {
				for (int j = 0; j < mux.length; j++)
					res[i * mux.length + j] = Math.pow(mux[j]
							- clusters.get(c)[1], 2)
							+ Math.pow(muy[j] - clusters.get(c)[2], 2);
				q[i++] = clusters.get(c)[0] / (double) n;
			}

		return getOptResidual(res, q, p);
	}

	public static double getOptResidual(final double[] res,
			final double[] cols, final double[] rows) {
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
		truePlot.setLocation(0, 30);
		truePlot.setVisible(true);

		currentPlot = new ScatterPlot(data, glabels, colors, "DP Model");
		currentPlot.pack();
		RefineryUtilities.centerFrameOnScreen(currentPlot);
		currentPlot.setSize(685, 650);
		currentPlot.setLocation(685, 30);
		currentPlot.setVisible(true);
		return currentPlot;
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) throws InterruptedException {
		final int n = 50000;
		final int maxIters = 1000000;
		final double alpha = 1;
		final double theta = 100;
		final double beta = 20;
		final double xi = 20;
        final int initClusters = 500;
		final int maxNumClusters;
		final int visual = 1;
		final int eval = 0;

		// Generating data ////////////////////////////////////////////////
		System.out.print("Generating data...");
		final ArrayList<Point2D> data = new ArrayList<Point2D>();
		int[] labels = new int[n];
		CRP crp = new CRP(alpha, theta, beta, xi);
		for (int i = 0; i < n; i++)
			data.add(crp.next());

		Collections.sort(data);
		for (int i = 0; i < n; i++)
			labels[i] = data.get(i).cluster;
		maxNumClusters = 50;// n; // crp.mux.size();

		double[] proportion = new double[crp.size.size()];
		double[] centroidx = new double[crp.mux.size()];
		double[] centroidy = new double[crp.mux.size()];

		for (int c = 0; c < proportion.length; c++) {
			proportion[c] = crp.size.get(c) / (double) n;
			centroidx[c] = crp.mux.get(c);
			centroidy[c] = crp.mux.get(c);
		}

		System.out.println("done.");
		System.out.println(n + " data points and " + crp.mux.size()
				+ " clusters with size " + crp.size.toString() + " generated.");
		crp = null;

		// Simulation ///////////////////////////////////////////////////////
		GibbsSampler gibbs = new GibbsSampler(alpha, theta, beta, xi, initClusters,
				maxNumClusters, data);

		ScatterPlot truePlot = null;
		ScatterPlot currentPlot = null;
		if (visual > 0)
			currentPlot = initPlots(truePlot, currentPlot, data, labels,
					gibbs.labels);

		long startTime = System.nanoTime();
		for (int i = 0; i < maxIters; i++) {
			System.out.format("Iteration %d" + "; %d" + " cluster(s); %f"
					+ " milliseconds per iteration\n", i,
					(gibbs.clusters.size() - 1),
					(System.nanoTime() - startTime) / 1000000.0 / (i + 1));
			gibbs.next();
			if (visual > 0 && i % visual == 0)
				updateColors(currentPlot.plot, gibbs.labels);

			if (eval > 0 && i % eval == 0)
				System.out.println("residual = "
						+ getResidual(n, gibbs.clusters, proportion, centroidx,
								centroidy));
		}
	}
}
