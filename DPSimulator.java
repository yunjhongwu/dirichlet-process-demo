/*
 * Created on Mon Nov 12 2014
 * JAVA SE 8
 * Author: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu
 */

import java.awt.Color;
import java.awt.Paint;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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
			Color.cyan, Color.darkGray, Color.lightGray, Color.gray };

	public static class Point2D implements Comparable<Point2D> {
		int cluster;
		double x;
		double y;

		public Point2D(int s, double t, double u) {
			cluster = s;
			x = t;
			y = u;
		}

		public int compareTo(Point2D that) {
			if (this.x != that.x)
				return Double.compare(this.x, that.x);
			else
				return Double.compare(this.y, that.y);
		}

		public String toString() {
			return cluster + " (" + x + ", " + y + ")";
		}
	}

	public static class CRP {
		// X \sim N(\mu_c, \sigma^2)
		// \mu_c \sim N(0, \xi\sigma^2)
		// \sigma^2 \sim Gamma(\theta, \beta)

		double alpha, theta, beta, xi;
		int n = 0;

		RandomDataGenerator sampler = new RandomDataGenerator();
		ArrayList<Integer> customers = new ArrayList<Integer>();
		ArrayList<Double> mux = new ArrayList<Double>();
		ArrayList<Double> muy = new ArrayList<Double>();
		ArrayList<Double> varx = new ArrayList<Double>();
		ArrayList<Double> vary = new ArrayList<Double>();

		public CRP(double a, double t, double b, double x) {
			alpha = a;
			theta = t;
			beta = 1 / b;
			xi = x;

		}

		public Point2D getPosition(int cluster) {
			return new Point2D(cluster, sampler.nextGaussian(mux.get(cluster),
					varx.get(cluster)), sampler.nextGaussian(muy.get(cluster),
					vary.get(cluster)));
		}

		public Point2D next() {
			int cluster = -1;
			double sample = sampler.nextUniform(0, 1) * (n + alpha);
			if (sample < n) {
				while (sample > 0)
					sample -= customers.get(++cluster);
				customers.set(cluster, customers.get(cluster) + 1);
			} else {
				customers.add(1);
				double s1 = sampler.nextGamma(theta, beta);
				double s2 = sampler.nextGamma(theta, beta);
				varx.add(s1);
				vary.add(s2);
				mux.add(sampler.nextGaussian(0, Math.sqrt(xi * s1)));
				muy.add(sampler.nextGaussian(0, Math.sqrt(xi * s2)));
				cluster++;
			}
			n++;

			return getPosition(cluster);
		}

	}

	public static class GibbsSampler {

		int n;
		int Kmax = 10000;
		double alpha, theta, beta, xi, xsigma2;
		int[] labels;

		// size, mux, muy, s2x, s2y
		HashMap<Integer, double[]> clusters = new HashMap<Integer, double[]>();
		Stack<Integer> emptyClusters = new Stack<Integer>();
		final ArrayList<Point2D> data;
		RandomDataGenerator sampler = new RandomDataGenerator();
		double[] p;
		double[] mux;
		double[] muy;
		int eval;
		int maxNumCluster;

		public GibbsSampler(double a, double t, double b, double x,
				ArrayList<Point2D> d, double[] p, double[] mux, double[] muy,
				int maxNumCluster, int eval) {
			data = d;
			n = data.size();
			labels = new int[n];
			alpha = a;
			theta = t;
			beta = b;
			xi = x;
			this.p = p;
			this.mux = mux;
			this.muy = muy;
			this.eval = eval;
			this.maxNumCluster = maxNumCluster;
			for (int i = n - 1; i > 0; i--)
				emptyClusters.push(i);
			double xsigma2 = (xi + 1) * theta / beta;

			clusters.put(-1, new double[5]);
			clusters.get(-1)[3] = xsigma2;
			clusters.get(-1)[4] = xsigma2;
			clusters.put(0, new double[5]);
			clusters.get(0)[0] = n;
			HashSet<Integer> s = new HashSet<Integer>();
			s.add(0);
			updateMoments(s);

		}

		public double posteriorGammaSampler(double mu, double m2, double size) {
			return (2 * theta + size)
					/ (2 * beta + (m2 - Math.pow(mu, 2) / (size + xi)));
			// return sampler.nextGamma(theta + size / 2.0,
			// 2 / (2 * beta + (m2 - Math.pow(mu, 2) / (size + xi))));
		}

		public double posteriorNormalSampler(double mu, double size) {
			return mu / (size + xi);
			// return sampler.nextGaussian(mu / (size + xi), Math.sqrt((1 +
			// xi/ size) * sigma2 ));
		}

		public double logNormalLikelihood(int i, int j) {
			return -Math.log(clusters.get(j)[3]) / 2
					- Math.pow(data.get(i).x - clusters.get(j)[1], 2)
					/ (2 * clusters.get(j)[3]) - Math.log(clusters.get(j)[4])
					/ 2 - Math.pow(data.get(i).y - clusters.get(j)[2], 2)
					/ (2 * clusters.get(j)[4]);
		}

		public double MHthreshold(int i, int j, int k) {
			double t = Math.log(clusters.get(j)[0] - 1)
					+ logNormalLikelihood(i, j);
			t -= Math.log((k == -1) ? alpha / clusters.size()
					: clusters.get(k)[0] - 1);
			t -= logNormalLikelihood(i, k);
			return t;
		}

		public void updateMoments(HashSet<Integer> s) {
			for (Integer c : s) {
				clusters.get(c)[1] = 0;
				clusters.get(c)[2] = 0;
				clusters.get(c)[3] = 0;
				clusters.get(c)[4] = 0;
			}

			for (int i = 0; i < n; i++)
				if (s.contains(labels[i])) {
					clusters.get(labels[i])[1] += data.get(i).x;
					clusters.get(labels[i])[2] += data.get(i).y;
					clusters.get(labels[i])[3] += Math.pow(data.get(i).x, 2);
					clusters.get(labels[i])[4] += Math.pow(data.get(i).y, 2);
				}

			for (Integer c : s)
				if (c > -1) {
					clusters.get(c)[3] = posteriorGammaSampler(
							clusters.get(c)[1], clusters.get(c)[3],
							clusters.get(c)[0]);
					clusters.get(c)[4] = posteriorGammaSampler(
							clusters.get(c)[2], clusters.get(c)[4],
							clusters.get(c)[0]);
					clusters.get(c)[1] = posteriorNormalSampler(
							clusters.get(c)[1], clusters.get(c)[0]);
					clusters.get(c)[2] = posteriorNormalSampler(
							clusters.get(c)[2], clusters.get(c)[0]);
				}

		}

		public void next(int iter) {
			int count = 0;
			ArrayList<Integer> lb = new ArrayList<Integer>(clusters.keySet());
			HashSet<Integer> s = new HashSet<Integer>();

			for (int i = 0; i < this.n; i++) {
				int r = sampler.nextInt(0, lb.size() - 1);
				if (lb.get(r) != -1 || lb.size() <= maxNumCluster)
					if (lb.get(r) != labels[i]
							&& sampler.nextExponential(1) > MHthreshold(i,
									labels[i], lb.get(r))) {
						int c = (lb.get(r) == -1) ? emptyClusters.peek() : lb
								.get(r);

						if (!clusters.containsKey(c))
							clusters.put(c, new double[5]);
						clusters.get(labels[i])[0]--;
						clusters.get(c)[0]++;
						s.add(labels[i]);
						s.add(c);
						labels[i] = c;
						count++;
					}
			}
			if (clusters.size() > lb.size())
				emptyClusters.pop();
			for (Integer c : lb)
				if (clusters.get(c)[0] == 0 && c > -1) {
					clusters.remove(c);
					emptyClusters.push(c);
					s.remove(c);
				}
			updateMoments(s);
			if (iter % eval == 0)
				System.out.println("Acc. rate = " + count / (double) n
						+ "; residual = " + getResidual(p, mux, muy));

		}

		public double getResidual(double[] p, double[] mux, double[] muy) {
			double[] res = new double[(clusters.size() - 1) * mux.length];
			double[] q = new double[clusters.size() - 1];
			int i = 0;
			for (Integer c : clusters.keySet()) {

				if (c > -1) {
					q[i] = clusters.get(c)[0] / (double) n;

					for (int j = 0; j < mux.length; j++)
						res[i * mux.length + j] = Math.pow(
								mux[j] - clusters.get(c)[1], 2)
								+ Math.pow(muy[j] - clusters.get(c)[2], 2);
					i++;
				}

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
				return colors[labels[col] % 12];
			}
		});

		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot
				.getRenderer();
		renderer.setSeriesShape(0, new Ellipse2D.Float(0f, 0f, 1f, 1f), false);

	}

	public static void main(String[] args) throws InterruptedException {
		int n = 1000;
		int iters = 1000000;
		double alpha = 1;
		double theta = 200;
		double beta = 50;
		double xi = 30;

		int maxNumCluster = n;
		int visual = 1;
		int eval = 1000000;

		// / Generating data ////////////////////////////////////////////////
		ArrayList<Point2D> data = new ArrayList<Point2D>();
		int[] labels = new int[n];
		CRP crp = new CRP(alpha, theta, beta, xi);
		for (int i = 0; i < n; i++)
			data.add(crp.next());
		Collections.sort(data);

		for (int i = 0; i < n; i++)
			labels[i] = data.get(i).cluster;

		double[] p = new double[crp.customers.size()];
		for (int i = 0; i < p.length; i++)
			p[i] = crp.customers.get(i) / (double) n;
		double[] centroidx = new double[crp.mux.size()];
		double[] centroidy = new double[crp.mux.size()];
		for (int i = 0; i < centroidx.length; i++) {
			centroidx[i] = crp.mux.get(i);
			centroidy[i] = crp.mux.get(i);
		}
		// //////////////////////////////////////////////////////////////////
		maxNumCluster = n; // crp.mux.size();

		GibbsSampler gibbs = new GibbsSampler(alpha, theta, beta, xi, data, p,
				centroidx, centroidy, maxNumCluster, eval);
		ScatterPlot truePlot = null;
		ScatterPlot currentPlot = null;
		if (visual > 0) {
			truePlot = new ScatterPlot(data, labels, colors, "Data");
			truePlot.pack();
			RefineryUtilities.centerFrameOnScreen(truePlot);
			truePlot.setVisible(true);
			truePlot.setSize(685, 650);
			truePlot.setLocation(0, 30);
			// truePlot.setExtendedState(truePlot.getExtendedState()
			// | Frame.MAXIMIZED_BOTH);

			currentPlot = new ScatterPlot(data, gibbs.labels, colors,
					"DP Model");
			currentPlot.pack();
			RefineryUtilities.centerFrameOnScreen(currentPlot);
			currentPlot.setVisible(true);
			currentPlot.setSize(685, 650);
			currentPlot.setLocation(685, 30);
			// currentPlot.setExtendedState(currentPlot.getExtendedState()
			// | Frame.MAXIMIZED_BOTH);
		}
		long startTime = System.nanoTime();
		for (int i = 0; i < iters; i++) {
			System.out.println("Iteration " + i + "; "
					+ (gibbs.clusters.size() - 1) + " cluster(s)");
			gibbs.next(i);
			if (visual > 0 && i % visual == 0)
				updateColors(currentPlot.plot, gibbs.labels);
			Thread.sleep(20);
		}
		System.out.println((System.nanoTime() - startTime) / 1000000.0 / iters
				+ " milliseconds per iteration");
	}
}
