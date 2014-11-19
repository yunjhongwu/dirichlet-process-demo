/*
 * Created on Nov 18 2014
 * JAVA SE 8
 * Author: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu
 */

import java.awt.Color;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;

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

public class ScatterPlot extends ApplicationFrame {
	private static final long serialVersionUID = 1L;
	public XYPlot plot;
	final static Color[] colors = { Color.red, Color.yellow, Color.blue,
			Color.green, Color.white, Color.magenta, Color.orange, Color.pink,
			Color.cyan, Color.lightGray };

	public ScatterPlot(final ArrayList<Point2D> data, int[] labels, String title) {
		super(title);
		XYSeriesCollection xyseries = getData(data);

		final JFreeChart chart = ChartFactory.createScatterPlot("DP Model",
				"X", "Y", xyseries, PlotOrientation.VERTICAL, false, false,
				false);

		final ChartPanel panel = new ChartPanel(chart);
		plot = (XYPlot) chart.getPlot();
		updateColors(labels);
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

	protected static ScatterPlot initPlots(ArrayList<Point2D> data,
			final int[] labels, final int[] glabels) {
		ScatterPlot truePlot = new ScatterPlot(data, labels, "Data");
		truePlot.pack();
		RefineryUtilities.centerFrameOnScreen(truePlot);
		truePlot.setSize(685, 650);
		truePlot.setLocation(0, 20);
		truePlot.setVisible(true);

		ScatterPlot currentPlot = new ScatterPlot(data, glabels, "DP Model");
		currentPlot.pack();
		RefineryUtilities.centerFrameOnScreen(currentPlot);
		currentPlot.setSize(685, 650);
		currentPlot.setLocation(685, 20);
		currentPlot.setVisible(true);
		return currentPlot;
	}

	private XYSeriesCollection getData(final ArrayList<Point2D> data) {
		XYSeriesCollection xyseries = new XYSeriesCollection();
		XYSeries series = new XYSeries("data");
		for (int i = 0; i < data.size(); i++)
			series.add(data.get(i).x, data.get(i).y);
		xyseries.addSeries(series);
		return xyseries;
	}

	protected void updateColors(int[] labels) {
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
}