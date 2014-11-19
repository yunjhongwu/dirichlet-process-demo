/*
 * Created on Nov 19, 2014
 * JAVA SE 8
 * Auther: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu 
 */

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

public class TrackNumClusters extends ApplicationFrame {
	private static final long serialVersionUID = 1L;
	private XYSeries series;
	private XYSeries trueValue;
	private int trueNum;
	protected final JFreeChart chart;

	public TrackNumClusters(int trueNum, int initNum) {
		super("Number of Clusters");
		this.series = new XYSeries("Estimates");
		this.trueValue = new XYSeries("Parameter");
		this.trueNum = trueNum;
		final XYSeriesCollection nums = new XYSeriesCollection(series);
		nums.addSeries(trueValue);
		chart = ChartFactory.createXYLineChart("DP Model", "Iterations",
				"Number of Clusters", nums, PlotOrientation.VERTICAL, false,
				false, false);
		final XYPlot plot = chart.getXYPlot();
		final JPanel content = new JPanel(new BorderLayout());
		final ChartPanel panel = new ChartPanel(chart);

		content.add(panel);
		setContentPane(content);
		plot.getRenderer().setSeriesPaint(0, Color.yellow);
		plot.getRenderer().setSeriesPaint(1, Color.red);
		plot.setDomainCrosshairVisible(false);
		plot.setRangeCrosshairVisible(false);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
		plot.setBackgroundPaint(Color.black);
		chart.setBackgroundPaint(Color.black);
		chart.setBorderVisible(false);
		panel.setPreferredSize(new java.awt.Dimension(650, 400));
		panel.setBackground(Color.black);
		series.add(0, initNum);
		trueValue.add(0, trueNum);
		initPlot();

	}

	private void initPlot() {
		this.pack();
		RefineryUtilities.centerFrameOnScreen(this);
		this.setVisible(true);
	};

	public void updateSeries(int iters, int num) {
		series.add(iters, num);
		trueValue.add(iters, trueNum);
	}

}
