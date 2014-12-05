/*
 * Created on Nov 19, 2014
 * JAVA SE 8
 * Auther: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu 
 */

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

public class Evaluation extends ApplicationFrame {
	private static final long serialVersionUID = 1L;
	private XYSeries series;
	protected final JFreeChart chart;

	public Evaluation(int initNum, String title) {
		super(title);
		this.series = new XYSeries("Estimates");
		final XYSeriesCollection nums = new XYSeriesCollection(series);
		chart = ChartFactory.createXYLineChart("DP Model", "Iterations", title,
				nums, PlotOrientation.VERTICAL, false, false, false);
		final XYPlot plot = chart.getXYPlot();
		final JPanel content = new JPanel(new BorderLayout());
		final ChartPanel panel = new ChartPanel(chart);

		content.add(panel);
		setContentPane(content);
		plot.getDomainAxis().setStandardTickUnits(
				NumberAxis.createIntegerTickUnits());
		plot.getRenderer().setSeriesPaint(0, Color.yellow);
		plot.setDomainCrosshairVisible(false);
		plot.setRangeCrosshairVisible(false);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
		plot.getDomainAxis().setTickLabelPaint(Color.lightGray);
		plot.getRangeAxis().setTickLabelPaint(Color.lightGray);
		plot.getDomainAxis().setAxisLinePaint(Color.lightGray);
		plot.getRangeAxis().setAxisLinePaint(Color.lightGray);
		plot.getDomainAxis().setLabelPaint(Color.lightGray);
		plot.getRangeAxis().setLabelPaint(Color.lightGray);
		plot.setBackgroundPaint(Color.black);
		chart.setBackgroundPaint(Color.black);
		chart.setBorderVisible(false);
		panel.setBackground(Color.black);
		panel.setPreferredSize(new Dimension(685, 220));
		if (initNum > 0)
			series.add(0, initNum);
		this.setLocation(685, 520);
		this.pack();
		this.setVisible(true);
	}

	public void updateSeries(int iters, int num) {
		series.add(iters, num);
	}

	public void updateSeries(int iters, double res) {
		series.add(iters, res);
	}

}
