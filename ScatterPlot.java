import java.awt.Color;
import java.util.ArrayList;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

public class ScatterPlot extends ApplicationFrame {
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
		DPSimulator.updateColors(plot, labels);
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