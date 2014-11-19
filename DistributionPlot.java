/*
 * Created on Nov 19, 2014
 * JAVA SE 8
 * Auther: Yun-Jhong Wu
 * E-mail: yjwu@umich.edu 
 */

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.util.ArrayList;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

public class DistributionPlot extends ApplicationFrame {
	private static final long serialVersionUID = 1L;

	public DistributionPlot(ArrayList<Integer> freq, int num) {
		super("Size of Clusters");
		double[][] freqDouble = new double[1][freq.size()];
		for (int i = 0; i < freq.size(); i++)
			freqDouble[0][i] = freq.get(i);

		final CategoryDataset data = DatasetUtilities.createCategoryDataset("",
				"", freqDouble);
		final JFreeChart chart = ChartFactory.createBarChart(
				"Size of Clusters", "Clusters", "Size", data,
				PlotOrientation.VERTICAL, false, false, false);

		final CategoryPlot plot = chart.getCategoryPlot();
		final JPanel content = new JPanel(new BorderLayout());
		final ChartPanel panel = new ChartPanel(chart);

		content.add(panel);
		setContentPane(content);
		plot.getRenderer().setSeriesPaint(0, Color.green);
		plot.setDomainCrosshairVisible(false);
		plot.setRangeCrosshairVisible(false);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
		plot.setBackgroundPaint(Color.black);
		chart.setBackgroundPaint(Color.black);
		chart.setBorderVisible(false);
		panel.setBackground(Color.black);
		panel.setPreferredSize(new Dimension(650, 400));
		this.pack();
		RefineryUtilities.centerFrameOnScreen(this);
		this.setVisible(true);
	}

}
