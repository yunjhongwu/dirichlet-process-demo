public class Point2D implements Comparable<Point2D> {
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
