import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;

public class XYPlotExample extends JPanel {

    private static final int WIDTH = 600;
    private static final int HEIGHT = 400;
    private static final int marg = 20;
    private static final int[] xData = {30, 20, 10, -10, -60, -40};
    private static final int[] yData = {1000, 850, 700, 500, 200, 100};

    static int C2K = 273;

    static int[] range = {203, 333};
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;

        //set color for points
        g2.setPaint(Color.BLACK);

        int w = getWidth();
        int h = getHeight();

        // Draw x and y axes
        g2.drawLine(marg, h - marg, w - marg, h - marg);
        g2.drawLine(marg, marg, marg, h - marg);

        // Draw x and y axis labels
        g2.drawString("X", w - marg/2, h - marg/2);
        g2.drawString("Y", marg - 10, marg);

        // Set the stroke to draw dotted lines
        float[] dashPattern = {5, 5}; // 5 pixels drawn, 5 pixels skipped
        g2.setStroke(new BasicStroke(1, BasicStroke.CAP_ROUND, // cap-the decoration of the ends of a BasicStroke
                                      BasicStroke.JOIN_ROUND, // join-the decoration applied where path segments meet
                                      0, // miterlimit - the limit to trim the miter join.
                                      dashPattern, // dash - the array representing the dashing pattern
                            0)); // dash_phase - the offset to start the dashing pattern
        //set color for points
        g2.setPaint(Color.GRAY);
        int max_x = 350;
        for (int i = C2K-70; i <= C2K+50; i+=50) {
            System.out.println(transformValue2X(i));
            int x = marg + (int) (transformValue2X(i) * (w - 2 * marg));
            g2.drawLine(x, marg, x, h - marg);
            g2.drawString(String.format("%d", i-C2K), x, h - marg/4);

        }

        //set color for points
        g2.setPaint(Color.BLUE);

        // Set the stroke to draw solid lines
        g2.setStroke(new BasicStroke());

        // Draw data points
        int max_y = getMax(yData);
        for (int i = 0; i < xData.length; i++) {
            int x = marg + (int) (transformValue2X(xData[i]+C2K) * (w - 2 * marg));
            //int y = h - marg - (yData[i] * (h - 2 * marg))/max_y;
            int y = marg + (yData[i] * (h - 2 * marg))/max_y; // Inverted Y coordinate
            //g2.fillOval(x - 2, y - 2, 4, 4);
            g2.fill(new Ellipse2D.Double(x-2, y-2, 4, 4));

            // Draw line connecting consecutive points
            if (i > 0) {
                int prevX = marg + (int) (transformValue2X(xData[i - 1]+C2K) * (w - 2 * marg));
                int prevY = marg + (yData[i - 1] * (h - 2 * marg))/max_y;
                g2.drawLine(prevX, prevY, x, y);
            }
        }
    }

    private float transformValue2X(int value) {
        int distance = range[1]-range[0];
        float XCor;
        if ((range[0] <= value) && (range[1] >= value)) {
            XCor = (float) (value - range[0])/distance;
        } else if (range[0] > value) {
            XCor = (float) value/range[1];
        } else {
            XCor = (float) value/range[1];
        }
        return XCor;
    }

    private int getMax(int[] data) {
        int max = Integer.MIN_VALUE;
        for (int value : data) {
            value += C2K;
            if (value > max) {
                max = value;
            }
        }
        return max;
    }

    public static void main(String[] args) {
        JFrame frame = new JFrame("XY Plot Example");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        XYPlotExample plot = new XYPlotExample();
        frame.add(plot);
        frame.setSize(WIDTH, HEIGHT);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}