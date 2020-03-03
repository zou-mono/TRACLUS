package mono;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class MainFrame2 extends JFrame {
    private Canvas p = null;

    public MainFrame2() {

    }

    public MainFrame2(ArrayList<Trajectory>  trajectoryAL,  ArrayList<CMDPoint> lineSegmentPointArray,
                      ArrayList<Cluster> clusterRepresentativeTrajectoryAL) {
        p = new Canvas();
        p.addTrajectorys(trajectoryAL);
        p.addClusterComponent(lineSegmentPointArray);
        p.addrTrajectorys(clusterRepresentativeTrajectoryAL);

        JFrame f = new JFrame("效果图");
        // 取屏幕分辨率的大小
        p.setPreferredSize(Toolkit.getDefaultToolkit().getScreenSize());
        f.add(p, BorderLayout.CENTER);
        f.pack();
        f.setVisible(true);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        p.repaint();

        try{
            saveImageToFile(p);
        }catch (Exception e){
            e.printStackTrace();
        }
    }


    private void saveImageToFile(Canvas p) throws IOException {
        Dimension imagesize = Toolkit.getDefaultToolkit().getScreenSize();
        BufferedImage image = new BufferedImage(imagesize.width,imagesize.height,BufferedImage.TYPE_INT_RGB);
        Graphics2D g = image.createGraphics();
        g.setColor(Color.WHITE);//设置笔刷白色
        g.fillRect(0,0,imagesize.width,imagesize.height);//填充整个屏幕
        g.setColor(Color.BLACK);
        p.paint(g);
        g.dispose();

        File f = new File("res/savePic.jpg");
        if(!f.exists()) {
            f.createNewFile();
        }
        ImageIO.write(image, "jpg", f);
    }
}

class Canvas extends JPanel {
    private static final long serialVersionUID = 1L;
    private ArrayList<Trajectory> trajectoryAL;
    private ArrayList<Cluster> clusterRepresentativeTrajectoryAL;
    private ArrayList<CMDPoint> lineSegmentPointArray;
    public Graphics g;

    public void paint(Graphics g) {
//		drawTrajectory(trajecotrys, rTrajectory);
        this.g = g;
        drawTrajectory(trajectoryAL);
//        drawClusterComponentLineSegments(lineSegmentPointArray);
        drawrTrajectory(clusterRepresentativeTrajectoryAL);
    }

    public void addClusterComponent(ArrayList<CMDPoint> lineSegmentPointArray) {
        this.lineSegmentPointArray = lineSegmentPointArray;
    }

    public void addTrajectorys(ArrayList<Trajectory>  trajectoryAL) {
        this.trajectoryAL = trajectoryAL;
    }

    public void addrTrajectorys(ArrayList<Cluster> clusterRepresentativeTrajectoryAL) {
        this.clusterRepresentativeTrajectoryAL = clusterRepresentativeTrajectoryAL;
    }

    private void drawTrajectory(ArrayList<Trajectory> trajectoryAL) {
        g.setColor(Color.gray);
        //java绘图左上角是原点
//        int icount = 0;
        for (int i = 0; i < trajectoryAL.size();i++) {
            for (int m = 0; m < trajectoryAL.get(i).getM_pointArray().size() - 2 ;m++) {
                double startX = trajectoryAL.get(i).getM_pointArray().get(m).getM_coordinate(0);
                double startY = trajectoryAL.get(i).getM_pointArray().get(m).getM_coordinate(1);
                double endX = trajectoryAL.get(i).getM_pointArray().get(m+1).getM_coordinate(0);
                double endY = trajectoryAL.get(i).getM_pointArray().get(m+1).getM_coordinate(1);

                g.drawLine((int)startX, (int)startY, (int)endX, (int)endY);
//                g.drawLine((int) (startX - 82154) * 1200 / 91202,
//                        600 - (int) (startY - 8145) * 600 / 47087,
//                        (int) (endX - 82154) * 1200 / 91202,
//                        600 - (int) (endY - 8145) * 600 / 47087);

//                icount++;
            }
        }
//        System.out.println(icount);
    }

    private void drawClusterComponentLineSegments(ArrayList<CMDPoint> lineSegmentPointArray) {
        g.setColor(Color.gray);

//        int icount = 0;

        for (int i = 0; i < lineSegmentPointArray.size(); i++) {
            CMDPoint lineSegmentEntry = lineSegmentPointArray.get(i);
            double startX = lineSegmentEntry.getM_coordinate(0);
            double startY = lineSegmentEntry.getM_coordinate(1);
            double endX = lineSegmentEntry.getM_coordinate(2);
            double endY = lineSegmentEntry.getM_coordinate(3);

            g.drawLine((int)startX, (int)startY, (int)endX, (int)endY);
//            g.drawLine((int) (startX - 82154) * 1200 / 91202,
//                    600 - (int) (startY - 8145) * 600 / 47087,
//                    (int) (endX - 82154) * 1200 / 91202,
//                    600 - (int) (endY - 8145) * 600 / 47087);
//            icount++;
        }
//        System.out.println(icount);
    }

    private void drawrTrajectory(ArrayList<Cluster> clusterRepresentativeTrajectoryAL) {
        g.setColor(Color.RED);

//        int icount = 0;
        for (int i = 0; i < clusterRepresentativeTrajectoryAL.size();i++) {
            for (int j = 0; j < clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().size() - 2; j++) {
                double startX = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j).getM_coordinate(0);
                double startY = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j).getM_coordinate(1);
                double endX = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j + 1).getM_coordinate(0);
                double endY = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j + 1).getM_coordinate(1);

                g.drawLine((int)startX, (int)startY, (int)endX, (int)endY);
//                g.drawLine((int) (startX - 82154) * 1200 / 91202,
//                        600 - (int) (startY - 8145) * 600 / 47087,
//                        (int) (endX - 82154) * 1200 / 91202,
//                        600 - (int) (endY - 8145) * 600 / 47087);
//                icount++;
            }
        }
//        System.out.println(icount);
    }
}

