package mono;

import com.sun.media.jfxmedia.logging.Logger;
import org.apache.logging.log4j.LogManager;
import org.geotools.data.DataUtilities;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.Transaction;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.SchemaException;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.MultiLineString;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.util.*;
import java.util.List;
import mono.ClusterGen.LineSegmentId;

public class MainFrame extends JFrame {
    private static final org.apache.logging.log4j.Logger logger = LogManager.getLogger(ClusterGen.class.getName());

    private Canvas p = null;

    public MainFrame() {

    }

    public MainFrame(ArrayList<Trajectory> trajectoryAL, ArrayList<CMDPoint> lineSegmentPointArray,
                     ArrayList<Cluster> clusterRepresentativeTrajectoryAL, ArrayList<Integer> componentIdArray,
                     ArrayList<LineSegmentId> idArray) {
        p = new Canvas();
        p.addTrajectorys(trajectoryAL);
//        p.addClusterComponent(lineSegmentPointArray);
        p.addrTrajectorys(clusterRepresentativeTrajectoryAL);

        JFrame f = new JFrame("效果图");
        // 取屏幕分辨率的大小
        p.setPreferredSize(Toolkit.getDefaultToolkit().getScreenSize());
        f.add(p, BorderLayout.CENTER);
        f.pack();
        f.setVisible(true);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        p.repaint();

        try {
//            saveImageToFile(p);
            logger.info("开始输出聚类后的线段集...");
            saveClustedLineSegmentToSHP(componentIdArray, lineSegmentPointArray, idArray);
            logger.info("聚类后线段集成功输出至shapefile.");

            logger.info("开始输出表达线...");
            saveRepresentiveTrajectoryToSHP(clusterRepresentativeTrajectoryAL);
            logger.info("表达线成功至shapefile.");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void saveRepresentiveTrajectoryToSHP(ArrayList<Cluster> clusterRepresentativeTrajectoryAL){
        String shpfile = "res/representiveTrajectory.shp";
        File newFile = new File(shpfile);
        //创建shapefileDataStore工厂
        ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();

        try {
            //参数设置
            Map<String, Serializable> params = new HashMap<>();
            params.put("url", newFile.toURI().toURL());
            params.put("create spatial index", Boolean.TRUE);

            //根据关键字创建shapefileDataStore
            ShapefileDataStore newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);

            //设置编码，防止中文乱码
            Charset charset = Charset.forName("GBK");
            newDataStore.setCharset(charset);

            //设置要素的字段名称及其类型
            final SimpleFeatureType TYPE =
                    DataUtilities.createType(
                            "representiveTrajectory",
                            "the_geom:LineString," + // geometry属性设置
                                    "rtraID:Integer," + //表达线ID
                                    "classID:Integer"// 存储线段所属的类别
                    );
            newDataStore.createSchema(TYPE);

            //创建要素集合
            List<SimpleFeature> features = new ArrayList<>();

            //创建要素模板
            GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();
            SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(TYPE);

            for (int i = 0; i < clusterRepresentativeTrajectoryAL.size(); i++) {
                Cluster clsEntry = clusterRepresentativeTrajectoryAL.get(i);
                Coordinate[] coords = new Coordinate[clsEntry.getM_PointArray().size()];

                for (int j = 0; j < clsEntry.getM_PointArray().size(); j++) {
                    double X = clsEntry.getM_PointArray().get(j).getM_coordinate(0);
                    double Y = clsEntry.getM_PointArray().get(j).getM_coordinate(1);

                    coords[j] = new Coordinate(X, Y);
                }
                //创建一个lineString
                LineString line = geometryFactory.createLineString(coords);
                //添加geometry属性
                featureBuilder.add(line);
                //添加rtraID
                featureBuilder.add(i);
                //添加classID属性
                featureBuilder.add(clsEntry.getM_clusterId());
                //构建要素
                SimpleFeature feature = featureBuilder.buildFeature(null);
                //将要素添加到要素几何中
                features.add(feature);
            }

            writeFeatureToSHP(newDataStore, TYPE, features);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void saveClustedLineSegmentToSHP(ArrayList<Integer> componentIdArray, ArrayList<CMDPoint> lineSegmentPointArray,
                                             ArrayList<LineSegmentId> idArray) {

        String shpfile = "res/clustedLineSegments.shp";
        File newFile = new File(shpfile);
        //创建shapefileDataStore工厂
        ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();

        try {
            //参数设置
            Map<String, Serializable> params = new HashMap<>();
            params.put("url", newFile.toURI().toURL());
            params.put("create spatial index", Boolean.TRUE);

            //根据关键字创建shapefileDataStore
            ShapefileDataStore newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);

            //设置编码，防止中文乱码
            Charset charset = Charset.forName("GBK");
            newDataStore.setCharset(charset);

            //设置要素的字段名称及其类型
            final SimpleFeatureType TYPE =
                    DataUtilities.createType(
                            "clustedLineSegment",
                            "the_geom:LineString," + // geometry属性设置
                                    "lineID:Integer," + //线段ID
                                    "classID:Integer," + // 存储线段所属的类别ID
                                    "trajID:Integer" //存储线段所属的轨迹ID
                    );
            newDataStore.createSchema(TYPE);

            //创建要素集合
            List<SimpleFeature> features = new ArrayList<>();

            //创建要素模板
            GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();
            SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(TYPE);

            for (int i = 0; i < lineSegmentPointArray.size(); i++) {
                CMDPoint lineSegmentEntry = lineSegmentPointArray.get(i);
                Integer clusterID = componentIdArray.get(i);
                Integer trajectoryID = idArray.get(i).trajectoryId;

                double startX = lineSegmentEntry.getM_coordinate(0);
                double startY = lineSegmentEntry.getM_coordinate(1);
                double endX = lineSegmentEntry.getM_coordinate(2);
                double endY = lineSegmentEntry.getM_coordinate(3);

                //创建一个lineString
                LineString line = geometryFactory.createLineString(
                        new Coordinate[]{new Coordinate(startX, startY), new Coordinate(endX, endY)});

                //添加属性
                featureBuilder.add(line);
                featureBuilder.add(i);
                featureBuilder.add(clusterID);
                featureBuilder.add(trajectoryID);
                //构建要素
                SimpleFeature feature = featureBuilder.buildFeature(null);
                //将要素添加到要素几何中
                features.add(feature);
            }

            writeFeatureToSHP(newDataStore, TYPE, features);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void writeFeatureToSHP(ShapefileDataStore newDataStore, SimpleFeatureType TYPE, List<SimpleFeature> features) {
        try {
            Transaction transaction = new DefaultTransaction("create");

            String typeName = newDataStore.getTypeNames()[0];
            SimpleFeatureSource featureSource = newDataStore.getFeatureSource(typeName);

//        SimpleFeatureType SHAPE_TYPE = featureSource.getSchema();
//            System.out.println("SHAPE:" + SHAPE_TYPE);

            if (featureSource instanceof SimpleFeatureStore) {
                SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;
                SimpleFeatureCollection collection = new ListFeatureCollection(TYPE, features);
                featureStore.setTransaction(transaction);
                try {
                    featureStore.addFeatures(collection);
                    transaction.commit();
                } catch (Exception problem) {
                    problem.printStackTrace();
                    transaction.rollback();
                } finally {
                    transaction.close();
                }
//                System.exit(0); // success!
            } else {
                logger.error(typeName + " does not support read/write access");
//                System.exit(1);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void saveImageToFile(Canvas p) throws IOException {
        Dimension imagesize = Toolkit.getDefaultToolkit().getScreenSize();
        BufferedImage image = new BufferedImage(imagesize.width, imagesize.height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = image.createGraphics();
        g.setColor(Color.WHITE);//设置笔刷白色
        g.fillRect(0, 0, imagesize.width, imagesize.height);//填充整个屏幕
        g.setColor(Color.BLACK);
        p.paint(g);
        g.dispose();

        File f = new File("res/savePic.jpg");
        if (!f.exists()) {
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

    public void addTrajectorys(ArrayList<Trajectory> trajectoryAL) {
        this.trajectoryAL = trajectoryAL;
    }

    public void addrTrajectorys(ArrayList<Cluster> clusterRepresentativeTrajectoryAL) {
        this.clusterRepresentativeTrajectoryAL = clusterRepresentativeTrajectoryAL;
    }

    private void drawTrajectory(ArrayList<Trajectory> trajectoryAL) {
        g.setColor(Color.gray);
        //java绘图左上角是原点
//        int icount = 0;
        for (int i = 0; i < trajectoryAL.size(); i++) {
            for (int m = 0; m < trajectoryAL.get(i).getM_pointArray().size() - 2; m++) {
                double startX = trajectoryAL.get(i).getM_pointArray().get(m).getM_coordinate(0);
                double startY = trajectoryAL.get(i).getM_pointArray().get(m).getM_coordinate(1);
                double endX = trajectoryAL.get(i).getM_pointArray().get(m + 1).getM_coordinate(0);
                double endY = trajectoryAL.get(i).getM_pointArray().get(m + 1).getM_coordinate(1);

//                g.drawLine((int)startX, (int)startY, (int)endX, (int)endY);
                g.drawLine((int) (startX - 82154) * 1200 / 91202,
                        600 - (int) (startY - 8145) * 600 / 47087,
                        (int) (endX - 82154) * 1200 / 91202,
                        600 - (int) (endY - 8145) * 600 / 47087);

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

            g.drawLine((int) startX, (int) startY, (int) endX, (int) endY);
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
        for (int i = 0; i < clusterRepresentativeTrajectoryAL.size(); i++) {
            for (int j = 0; j < clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().size() - 2; j++) {
                double startX = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j).getM_coordinate(0);
                double startY = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j).getM_coordinate(1);
                double endX = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j + 1).getM_coordinate(0);
                double endY = clusterRepresentativeTrajectoryAL.get(i).getM_PointArray().get(j + 1).getM_coordinate(1);

//                g.drawLine((int)startX, (int)startY, (int)endX, (int)endY);
                g.drawLine((int) (startX - 82154) * 1200 / 91202,
                        600 - (int) (startY - 8145) * 600 / 47087,
                        (int) (endX - 82154) * 1200 / 91202,
                        600 - (int) (endY - 8145) * 600 / 47087);
//                icount++;
            }
        }
//        System.out.println(icount);
    }
}

