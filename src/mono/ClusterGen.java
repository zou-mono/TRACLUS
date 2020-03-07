package mono;

import java.awt.*;
import java.util.*;
import java.util.List;

import mono.TraClusterDoc.Parameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.strtree.STRtree;

public class ClusterGen {
    private static final Logger logger = LogManager.getLogger(ClusterGen.class.getName());

    public TraClusterDoc m_document;

    private double m_epsParam;
    private int m_minLnsParam;
    private int m_nTotalLineSegments; //存储分段后的线段数量
    private int m_currComponentId;
    // the number of dense components discovered until now
    private ArrayList<Integer> m_componentIdArray = new ArrayList<Integer>();

    // the list of line segment clusters
    private LineSegmentCluster[] m_lineSegmentClusters;

    // programming trick: avoid frequent execution of the new and delete operations
    private CMDPoint m_startPoint1, m_endPoint1, m_startPoint2, m_endPoint2;
    private CMDPoint m_vector1; //  = new CMDPoint(m_document.m_nDimensions);
    private CMDPoint m_vector2; // = new CMDPoint(m_document.m_nDimensions);;
    private CMDPoint m_projectionPoint; // = new CMDPoint( m_document.m_nDimensions);;
    double m_coefficient;
    private STRtree m_tree; //空间索引
//    private HashMap<Integer, Deque<Integer>> m_searchRegion = new HashMap<>(); //存储所有点的空间索引范围

    private ArrayList<LineSegmentId> m_idArray = new ArrayList<ClusterGen.LineSegmentId>();

    public ArrayList<CMDPoint> get_lineSegmentPointArray() {
        return m_lineSegmentPointArray;
    }

    private ArrayList<CMDPoint> m_lineSegmentPointArray = new ArrayList<CMDPoint>();

    // used for performing the DBSCAN algorithm
    public static final int UNCLASSIFIED = -2;
    public static final int NOISE = -1;

    private double MIN_LINESEGMENT_LENGTH = 100; //100 500
    private final double MIN_DIFFERENT_LENGTH = MIN_LINESEGMENT_LENGTH / 1.414;

    public void setMdlCostAdwantage(int mdlCostAdwantage) {
        MDL_COST_ADWANTAGE = mdlCostAdwantage;
    }

    private int MDL_COST_ADWANTAGE = 25; //25 50
    private static final int INT_MAX = Integer.MAX_VALUE;

    public void setMinLinesegmentLength(double minLinesegmentLength) {
        MIN_LINESEGMENT_LENGTH = minLinesegmentLength;
    }

    // used for InsertClusterPoint() and ReplaceClusterPoint()
    enum PointLocation {
        HEAD, TAIL
    }

    class LineSegmentId {

        int trajectoryId;
        int order;
    }

    class CandidateClusterPoint {

        double orderingValue;
        int lineSegmentId;
        boolean startPointFlag;

    }

    class LineSegmentCluster {

        int lineSegmentClusterId;
        int nLineSegments;
        CMDPoint avgDirectionVector;
        double cosTheta, sinTheta;
        ArrayList<CandidateClusterPoint> candidatePointList = new ArrayList<ClusterGen.CandidateClusterPoint>();
//        Queue<CandidateClusterPoint> candidatePointHeap= new PriorityQueue<>(
//                (e1, e2) -> (int) (e1.orderingValue - e2.orderingValue));
        int nClusterPoints;
        ArrayList<CMDPoint> clusterPointArray = new ArrayList<CMDPoint>();
        int nTrajectories;
        ArrayList<Integer> trajectoryIdList = new ArrayList<Integer>();
        boolean enabled;
    }

    // this default constructor should be never used
    public ClusterGen() {

    }

    // use the following constructor instead
    public ClusterGen(TraClusterDoc document) {
        m_document = document;

        m_startPoint1 = new CMDPoint(m_document.m_nDimensions);
        m_startPoint2 = new CMDPoint(m_document.m_nDimensions);
        m_endPoint1 = new CMDPoint(m_document.m_nDimensions);
        m_endPoint2 = new CMDPoint(m_document.m_nDimensions);

        m_vector1 = new CMDPoint(m_document.m_nDimensions);
        m_vector2 = new CMDPoint(m_document.m_nDimensions);
        m_projectionPoint = new CMDPoint(m_document.m_nDimensions);


        m_idArray.clear();
        m_lineSegmentPointArray.clear();
        m_tree = new STRtree(); //R树索引
    }

    public boolean constructCluster() {
        // this step consists of two sub-steps
        // notice that the result of the previous sub-step is used in the following sub-steps
        logger.info("开始构建新的轨迹...\n");
        if (!constructLineSegmentCluster()) {
            return false;
        }
        if (!storeLineSegmentCluster()) {
            return false;
        }
        logger.info("新轨迹构建完成.\n");
        return true;
    }

    public boolean partitionTrajectory() {

        for (int i = 0; i < m_document.m_trajectoryList.size(); i++) {
            Trajectory pTrajectory = m_document.m_trajectoryList.get(i);

            findOptimalPartition(pTrajectory);

            m_document.m_trajectoryList.set(i, pTrajectory);
        }
        if (!storeClusterComponentIntoIndex()) {
            return false;
        }
        return true;
    }

    public boolean performDBSCAN(double eps, int minLns) {

        m_epsParam = eps;
        m_minLnsParam = minLns;

        m_currComponentId = 0;

//		createRtree();

        for (int i = 0; i < m_nTotalLineSegments; i++) {
            m_componentIdArray.add(UNCLASSIFIED);
        }

        for (int i = 0; i < m_nTotalLineSegments; i++) {
            if (i % 100 == 0)
                logger.debug("processed(cluster) " + i + "...");
            if (m_componentIdArray.get(i) == UNCLASSIFIED && expandDenseComponent2(i, m_currComponentId, eps, minLns)) {
                m_currComponentId++;
            }
        }
        logger.info("聚类后的类别数目是:"+m_currComponentId);
        return true;
    }


    private boolean storeClusterComponentIntoIndex() {

        int nDimensions = m_document.m_nDimensions;
        CMDPoint startPoint;
        CMDPoint endPoint;

//        ArrayList<CMDPoint> cmdpt = new ArrayList<CMDPoint>();
        m_nTotalLineSegments = 0;
        for (int i = 0; i < m_document.m_trajectoryList.size(); i++) {
            Trajectory pTrajectory = m_document.m_trajectoryList.get(i);
//            if(m_nTotalLineSegments==115){
//                System.out.println("test");
//            }

            for (int j = 0; j < pTrajectory.getM_nPartitionPoints() - 1; j++) {
                // convert an n-dimensional line segment into a 2n-dimensional point
                // i.e., the first n-dimension: the start point
                //       the last n-dimension: the end point
                startPoint = pTrajectory.getM_partitionPointArray().get(j);
                endPoint = pTrajectory.getM_partitionPointArray().get(j + 1);

                // 太短的线段过滤掉
                if (measureDistanceFromPointToPoint(startPoint, endPoint) < MIN_LINESEGMENT_LENGTH) {
                    continue;
                }

                CMDPoint lineSegmentPoint = new CMDPoint(nDimensions * 2);
                for (int m = 0; m < nDimensions; m++) {
                    lineSegmentPoint.setM_coordinate(m, startPoint.getM_coordinate(m));
                    lineSegmentPoint.setM_coordinate(nDimensions + m, endPoint.getM_coordinate(m));
                }

                LineSegmentId id = new LineSegmentId();
                id.trajectoryId = pTrajectory.getM_trajectoryId();
                id.order = j;

                m_idArray.add(id);
                m_lineSegmentPointArray.add(lineSegmentPoint);

//                if(cmdpt.size() ==0) cmdpt.add(lineSegmentPoint);
//                for(int k = 0; k < cmdpt.size(); k++){
//                    if(cmdpt.get(k).getM_coordinate(0) != startPoint.getM_coordinate(0) ||
//                            cmdpt.get(k).getM_coordinate(1) != startPoint.getM_coordinate(1) ||
//                            cmdpt.get(k).getM_coordinate(2) != endPoint.getM_coordinate(0) ||
//                            cmdpt.get(k).getM_coordinate(3) != endPoint.getM_coordinate(1)){
//                        cmdpt.add(lineSegmentPoint);
//                        break;
//                    }else{
//                        System.out.println(m_nTotalLineSegments);
//                        break;
//                    }
//                }

                insertRtree(lineSegmentPoint, m_nTotalLineSegments); //加入索引树

                m_nTotalLineSegments++;
            }
        }

        logger.info("分段后的线段数目是:" + m_nTotalLineSegments);
        m_tree.build();

//        long startTime=System.currentTimeMillis();
//        for(int i = 0; i < m_nTotalLineSegments; i++){
//            Deque<Integer> seeds = new ArrayDeque<>();
//            extractStartAndEndPoints(i, m_startPoint1, m_endPoint1);
//            computeEPSNeighborhoodByRtree2(m_startPoint1, m_endPoint1, eps, seeds);
//            m_searchRegion.put(i, seeds);
//            System.out.println(i);
//        }
//        long endTime=System.currentTimeMillis();
//        System.out.println("建立索引花费时间" + (endTime-startTime) / 1000);
        logger.info("R树索引创建成功！");

        return true;
    }

    private void findOptimalPartition(Trajectory pTrajectory) {

        int nPoints = pTrajectory.getM_nPoints();
        int startIndex = 0, length;
        int fullPartitionMDLCost, partialPartitionMDLCost;

        // add the start point of a trajectory
        CMDPoint startP = pTrajectory.getM_pointArray().get(0);
        pTrajectory.addPartitionPointToArray(startP);

        for (; ; ) {
            fullPartitionMDLCost = partialPartitionMDLCost = 0;

            for (length = 1; startIndex + length < nPoints; length++) {
                // compute the total length of a trajectory
                fullPartitionMDLCost += computeModelCost(pTrajectory, startIndex + length - 1, startIndex + length);

//				fullPartitionMDLCost += computeModelCost(pTrajectory, startIndex, startIndex + length);

                // compute the sum of (1) the length of a cluster component and
                // 					 (2) the perpendicular and angle distances
                partialPartitionMDLCost = computeModelCost(pTrajectory, startIndex, startIndex + length) +
                        computeEncodingCost(pTrajectory, startIndex, startIndex + length);

                if (fullPartitionMDLCost + MDL_COST_ADWANTAGE < partialPartitionMDLCost) {

                    pTrajectory.addPartitionPointToArray(pTrajectory.getM_pointArray().get(startIndex + length - 1));
                    startIndex = startIndex + length - 1;
                    length = 0;
                    break;
                }
            }
            // if we reach at the end of a trajectory
            if (startIndex + length >= nPoints) {
                break;
            }
        }

        // add the end point of a trajectory
        pTrajectory.addPartitionPointToArray(pTrajectory.getM_pointArray().get(nPoints - 1));

        return;
    }

    private double LOG2(double x) {
        return Math.log(x) / Math.log(2);
    }

    private int computeModelCost(Trajectory pTrajectory, int startPIndex, int endPIndex) {

        CMDPoint lineSegmentStart = pTrajectory.getM_pointArray().get(startPIndex);
        CMDPoint lineSegmentEnd = pTrajectory.getM_pointArray().get(endPIndex);

        double distance = measureDistanceFromPointToPoint(lineSegmentStart, lineSegmentEnd);

        if (distance < 1.0) {
            distance = 1.0;        // to take logarithm
        }

        return (int) Math.ceil(LOG2(distance));

    }

    private int computeEncodingCost(Trajectory pTrajectory, int startPIndex, int endPIndex) {

        CMDPoint clusterComponentStart;
        CMDPoint clusterComponentEnd;
        CMDPoint lineSegmentStart;
        CMDPoint lineSegmentEnd;
        double perpendicularDistance;
        double angleDistance;
        int encodingCost = 0;

        clusterComponentStart = pTrajectory.getM_pointArray().get(startPIndex);
        clusterComponentEnd = pTrajectory.getM_pointArray().get(endPIndex);

        for (int i = startPIndex; i < endPIndex; i++) {
            lineSegmentStart = pTrajectory.getM_pointArray().get(i);
            lineSegmentEnd = pTrajectory.getM_pointArray().get(i + 1);

            perpendicularDistance = measurePerpendicularDistance(clusterComponentStart,
                    clusterComponentEnd, lineSegmentStart, lineSegmentEnd);
            angleDistance = measureAngleDisntance(clusterComponentStart,
                    clusterComponentEnd, lineSegmentStart, lineSegmentEnd);

            if (perpendicularDistance < 1.0) perpendicularDistance = 1.0;    //  to take logarithm
            if (angleDistance < 1.0) angleDistance = 1.0;                    //  to take logarithm

            encodingCost += ((int) Math.ceil(LOG2(perpendicularDistance)) + (int) Math.ceil(LOG2(angleDistance)));
        }
        return encodingCost;

    }

    private double measurePerpendicularDistance(CMDPoint s1, CMDPoint e1, CMDPoint s2, CMDPoint e2) {

        //  we assume that the first line segment is longer than the second one
        double distance1;    //  the distance from a start point to the cluster component
        double distance2;    //  the distance from an end point to the cluster component

        distance1 = measureDistanceFromPointToLineSegment(s1, e1, s2);
        distance2 = measureDistanceFromPointToLineSegment(s1, e1, e2);

        //  if the first line segment is exactly the same as the second one,
        //  the perpendicular distance should be zero
        if (distance1 == 0.0 && distance2 == 0.0) return 0.0;

        //  return (d1^2 + d2^2) / (d1 + d2) as the perpendicular distance
        return ((Math.pow(distance1, 2) + Math.pow(distance2, 2)) / (distance1 + distance2));

    }

    private double measureDistanceFromPointToLineSegment(CMDPoint s, CMDPoint e, CMDPoint p) {

        int nDimensions = p.getM_nDimensions();

        //  NOTE: the variables m_vector1 and m_vector2 are declared as member variables

        //  construct two vectors as follows
        //  1. the vector connecting the start point of the cluster component and a given point
        //  2. the vector representing the cluster component
        for (int i = 0; i < nDimensions; i++) {
            m_vector1.setM_coordinate(i, p.getM_coordinate(i) - s.getM_coordinate(i));
            m_vector2.setM_coordinate(i, e.getM_coordinate(i) - s.getM_coordinate(i));
        }

        //  a coefficient (0 <= b <= 1)
        m_coefficient = computeInnerProduct(m_vector1, m_vector2) / computeInnerProduct(m_vector2, m_vector2);

        //  the projection on the cluster component from a given point
        //  NOTE: the variable m_projectionPoint is declared as a member variable

        for (int i = 0; i < nDimensions; i++) {
            m_projectionPoint.setM_coordinate(i, s.getM_coordinate(i) + m_coefficient * m_vector2.getM_coordinate(i));
        }

        //  return the distance between the projection point and the given point
        return measureDistanceFromPointToPoint(p, m_projectionPoint);

    }

    private double measureDistanceFromPointToPoint(CMDPoint point1, CMDPoint point2) {

        int nDimensions = point1.getM_nDimensions();
        double squareSum = 0.0;

        for (int i = 0; i < nDimensions; i++) {
            squareSum += Math.pow((point2.getM_coordinate(i) - point1.getM_coordinate(i)), 2);
        }
        return Math.sqrt(squareSum);

    }

    private double computeVectorLength(CMDPoint vector) {

        int nDimensions = vector.getM_nDimensions();
        double squareSum = 0.0;

        for (int i = 0; i < nDimensions; i++) {
            squareSum += Math.pow(vector.getM_coordinate(i), 2);
        }

        return Math.sqrt(squareSum);
    }

    private double computeInnerProduct(CMDPoint vector1, CMDPoint vector2) {
        int nDimensions = vector1.getM_nDimensions();
        double innerProduct = 0.0;

        for (int i = 0; i < nDimensions; i++) {
            innerProduct += (vector1.getM_coordinate(i) * vector2.getM_coordinate(i));
        }

        return innerProduct;
    }

    private double measureAngleDisntance(CMDPoint s1, CMDPoint e1, CMDPoint s2, CMDPoint e2) {

        int nDimensions = s1.getM_nDimensions();

        //  NOTE: the variables m_vector1 and m_vector2 are declared as member variables
        //  construct two vectors representing the cluster component and a line segment, respectively
        for (int i = 0; i < nDimensions; i++) {
            m_vector1.setM_coordinate(i, e1.getM_coordinate(i) - s1.getM_coordinate(i));
            m_vector2.setM_coordinate(i, e2.getM_coordinate(i) - s2.getM_coordinate(i));
        }

        //  we assume that the first line segment is longer than the second one
        //  i.e., vectorLength1 >= vectorLength2
        double vectorLength1 = computeVectorLength(m_vector1);
        double vectorLength2 = computeVectorLength(m_vector2);

        //  if one of two vectors is a point, the angle distance becomes zero
        if (vectorLength1 == 0.0 || vectorLength2 == 0.0) return 0.0;

        //  compute the inner product of the two vectors
        double innerProduct = computeInnerProduct(m_vector1, m_vector2);

        //  compute the angle between two vectors by using the inner product
        double cosTheta = innerProduct / (vectorLength1 * vectorLength2);
        //  compensate the computation error (e.g., 1.00001)
        //  cos(theta) should be in the range [-1.0, 1.0]
        //  START ...
        if (cosTheta > 1.0) cosTheta = 1.0;
        if (cosTheta < -1.0) cosTheta = -1.0;
        //  ... END
        double sinTheta = Math.sqrt(1 - Math.pow(cosTheta, 2));
        //  if 90 <= theta <= 270, the angle distance becomes the length of the line segment
        //  if (cosTheta < -1.0) sinTheta = 1.0;

        return (vectorLength2 * sinTheta);
    }

    private boolean expandDenseComponent(int index, int componentId, double eps, int minDensity) {
//        if (index % 100 == 0)
        System.out.println("processed(cluster) " + index + "...");

//        if(index==43){
//            System.out.println("slow");
//        }

        Set<Integer> seeds = new HashSet<Integer>();
        Set<Integer> seedResult = new HashSet<Integer>();

        int currIndex;

        extractStartAndEndPoints(index, m_startPoint1, m_endPoint1);
        computeEPSNeighborhood(m_startPoint1, m_endPoint1, eps, seeds); //computeEPSNeighborhoodByRtree computeEPSNeighborhood
//        seeds = m_searchRegion.get(index);

        if ((int) seeds.size() < minDensity) { //  not a core line segment
            m_componentIdArray.set(index, NOISE);
            return false;
        }
        // else...
        for (int i = 0; i < seeds.size(); i++) {
            m_componentIdArray.set((Integer) (seeds.toArray()[i]), componentId);
        }
        seeds.remove(index);
        while (!seeds.isEmpty()) {
            currIndex = (Integer) seeds.toArray()[0];

            extractStartAndEndPoints(currIndex, m_startPoint1, m_endPoint1);
            computeEPSNeighborhood(m_startPoint1, m_endPoint1, eps, seedResult);
//            seedResult = m_searchRegion.get(currIndex);

            if ((int) seedResult.size() >= minDensity) {
                for (int iter = 0; iter < seedResult.size(); iter++) {
                    if (m_componentIdArray.get((Integer) (seedResult.toArray()[iter])) == UNCLASSIFIED ||
                            m_componentIdArray.get((Integer) (seedResult.toArray()[iter])) == NOISE) {
                        if (m_componentIdArray.get((Integer) (seedResult.toArray()[iter])) == UNCLASSIFIED) {
                            seeds.add((Integer) (seedResult.toArray()[iter]));
                        }
                        m_componentIdArray.set((Integer) (seedResult.toArray()[iter]), componentId);
                    }
                }
            }

            seeds.remove(currIndex);
        }

        return true;
    }

    private boolean expandDenseComponent2(int index, int componentId, double eps, int minDensity) {

//        if(index==43){
//            System.out.println("slow");
//        }

        Deque<Integer> seeds = new ArrayDeque<Integer>();
        Deque<Integer> seedResult = new ArrayDeque<Integer>();

        int currIndex;

        extractStartAndEndPoints(index, m_startPoint1, m_endPoint1);
        computeEPSNeighborhoodByRtree(m_startPoint1, m_endPoint1, eps, seeds); //computeEPSNeighborhoodByRtree computeEPSNeighborhood
//        seeds = m_searchRegion.get(index);

        if (seeds.size() < minDensity) { //  not a core line segment
            m_componentIdArray.set(index, NOISE);
            return false;
        }
        // else...
        for (int i = 0; i < seeds.size(); i++) {
            m_componentIdArray.set((Integer) (seeds.toArray()[i]), componentId);
        }
        seeds.remove(index);
        while (!seeds.isEmpty()) {
            currIndex = seeds.getFirst();

            extractStartAndEndPoints(currIndex, m_startPoint1, m_endPoint1);
            computeEPSNeighborhoodByRtree(m_startPoint1, m_endPoint1, eps, seedResult);
//            seedResult = m_searchRegion.get(currIndex);

            if (seedResult.size() >= minDensity) {
                for(Integer iter: seedResult){
                    if (m_componentIdArray.get(iter) == UNCLASSIFIED ||
                            m_componentIdArray.get(iter) == NOISE) {
                        if (m_componentIdArray.get(iter) == UNCLASSIFIED) {
                            seeds.addLast(iter);
                        }
                        m_componentIdArray.set(iter, componentId);
                    }
                }
            }

            seeds.removeFirst();
        }

        return true;
    }

    boolean constructLineSegmentCluster() {
        int nDimensions = m_document.m_nDimensions;
        m_lineSegmentClusters = new LineSegmentCluster[m_currComponentId];

        //  initialize the list of line segment clusters
        //  START ...
        logger.info("开始初始化聚类线段集...");
        for (int i = 0; i < m_currComponentId; i++) {
            m_lineSegmentClusters[i] = new LineSegmentCluster();
            m_lineSegmentClusters[i].avgDirectionVector = new CMDPoint(nDimensions);
            m_lineSegmentClusters[i].lineSegmentClusterId = i;
            m_lineSegmentClusters[i].nLineSegments = 0;
            m_lineSegmentClusters[i].nClusterPoints = 0;
            m_lineSegmentClusters[i].nTrajectories = 0;
            m_lineSegmentClusters[i].enabled = false;
        }
        //  ... END
        logger.info("完成初始化聚类线段集.");

        //  accumulate the direction vector of a line segment
        logger.info("开始计算线段的方向向量...");
        for (int i = 0; i < m_nTotalLineSegments; i++) {
            int componentId = m_componentIdArray.get(i);
            if (componentId >= 0) {
//                if(componentId == 0){
//                    System.out.println("test");
//                }
                for (int j = 0; j < nDimensions; j++) {
                    double difference = m_lineSegmentPointArray.get(i).getM_coordinate(nDimensions + j)
                            - m_lineSegmentPointArray.get(i).getM_coordinate(j);
                    double currSum = m_lineSegmentClusters[componentId].avgDirectionVector.getM_coordinate(j)
                            + difference;
                    m_lineSegmentClusters[componentId].avgDirectionVector.setM_coordinate(j, currSum);
                }
                m_lineSegmentClusters[componentId].nLineSegments++;
            }
        }
        logger.info("完成计算线段的方向向量.");

        //  compute the average direction vector of a line segment cluster
        //  START ...
        logger.info("开始计算聚类线段集的平均方向向量...");
        double vectorLength1, vectorLength2, innerProduct;
        double cosTheta, sinTheta;

        m_vector2.setM_coordinate(0, 1.0);
        m_vector2.setM_coordinate(1, 0.0);

        for (int i = 0; i < m_currComponentId; i++) {
            LineSegmentCluster clusterEntry = m_lineSegmentClusters[i];

            for (int j = 0; j < nDimensions; j++) {
                clusterEntry.avgDirectionVector.setM_coordinate(j, clusterEntry.avgDirectionVector.getM_coordinate(j) / (double) clusterEntry.nLineSegments);
            }
            vectorLength1 = computeVectorLength(clusterEntry.avgDirectionVector);
            vectorLength2 = 1.0;

            innerProduct = computeInnerProduct(clusterEntry.avgDirectionVector, m_vector2);
            cosTheta = innerProduct / (vectorLength1 * vectorLength2);
            if (cosTheta > 1.0) cosTheta = 1.0;
            if (cosTheta < -1.0) cosTheta = -1.0;
            sinTheta = Math.sqrt(1 - Math.pow(cosTheta, 2));

            if (clusterEntry.avgDirectionVector.getM_coordinate(1) < 0) {
                sinTheta = -sinTheta;
            }

            clusterEntry.cosTheta = cosTheta;
            clusterEntry.sinTheta = sinTheta;

        }
        logger.info("完成计算聚类线段集的平均方向向量.");
        //  ... END

        //  summarize the information about line segment clusters
        //  the structure for summarization is as follows
        //  [lineSegmentClusterId, nClusterPoints, clusterPointArray, nTrajectories, { trajectoryId, ... }]
        logger.info("开始汇总聚类线段集的信息...");
        for (int i = 0; i < m_nTotalLineSegments; i++) {
            if (m_componentIdArray.get(i) >= 0)        //  if the componentId < 0, it is a noise
                RegisterAndUpdateLineSegmentCluster2(m_componentIdArray.get(i), i);
        }
        logger.info("完成汇总聚类线段集的信息.");

        logger.info("开始计算表达线段...");
        Set<Integer> trajectories = new HashSet<Integer>();
        for (int i = 0; i < m_currComponentId; i++) {
            LineSegmentCluster clusterEntry = (m_lineSegmentClusters[i]);

            //  a line segment cluster must have trajectories more than the minimum threshold
            if (clusterEntry.nTrajectories >= m_minLnsParam) {
                clusterEntry.enabled = true;
                // m_lineSegmentClusters[i].enabled = true;
                //  DEBUG: count the number of trajectories that belong to clusters
//                for (int j = 0; j < (int) clusterEntry.trajectoryIdList.size(); j++) {
//                    trajectories.add(clusterEntry.trajectoryIdList.get(j));
//                }

                computeRepresentativeLines2(clusterEntry);
                // computeRepresentativeLines(m_lineSegmentClusters[i]);
            } else {
                clusterEntry.candidatePointList.clear();
                clusterEntry.clusterPointArray.clear();
                clusterEntry.trajectoryIdList.clear();
            }
            System.out.println(i);
//            if(i == 19){
//                System.out.println("slow");
//            }
        }
        logger.info("完成计算表达线段.");
        //  DEBUG: compute the ratio of trajectories that belong to clusters
//        m_document.m_clusterRatio = (double) trajectories.size() / (double) m_document.m_nTrajectories;

        return true;
    }


    private void computeRepresentativeLines(LineSegmentCluster clusterEntry) {

        Set<Integer> lineSegments = new LinkedHashSet<Integer>();
        Set<Integer> insertionList = new LinkedHashSet<Integer>();
        Set<Integer> deletionList = new LinkedHashSet<Integer>();

        int iter = 0;
        CandidateClusterPoint candidatePoint, nextCandidatePoint;
        double prevOrderingValue = 0.0;

        int nClusterPoints = 0;
        lineSegments.clear();

        //  sweep the line segments in a line segment cluster

        while (iter != (clusterEntry.candidatePointList.size() - 1) && clusterEntry.candidatePointList.size() > 0) {
            insertionList.clear();
            deletionList.clear();

            do {
                candidatePoint = clusterEntry.candidatePointList.get(iter);
                iter++;
                //  check whether this line segment has begun or not
                if (!lineSegments.contains(candidatePoint.lineSegmentId)) {
                    // iter1 = lineSegments.find(candidatePoint.lineSegmentId);
                    // if (iter1 == lineSegments.end())	{				//  if there is no matched element,
                    insertionList.add(candidatePoint.lineSegmentId);        //  this line segment begins at this point
                    lineSegments.add(candidatePoint.lineSegmentId);
                } else {                        //  if there is a matched element,
                    deletionList.add(candidatePoint.lineSegmentId);        //  this line segment ends at this point
                }
                //  check whether the next line segment begins or ends at the same point
                if (iter != (clusterEntry.candidatePointList.size() - 1)) {
                    nextCandidatePoint = clusterEntry.candidatePointList.get(iter);
                } else {
                    break;
                }
            } while (candidatePoint.orderingValue == nextCandidatePoint.orderingValue);

            //  check if a line segment is connected to another line segment in the same trajectory
            //  if so, delete one of the line segments to remove duplicates
            // for (iter2 = insertionList.begin(); iter2 != insertionList.end(); iter2++)
            for (int iter2 = 0; iter2 < insertionList.size(); iter2++) {
                for (int iter3 = 0; iter3 < deletionList.size(); iter3++) {
                    int a = (Integer) (insertionList.toArray()[iter2]);
                    int b = (Integer) (deletionList.toArray()[iter3]);
                    if (a == b) {
                        lineSegments.remove((Integer) (deletionList.toArray()[iter3]));
                        deletionList.remove((Integer) (deletionList.toArray()[iter3]));
                        break;
                    }
                }

                for (int iter3 = 0; iter3 < deletionList.size(); iter3++) {
                    if (m_idArray.get((Integer) (insertionList.toArray()[iter2])).trajectoryId
                            == m_idArray.get((Integer) (deletionList.toArray()[iter3])).trajectoryId) {
                        lineSegments.remove((Integer) (deletionList.toArray()[iter3]));
                        deletionList.remove((Integer) (deletionList.toArray()[iter3]));
                        break;
                    }
                }
            }

            // if the current density exceeds a given threshold
            if ((int) (lineSegments.size()) >= m_minLnsParam) {
                if (Math.abs(candidatePoint.orderingValue - prevOrderingValue) > MIN_DIFFERENT_LENGTH) {
                    computeAndRegisterClusterPoint(clusterEntry, candidatePoint.orderingValue, lineSegments);
                    prevOrderingValue = candidatePoint.orderingValue;
                    nClusterPoints++;
                }
            }

            //  delete the line segment that is not connected to another line segment
            for (int iter3 = 0; iter3 < deletionList.size(); iter3++) {
                lineSegments.remove((Integer) (deletionList.toArray()[iter3]));
            }
        }

        if (nClusterPoints >= 2) {
            clusterEntry.nClusterPoints = nClusterPoints;
        } else {
            //  there is no representative trend in this line segment cluster
            clusterEntry.enabled = false;
            clusterEntry.candidatePointList.clear();
            clusterEntry.clusterPointArray.clear();
            clusterEntry.trajectoryIdList.clear();
        }
        return;
    }

    private void computeRepresentativeLines2(LineSegmentCluster clusterEntry) {

        Set<Integer> lineSegments = new LinkedHashSet<Integer>();
        Set<Integer> insertionList = new LinkedHashSet<Integer>();
        Set<Integer> deletionList = new LinkedHashSet<Integer>();

        int iter = 0;
        CandidateClusterPoint candidatePoint, nextCandidatePoint;
        double prevOrderingValue = 0.0;

        int nClusterPoints = 0;
        lineSegments.clear();

        //  sweep the line segments in a line segment cluster
        clusterEntry.candidatePointList.sort((o1, o2) -> {
            if (Math.abs(o1.orderingValue - o2.orderingValue) < .00001) //排序时候避免浮点型数导致的Bug
                return 0;
            else return Double.compare(o1.orderingValue, o2.orderingValue);
        });
//        clusterEntry.candidatePointList.sort((o1,o2)->(int)(o1.orderingValue-o2.orderingValue));

        while (iter != (clusterEntry.candidatePointList.size() - 1) && clusterEntry.candidatePointList.size() > 0) {
            insertionList.clear();
            deletionList.clear();

            do {
                candidatePoint = clusterEntry.candidatePointList.get(iter);
                iter++;
                //  check whether this line segment has begun or not
                if (!lineSegments.contains(candidatePoint.lineSegmentId)) {
                    // iter1 = lineSegments.find(candidatePoint.lineSegmentId);
                    // if (iter1 == lineSegments.end())	{				//  if there is no matched element,
                    insertionList.add(candidatePoint.lineSegmentId);        //  this line segment begins at this point
                    lineSegments.add(candidatePoint.lineSegmentId);
                } else {                        //  if there is a matched element,
                    deletionList.add(candidatePoint.lineSegmentId);        //  this line segment ends at this point
                }
                //  check whether the next line segment begins or ends at the same point
                if (iter != (clusterEntry.candidatePointList.size() - 1)) {
                    nextCandidatePoint = clusterEntry.candidatePointList.get(iter);
                } else {
                    break;
                }
            } while (candidatePoint.orderingValue == nextCandidatePoint.orderingValue);

            //  check if a line segment is connected to another line segment in the same trajectory
            //  if so, delete one of the line segments to remove duplicates
            // for (iter2 = insertionList.begin(); iter2 != insertionList.end(); iter2++)
            for (int iter2 = 0; iter2 < insertionList.size(); iter2++) {
                for (int iter3 = 0; iter3 < deletionList.size(); iter3++) {
                    int a = (Integer) (insertionList.toArray()[iter2]);
                    int b = (Integer) (deletionList.toArray()[iter3]);
                    if (a == b) {
                        lineSegments.remove((Integer) (deletionList.toArray()[iter3]));
                        deletionList.remove((Integer) (deletionList.toArray()[iter3]));
                        break;
                    }
                }

                for (int iter3 = 0; iter3 < deletionList.size(); iter3++) {
                    if (m_idArray.get((Integer) (insertionList.toArray()[iter2])).trajectoryId
                            == m_idArray.get((Integer) (deletionList.toArray()[iter3])).trajectoryId) {
                        lineSegments.remove((Integer) (deletionList.toArray()[iter3]));
                        deletionList.remove((Integer) (deletionList.toArray()[iter3]));
                        break;
                    }
                }
            }

            // if the current density exceeds a given threshold
            if (lineSegments.size() >= m_minLnsParam) {
                if (Math.abs(candidatePoint.orderingValue - prevOrderingValue) > MIN_DIFFERENT_LENGTH) {
                    computeAndRegisterClusterPoint(clusterEntry, candidatePoint.orderingValue, lineSegments);
                    prevOrderingValue = candidatePoint.orderingValue;
                    nClusterPoints++;
                }
            }

            //  delete the line segment that is not connected to another line segment
            for (int iter3 = 0; iter3 < deletionList.size(); iter3++) {
                lineSegments.remove((Integer) (deletionList.toArray()[iter3]));
            }
        }

        if (nClusterPoints >= 2) {
            clusterEntry.nClusterPoints = nClusterPoints;
        } else {
            //  there is no representative trend in this line segment cluster
            clusterEntry.enabled = false;
            clusterEntry.candidatePointList.clear();
            clusterEntry.clusterPointArray.clear();
            clusterEntry.trajectoryIdList.clear();
        }
        return;
    }

    private void computeAndRegisterClusterPoint(
            LineSegmentCluster clusterEntry,
            double currValue,
            Set<Integer> lineSegments) {
        int nDimensions = m_document.m_nDimensions;
        int nLineSegmentsInSet = (int) (lineSegments.size());
        CMDPoint clusterPoint = new CMDPoint(nDimensions);
        CMDPoint sweepPoint = new CMDPoint(nDimensions);

        for (int iter = 0; iter < lineSegments.size(); iter++) {
            // get the sweep point of each line segment
            // this point is parallel to the current value of the sweeping direction
            getSweepPointOfLineSegment(clusterEntry, currValue,
                    (Integer) (lineSegments.toArray()[iter]), sweepPoint);
            for (int i = 0; i < nDimensions; i++) {
                clusterPoint.setM_coordinate(i, clusterPoint.getM_coordinate(i) +
                        (sweepPoint.getM_coordinate(i) / (double) nLineSegmentsInSet));
            }
        }

        // NOTE: this program code works only for the 2-dimensional data
        double origX, origY;
        origX = GET_X_REV_ROTATION(clusterPoint.getM_coordinate(0), clusterPoint.getM_coordinate(1), clusterEntry.cosTheta, clusterEntry.sinTheta);
        origY = GET_Y_REV_ROTATION(clusterPoint.getM_coordinate(0), clusterPoint.getM_coordinate(1), clusterEntry.cosTheta, clusterEntry.sinTheta);
        clusterPoint.setM_coordinate(0, origX);
        clusterPoint.setM_coordinate(1, origY);

        // register the obtained cluster point (i.e., the average of all the sweep points)
        clusterEntry.clusterPointArray.add(clusterPoint);

        return;
    }

    private void getSweepPointOfLineSegment(LineSegmentCluster clusterEntry,
                                            double currValue, int lineSegmentId, CMDPoint sweepPoint) {

        CMDPoint lineSegmentPoint = m_lineSegmentPointArray.get(lineSegmentId);        //  2n-dimensional point
        double coefficient;

        //  NOTE: this program code works only for the 2-dimensional data
        double newStartX, newEndX, newStartY, newEndY;
        newStartX = GET_X_ROTATION(lineSegmentPoint.getM_coordinate(0), lineSegmentPoint.getM_coordinate(1), clusterEntry.cosTheta, clusterEntry.sinTheta);
        newEndX = GET_X_ROTATION(lineSegmentPoint.getM_coordinate(2), lineSegmentPoint.getM_coordinate(3), clusterEntry.cosTheta, clusterEntry.sinTheta);
        newStartY = GET_Y_ROTATION(lineSegmentPoint.getM_coordinate(0), lineSegmentPoint.getM_coordinate(1), clusterEntry.cosTheta, clusterEntry.sinTheta);
        newEndY = GET_Y_ROTATION(lineSegmentPoint.getM_coordinate(2), lineSegmentPoint.getM_coordinate(3), clusterEntry.cosTheta, clusterEntry.sinTheta);

        coefficient = (currValue - newStartX) / (newEndX - newStartX);
        sweepPoint.setM_coordinate(0, currValue);
        sweepPoint.setM_coordinate(1, newStartY + coefficient * (newEndY - newStartY));

        return;
    }


    private double GET_X_ROTATION(double _x, double _y, double _cos, double _sin) {
        return ((_x) * (_cos) + (_y) * (_sin));
    }

    private double GET_Y_ROTATION(double _x, double _y, double _cos, double _sin) {
        return (-(_x) * (_sin) + (_y) * (_cos));
    }

    private double GET_X_REV_ROTATION(double _x, double _y, double _cos, double _sin) {
        return ((_x) * (_cos) - (_y) * (_sin));
    }

    private double GET_Y_REV_ROTATION(double _x, double _y, double _cos, double _sin) {
        return ((_x) * (_sin) + (_y) * (_cos));
    }

    private void RegisterAndUpdateLineSegmentCluster(int componentId, int lineSegmentId) {
        LineSegmentCluster clusterEntry = m_lineSegmentClusters[componentId];

        //  the start and end values of the first dimension (e.g., the x value in the 2-dimension)
        //  NOTE: this program code works only for the 2-dimensional data

        CMDPoint aLineSegment = m_lineSegmentPointArray.get(lineSegmentId);
        double orderingValue1 = GET_X_ROTATION(aLineSegment.getM_coordinate(0),
                aLineSegment.getM_coordinate(1), clusterEntry.cosTheta, clusterEntry.sinTheta);
        double orderingValue2 = GET_X_ROTATION(aLineSegment.getM_coordinate(2),
                aLineSegment.getM_coordinate(3), clusterEntry.cosTheta, clusterEntry.sinTheta);

        CandidateClusterPoint existingCandidatePoint, newCandidatePoint1, newCandidatePoint2;
        int i, j;
        //  sort the line segment points by the coordinate of the first dimension
        //  simply use the insertion sort algorithm
        //  START ...
        int iter1 = 0;
        for (i = 0; i < (int) clusterEntry.candidatePointList.size(); i++) {
            existingCandidatePoint = clusterEntry.candidatePointList.get(iter1);
            if (existingCandidatePoint.orderingValue >= orderingValue1) {
                break;
            }
            iter1++;
        }
        newCandidatePoint1 = new CandidateClusterPoint();

        newCandidatePoint1.orderingValue = orderingValue1;
        newCandidatePoint1.lineSegmentId = lineSegmentId;
        newCandidatePoint1.startPointFlag = true;
        if (i == 0) {
            clusterEntry.candidatePointList.add(0, newCandidatePoint1);
        } else if (i >= clusterEntry.candidatePointList.size()) {
            clusterEntry.candidatePointList.add(newCandidatePoint1);
        } else {
            clusterEntry.candidatePointList.add(iter1, newCandidatePoint1);
        }
        int iter2 = 0;
        for (j = 0; j < clusterEntry.candidatePointList.size(); j++) {
            existingCandidatePoint = clusterEntry.candidatePointList.get(iter2);
            if (existingCandidatePoint.orderingValue >= orderingValue2) {
                break;
            }
            iter2++;
        }

        newCandidatePoint2 = new CandidateClusterPoint();
        newCandidatePoint2.orderingValue = orderingValue2;
        newCandidatePoint2.lineSegmentId = lineSegmentId;
        newCandidatePoint2.startPointFlag = false;

        if (j == 0) {
            clusterEntry.candidatePointList.add(0, newCandidatePoint2);
        } else if (j >= (int) clusterEntry.candidatePointList.size()) {
            clusterEntry.candidatePointList.add(newCandidatePoint2);
        } else {
            clusterEntry.candidatePointList.add(iter2, newCandidatePoint2);
        }
        //  ... END

        int trajectoryId = m_idArray.get(lineSegmentId).trajectoryId;

        //  store the identifier of the trajectories that belong to this line segment cluster
        if (!clusterEntry.trajectoryIdList.contains(trajectoryId)) {
            clusterEntry.trajectoryIdList.add(trajectoryId);
            clusterEntry.nTrajectories++;
        }
    }

    private void RegisterAndUpdateLineSegmentCluster2(int componentId, int lineSegmentId) {
        LineSegmentCluster clusterEntry = m_lineSegmentClusters[componentId];

        //  the start and end values of the first dimension (e.g., the x value in the 2-dimension)
        //  NOTE: this program code works only for the 2-dimensional data

        CMDPoint aLineSegment = m_lineSegmentPointArray.get(lineSegmentId);
        double orderingValue1 = GET_X_ROTATION(aLineSegment.getM_coordinate(0),
                aLineSegment.getM_coordinate(1), clusterEntry.cosTheta, clusterEntry.sinTheta);
        double orderingValue2 = GET_X_ROTATION(aLineSegment.getM_coordinate(2),
                aLineSegment.getM_coordinate(3), clusterEntry.cosTheta, clusterEntry.sinTheta);

        CandidateClusterPoint existingCandidatePoint, newCandidatePoint1, newCandidatePoint2;
        int i, j;
        //  插入的时候不排序，后面用的时候再排序
        //  START ...

//        if(componentId == 0){
//            System.out.println("0");
//        }

        newCandidatePoint1 = new CandidateClusterPoint(); //线段起始点
        newCandidatePoint1.orderingValue = orderingValue1;
        newCandidatePoint1.lineSegmentId = lineSegmentId;
        newCandidatePoint1.startPointFlag = true;

        clusterEntry.candidatePointList.add(newCandidatePoint1);

        newCandidatePoint2 = new CandidateClusterPoint();
        newCandidatePoint2.orderingValue = orderingValue2;
        newCandidatePoint2.lineSegmentId = lineSegmentId;
        newCandidatePoint2.startPointFlag = false;

        clusterEntry.candidatePointList.add(newCandidatePoint2);
        //  ... END

        int trajectoryId = m_idArray.get(lineSegmentId).trajectoryId;

        //  store the identifier of the trajectories that belong to this line segment cluster
        if (!clusterEntry.trajectoryIdList.contains(trajectoryId)) {
            clusterEntry.trajectoryIdList.add(trajectoryId);
            clusterEntry.nTrajectories++;
        }
    }

    private void computeEPSNeighborhood(CMDPoint startPoint, CMDPoint endPoint, double eps, Set<Integer> result) {
        result.clear();
        for (int j = 0; j < m_nTotalLineSegments; j++) {
            extractStartAndEndPoints(j, m_startPoint2, m_endPoint2);
            double distance = computeDistanceBetweenTwoLineSegments(startPoint, endPoint, m_startPoint2, m_endPoint2);
            //  if the distance is below the threshold, this line segment belongs to the eps-neighborhood
            if (distance <= eps) result.add(j);
        }
    }

	private void computeEPSNeighborhoodByRtree(CMDPoint startPoint, CMDPoint endPoint, double eps, Set<Integer> result) {
		result.clear();
		Geometry searchRegion = get_bounding_box_of_line_segment(startPoint, endPoint, eps); //直接最大外接矩形
        List lst = m_tree.query(searchRegion.getEnvelopeInternal());

        for (int i = 0; i < lst.size(); i++) {
            int segmentID = (int)lst.get(i);
//            System.out.println(segmentID);
            extractStartAndEndPoints(segmentID, m_startPoint2, m_endPoint2);
            double distance = computeDistanceBetweenTwoLineSegments(startPoint, endPoint, m_startPoint2, m_endPoint2);
//            double distance = measureDistanceFromPointToLineSegment(m_startPoint2, m_endPoint2, startPoint);

            if (distance <= eps)
                result.add(segmentID);
//            result.add(segmentID);
        }
	}

    private void computeEPSNeighborhoodByRtree(CMDPoint startPoint, CMDPoint endPoint, double eps, Deque<Integer> result) {
        result.clear();
        Geometry searchRegion = get_bounding_box_of_line_segment(startPoint, endPoint, eps); //直接最大外接矩形
        List lst = m_tree.query(searchRegion.getEnvelopeInternal());

        GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();

        for (int i = 0; i < lst.size(); i++) {
            int segmentID = (int)lst.get(i);
//            if(segmentID == 0){
//                System.out.println("0");
//            }
            extractStartAndEndPoints(segmentID, m_startPoint2, m_endPoint2);
            double distance = computeDistanceBetweenTwoLineSegments(startPoint, endPoint, m_startPoint2, m_endPoint2);

            if (distance <= eps)
                result.addFirst(segmentID);

//            Coordinate[] coords  =
//                    new Coordinate[] {new Coordinate(m_startPoint2.getM_coordinate(0),m_startPoint2.getM_coordinate(1)),
//                            new Coordinate(m_endPoint2.getM_coordinate(0),m_endPoint2.getM_coordinate(1))};
//
//            LineString line = geometryFactory.createLineString(coords);
//            if(line.intersects(searchRegion)) {
//                double distance = computeDistanceBetweenTwoLineSegments(startPoint, endPoint, m_startPoint2, m_endPoint2);
//
//                if (distance <= eps)
//                    result.addFirst(segmentID);
//            }
        }
//        logger.debug(result.size());
    }

	private Geometry get_bounding_box_of_line_segment(CMDPoint startPoint, CMDPoint endPoint, double radius) {
		double btm = Math.min(startPoint.getM_coordinate(1), endPoint.getM_coordinate(1)) - radius;
		double top = Math.max(startPoint.getM_coordinate(1), endPoint.getM_coordinate(1)) + radius;
		double left = Math.min(startPoint.getM_coordinate(0), endPoint.getM_coordinate(0)) - radius;
		double right = Math.max(startPoint.getM_coordinate(0), endPoint.getM_coordinate(0)) + radius;

		Coordinate[] coords = new Coordinate[5];
		coords[0] = new Coordinate(left, top);
		coords[1] = new Coordinate(right, top);
		coords[2] = new Coordinate(right, btm);
		coords[3] = new Coordinate(left, btm);
		coords[4] = new Coordinate(left, top);
		GeometryFactory factory = new GeometryFactory();
//        LinearRing ring= new LinearRing(new CoordinateArraySequence(coors),new GeometryFactory());
        Polygon polygon = factory.createPolygon(coords);
		return polygon;
	}

    /**
     * 建立R树索引
     *
     * @return
     */
    public STRtree insertRtree(CMDPoint lineSegmentPoint, int nTotalLineSegment) {
//        for (int index = 0; index < m_lineSegmentPointArray.size(); index++) {
        GeometryFactory factory = new GeometryFactory();
//        CMDPoint lineSegmentPoint = m_lineSegmentPointArray.get(index);

        Coordinate[] coords =
                new Coordinate[]{
                        new Coordinate(lineSegmentPoint.getM_coordinate(0), lineSegmentPoint.getM_coordinate(1)),
                        new Coordinate(lineSegmentPoint.getM_coordinate(2), lineSegmentPoint.getM_coordinate(3))};
        LineString lineStr = factory.createLineString(coords);
        Envelope env = lineStr.getEnvelopeInternal();

        m_tree.insert(env, nTotalLineSegment);
//        }
//		m_tree.build();
//		System.out.println("R树索引创建成功！");
//		for (Line line : lines) {
//			GeometryFactory factory = new GeometryFactory();
//			Coordinate[] coords =
//					new Coordinate[]{new Coordinate(line.getS().x, line.getS().y), new Coordinate(line.getE().x, line.getE().y)};
//			LineString lineStr = factory.createLineString(coords);
//			Envelope env = lineStr.getEnvelopeInternal();
//			m_tree.insert(env, lineStr);
//		}
//
//		m_tree.build();
//		System.out.println("R树索引创建成功！");
        return m_tree;
    }

    private double computeDistanceBetweenTwoLineSegments(CMDPoint startPoint1,
                                                         CMDPoint endPoint1, CMDPoint startPoint2, CMDPoint endPoint2) {
        double perpendicularDistance = 0;
        double parallelDistance = 0;
        double angleDistance = 0;

        return subComputeDistanceBetweenTwoLineSegments(startPoint1, endPoint1, startPoint2, endPoint2, perpendicularDistance, parallelDistance, angleDistance);
    }


    private boolean storeLineSegmentCluster() {

        int currClusterId = 0;

        for (int i = 0; i < m_currComponentId; i++) {

            if (!m_lineSegmentClusters[i].enabled) {
                continue;
            }

            //  store the clusters finally identified
            //  START ...
            Cluster pClusterItem = new Cluster(currClusterId, m_document.m_nDimensions);

            for (int j = 0; j < m_lineSegmentClusters[i].nClusterPoints; j++) {
                pClusterItem.addPointToArray(m_lineSegmentClusters[i].clusterPointArray.get(j));
            }

            pClusterItem.setDensity(m_lineSegmentClusters[i].nTrajectories);

            m_document.m_clusterList.add(pClusterItem);

            currClusterId++;    //  increase the number of final clusters
            //  ... END
        }

        m_document.m_nClusters = currClusterId;
        return true;
    }


    private double subComputeDistanceBetweenTwoLineSegments(CMDPoint startPoint1,
                                                            CMDPoint endPoint1, CMDPoint startPoint2, CMDPoint endPoint2,
                                                            double perpendicularDistance, double parallelDistance,
                                                            double angleDistance) {

        double perDistance1, perDistance2;
        double parDistance1, parDistance2;
        double length1, length2;

        //  the length of the first line segment
        length1 = measureDistanceFromPointToPoint(startPoint1, endPoint1);
        //  the length of the second line segment
        length2 = measureDistanceFromPointToPoint(startPoint2, endPoint2);

        //  compute the perpendicular distance and the parallel distance
        //  START ...
        if (length1 > length2) {
            perDistance1 = measureDistanceFromPointToLineSegment(startPoint1, endPoint1, startPoint2);
            if (m_coefficient < 0.5) parDistance1 = measureDistanceFromPointToPoint(startPoint1, m_projectionPoint);
            else parDistance1 = measureDistanceFromPointToPoint(endPoint1, m_projectionPoint);

            perDistance2 = measureDistanceFromPointToLineSegment(startPoint1, endPoint1, endPoint2);
            if (m_coefficient < 0.5) parDistance2 = measureDistanceFromPointToPoint(startPoint1, m_projectionPoint);
            else parDistance2 = measureDistanceFromPointToPoint(endPoint1, m_projectionPoint);
        } else {
            perDistance1 = measureDistanceFromPointToLineSegment(startPoint2, endPoint2, startPoint1);
            if (m_coefficient < 0.5) parDistance1 = measureDistanceFromPointToPoint(startPoint2, m_projectionPoint);
            else parDistance1 = measureDistanceFromPointToPoint(endPoint2, m_projectionPoint);

            perDistance2 = measureDistanceFromPointToLineSegment(startPoint2, endPoint2, endPoint1);
            if (m_coefficient < 0.5) parDistance2 = measureDistanceFromPointToPoint(startPoint2, m_projectionPoint);
            else parDistance2 = measureDistanceFromPointToPoint(endPoint2, m_projectionPoint);

        }

        //  compute the perpendicular distance; take (d1^2 + d2^2) / (d1 + d2)
        if (!(perDistance1 == 0.0 && perDistance2 == 0.0)) {
            perpendicularDistance = ((Math.pow(perDistance1, 2) + Math.pow(perDistance2, 2)) / (perDistance1 + perDistance2));
        } else {
            perpendicularDistance = 0.0;
        }
        //  compute the parallel distance; take the minimum
        parallelDistance = (parDistance1 < parDistance2) ? parDistance1 : parDistance2;
        //  ... END

        //  compute the angle distance
        //  START ...
        //  MeasureAngleDisntance() assumes that the first line segment is longer than the second one
        if (length1 > length2) {
            angleDistance = measureAngleDisntance(startPoint1, endPoint1, startPoint2, endPoint2);
        } else {
            angleDistance = measureAngleDisntance(startPoint2, endPoint2, startPoint1, endPoint1);
        }
        //  ... END

        return (perpendicularDistance + parallelDistance + angleDistance);


    }

    private void extractStartAndEndPoints(int index, CMDPoint startPoint, CMDPoint endPoint) {//  for speedup
        //  compose the start and end points of the line segment
        for (int i = 0; i < m_document.m_nDimensions; i++) {
            startPoint.setM_coordinate(i, m_lineSegmentPointArray.get(index).getM_coordinate(i));
            endPoint.setM_coordinate(i, m_lineSegmentPointArray.get(index).getM_coordinate(m_document.m_nDimensions + i));
            ;
        }
    }

    public boolean estimateParameterValue(Parameter p) {

        double entropy, minEntropy = (double) INT_MAX;
        double eps, minEps = (double) INT_MAX;
        int totalSize, minTotalSize = INT_MAX;
        Set<Integer> seeds = new HashSet<Integer>();

        int[] EpsNeighborhoodSize = new int[m_nTotalLineSegments];

        for (eps = 50; eps <= 1000; eps += 1.0) {
            entropy = 0.0;
            totalSize = 0;
            seeds.clear();
            for (int i = 0; i < m_nTotalLineSegments; i++) {
                extractStartAndEndPoints(i, m_startPoint1, m_endPoint1);
                computeEPSNeighborhoodByRtree(m_startPoint1, m_endPoint1, eps, seeds);
                EpsNeighborhoodSize[i] = (int) seeds.size();
                totalSize += (int) seeds.size();
                seeds.clear();
            }
            for (int i = 0; i < m_nTotalLineSegments; i++) {
                entropy += ((double) EpsNeighborhoodSize[i] / (double) totalSize) * LOG2((double) EpsNeighborhoodSize[i] / (double) totalSize);
            }
            entropy = -entropy;

            if (entropy < minEntropy) {
                minEntropy = entropy;
                minTotalSize = totalSize;
                minEps = eps;
            }
        }
        // setup output arguments
        p.epsParam = minEps;
        p.minLnsParam = (int) Math.ceil((double) minTotalSize / (double) m_nTotalLineSegments);
        return true;
    }

}
