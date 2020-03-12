package mono;

import com.sun.media.jai.mlib.MlibDCTRIF;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Locale;
import mono.ClusterGen.LineSegmentId;

public class TraClusterDoc {
	private static final Logger logger = LogManager.getLogger(TraClusterDoc.class.getName());

	public int m_nDimensions;
	public int m_nTrajectories;
	public int m_nClusters;
	public double m_clusterRatio;
	public int m_maxNPoints;
	public ArrayList<Trajectory> m_trajectoryList;
	public ArrayList<CMDPoint> m_lineSegmentPointArray;
	public ArrayList<Integer> m_componentIdArray;
	public ArrayList<Cluster> m_clusterList;
	public ArrayList<LineSegmentId> m_idArray;
	public double epsParam;
	public int minLnsParam;

	public TraClusterDoc() {
			
		m_nTrajectories = 0;
		m_nClusters = 0;
		m_clusterRatio = 0.0;	
		m_trajectoryList = new ArrayList<Trajectory>();
		m_clusterList = new ArrayList<Cluster>();
	}
	
	public class Parameter {
		double epsParam;
		int minLnsParam;
	}
	
	boolean onOpenDocument(String inputFileName) {
		
		int nDimensions = 2;		// default dimension = 2
		int nTrajectories = 0;
		int nTotalPoints = 0;		//no use
		int trajectoryId;
		int nPoints;
		double value;

		DataInputStream in;
		BufferedReader inBuffer = null;
		try {
			in = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFileName)));
			
			inBuffer = new BufferedReader(new InputStreamReader(in));			
			nDimensions = Integer.parseInt(inBuffer.readLine()); // the number of dimensions
			m_nDimensions = nDimensions;
			nTrajectories = Integer.parseInt(inBuffer.readLine()); // the number of trajectories
			m_nTrajectories = nTrajectories;
			
			m_maxNPoints = -1; // initialize for comparison
			
			// the trajectory Id, the number of points, the coordinate of a point ...
			for (int i = 0; i < nTrajectories; i++) {				
	
				String str = inBuffer.readLine();
				
				Scanner sc = new Scanner(str); 
				sc.useLocale(Locale.US);
				
				trajectoryId = sc.nextInt(); //trajectoryID
				nPoints = sc.nextInt(); // number of points in the trajectory
				
				if (nPoints > m_maxNPoints) {
					m_maxNPoints = nPoints;
				}
				nTotalPoints += nPoints;
				Trajectory pTrajectoryItem = new Trajectory(trajectoryId, nDimensions);		
				for (int j = 0; j < nPoints; j++) {
					CMDPoint point = new CMDPoint(nDimensions);   // initialize the CMDPoint class for each point
					
					for (int k = 0; k < nDimensions; k++) {						
						value = sc.nextDouble();
						point.setM_coordinate(k, value);						
					}
					pTrajectoryItem.addPointToArray(point);				
				}
				
				m_trajectoryList.add(pTrajectoryItem);
			}					
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			logger.error("Unable to open input file");
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}  finally {
			try {
				inBuffer.close();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
        		
		return true;
	}

	void loadTrajectoryFromShp(String srcfilepath) {
		int nDimensions = 2;		// default dimension = 2
		int nTrajectories = 0;
		int trajectoryId; //轨迹序号
		int order = 0; //所有点的序号
		int nPoints; //每条轨迹点的数量
		double value;
		SimpleFeatureIterator iterator = null;

		m_nDimensions = nDimensions;

		File srcfile = new File(srcfilepath);
		FileDataStore store = null;
		//创建目标shape文件对象
		try {
			store= FileDataStoreFinder.getDataStore(srcfile);
			SimpleFeatureSource featureSource = store.getFeatureSource();
			iterator = featureSource.getFeatures().features();

			trajectoryId = 0;
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
//                System.out.println(feature.getAttribute(4).toString());
				Geometry geom = (Geometry) feature.getDefaultGeometry();
				Coordinate[] coords = geom.getCoordinates();

				nPoints = coords.length;
				if (nPoints > m_maxNPoints) {
					m_maxNPoints = nPoints;
				}

				Trajectory pTrajectoryItem = new Trajectory(trajectoryId, nDimensions);
				for (int i = 0; i < coords.length; i++) {
					double x = coords[i].getX();
					double y = coords[i].getY();

					CMDPoint point = new CMDPoint(nDimensions);   // initialize the CMDPoint class for each point

					point.setM_coordinate(0, x);
					point.setM_coordinate(1, y);

					pTrajectoryItem.addPointToArray(point);
				}

				m_trajectoryList.add(pTrajectoryItem);
				trajectoryId++;
			}
			m_nTrajectories = trajectoryId;
//			iterator.close();

		} catch (IOException e) {
			logger.error("打开文件错误.");
			e.printStackTrace();
		} catch(Exception e){
			e.printStackTrace();
		} finally {
			try{
				iterator.close();
				store.dispose();
			} catch (Exception e){
				e.printStackTrace();
			}
		}
	}

	boolean onClusterGenerate(String clusterFileName, double epsParam, int minLnsParam, double minLnLen, int mdlcost) {
//////////////////////////////////////////////////still to be written
		this.epsParam = epsParam;
		this.minLnsParam = minLnsParam;
		ClusterGen generator = new ClusterGen(this);
		generator.setMinLinesegmentLength(minLnLen);
		generator.setMdlCostAdwantage(mdlcost);

		if(m_nTrajectories == 0) {
			logger.error("Load a trajectory data set first");
		}
		
		// FIRST STEP: Trajectory Partitioning
		logger.info("开始分段...");
		if (!generator.partitionTrajectory()) {
			logger.error("Unable to partition a trajectory\n");
			return false;
		}
		logger.info("分段完成.");

		// SECOND STEP: Density-based Clustering
		logger.info("开始聚类...");
		if (!generator.performDBSCAN())
		{
			logger.error("Unable to perform the DBSCAN algorithm\n");
			return false;
		}
		logger.info("聚类完成.");

		logger.info("开始生成轨迹结果...");
		// THIRD STEP: Cluster Construction
		if (!generator.constructCluster())
		{
			logger.error( "Unable to construct a cluster\n");
			return false;
		}
		logger.info("轨迹结果生成完成.");

		m_lineSegmentPointArray = generator.get_lineSegmentPointArray();
		m_componentIdArray = generator.getM_componentIdArray();
		m_idArray = generator.getM_idArray();

//		for (int i = 0; i <m_clusterList.size(); i++) {
//			//m_clusterList.
//			System.out.println(m_clusterList.get(i).getM_clusterId());
//			for (int j = 0; j<m_clusterList.get(i).getM_PointArray().size(); j++) {
//
//				double x = m_clusterList.get(i).getM_PointArray().get(j).getM_coordinate(0);
//				double y = m_clusterList.get(i).getM_PointArray().get(j).getM_coordinate(1);
//			System.out.print("   "+ x +" "+ y +"   ");
//			}
//			System.out.println();
//		}
		FileOutputStream fos = null;
		BufferedWriter bw = null;
		OutputStreamWriter osw = null;
		try {
			fos = new FileOutputStream(clusterFileName);
			osw = new OutputStreamWriter(fos);
			bw = new BufferedWriter(osw);
			
			bw.write("epsParam:"+epsParam +"   minLnsParam:"+minLnsParam);
			
			for (int i = 0; i < m_clusterList.size(); i++) {
				// m_clusterList.
				bw.write("\nclusterID: "+ m_clusterList.get(i).getM_clusterId() + "  Points Number:  " + m_clusterList.get(i).getM_PointArray().size() + "\n");
				for (int j = 0; j < m_clusterList.get(i).getM_PointArray().size(); j++) {
					
					double x = m_clusterList.get(i).getM_PointArray().get(j).getM_coordinate(0);
					double y = m_clusterList.get(i).getM_PointArray().get(j).getM_coordinate(1);
					bw.write(x+" "+y+"   ");
				}
			}						
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return true;		
	}
	
	Parameter onEstimateParameter() {
		Parameter p = new Parameter();
		ClusterGen generator = new ClusterGen(this);
		if (!generator.partitionTrajectory()) {
			logger.error("Unable to partition a trajectory\n");
			return null;
		}
		if (!generator.estimateParameterValue(p)) {
			logger.error("Unable to calculate the entropy\n");
			return null;
		}
		return p;
	}

}
