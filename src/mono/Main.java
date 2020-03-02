package mono;

import org.apache.commons.cli.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;

public class Main {
	private static CommandLine cl;
	private static String HELP_STRING = null;
	private static Options OPTIONS = new Options();

	public static void main(String[] args) {
		String input_file="";

//		if (args.length == 4) {
//			TraClusterDoc tcd = new TraClusterDoc();
//			tcd.onOpenDocument(args[0]);
//			tcd.onClusterGenerate(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3])); // 25, 5~7
//		} else if (args.length == 2) {
//			TraClusterDoc tcd = new TraClusterDoc();
//			tcd.onOpenDocument(args[0]);
//
//			Parameter p = tcd.onEstimateParameter();
//			if (p != null) {
//				System.out.println("Based on the algorithm, the suggested parameters are:\n"+"eps:" + p.epsParam + "  minLns:" + p.minLnsParam);
//			}
//			tcd.onClusterGenerate(args[1], p.epsParam, p.minLnsParam);
//		} else {
//			System.out.println("Please give me 2 or 4 input parameters! \n "
//					+ "If you have no idea how to decide eps and minLns, just feed in 2 parameters (inputFilePath, outputFilePath):\n"
//					+ "--e.g. java mono.Main deer_1995.tra testOut.txt \n"
//					+ "If you know the two parameters, just feed in all the 4 parameters (inputFilePath, outputFilePath, eps, minLns)"
//					+ "--e.g. java mono.Main deer_1995.tra testOut.txt 29 8 \n");
//		}

/**
 * To use the GUI, Remove the below comment and comment out the above section of code
 * An adjustable GUI is to be added.
 */
		CommandLineParser commandLineParser = new DefaultParser();
//		Options OPTIONS = new Options();
		OPTIONS.addOption(Option.builder("i").
				argName("input file").
				longOpt("input").required().
				hasArg(true).type(String.class).build());

		try {
			cl = commandLineParser.parse(OPTIONS, args);

			if(cl.hasOption("i")){
				input_file = cl.getOptionValue("i");
			}
		} catch (ParseException e) {
			System.out.println(e.getMessage() + "\n" + getHelpString());
			System.exit(0);
		}

		long startTime=System.currentTimeMillis();

		TraClusterDoc tcd = new TraClusterDoc();

		//tcd.onOpenDocument("src/deer_1995.tra");
		//tcd.onClusterGenerate("testDeerResult.txt", 29, 8);

		//tcd.onOpenDocument("src/hurricane1950_2006.tra");
		//tcd.onClusterGenerate("testHurricaneResult.txt", 32, 6);

//		tcd.onOpenDocument("data/deer_1995.tra");
//		tcd.onClusterGenerate("testDeerResult.txt", 29, 8);// 25, 5~7

		tcd.loadTrajectoryFromShp(input_file);   //经过22和29节点的线.shp res_sp_s.shp
		tcd.onClusterGenerate("res/testDeerResult.txt", 200, 5);// 25, 5~7

		MainFrame2 mf = new MainFrame2(tcd.m_trajectoryList, tcd.m_lineSegmentPointArray, tcd.m_clusterList);


//		Parameter p = tcd.onEstimateParameter();
//		if (p != null) {
//			System.out.println("Based on the algorithm, the suggested parameters are:\n" + "eps:" + p.epsParam + "  minLns:" + p.minLnsParam);
//		}
		long endTime=System.currentTimeMillis();
		System.out.println("over!" + (endTime-startTime) / 1000);
	}

	private static String getHelpString() {
		if (HELP_STRING == null) {
			HelpFormatter helpFormatter = new HelpFormatter();

			ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
			PrintWriter printWriter = new PrintWriter(byteArrayOutputStream);
			helpFormatter.printHelp(printWriter, HelpFormatter.DEFAULT_WIDTH, "traclus -help", null,
					OPTIONS, HelpFormatter.DEFAULT_LEFT_PAD, HelpFormatter.DEFAULT_DESC_PAD, null);
			printWriter.flush();
			HELP_STRING = new String(byteArrayOutputStream.toByteArray());
			printWriter.close();
		}
		return HELP_STRING;
	}
	
}
