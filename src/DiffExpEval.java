import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;

public class DiffExpEval {

	public static void main(String[] args) {
		
		String countsfilesPath ="";
		String labelsPath ="";
		String outputPath ="";
		String configPath = "";
		for(int i =0; i < args.length-1; i++){
			if(args[i].equals("-countfiles")){
				countsfilesPath = args[i+1];
				i++;
			}
			else if(args[i].equals("-labels")){
				labelsPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-outdir")){
				outputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-config")){
				configPath = args[i+1];
				i++; 
			}
		}
		
		if(countsfilesPath.equals("") || labelsPath.equals("") || outputPath.equals("") || configPath.equals("")){
			System.out.println("Usage Info:\n-countfiles <list of countfiles>\n-labels <label file>\n-outdir <output directory>\n-config <path to config file>");
		}
		else{
	
			Path countsFilePath = Paths.get(countsfilesPath);
			Path labelsFilePath = Paths.get(labelsPath);
			Path configFilePath = Paths.get(configPath);
			
			
			
			File countsFile = countsFilePath.toFile();
			File labelFile = labelsFilePath.toFile();
			File configFile = configFilePath.toFile();

			try {
				DiffExpEval divexpeval = new DiffExpEval(configFile);
				
				divexpeval.getEnrichmentBrowser(countsFile,labelFile,outputPath);
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private void getEnrichmentBrowser(File countFiles, File labelFile, String outputPath) throws IOException {
		
		String fPath = outputPath.concat("/f_data.txt");
		String pPath = outputPath.concat("/p_data.txt");
		String exPath = outputPath.concat("/exprs.txt");
		
		Path fFilePath = Paths.get(fPath);
		Path pFilePath = Paths.get(pPath);
		Path exFilePath = Paths.get(exPath);
		
		String dsPath = outputPath.concat("/DEseq.out");
		String erPath = outputPath.concat("/edgeR.out");
		String liPath = outputPath.concat("/limma.out");
		
		ArrayList<String> geneList = new ArrayList<String>();
		
		HashMap<String,ArrayList<Integer>> counts = new HashMap<String,ArrayList<Integer>>();

		BufferedReader lbr = new BufferedReader (new FileReader(labelFile));
        String lline;

        while ((lline = lbr.readLine()) != null){
        	String[] lineSplit = lline.split("\t");
        	if(lineSplit.length > 1){
        		geneList.add(lineSplit[0]);
        		counts.put(lineSplit[0], new ArrayList<Integer>());
        	}
        }
		
        lbr.close();
        
        HashSet<String> genes = new HashSet<String>(geneList);
        
        StringBuilder fBuilder = new StringBuilder();
        StringBuilder pBuilder = new StringBuilder();
        StringBuilder eBuilder = new StringBuilder();
        
        String tab = "\t";
        String brk = "\n";
        char dot = '.';
        char sp = ' ';
        char com = ',';
        
        BufferedReader clbr = new BufferedReader (new FileReader(countFiles));
        String clline;
       
        
        String curCondition = "";
        int i = -1;
        while ((clline = clbr.readLine()) != null){
        	String[] lineSplit = clline.split("\t");
        	if(lineSplit.length > 1){
        		String condition = lineSplit[0];
        		if(!condition.equals(curCondition)){
        			i ++;
        			curCondition = condition;
        		}
        		String replicate = lineSplit[2];
      
        		
        		pBuilder.append(condition);
        		pBuilder.append(dot);
        		pBuilder.append(replicate);
        		pBuilder.append(tab);
        		pBuilder.append(i);
        		pBuilder.append(brk);
        		
        		String countPath = lineSplit[1];
        		Path countFilePath = Paths.get(countPath);
        		File countFile = countFilePath.toFile();
        		
                BufferedReader cfbr = new BufferedReader (new FileReader(countFile));
                
                HashSet<String> foundGenes = new HashSet<String>();
                String cfline;
                cfbr.readLine();
                while ((cfline = cfbr.readLine()) != null){
                	String[] clineSplit = cfline.split("\t");
                	if(clineSplit.length > 1){
                		String gene = clineSplit[0];
                		if(counts.containsKey(gene)){
                			counts.get(gene).add(Integer.parseInt(clineSplit[8]));
                			foundGenes.add(gene);
                		}
                	}
                }
                cfbr.close();
                HashSet<String> cgenes = new HashSet<String>(genes);
        		foundGenes.retainAll(cgenes);
        		cgenes.removeAll(foundGenes);
        		for(String gene: cgenes){
        			counts.get(gene).add(0);
        		}
        	}
        }
		
        clbr.close();
        
        for(String gene: geneList){
        	fBuilder.append(gene);
        	for(int j = 1; j <= i; j++){
        		fBuilder.append(tab);
        		fBuilder.append(gene);
        	}
        	fBuilder.append(brk);
        	
        	ArrayList<Integer> geneCount = counts.get(gene);
        	eBuilder.append(geneCount.get(0));
        	for(int k = 1; k< geneCount.size(); k++){
        		eBuilder.append(tab);
        		eBuilder.append(geneCount.get(k));
        	}
        	eBuilder.append(brk);
        }
        
        ArrayList<String> fList = new ArrayList<String>();
        fList.add(fBuilder.toString());
        Files.write(fFilePath, fList, Charset.forName("UTF-8"));
        
        ArrayList<String> pList = new ArrayList<String>();
        pList.add(pBuilder.toString());
        Files.write(pFilePath, pList, Charset.forName("UTF-8"));

        ArrayList<String> exList = new ArrayList<String>();
        exList.add(eBuilder.toString());
        Files.write(exFilePath, exList, Charset.forName("UTF-8"));
        
        String rInput = rPath.concat(" ").concat(diffScriptPath);
        String enrichInput = exPath.concat(" ").concat(pPath).concat(" ").concat(fPath);
        
        String dsOutput = "DEseq ".concat(dsPath);
        String erOutput = "edgeR ".concat(erPath);
        String liOutput = "limma ".concat(liPath);
        
        Runtime.getRuntime().exec(rInput.concat(" ").concat(enrichInput).concat(" ").concat(dsOutput)); 
        Runtime.getRuntime().exec(rInput.concat(" ").concat(enrichInput).concat(" ").concat(erOutput)); 
        Runtime.getRuntime().exec(rInput.concat(" ").concat(enrichInput).concat(" ").concat(liOutput)); 
        
        //Funktioniert nicht mangels R Version
        //TODO Script Outputs einlesen, mit DoubleComparator sortieren, in BenjaminiHochberg als Parameter und wieder ausgeben.
	}

	private String rPath;
	private String diffScriptPath;
	
	public DiffExpEval(File configFile) throws IOException{
		BufferedReader br = new BufferedReader (new FileReader(configFile));
        String line;
        
        while ((line = br.readLine()) != null){
        	String[] lineSplit = line.split("\t");
        	if(lineSplit.length > 1){
        		if(lineSplit[0].equals("R")){
        			this.rPath = lineSplit[1];
        		}
        		else if(lineSplit[0].equals("diffscript")){
        			this.diffScriptPath = lineSplit[1];
        		}
        	}
        }
        br.close();
        
//		System.out.println(rPath);
//		System.out.println(diffScriptPath);
	}
	
	public double[] BenjaminiHochberg(ArrayList<Double> pValues) {
		
		int length = pValues.size();
		
        double[] apValues = new double[length];

        for (int i = length - 1; i >= 0; i--) {
            if (i == length - 1) {
                apValues[i] = pValues.get(i);
            } else {
                double l = apValues[i + 1];
                double r = (length / (double) (i+1)) * pValues.get(i);
                apValues[i] = Math.min(l, r);
            }
        }
        
        return apValues;
	}
	
	 class DoubleComparator implements Comparator<Double> {
		  @Override
		  public int compare(Double x1, Double x2) {
		    return x1.compareTo(x2);
		  }
	 }

}
