package phylogenetic_tree;

import java.io.*;
import java.util.*;

public class Main {
    public static void main(String[] args) {


        String clustalOutput = "clustalOutput.txt"; 
        String matrixFile = "matrix.txt"; 
        processDists(clustalOutput, matrixFile); 
        System.exit(1); 

    }
    

    static void credit() {
        System.out.println("Clustal Matrix Generation script, Brian Chen, Lehigh University, 2011. Translated to java by Luke Coolbaugh");
    }
    // ... The remainder of the functions would follow the same pattern ...

    // I will continue to translate the remaining functions in the next message.

    static double getAlnScore(String line) {
        String[] lineList = line.trim().split("\\s+");
        if (lineList.length != 5 || !lineList[0].equals("Sequences")) {
            System.out.println("ERROR: this line [" + line + "] is not a Sequence Pair Line");
            System.out.println("Which should look like this:");
            System.out.println("Sequences (1:2) Aligned. Score:  100");
            System.exit(1);
        }
        return 1.0 - Double.parseDouble(lineList[lineList.length - 1]) / 100.0;
    }
    
    static int getSecondSeqId(String line) {
        String[] lineList = line.trim().split("\\s+");
        if (lineList.length != 5 || !lineList[0].equals("Sequences")) {
            System.out.println("ERROR: this line [" + line + "] is not a Sequence Pair Line");
            System.out.println("Which should look like this:");
            System.out.println("Sequences (1:2) Aligned. Score:  100");
            System.exit(1);
        }
        return Integer.parseInt(lineList[1].split(":")[1].substring(0, lineList[1].split(":")[1].length() - 1));
    }
    
    static int getFirstSeqId(String line) {
        String[] lineList = line.trim().split("\\s+");
        if (lineList.length != 5 || !lineList[0].equals("Sequences")) {
            System.out.println("ERROR: this line [" + line + "] is not a Sequence Pair Line");
            System.out.println("Which should look like this:");
            System.out.println("Sequences (1:2) Aligned. Score:  100");
            System.exit(1);
        }
        return Integer.parseInt(lineList[1].split(":")[0].substring(1));
    }
    
    
    static String getSequenceName(String line) {
        String[] lineList = line.split(" ");
        if (lineList.length < 5 || !lineList[lineList.length - 1].equals("aa")) {
            System.out.println("ERROR: this line [" + line + "] is not a sequence Number Line");
            System.exit(1);
        }
        // Assuming the sequence name is the 0th element in the split line.
        return lineList[2];
    }
       
    
    int getSequenceNumber(String line) {
        String[] lineList = line.split(" ");
        if (lineList.length < 5 || !lineList[lineList.length - 1].equals("aa")) {
            System.out.println("ERROR: this line [" + line + "] is not a sequence Number Line");
            System.exit(1);
        }
        return Integer.parseInt(lineList[1].substring(0, lineList[1].length() - 1));
    }    
    
    static String taxaNameProcessing(String name) {
        if (name.length() == 0) {
            System.out.println("ERR: no name found! exiting!");
            System.exit(1);
        }
        if (name.length() > 10) {
            return name.substring(0, 10) + " ";
        } else {
            return String.format("%-10s", name) + " ";
        }
    }
    
    // ... The processDists function is too long to fit into this message, so I will continue to translate it in the next message ...
    

    static void processDists(String clustalOutput, String matrixFile) {
        try {
            BufferedReader scoresFile = new BufferedReader(new FileReader(clustalOutput));
            String line = scoresFile.readLine();
    
            List<String> nameList = new ArrayList<>();
    
            System.out.println("Parsing Data File");
    
            // Skip to the sequence names and indices
            while (line != null) {
                if (line.startsWith("Sequence format is Pearson")) {
                    line = scoresFile.readLine();
                    break;
                }
                line = scoresFile.readLine();
            }
    
// Read and store the sequence names and indices
while (line != null) {
    if (line.startsWith("Start of Pairwise alignments")) {
        break;
    }

    if (line.startsWith("Sequence")) {
        String seqName = getSequenceName(line);
        nameList.add(seqName);
    }

    line = scoresFile.readLine();
}

    
            // Set up the storage datastructure
            Double[][] matrix = new Double[nameList.size()][nameList.size()];
            for (int i = 0; i < nameList.size(); i++) {
                Arrays.fill(matrix[i], 0.0);  // fill each row with 0.0 by default
            }
            // Now parse out the scores and put them in the matrix
            while (line != null) {
                String[] lineList = line.split(" ");
                if (lineList.length == 0) {
                    line = scoresFile.readLine();
                    continue;
                }
            
                if (!line.startsWith("Sequences")) {
                    line = scoresFile.readLine();
                    continue;
                }
            
                int seqNum1 = getFirstSeqId(line);
                int seqNum2 = getSecondSeqId(line);
                double score = getAlnScore(line);
            
                System.out.println("Got Seq1: [" + seqNum1 + "] Got Seq2: [" + seqNum2 + "] got Score: [" + score + "]");
            
                // Symmetry
                matrix[seqNum1 - 1][seqNum2 - 1] = score;
                matrix[seqNum2 - 1][seqNum1 - 1] = score;
                line = scoresFile.readLine();
            }
            
    
            // Finally output the score matrix.
            PrintWriter outputFile = new PrintWriter(new FileWriter(matrixFile));
            outputFile.println("   " + nameList.size());
    
// Now the output values
for (int i = 0; i < nameList.size(); i++) {
    StringBuilder lineOutput = new StringBuilder(taxaNameProcessing(nameList.get(i)));
    for (int j = 0; j < nameList.size(); j++) {
        if (i == j) {
            lineOutput.append("1.0  ");  // Assuming sequence similarity with itself is 1.0
        } else if (matrix[i][j] == null) {
            lineOutput.append("0.0  ");  // Assuming sequence similarity is 0.0 if not defined
        } else {
            lineOutput.append(matrix[i][j]).append("  ");
        }
    }

    lineOutput.append("\n");
    outputFile.write(lineOutput.toString());
}

    
            // Done
            outputFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
