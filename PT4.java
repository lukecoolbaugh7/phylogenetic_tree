package phylogenetic_tree;

import java.util.*;
import java.io.*;

public class PT4 {
    public static void main(String[] args) throws FileNotFoundException {
        Map.Entry<List<String>, double[][]> result = readMatrixFromFile("matrix.txt");
        List<String> sequences = result.getKey();
        System.out.println("");
        double[][] matrix = result.getValue();
        System.out.println("");
        System.out.println("UPGMA Newick output: ");

        createUPGMATree(sequences, matrix);
        System.out.println("\n");

        double[][] matrixNJ = result.getValue();
        List<String> sequencesNJ = result.getKey();
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < sequencesNJ.size(); i++) {
            clusters.add(new Cluster(sequencesNJ.get(i), i));
        }

        System.out.println("NJ Newick output: ");

        createNJTree(clusters, matrixNJ);

    }

    public static void printArr(double[][] array) {
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[i].length; j++) {
                System.out.print(array[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static Map.Entry<List<String>, double[][]> readMatrixFromFile(String filename) throws FileNotFoundException {
        Scanner scanner = new Scanner(new File(filename));
        int n = scanner.nextInt();
        scanner.nextLine(); // skip to next line
        double[][] input = new double[n][n];
        List<String> sequences = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            String[] line = scanner.nextLine().split("\\s+");
            sequences.add(line[0]); // store sequence name
            for (int j = 1; j <= n; j++) { // start from 1 to skip sequence name
                input[i][j - 1] = Double.parseDouble(line[j]); // subtract 1 to align with matrix array
            }
        }
        for (int i = 0; i < n; i++) {
            input[i][i] = 0.0;
        }
        scanner.close();
        return new AbstractMap.SimpleEntry<>(sequences, input);
    }

    public static void createUPGMATree(List<String> sequences, double[][] matrix) {
        List<Cluster> clusters = new ArrayList<>();
        for (int i = 0; i < sequences.size(); i++) {
            clusters.add(new Cluster(sequences.get(i)));
        }

        while (clusters.size() > 1) {
            // Find the pair of clusters with the smallest distance
            int[] minIndex = findMinDistance(clusters, matrix);
            Cluster c1 = clusters.get(minIndex[0]);
            Cluster c2 = clusters.get(minIndex[1]);

            // Merge the two clusters and update the distance matrix
            Cluster newCluster = new Cluster(c1, c2, matrix[minIndex[0]][minIndex[1]] / 2);
            clusters.remove(c1);
            clusters.remove(c2);
            clusters.add(newCluster);
            matrix = updateMatrix(matrix, clusters, newCluster, minIndex);
        }

        // Print the UPGMA tree in Newick format
        String root = clusters.get(0).toString() + ";";
        System.out.println(clusters.get(0).toString() + ";");
        System.out.println("");
        NewickParser parser = new NewickParser();
        TreeNode root1 = parser.parseNewick(root);
        root1.printTree("", true);

    }

    private static int[] findMinDistance(List<Cluster> clusters, double[][] matrix) {
        int[] minIndex = new int[2];
        double minDistance = Double.MAX_VALUE;

        for (int i = 0; i < clusters.size(); i++) {
            for (int j = i + 1; j < clusters.size(); j++) {
                if (matrix[i][j] < minDistance) {
                    minDistance = matrix[i][j];
                    minIndex[0] = i;
                    minIndex[1] = j;
                }
            }
        }

        return minIndex;
    }

    private static double[][] updateMatrix(double[][] matrix, List<Cluster> clusters, Cluster newCluster, int[] index) {
        int size = clusters.size();
        double[][] newMatrix = new double[size][size];

        // Get the indices of the new cluster and the remaining clusters in the original matrix
        
        int newIndex = clusters.indexOf(newCluster);
        List<Integer> oldIndices = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            if (i != index[0] && i != index[1]) {
                oldIndices.add(i);
            }
        }

        // Calculate the distances between the new cluster and the remaining clusters
        for (int i = 0; i < oldIndices.size(); i++) {
            double newDistance = (matrix[oldIndices.get(i)][index[0]] + matrix[oldIndices.get(i)][index[1]]) / 2.0;
            newMatrix[newIndex][i] = newDistance;
            newMatrix[i][newIndex] = newDistance;
        }

        // Copy the distances between the remaining clusters from the original matrix
        for (int i = 0; i < oldIndices.size(); i++) {
            for (int j = 0; j < oldIndices.size(); j++) {
                newMatrix[i][j] = matrix[oldIndices.get(i)][oldIndices.get(j)];
            }
        }

        return newMatrix;
    }

    public static Cluster createNJTree(List<Cluster> clusters, double[][] distanceMatrix) {

        // While there is more than one cluster
        while (clusters.size() > 1) {
            int n = clusters.size(); // Get the current size of clusters and distanceMatrix
            double[][] totalDistances = new double[n][n];
            double smallestDistance = Double.MAX_VALUE;
            int[] smallestPair = new int[2];

            // Calculate total distances
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < n; k++) {
                        sum += distanceMatrix[i][k] + distanceMatrix[j][k];
                    }
                    totalDistances[i][j] = 0.5 * n * distanceMatrix[i][j] - sum;
                }
            }

            // Find the pair with the smallest total distance
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    if (totalDistances[i][j] < smallestDistance) {
                        smallestDistance = totalDistances[i][j];
                        smallestPair[0] = i;
                        smallestPair[1] = j;
                    }
                }
            }

            // Create a new cluster from the smallest pair
            Cluster cluster1 = clusters.get(smallestPair[0]);
            Cluster cluster2 = clusters.get(smallestPair[1]);
            double distance = distanceMatrix[smallestPair[0]][smallestPair[1]];
            Cluster newCluster = new Cluster(cluster1, cluster2, distance);

            // Remove the original clusters and add the new cluster
            clusters.remove(cluster2);
            clusters.remove(cluster1);
            clusters.add(newCluster);

            // Update the distance matrix
            distanceMatrix = updateMatrixNJ(distanceMatrix, clusters, newCluster, smallestPair);
        }

        // Return the final cluster
        Cluster finalCluster = clusters.get(0);
        System.out.println(finalCluster.toString() + ";");
        NewickParser parser = new NewickParser();
        TreeNode root = parser.parseNewick(finalCluster.toString());
        root.printTree("", true);
        return finalCluster;
    }

    private static double[][] updateMatrixNJ(double[][] matrix, List<Cluster> clusters, Cluster newCluster,
            int[] index) {
        int size = clusters.size();
        double[][] newMatrix = new double[size][size];

        // Get the indices of the new cluster and the remaining clusters in the original matrix
        
        int newIndex = clusters.indexOf(newCluster);
        List<Integer> oldIndices = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            if (i != index[0] && i != index[1]) {
                oldIndices.add(i);
            }
        }

        // Calculate the new distances
        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                if (i == newIndex || j == newIndex) {
                    int oldIndex = oldIndices.get(i == newIndex ? j : i);
                    double distance = (matrix[oldIndex][index[0]] + matrix[oldIndex][index[1]]
                            - matrix[index[0]][index[1]]) / 2;
                    newMatrix[i][j] = distance;
                    newMatrix[j][i] = distance;
                } else {
                    int oldIndex1 = oldIndices.get(i);
                    int oldIndex2 = oldIndices.get(j);
                    double distance = matrix[oldIndex1][oldIndex2];
                    newMatrix[i][j] = distance;
                    newMatrix[j][i] = distance;
                }
            }
        }

        // Return the updated matrix
        return newMatrix;
    }

}

class Cluster {
    int index;
    String name;
    double height = 0.0;
    Cluster left;
    Cluster right;
    int size = 1;

    public int getIndex() {
        return this.index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public String getName() {
        return this.name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public double getHeight() {
        return this.height;
    }

    public void setHeight(double height) {
        this.height = height;
    }

    public Cluster getLeft() {
        return this.left;
    }

    public void setLeft(Cluster left) {
        this.left = left;
    }

    public Cluster getRight() {
        return this.right;
    }

    public void setRight(Cluster right) {
        this.right = right;
    }

    public int getSize() {
        return this.size;
    }

    public void setSize(int size) {
        this.size = size;
    }

    Cluster(String name, int index) {
        this.name = name;
        this.index = index;
    }

    // leaf node
    Cluster(String name) {
        this.name = name;
        this.size = 1;
    }

    // merge clusters
    Cluster(Cluster left, Cluster right, double height) {
        this.name = "(" + left.name + ":" + (height - left.height) + "," + right.name + ":" + (height - right.height)
                + ")";
        this.height = height;
        this.left = left;
        this.right = right;
        this.size = left.size + right.size;
    }

    @Override
    public String toString() {
        if (left == null && right == null) {
            return name;
        } else {
            return "(" + left.toString() + ":" + height + "," + right.toString() + ":" + height + ")";
        }
    }
}

class Distance {
    double value;
    int row;
    int column;

    public double getValue() {
        return this.value;
    }

    public void setValue(double value) {
        this.value = value;
    }

    public int getRow() {
        return this.row;
    }

    public void setRow(int row) {
        this.row = row;
    }

    public int getColumn() {
        return this.column;
    }

    public void setColumn(int column) {
        this.column = column;
    }

    Distance(double value, int row, int column) {
        this.value = value;
        this.row = row;
        this.column = column;
    }
}

class TreeNode {
    String name;
    TreeNode left;
    TreeNode right;
    double length;

    TreeNode() {
        this.name = "";
        this.left = null;
        this.right = null;
    }

    TreeNode(String name) {
        this.name = name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setLeft(TreeNode left) {
        this.left = left;
    }

    public void setRight(TreeNode right) {
        this.right = right;
    }

    public void setLength(double length) {
        this.length = length;
    }

    public TreeNode getRight() {
        return this.right;
    }

    public TreeNode getLeft() {
        return this.left;
    }

    public String getName() {
        return this.name;
    }

    public double getLength() {
        return this.length;
    }



    public void printTree(String prefix, boolean isTail) {
        System.out.println(prefix + (isTail ? "└── " : "├── ") + this.getName() + ":" + this.getLength());
        List<TreeNode> children = new ArrayList<>();
        if (this.getLeft() != null)
            children.add(this.getLeft());
        if (this.getRight() != null)
            children.add(this.getRight());
        for (int i = 0; i < children.size() - 1; i++) {
            children.get(i).printTree(prefix + (isTail ? "    " : "│   "), false);
        }
        if (children.size() > 0) {
            children.get(children.size() - 1).printTree(prefix + (isTail ? "    " : "│   "), true);
        }
    }
}

class NewickParser {
    public TreeNode parseNewick(String newick) {
        return parseHelper(newick.substring(1, newick.length() - 1)); // remove the outer parentheses and the ending semicolon
                                                                      
    }

    public TreeNode parseHelper(String newick) {
        TreeNode tree = new TreeNode();
        if (newick.indexOf('(') == -1) { // leaf node
            String[] split = newick.split(":");
            tree.setName(split[0]);
            tree.setLength(Double.parseDouble(split[1].replace(")", "")));
        } else { // internal node
            int start = newick.indexOf('(') + 1;
            int count = 1;
            int commaIndex = -1;
            for (int i = start; i < newick.length(); i++) {
                if (newick.charAt(i) == '(') {
                    count++;
                } else if (newick.charAt(i) == ')') {
                    count--;
                } else if (newick.charAt(i) == ',') {
                    if (count == 1) {
                        commaIndex = i;
                        break;
                    }
                }
            }
            if (commaIndex == -1) {
                throw new RuntimeException("Invalid Newick format");
            }
            tree.setLeft(parseHelper(newick.substring(start, commaIndex)));
            tree.setRight(parseHelper(newick.substring(commaIndex + 1, newick.lastIndexOf(')'))));

            String lengthString = newick.substring(newick.lastIndexOf(':') + 1);
            if (lengthString.endsWith(")")) {
                lengthString = lengthString.substring(0, lengthString.length() - 1);
            }
            tree.setLength(Double.parseDouble(lengthString));
        }
        return tree;
    }

}