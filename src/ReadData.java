import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class ReadData {
    public static void main(String[] args) {
        try {
            // Open the file containing the text
            File file = new File("src/data.txt");
            Scanner scanner = new Scanner(file);

            // Skip the first few lines (header)
            for (int i = 0; i < 4; i++) {
                scanner.nextLine();
            }

            double[] presArray = new double[114];
            double[] highArray = new double[114];
            double[] tempArray = new double[114];
            double[] dwptArray = new double[114];
            double[] mrArray = new double[114];

            // Read and process the data lines
            int index=0;
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                String[] tokens = line.split("\\s+"); // Split the line by whitespace
                // Process the tokens (e.g., convert them to appropriate data types)
                for (String token : tokens) {
                    System.out.print(token + " ");
                }
                System.out.println(); // Print a new line after each row
                presArray[index] = Double.parseDouble(tokens[1]);
                highArray[index] = Double.parseDouble(tokens[2]);
                tempArray[index] = Double.parseDouble(tokens[3]);
                dwptArray[index] = Double.parseDouble(tokens[4]);
                mrArray[index] = Double.parseDouble(tokens[6]);

                // Increment index
                index++;

                
            }
            System.out.println(index);

            // Close the scanner
            scanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}

