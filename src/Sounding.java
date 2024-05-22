import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.lang.Math;

public class Sounding {
    static double C2K = 273.15;
    static double ROCP = 0.2855721; //0.28571426
    private String stationID;

    private int index=0;

    // Data arrays
    private ArrayList<Double> presArray = new ArrayList<>();
    private ArrayList<Double> highArray = new ArrayList<>();
    private ArrayList<Double> tempArray = new ArrayList<>();
    private ArrayList<Double> dwptArray = new ArrayList<>();
    private ArrayList<Double> mrArray = new ArrayList<>();
    private ArrayList<Double> potArray;

    // Constructor
    public Sounding (String filename, String station){

        this.stationID = station;

        try {
            // Open the file containing the text
            File file = new File("src/data.txt");
            Scanner scanner = new Scanner(file);

            // Skip the first few lines (header)
            for (int i = 0; i < 4; i++) {
                scanner.nextLine();
            }
            // Read and process the data lines
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                String[] tokens = line.split("\\s+"); // Split the line by whitespace
                // Process the tokens (e.g., convert them to appropriate data types)

                presArray.add(Double.parseDouble(tokens[1]));
                highArray.add(Double.parseDouble(tokens[2]));
                tempArray.add(Double.parseDouble(tokens[3]));
                dwptArray.add(Double.parseDouble(tokens[4]));
                mrArray.add(Double.parseDouble(tokens[6]));
                // Increment index
                index++;
            }
            // Close the scanner
            scanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        potArray = potential_temp(presArray, tempArray, 1000.0 );

    }

    public void MetaInfo(){
        System.out.println("Station: " + this.stationID);
    }

    public int getDatapointCount(){
        return this.index;
    }

    public double getaverage_pottemp(){
        return this.potArray.stream()
                .mapToDouble(Double::doubleValue)
                .average()
                .orElse(0.0);
    }
    // potential temperature
    private ArrayList<Double> potential_temp (ArrayList<Double> presArray, ArrayList<Double> tempArray, double p0){
        ArrayList<Double> potArray = new ArrayList<>();
        for (int i = 0; i < presArray.size(); i++) {
            potArray.add((tempArray.get(i)+C2K) * Math.pow((p0/presArray.get(i)), ROCP));
        }
        return potArray;
    }

}
