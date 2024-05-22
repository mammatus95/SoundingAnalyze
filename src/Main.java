public class Main {

    public static void main(String[] args) {
        String filename = "src/data.txt";
        String station = "10184";
        Sounding Data = new Sounding(filename,station);
        Data.MetaInfo();
        System.out.println(Data.getDatapointCount());
        System.out.println(Data.getaverage_pottemp());
    }

}
