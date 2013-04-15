within BrineProp.Examples;
model DensityValidationFluMo "mit Messwerten von Elvira"
//Unzahl an Werten wurde verdichtet mit Guidos Prog (BestFitPlane), dann Stützpunkte berechnet
package Medium = Brine_5salts_noGas;
//  Modelica.SIunits.Density d= props.d;  /**/

  constant Real data[:,:]=DataFiles.readCSVmatrix("FluMoFit.csv");
  constant Integer n=size(data,1);
//  Integer nd=size(data,1);
  Medium.BaseProperties[n] props;
//  String csvFilename = "PowerPlant/pT_profil_static.csv";
  String csvFilename = "pT_profil.csv";
  Real depth= time;
  Modelica.SIunits.Density[n] d;
equation

//calculate VLE at in-situ conditions
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), time);
  for i in 1:n loop
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth);
    props[i].p=data[i,1]*1e5;
    props[i].T=data[i,2]+273.15;
//    props[i].Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842} "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//    props[i].Xi = {} "Elvira 2-2013 STP";
//    props[i].Xi = {0.0836450885895266,0.00252572141557871,0.122402440677945,0.000609924899400504,0.0021327472452556}       "Elvira 2-2013 23°C";
    props[i].Xi = {0.083945671051201,0.00253479771131107,0.122842299461699,0.000612116692496665,0.00214041137028575}
      "Elvira 2-2013 1.1775g/ml";

    d[i]=props[i].d;
  end for;
algorithm
//  PowerPlant.saveCSV("FluMoCalc",{"rho"},transpose({d}));
  DataFiles.writeCSVmatrix("FluMoCalc.csv", {"rho"}, transpose({d}), ";");
//  Modelica.Utilities.Streams.print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  Modelica.Utilities.Streams.print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end DensityValidationFluMo;
