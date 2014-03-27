within BrineProp.Examples;
model DensityCheckProbeMeasurement "Vergleich mit p-T-Log"
//package Medium = Brine_5salts_noGas;
package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;

  constant Real data[:,:]=DataFiles.readCSVmatrix("e:/francke/Eigene Dateien/GFZ/Promotion/Material/Loggingdaten/2011-09/pT_Profil_Einfahrt.csv");
  constant Integer n=size(data,1);
//  Integer nd=size(data,1);
  Medium.BaseProperties[n] props;
//  Real depth= time;
  SI.Density[n] d;
/*  SI.DynamicViscosity[n] eta=Medium.dynamicViscosity(props.state);
  SI.ThermalConductivity[n] lambda=Medium.thermalConductivity(props.state);
*/
parameter Real s=0.9365 "Vervielfachung des Salzgehalts";
equation

//calculate VLE at in-situ conditions
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), time);
  for i in 1:n loop
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth);
    props[i].p=data[i,1]*1e5;
    props[i].T=data[i,2]+273.15;
  props[i].Xi[:]={s*0.0830158933516516,s*0.00536930693677626,s*0.111697141833139,s*0.00112120951982728,s*0.00198855862851682, g*0.00018083,  g*0.00074452,  g*6.661e-005}
      "Elvira 9/11";
//        props[i].Xi = {0.083945671051201,0.00253479771131107,0.122842299461699,0*0.000612116692496665,0*0.00214041137028575,  0.00016883,  0.00073459, 6.5652e-005} "Elvira 2-2013 1.1775g/ml";
//    props[i].Xi = {0.083945671051201,0.00253479771131107,0.122842299461699,0*0.000612116692496665,0*0.00214041137028575}       "Elvira 2-2013 1.1775g/ml";
    d[i]=props[i].d;
//    eta[i]=Medium.dynamicViscosity(props[i].state);
  end for;
algorithm
//  PowerPlant.saveCSV("FluMoCalc",{"rho"},transpose({d}));
//  DataFiles.writeCSVmatrix("DensityCheckProbeMeasurement.csv", {"p","T","rho"}, cat(2,data[:,1:2],transpose({d})), ";");
  DataFiles.writeCSVmatrix("DensityCheckProbeMeasurement.csv", {"rho"}, transpose({d}), ";");

//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end DensityCheckProbeMeasurement;
