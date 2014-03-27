within BrineProp.Examples;
model DensityValidationLabor2 "mit Messungen von Ulrike"
package Medium = Brine_5salts_noGas;
//  SI.Density d= props.d;  /**/

  constant Real data[:,:]=DataFiles.readCSVmatrix("e:\\francke\\Eigene Dateien\\GFZ\\Promotion\\Material\\BrineProperties\\Density\\Messungen\\Validierung Labor Ulrike\\Zusammenfassung_final.csv");
//  constant Real data[:,:]=DataFiles.readCSVmatrix("e:\\francke\\Eigene Dateien\\GFZ\\Promotion\\Material\\BrineProperties\\Density\\Messungen\\Validierung Labor Ulrike\\Zhang1996.csv");
  constant Integer n=size(data,1);
//  Integer nd=size(data,1);
  Medium.BaseProperties[n] props;
//  String csvFilename = "PowerPlant/pT_profil_static.csv";
  String csvFilename = "pT_profil.csv";
  Real depth= time;
  SI.Density[n] d;
  SI.DynamicViscosity[n] eta=Medium.dynamicViscosity(props.state);
  SI.ThermalConductivity[n] lambda=Medium.thermalConductivity(props.state);
equation

//calculate VLE at in-situ conditions
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), time);
  for i in 1:n loop
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth);
    props[i].p=data[i,4]*1e5;
    props[i].T=data[i,5];
//    props[i].Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842} "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//    props[i].Xi = {} "Elvira 2-2013 STP";
//    props[i].Xi = {0.0836450885895266,0.00252572141557871,0.122402440677945,0.000609924899400504,0.0021327472452556}       "Elvira 2-2013 23°C";
    props[i].Xi = {data[i,1]*SaltData.M_NaCl,data[i,2]*SaltData.M_KCl,data[i,3]*SaltData.M_CaCl2,0,0} /
        (1+data[i,1]*SaltData.M_NaCl+data[i,2]*SaltData.M_KCl+data[i,3]*SaltData.M_CaCl2);

    d[i]=props[i].d;
//    eta[i]=Medium.dynamicViscosity(props[i].state);
  end for;
algorithm
//  PowerPlant.saveCSV("FluMoCalc",{"rho"},transpose({d}));
  DataFiles.writeCSVmatrix("LabCalc.csv", {"b_NaCl","b_KCl2","b_CaCl2","p","T","eta","rho","lambda"}, cat(2,data[:,1:5],transpose({eta,d,lambda})), ";");
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end DensityValidationLabor2;
