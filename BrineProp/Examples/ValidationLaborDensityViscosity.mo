within BrineProp.Examples;
model ValidationLaborDensityViscosity "with lab measurements"

  //see 6.1.1 and 6.2 in PhD-Thesis (http://nbn-resolving.de/urn:nbn:de:kobv:83-opus4-47126)
package Medium = Brine5salts(outOfRangeMode= 1) "TO AVOID ERROR MESSAGES";
  constant Real data[:,:]=DataFiles.readCSVmatrix(DataDir + "/BrineLabMeasurements_Hoffert.csv");
//  constant Real data[:,:]=DataFiles.readCSVmatrix(DataDir + "Zhang1996.csv");
  constant Integer n=size(data,1);
  Medium.BaseProperties[n] props;
  String csvFilename = "pT_profil.csv";
  Real depth= time;
  SI.Density[n] d;
  SI.DynamicViscosity[n] eta=Medium.dynamicViscosity(props.state);
  SI.ThermalConductivity[n] lambda=Medium.thermalConductivity(props.state);
equation

  for i in 1:n loop
    props[i].p=data[i,4]*1e5;
    props[i].T=data[i,5];
    props[i].Xi = {data[i,1]*SaltData.M_NaCl,data[i,2]*SaltData.M_KCl,data[i,3]*SaltData.M_CaCl2,0,0} /
        (1+data[i,1]*SaltData.M_NaCl+data[i,2]*SaltData.M_KCl+data[i,3]*SaltData.M_CaCl2);

    d[i]=props[i].d;
  end for;
algorithm
  Modelica.Utilities.Files.createDirectory(BrineProp.OutputDir);
  DataFiles.writeCSVmatrix(BrineProp.OutputDir + "/DensityValidationLabor_out.csv", {"b_NaCl","b_KCl2","b_CaCl2","p","T","eta","rho","lambda"}, cat(2,data[:,1:5],transpose({eta,d,lambda})), ";");

  annotation (experiment(__Dymola_NumberOfIntervals=1),
      __Dymola_experimentSetupOutput);
end ValidationLaborDensityViscosity;
