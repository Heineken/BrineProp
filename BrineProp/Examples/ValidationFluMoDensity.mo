within BrineProp.Examples;
model ValidationFluMoDensity "Validation with density measured by Flumo"
  //working directory must be parent folder of \Data in BrineProp package
  //multitude of values from online measurement have been condensed by finding "BestFitPlane" and extracting the interpolation points
  //see 6.1.2 in PhD-Thesis (http://nbn-resolving.de/urn:nbn:de:kobv:83-opus4-47126)
package Medium = Brine_5salts;

  constant Real data[:,:]=DataFiles.readCSVmatrix("Data\FluMoFit.csv");
  constant Integer n=size(data,1);
  Medium.BaseProperties[n] props;
  Real depth= time;
  SI.Density[n] d;
equation

//calculate VLE at in-situ conditions
  for i in 1:n loop
    props[i].p=data[i,3]*1e5;
    props[i].T=data[i,2]+273.15;
    props[i].Xi = {0.083945671051201,0.00253479771131107,0.122842299461699,0*0.000612116692496665,0*0.00214041137028575}
      "Feldbusch 2-2013 1.1775g/ml";
    d[i]=props[i].d;
  end for;
algorithm
  DataFiles.writeCSVmatrix("DensityValidationFluMo_out.csv", {"rho_meas/10kg/m^3","p/bar","T/C","rho_calc"}, cat(2,data,transpose({d})), ";");

  annotation (experiment(__Dymola_NumberOfIntervals=1),
      __Dymola_experimentSetupOutput);
end ValidationFluMoDensity;
