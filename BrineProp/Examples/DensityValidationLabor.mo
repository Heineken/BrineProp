within BrineProp.Examples;
model DensityValidationLabor "mit Messungen von Ulrike"
  //p:\GEOFLUIDS\Fluidphysics\Daten\Dichte Reihe 8 20°C-80°C.xlsx

package Medium = Brine_5salts;
//  SI.Density d= props.d;  /**/

  Medium.BaseProperties[n_b,n_T] props;
//  String csvFilename = "PowerPlant/pT_profil_static.csv";
  String csvFilename = "pT_profil.csv";
  SI.Density[n_b,n_T] d;
  constant SI.Temp_C T[:]=273.16.+{20,25,30,40,50,60,70,80};
  constant Integer n_T = size(T,1);
  constant Real b_NaCl[:] = {0,0.381608,0.755841,1.120155,1.471588,1.822351,2.168757};
  constant Integer n_b = size(b_NaCl,1);
  constant Real b_CaCl2[n_b] = {0,0.314340,0.622386,0.922500,1.211921,1.500710,1.786067};

equation
//calculate VLE at in-situ conditions
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), time);
for i in 1:n_b loop
  for j in 1:n_T loop
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth);
    props[i,j].p=1.01325e5;
    props[i,j].T=T[j]+273.16;
//    props[i].Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842} "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//    props[i].Xi = {} "Elvira 2-2013 STP";
//    props[i].Xi = {0.0836450885895266,0.00252572141557871,0.122402440677945,0.000609924899400504,0.0021327472452556}       "Elvira 2-2013 23°C";
  props[i,j].Xi = {b_NaCl[i]*SaltData.M_NaCl,0,b_CaCl2[i]*SaltData.M_CaCl2,0,0} /
        (1+b_NaCl[i]*SaltData.M_NaCl+b_CaCl2[i]*SaltData.M_CaCl2)
        "NaCl, KCl, CaCl2, MgCl2, SrCl2";

    d[i,j]=props[i,j].d;
  end for;
end for;

algorithm
  DataFiles.writeCSVmatrix("LabCalc.csv", cat(1,{"b_NaCl","b_CaCl2"},String(T)), cat(2,transpose({b_NaCl,b_CaCl2}),d), ";");
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end DensityValidationLabor;
