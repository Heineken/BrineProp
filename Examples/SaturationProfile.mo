within BrineProp.Examples;
model SaturationProfile
  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Medium.BaseProperties props_insitu(phase=0);
  Medium.BaseProperties props_static(phase=0);
  Medium.BaseProperties props_STP(phase=0);

//  SI.Density d= props.d;  /**/

//  Real[:] y_l=if not max(props.X_l[6:8])>0 then fill(0,Medium.nX_gas) else props.X_l[6:8]./Medium.MM_gas / sum(props.X_l[6:8]./Medium.MM_gas) "mol fraction of dissolved gases";
  Real V_l = sum(props_insitu.X_l[6:8]./Medium.MM_gas)*22.4/props_insitu.X_l[end]
    "Liter of dissolved gas per kg_brine would have after complete degassing at standard conditions";
  Real ratioGasLiquid_1 = V_l*props_insitu.d_l/1000 "complete degassing";

  Real ratioGasLiquid = props_STP.GVF/(1-props_STP.GVF) "STP";
  SI.Density d_l = props_insitu.d_l;
  Real ratioGasLiquid_vec[:] = ratioGasLiquid * props_STP.p_gas/props_STP.p;

  Real val[2];
  Real val_static[2];
//  String csvFilename = "PowerPlant/pT_profil_static.csv";
  String csvFilename = "e:/francke/Eigene Dateien/Dymola/PowerPlant/pT_profil_idle.csv";
  Real depth= time;
//  Real depth=-4100;
  Real m_g=sum(props_insitu.X_l[6:8]);
equation

//calculate VLE at in-situ conditions
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), time);
    val=PowerPlant.getDataFromMatrix(DataFiles.readCSVmatrix(csvFilename), depth);
    props_insitu.p=val[1];
    props_insitu.T=val[2];
//    props_insitu.Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842,  2*0.00016889,  2*0.00073464, 2*6.5657e-005} "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//    props_insitu.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0,2*7.2426359111e-05,2*0.000689505657647,2*6.14906384726e-05}     "Elvira 2-2013 1.1775g/ml V2";
    props_insitu.Xi = {0.0839112,0.00253376,0.122792,0,0,3*5.43473e-005,3*0.000414069,3*7.97687e-005}
    "Elvira 3";

    //calculate new VLE at static conditions after degassing at in-situ conditions
    val_static=PowerPlant.getDataFromMatrix(DataFiles.readCSVmatrix(csvFilename), depth);
    props_static.p = val_static[1];
    props_static.T = val_static[2];
    props_static.Xi = props_insitu.X_l[1:end-1];

//calculate degassing at STP
    props_STP.p = 1.01325e5;
    props_STP.T = 273.16;
    props_STP.Xi = props_insitu.X_l[1:end-1];
//    props_STP.Xi = props_static.X_l[1:end-1];

algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
  annotation (experiment(
      StartTime=-4257,
      StopTime=-1174,
      __Dymola_NumberOfIntervals=100,
      Tolerance=0.001), __Dymola_experimentSetupOutput);
end SaturationProfile;
