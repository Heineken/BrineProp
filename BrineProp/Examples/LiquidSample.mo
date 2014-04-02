within BrineProp.Examples;
model LiquidSample
  "Wie sieht die Gaszusammensetzung aus, wenn bei 10 bar eine Flüssigprobe gezogen und entgast wird"
  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Medium.BaseProperties props_string(phase=0);
  Medium.BaseProperties props_sample(phase=0);

/*
//  Real[:] y_l=if not max(props.X_l[6:8])>0 then fill(0,Medium.nX_gas) else props.X_l[6:8]./Medium.MM_gas / sum(props.X_l[6:8]./Medium.MM_gas) "mol fraction of dissolved gases";
  Real V_l = sum(props_insitu.X_l[6:8]./Medium.MM_gas)*22.4/props_insitu.X_l[end] 
    "Liter of dissolved gas per kg_brine would have after complete degassing at standard conditions";
  Real ratioGasLiquid_1 = V_l*props_insitu.d_l/1000 "complete degassing";

  Real ratioGasLiquid = props_STP.state.GVF/(1-props_STP.state.GVF) "STP";
  SI.Density d_l = props_insitu.d_l;

  Real val[2];
  Real val_static[2];
//  String csvFilename = "PowerPlant/pT_profil_static.csv";
  String csvFilename = "pT_profil.csv";
  Real depth= time;
//  Real depth=-4100;
  Real m_g=sum(props_insitu.X_l[6:8]);
  */
  Real[:] y_string=cat(1,props_string.p_gas,{props_string.p_H2O})/props_string.p;
  Real[:] y_sample=cat(1,props_sample.p_gas,{props_sample.p_H2O})/props_sample.p;
equation

//calculate VLE at in-situ conditions
//    val=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), time);
    props_string.p=10e5;
    props_string.T=330;
    props_string.Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842,  0.00016889,  0.00073464, 6.5657e-005}
    "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";

//calculate new VLE at static conditions after degassing at in-situ conditions
    props_sample.p = 1e5;
    props_sample.T = 300;
    props_sample.Xi = props_string.X_l[1:end-1];

algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
  annotation (experiment(
      StartTime=-4257,
      StopTime=0,
      NumberOfIntervals=100,
      Tolerance=0.001), __Dymola_experimentSetupOutput);
end LiquidSample;
