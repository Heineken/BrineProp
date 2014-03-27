within BrineProp.Examples;
model SaturationProfile_Stripping
  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  constant Integer n=2;
  Medium.BaseProperties[n] props;
//  SI.Density d= props.d;  /**/
  Real ratio;
//  Real[n] ratio;
  SI.MassFraction[n,Medium.nX] X_g;
protected
  SI.MassFraction[n] x;
  SI.MassFraction[n,Medium.nX] X;
  SI.Pressure[n] p=linspace(400e5,10e5,n);
  SI.Temperature[n] T=linspace(400,300,n);
  String csvFilename = "pT_profil.csv";
  SI.Length depth[:] = linspace(-4100,0,n);
  Real val[n,2];
equation
//calculate VLE at in-situ conditions

    val[1,:]=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth[1]);
    props[1].p=val[1,1];
    props[1].T=val[1,2];
/*    props[1].p=p[1];
    props[1].T=T[1];*/
    props[1].Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842,  0.00016889,  0.00073464, 6.5657e-005};
//    props[1].Xi = {    0,0,0,0,0, 0,  5e-3, 5e-3} "NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4";
    X[1,:]=props[1].X;
    for i in 2:n loop
/*      props[i].p=p[i];
      props[i].T=T[i];*/

      val[i,:]=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth[i]);
      props[i].p=val[i,1];
      props[i].T=val[i,2];

//      props[i].Xi = (props[i-1].state.x*PowerPlant.max_vec(0,(props[i-1].Xi-props[i-1].X_l[1:end-1])) + props[i-1].Xi)/sum((-props[i-1].X_l.+1).+1);
//      props[i].Xi = (2*props[1].Xi-props[i-1].X_l[1:end-1]*(1-props[i-1].state.x))/(1+props[i-1].state.x)
      X[i,:] = (props[1].X+X_g[i-1,:]*props[i-1].state.x)/(1+props[i-1].state.x);
      props[i].Xi = X[i,1:end-1];
//      props[i].Xi = (props[1].Xi+X_g[i-1,1:end-1]*props[i-1].state.x)/(1+props[i-1].state.x)       "Basic composition + gas phase from below";

    end for;

    for i in 1:n loop
      x[i]=props[i].state.x;
      X_g[i,:]= if props[i].state.x>0 then PowerPlant.max_vec(0,(props[i].X -props[i].state.X_l * (1- props[i].state.x))/props[i].state.x) else fill(0,Medium.nX);
//      ratio[i]= X_g[i,7]/X_g[i,8];
//      assert(max(X_g[i,6:8])<=1 and max(X_g[i,6:8])>0, "X out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X_g[i,6:8]]))+"");
    end for;
      ratio= X_g[end,7]/X_g[end,8];

algorithm
//  Modelica.Utilitiprops[1].Xi "es.Streams.print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
//print(Modelica.Math.Matrices.toString(transpose([(props[1].Xi+X_g[n-1,1:end-1]*props[n-1].state.x)/(1+props[n-1].state.x)])));
  annotation (experiment(
      StartTime=-4257,
      StopTime=0,
      NumberOfIntervals=100,
      Tolerance=0.001), __Dymola_experimentSetupOutput);
end SaturationProfile_Stripping;
