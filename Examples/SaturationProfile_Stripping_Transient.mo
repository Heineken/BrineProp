within BrineProp.Examples;
model SaturationProfile_Stripping_Transient
  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  // parameter SI.MassFraction[Medium.nXi] Xi={0,0,0,0,0,0,5e-3,5e-3};
  parameter SI.MassFraction[Medium.nXi] Xi={    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842,  0.00016889,  0.00073464, 6.5657e-005};

  constant Integer n=2 "number of elements";
  constant SI.Velocity u_B=0.2 "bubble velocity";

  Medium.BaseProperties[n] props;
//  Real ratio= Xi_g[end,7]/max(0.01,Xi_g[end,8]);
  parameter Real w_dg=1 "degassing coefficient";
 SI.MassFraction[n] x(each start=0) "actual gas fraction";
  SI.MassFraction[n,Medium.nXi] Xi_g;

  SI.Pressure[n] p=linspace(400e5,10e5,n);
  SI.Temperature[n] T=linspace(400,300,n);
protected
  Real[n,Medium.nXi] der_X_rho;
  SI.MassFraction[n,Medium.nXi] Xi_g_VLE;
//  String csvFilename = "pT_profil.csv";
  constant SI.Length depth[:] = linspace(-4100,0,n);
  constant SI.Length Delta_s=(depth[end]-depth[1])/n "cell length";
//  Real val[n,2];
initial equation
      for i in 1:n loop
          //props[i].Xi = Xi;
          //Xi_act[i,:] = Xi;
          //Xi_l[i,6:8] = props[i].Xi[6:8];
          Xi_g[i,:] = fill(0,Medium.nXi);
          x[i]=0;
      end for;
equation
//calculate VLE at in-situ conditions

     //total composition change by bubble flow
      //der(props[1].Xi*props[1].d) = der_X_rho[1,:];
      der_X_rho[1,:] = -Xi_g[1,:]*x[1]*props[1].d*u_B/Delta_s;
     //gas fraction change by bubble flow
      der(props[1].d*x[1]) = props[1].d*(props[1].state.x-x[1])*w_dg-x[1]*props[1].d*u_B/Delta_s;

    for i in 2:n loop
      //total composition change by bubble flow
      //der(props[i].Xi*props[i].d) = der_X_rho[i,:];
      der_X_rho[i,:] =(Xi_g[i-1,:]*x[i-1]*props[i-1].d - Xi_g[i,:]*x[i]*props[i].d)*u_B/Delta_s;
      //gas fraction change by bubble flow
      der(props[i].d*x[i]) = props[i].d*(props[i].state.x-x[i])*w_dg+(x[i-1]*props[i-1].d - x[i]*props[i].d)*u_B/Delta_s;

    end for;

    for i in 1:n loop
/*      val[i,:]=PowerPlant.getCSVData(DataFiles.readCSVmatrix(csvFilename), depth[i]);
      props[i].p=val[i,1];
      props[i].T=val[i,2];*/
      props[i].p=p[i];
      props[i].T=T[i];

      //props[i].Xi = Xi;
      //Xi_l[i,:] = Xi;
      //x[i]=0.5;

      // VLE
      Xi_g_VLE[i,:]= if x[i]>0 then PowerPlant.max_vec(0,(props[i].Xi - props[i].state.X_l[1:end-1] * (1-x[i]))/x[i]) else fill(0,Medium.nXi);

      //gas phase composition change by degassing + bubble rise
      der(Xi_g[i,:]*x[i]*props[i].d) = props[i].d*(Xi_g_VLE[i,:]-Xi_g[i,:])*x[i]*w_dg + der_X_rho[i,:];
      der(props[i].Xi[:]*props[i].d)=der_X_rho[i,:];
//      der(Xi_act[i,:]*props[i].d)=der_X_rho[i,:];
//      assert(max(X_g[i,:])<=1 and max(X_g[i,:])>=0, "X out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X_g[i,6:8]]))+"");
    end for;

algorithm
      print("Simulation time: " + String(time) + " s");
  annotation (experiment(StopTime=100, NumberOfIntervals=100),
                        __Dymola_experimentSetupOutput);
end SaturationProfile_Stripping_Transient;
