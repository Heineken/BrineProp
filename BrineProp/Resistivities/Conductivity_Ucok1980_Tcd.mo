within BrineProp.Resistivities;
function Conductivity_Ucok1980_Tcd
  //returns conductivity instead of resistivity to allow easy handling of zero salinities
  //salinity validity can't be checked here because it's given in weight-%, which is unknown here
  //measured at 31 MPa
  input SI.Temperature T_K;
  input Types.Molarity[:] c_vec
    "molar concentration (molarity) in mol / litre solution";
  input SI.Density d;
  input Integer AssertLevel=2
    "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";
  output SI.Conductivity[3] gamma;
protected
  Types.Molarity c;
  SI.Temp_C T_C=T_K-273.15;
  Real B_NaCl[:,:]= transpose([3.47,-6.65,2.633;-59.21,198.1,-64.8;0.4551,-0.2058,0.005799;-9.35E-05,7.37E-05,6.74E-05;-1.77E-06,8.768e-7,-2.14E-07]);
  Real B_KCl[:,:] = transpose([5.783,-6.607,1.665;-59.23,149.7,-31.21;0.2051,0.1064,-0.03418;1.82E-04,-7.04E-04,1.54E-04;-1.09E-06,1.08E-06,-1.95E-07]);
  Real B_CaCl2[:,:]=transpose([-34.62,24.64,-3.907;780.3,-492.3,64.59;1.05,-0.5922,0.06735;-0.002459,0.001461,-1.22E-04;9.99E-07,-7.11E-07,-4.73E-09]);
  Real C_[1,3];
  Real T_vec[5]={1,1/T_C,T_C,T_C^2,T_C^3};
  Real D[1,5];
  Real BB[3,3,5];
  constant AssertionLevel aLevel = if AssertLevel==1 then AssertionLevel.warning else AssertionLevel.error;
algorithm
  if AssertLevel>0 then
    assert(T_C >= 22 and T_C<= 375,"Temperature  T=" + String(T_C) + " C is out of validity range.\nTo ignore validity range call with AssertLevel=1 or 0",aLevel);
  end if;  BB[1,:,:] :=B_NaCl;
  BB[2,:,:] :=B_KCl;
  BB[3,:,:] :=B_CaCl2;

//print("c_vec="+Modelica.Math.Matrices.toString(transpose([c_vec])));
//print("T="+Modelica.Math.Matrices.toString(transpose([T_vec])));

  for i in 1:3 loop
    c:=c_vec[i];
    if not c>0 then
      gamma[i] := 0;
    else
      C_:=[c,c*sqrt(c),c*c*log(c)];
//      print("C="+Modelica.Math.Matrices.toString(C_));
      D:=C_*BB[i, :, :];
      gamma[i] :=D*T_vec*{1};
//      print("D="+Modelica.Math.Matrices.toString(D));
    end if;
  end for;

//  print("gamma="+Modelica.Math.Matrices.toString(transpose([gamma])));
end Conductivity_Ucok1980_Tcd;
