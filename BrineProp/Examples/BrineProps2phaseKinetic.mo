within BrineProp.Examples;
model BrineProps2phaseKinetic "Slow degassing"
package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;

  Medium.BaseProperties props(phase=0,n_g_norm_start=fill(0.5,
                                                             Medium.nX_gas+1));
  parameter SI.MassFraction[Medium.nXi] Xi={0,0,0,0,0,0,5e-3,5e-3}
    "fill(0,Medium.nXi)";
  parameter SI.Pressure p=400e5;
  parameter Real w_dg=0.5;
//  SI.Density d= props.d;  /**/
  SI.MassFraction x(start=0) "actual gas fraction";
//  SI.MassFraction[Medium.nXi] Xi_l(start=Xi);
  SI.MassFraction[:] Xi_g(start=fill(0,Medium.nXi));
  SI.MassFraction[:] Xi_g_VLE(start=fill(0,Medium.nXi));
initial equation
// props.T = 273.16+60;
  x=0;
//  Xi_l=Xi;
  Xi_g=fill(0,Medium.nXi);
equation
 props.p = p;
 props.T = 400;
 props.Xi = Xi;
// der(Xi_l) = (props.state.X_l[1:end-1]-Xi_l)*w_dg;
// der(Xi_l*(1-x)) = (props.state.X_l[1:end-1]-Xi_l)*(1-x)*w_dg;

 Xi_g_VLE= if x>0 then PowerPlant.max_vec(0,(props.Xi - props.state.X_l[1:end-1] * (1-x))/x) else fill(0,Medium.nXi);
// der(Xi_g*x) = (props.state.X_l[1:end-1]-Xi_l)*(1-x)*w_dg
 der(Xi_g*x) =     (Xi_g_VLE-Xi_g)*x*w_dg;

// der(x) = -sum(der(Xi_l*(1-x)));
 der(x) = (props.state.x-x)*w_dg;
algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
//  print(Modelica.Math.Matrices.toString(transpose([props.Xi])));
end BrineProps2phaseKinetic;
