within BrineProp.Examples;
model BrineProps1Phase_Convect
package Medium = Brine_5salts;
  parameter Integer n=1;
  Medium.BaseProperties[n+1] props;

//  SI.Density d= props.d;
  SI.SpecificEnthalpy h[n];
//  SI.Temp_C T_C= props.T-273.15;
  parameter SI.Mass m=1;
  parameter SI.MassFlowRate m_dot=1;
  parameter SI.HeatFlowRate Q_dot=1e4;
equation
  props[1].p = 15e5;
  props[1].T = 273.16+99.61;
//  props.h = 379778;
//  props.p = (10+time)*1.01325e5 "STP";
/*
 props.T = 273.16+60+time*80;
  props.p = 9.13e5;
 */
   props[1].Xi = {    0.081109,   0.0047275,     0.12519,   0.0013225,  0.0023842}
    "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";

for i in 1:n loop
//  h[i]=(props[i].h+props[i+1].h)/2;
  h[i]=props[i+1].h;
  props[i+1].p=props[i].p;
  props[i+1].Xi=props[i].Xi;

    m*der(h[i])=Q_dot-m_dot*(props[i+1].h-props[i].h);

end for;
initial equation
for i in 1:n loop
  props[i+1].T = 273.16+60;
end for;/**/
equation
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
algorithm
end BrineProps1Phase_Convect;
