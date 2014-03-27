within BrineProp.PartialGasData;
function solubility_res
extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
  input SI.MolarMass MM_vec[:] "molar masses of components";
  /**/
  input partial_solubility_pTX solufun;
//  input Modelica.Icons.Function solufun;
  input SI.MassFraction c_gas;
//  input SI.Pressure p_gas=u;
protected
  SI.MassFraction solu;
algorithm
    y:=c_gas-solufun(p=p,T=T,X=X,MM_vec=MM_vec,p_gas=u) "*X[end]";
//    print("T="+String(T-273.16)+"°C");
end solubility_res;
