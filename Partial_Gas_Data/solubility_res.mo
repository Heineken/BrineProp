within BrineProp.Partial_Gas_Data;
function solubility_res
extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
  input Modelica.SIunits.MolarMass MM_vec[:] "molar masses of components";
  /**/
  input Partial_Gas_Data.partial_solubility_pTX solufun;
//  input Modelica.Icons.Function solufun;
  input Modelica.SIunits.MassFraction c_gas;
//  input Modelica.SIunits.Pressure p_gas=u;
protected
  Modelica.SIunits.MassFraction solu;
algorithm
    y:=c_gas-solufun(p=p,T=T,X=X,MM_vec=MM_vec,p_gas=u) "*X[end]";
end solubility_res;
