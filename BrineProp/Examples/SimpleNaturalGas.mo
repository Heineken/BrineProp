within BrineProp.Examples;
model SimpleNaturalGas
  "Comparison of ideal gasmixture model with SimpleNaturalGas model"
  parameter SI.Volume V=1 "Fixed size of volume 1 and volume 2";
/*  parameter SI.MassFlowRate m_flow_ext=0.01 
    "Fixed mass flow rate in to volume 1 and in to volume 2";
  parameter SI.EnthalpyFlowRate H_flow_ext=5000 
    "Fixed enthalpy flow rate in to volume and in to volume 2";
*/
  package Medium2 = Modelica.Media.IdealGases.MixtureGases.SimpleNaturalGas
    "Medium model";
  Medium2.BaseProperties medium2(p(start=1.e5, stateSelect=StateSelect.prefer),
     T(start=300, stateSelect=StateSelect.prefer),
     X(start={0.1,0.1,0.1,0.2,0.2,0.3})) "CH4,C2H6,C3H8,N-Butane,N2,CO2";
  Medium2.SpecificHeatCapacity cp2=Medium2.specificHeatCapacityCp(medium2.state);
  Medium2.DynamicViscosity eta2= Medium2.dynamicViscosity(medium2.state);
  Medium2.ThermalConductivity lambda2= Medium2.thermalConductivity(medium2.state);

  constant SI.MolarMass M_CO2 = Modelica.Media.IdealGases.SingleGases.CO2.data.MM
    "0.0440095 [kg/mol]";
  constant Integer nM_CO2 = 1 "number of ions per molecule";
   constant SI.MolarMass M_N2 = Modelica.Media.IdealGases.SingleGases.N2.data.MM
    "0.0280134 [kg/mol]";
  constant Integer nM_N2 = 1 "number of ions per molecule";
  constant SI.MolarMass M_CH4 = Modelica.Media.IdealGases.SingleGases.CH4.data.MM
    "0.01604246 [kg/mol]";
  constant Real[:] MM_gas = {M_CO2,M_N2,M_CH4};
  Real[:] X = {medium2.X[6],medium2.X[5],medium2.X[1]};
  Real[:] n_g = X./MM_gas;
  Real R=sum(Modelica.Constants.R*X./MM_gas);

  SI.Density d_g1= medium2.d;
  SI.Density d_g2 = medium2.p/(Modelica.Constants.R*medium2.T)*(n_g*MM_gas)/sum(n_g)
    "From PartialBrine_ngas_Newton.setState_pTX()";
  SI.Density d_g3 = medium2.p/(medium2.T*R);

//  medium2.p/(Modelica.Constants.R*medium2.T)*(X)/sum(X./MM_gas)

equation
  medium2.X={0.15,0,0,0,0.8,0.05};
  medium2.p=1e5 "*(1+time)";
  medium2.T=400 "*(1-time)";

end SimpleNaturalGas;
