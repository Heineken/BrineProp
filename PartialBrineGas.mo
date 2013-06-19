within BrineProp;
package PartialBrineGas "nX_gas gases and Water"
redeclare record extends ThermodynamicState
    "a selection of variables that uniquely defines the thermodynamic state"
/*  AbsolutePressure p "Absolute pressure of medium";
  Temperature T(unit="K") "Temperature of medium";
  SpecificEnthalpy h "Specific enthalpy";
  Modelica.SIunits.SpecificEnthalpy h_g "Specific enthalpy gas phase";
  Modelica.SIunits.SpecificEnthalpy h_l "Specific enthalpy liquid phase";
  SpecificEntropy s "Specific entropy";
  Density d(start=300) "density";
  Real GVF "Gas Void Fraction";
  Density d_l(start=300) "density liquid phase";
  Density d_g(start=300) "density gas phase";
  Real x(start=0) "vapor quality on a mass basis [mass vapor/total mass]";
  AbsolutePressure p_H2O;
  AbsolutePressure p_gas[PartialBrine_ngas_Newton.nX_gas];
  AbsolutePressure[PartialBrine_ngas_Newton.nX_gas + 1] p_degas 
    "should be in SatProp, but is calculated in setState which returns a state";*/
   annotation (Documentation(info="<html>

</html>"));
end ThermodynamicState;

  extends BrineProp.PartialGasData;
  extends Modelica.Media.Interfaces.PartialMixtureMedium(
  reference_X=cat(1,fill(0,nX-1),{1}),
  nX = size(gasNames, 1));

 constant Modelica.SIunits.MolarMass[:] MM_vec;
 constant Modelica.SIunits.MolarMass[:] nM_vec "number of ions per molecule";

constant String gasNames[:]={""};
/*

  redeclare replaceable function extends specificHeatCapacityCp 
    "calculation of gas specific heat capacity"
  algorithm 
    //der Aufruf hier funzt seltsamerweise nur mit "BrineGas_3Gas.", bei rho, eta und lambda gehts ohne !?
      cp := specificHeatCapacityCp_pTX(
        state.p,
        state.T,
        state.X[end - nX+1:end]);
      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
  end specificHeatCapacityCp;

  redeclare replaceable function extends density "Density from state"

  algorithm 
    d := density_pTX(
        state.p,
        state.T,
        state.X[end - nX+1:end]);
  //  assert(lambda>0,"lambda="+String(lambda));
  end density;

  redeclare replaceable function extends dynamicViscosity 
    "Thermal conductivity of water"
  //very little influence of salinity  
  algorithm 
    eta := dynamicViscosity_pTX(
        state.p,
        state.T,
        state.X[end - nX+1:end]);
  //  assert(lambda>0,"lambda="+String(lambda));
  end dynamicViscosity;

  redeclare replaceable function extends thermalConductivity 
    "Thermal conductivity of water"
  //very little influence of salinity  
  algorithm 
    lambda := thermalConductivity_pTX(
        state.p,
        state.T,
        state.X[end - nX+1:end]);
  //  assert(lambda>0,"lambda="+String(lambda));
  if lambda<0 then
    Modelica.Utilities.Streams.print("lambda = " + String(lambda) + "W/(m·K)");
  end if;

  end thermalConductivity;
*/
  function specificHeatCapacityCp_pTX
    "calculation of gas specific heat capacity"
  //  import SG = Modelica.Media.IdealGases.SingleGases;
    input SI.Pressure p;
    input SI.Temp_K T;
    input SI.MassFraction X[nX]=reference_X "Mass fractions";
    output SI.SpecificHeatCapacity cp
      "Specific heat capacity at constant pressure";
  end specificHeatCapacityCp_pTX;

  redeclare function density_pTX "Density of an mixture of gases"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[nX] "Mass fractions";
    output Modelica.SIunits.Density d;
  end density_pTX;

  replaceable function dynamicViscosity_pTX
    "calculation of gas dynamic Viscosity"
    import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
    input SI.Pressure p;
    input SI.Temperature T;
    input SI.MassFraction[nX] X "Mass fractions of mixture";
    output SI.DynamicViscosity eta;

  end dynamicViscosity_pTX;

  replaceable function thermalConductivity_pTX
    "calculation of gas thermalConductivity"
    import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
    input SI.Pressure p;
    input SI.Temperature T;
    input SI.MassFraction[nX] X "Mass fractions of mixture";
    output SI.ThermalConductivity lambda;

  end thermalConductivity_pTX;
end PartialBrineGas;
