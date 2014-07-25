within BrineProp;
partial package PartialBrineGas "Medium template for gas mixture of nX_gas gases and water based on PartialMixtureMedium"
  //TODO add inverse functions


  extends BrineProp.GasData;


  extends Modelica.Media.Interfaces.PartialMixtureMedium(reference_X=cat(1,fill(0,nX-1),{1}));

 constant Boolean ignoreNoCompositionInBrineGas=false;


redeclare record extends ThermodynamicState
end ThermodynamicState;


 redeclare model extends BaseProperties "Base properties of medium"
    import BrineProp;

 //  SI.Pressure p_H2O;
  //   BrineProp.Partial_Units.Molality y_vec[:]=BrineProp.massToMoleFractions();
   SI.MoleFraction y_vec[:]=Utilities.massToMoleFractions(X,MM_vec);
 equation
    //   assert(nX_gas==2,"Wrong number of gas mass fractions specified (2 needed - CO2,N2)");
 //  assert(max(X)<=1 and min(X)>=0, "X out of range [0...1] = "+PowerPlant.vector2string(X)+" (saturationPressure_H2O())");
     MM = y_vec*MM_vec;
 //  R  = Modelica.Constants.R/MM;
   u = h - p/d;

 //  (h,x,d,d_g,d_l) = specificEnthalpy_pTX(p,T,X) unfortunately, this is not invertable;

     h = specificEnthalpy_pTX(p,T,X);
 //    d=density_pTX(p,T,X);
     (d,R) = density_pTX(p,T,X);
     state=ThermodynamicState(p=p,T=T,X=X);
 end BaseProperties;
 constant SI.MolarMass[:] MM_vec;
 constant Integer[:] nM_vec "number of ions per molecule";

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
    print("lambda = " + String(lambda) + "W/(m.K)");
  end if;

  end thermalConductivity;
*/


  replaceable function specificHeatCapacityCp_pTX
  "calculation of gas specific heat capacity"
  //  import SG = Modelica.Media.IdealGases.SingleGases;
    input SI.Pressure p;
    input SI.Temp_K T;
    input SI.MassFraction X[nX]=reference_X "Mass fractions";
    output SI.SpecificHeatCapacity cp
    "Specific heat capacity at constant pressure";
  end specificHeatCapacityCp_pTX;


  redeclare replaceable function density_pTX "Density of a mixture of gases"
    input SI.Pressure p;
    input SI.Temp_K T;
    input MassFraction X[nX] "Mass fractions";
    output SI.Density d;
    output SpecificHeatCapacity R_gas;
  //algorithm

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


  annotation (Documentation(info="<html>
<p>Ideal mixture of gases.</p>
<h5>Usage</h5>
<p>This partial package cannot be used as is. See <a href=\"Modelica://BrineProp.Examples.BrineGas\">BrineProp.Examples.BrineGas</a> or info of <a href=\"Modelica://BrineProp.BrineGas_3Gas\">BrineProp.BrineGas_3Gas</a> for examples.</p>
</html>"));
end PartialBrineGas;
