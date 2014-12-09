within BrineProp;
partial package PartialBrineGas "Medium template for gas mixture of nX_gas gases and water based on PartialMixtureMedium"
  //TODO add inverse functions


 extends GasData;


 extends Modelica.Media.Interfaces.PartialMixtureMedium(reference_X=cat(1,fill(0,nX-1),{1}));


 redeclare record extends ThermodynamicState
 //Dummy for OM
 end ThermodynamicState;


 redeclare replaceable model extends BaseProperties "Base properties of medium"
    SI.MoleFraction y_vec[:]=Utilities.massToMoleFractions(X,MM_vec);
 equation
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
  //    import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
    input SI.Pressure p;
    input SI.Temperature T;
    input SI.MassFraction[nX] X "Mass fractions of mixture";
    output SI.DynamicViscosity eta;

  end dynamicViscosity_pTX;


  replaceable function thermalConductivity_pTX
  "calculation of gas thermalConductivity"
  //    import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
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
