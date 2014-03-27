within BrineProp;
package BrineGas_3Gas
  "BrineGas_3Gas - water saturated gas mixture CO2+N2+CH4+H2O"
  //returns properties for given composition when _pTX functions are called directly
  //returns properties for given gas composition + saturated water when called via state functions ()e.g. density)
  //  TODO: Why?
  //speedup: calculate water saturated composition externally once and pass on
  extends PartialBrineGas(
    final substanceNames={"carbondioxide","nitrogen","methane","water"},
    final MM_vec = {M_CO2,M_N2,M_CH4, M_H2O},
    final nM_vec = {nM_CO2,nM_N2,nM_CH4, nM_CH4});

  constant Boolean waterSaturated=false "activates water saturation";

  replaceable function waterSaturatedComposition_pTX
    "calculates the water saturated mass vector for a given Temperature"
  //saturates the mixture with water
    extends Modelica.Icons.Function;
    input SI.Pressure p;
    input SI.Temperature T;
    input SI.MassFraction[nX] X_in "Mass fractions of mixture";
    output SI.MassFraction X[nX] "Mass fractions of mixture";
  protected
      SI.Pressure y_H2O=Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T)/p;
      Real y[nX] "mole fractions";
  algorithm
    if debugmode then
      print("Running waterSaturatedComposition_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X_in]))+")");
  //    print("y_H2O"+String(y_H2O)+", X[end]="+String(X_in[end]));
    end if;

    y:= X_in./MM_vec;

    if y_H2O<1 and X_in[end]>0 then
      //print(""+String(y_H2O));
      y:=cat(1,y[1:nX-1]/(sum(y[1:nX-1]))*(1-y_H2O), {y_H2O})
        "gases + fixed water fraction";
    else
      y:=cat(1,fill(0,nX-1), {y_H2O}) "only water vapour";
    end if;
    X:=y.*MM_vec "convert to mass fractions";
    X:=X/sum(X) "normalize";

  end waterSaturatedComposition_pTX;

  redeclare function extends density "water-saturated density from state"

  algorithm
    d := density_pTX(
      p=state.p,
      T=state.T,
      X= if waterSaturated then
        waterSaturatedComposition_pTX(state.p,state.T,state.X)
    else state.X);
  //  else state.X[end - nX + 1:end]);
  //      waterSaturatedComposition_pTX(state.p,state.T,state.X[end - nX+1:end])
  //  assert(lambda>0,"lambda="+String(lambda));
  end density;

  redeclare function extends density_pTX
    "Density of an ideal mixture of ideal gases"
  //    SpecificHeatCapacity R_gas; //= sum(Modelica.Constants.R*X ./ MM_vec);
  //    MassFraction[:] X_=cat(1,fill(nX-1,0),{1});
  //    MassFraction[size(X,1)] X_=cat(1,fill(size(X,1),0),{1});
  algorithm
  /*  if not R_gas >0 then
    print("R_gas="+String(R_gas)+", (MM="+Modelica.Math.Matrices.toString({MM_vec})+", X="+Modelica.Math.Matrices.toString({X})+")");
  end if;*/
  //  print("size(X_,1)="+String(size(X_,1))+",size(X,1)="+String(size(X,1)));
    if debugmode then
      print("Running density_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;
  // assert(min(X)>0,"Cannot calculate with empty composition.");
    if not min(X)>0 then
      print("No gas composition, assuming water vapour.(BrineProp.BrineGas_3Gas.density_pTX)");
  //  else
  //    X_:=X;
    end if;
    R_gas :=Modelica.Constants.R*sum(
    cat(1,X[1:end-1],{if min(X)>0 then X[end] else 1})
     ./ MM_vec);
  /*  if waterSaturated then
    R_gas :=sum(Modelica.Constants.R*waterSaturatedComposition_pTX(
        p,
        T,
        X[end - nX + 1:end]) ./ MM_vec);
    d :=p/(T*R_gas);
  else*/
        d :=p/(T*R_gas);
  //  end if;
  //  print("d="+String(d)+" kg/m³");
  end density_pTX;

  redeclare function extends specificHeatCapacityCp
    "water-saturated heat capacity of gas phase"
  algorithm
       cp := specificHeatCapacityCp_pTX(
          p=state.p,
          T=state.T,
          X= if waterSaturated then
       waterSaturatedComposition_pTX(state.p,state.T,state.X)
    else state.X);
  //  else state.X[end - nX + 1:end]);

  end specificHeatCapacityCp;

  function extends specificHeatCapacityCp_pTX
    "calculation of specific heat capacities of gas mixture"
    import Modelica.Media.IdealGases.Common.SingleGasNasa;
    import Modelica.Media.IdealGases.Common.SingleGasesData;

    import Modelica.Media.Water;
  protected
    SI.SpecificHeatCapacity cp_vec[:]={
      SingleGasNasa.cp_T(data=SingleGasesData.CO2,T=T),
      SingleGasNasa.cp_T(data=SingleGasesData.N2,T=T),
      SingleGasNasa.cp_T(data=SingleGasesData.CH4,T=T),
      Water.IF97_Utilities.cp_pT(min(p,Water.IF97_Utilities.BaseIF97.Basic.psat(T)-1),T=T)}
      "below psat -> gaseous";
  algorithm
    if debugmode then
      print("Running specificHeatCapacityCp_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;

    if not min(X)>0 then
      print("No gas composition, assuming water vapour.(BrineProp.BrineGas_3Gas.specificHeatCapacityCp_pTX)");
    end if;

  /*  if waterSaturated then
    cp := cp_vec * waterSaturatedComposition_pTX(p,T,X[end - nX+1:end]);
  else */
  //    cp := cp_vec * X[end - nX+1:end];
      cp := cp_vec * cat(1,X[1:end-1],{if min(X)>0 then X[end] else 1});
      //  end if;

  /*  print("cp_CO2: "+String(cp_vec[1])+" J/kg");
  print("cp_N2: "+String(cp_vec[2])+" J/kg");
  print("cp_CH4: "+String(cp_vec[3])+" J/kg");
  print("cp_H2O: "+String(cp_vec[4])+" J/kg"); */

  end specificHeatCapacityCp_pTX;

  redeclare function extends dynamicViscosity
    "water-saturated  thermal conductivity of water"
  //very little influence of salinity
  algorithm
    eta := dynamicViscosity_pTX(
          p=state.p,
          T=state.T,
          X= if waterSaturated then
        waterSaturatedComposition_pTX(state.p,state.T,state.X)
    else state.X);
  //  else state.X[end - nX + 1:end]);
  //  assert(lambda>0,"lambda="+String(lambda));
  end dynamicViscosity;

  redeclare function extends dynamicViscosity_pTX
    "calculation of gas dynamic Viscosity"
  /*  import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
  input SI.Pressure p;
  input SI.Temperature T;
  input SI.MassFraction[nX] X "Mass fractions of mixture";
  output SI.DynamicViscosity eta;*/
  algorithm
    eta:=Modelica.Media.Air.MoistAir.dynamicViscosity(
      Modelica.Media.Air.MoistAir.ThermodynamicState(
      p=0,
      T=T,
      X={0,0}));
  end dynamicViscosity_pTX;

  redeclare function extends thermalConductivity
    "water-saturated  thermal conductivity of water"
  //very little influence of salinity
  algorithm
    lambda := thermalConductivity_pTX(
          p=state.p,
          T=state.T,
          X= if waterSaturated then
        waterSaturatedComposition_pTX(state.p,state.T,state.X)
    else state.X);
  //  else state.X[end - nX + 1:end]);
  //  assert(lambda>0,"lambda="+String(lambda));
  if lambda<0 then
    print("lambda = " + String(lambda) + "W/(m·K)");
  end if;

  end thermalConductivity;

  redeclare function extends thermalConductivity_pTX
    "calculation of gas thermal conductivity"
  /*  import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
  input SI.Pressure p;
  input SI.Temperature T;
  input SI.MassFraction[nX] X "Mass fractions of mixture";
  output SI.DynamicViscosity eta;*/
  algorithm
    lambda:=Modelica.Media.Air.MoistAir.thermalConductivity(
      Modelica.Media.Air.MoistAir.ThermodynamicState(
      p=0,
      T=T,
      X={0,0}));
  end thermalConductivity_pTX;

  redeclare function specificEnthalpy_pTX
    "calculation of specific enthalpy of gas mixture"
    import Modelica.Media.IdealGases.Common.SingleGasNasa;
    import Modelica.Media.IdealGases.SingleGases;
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output SpecificEnthalpy h "Specific enthalpy";
  protected
    SI.SpecificEnthalpy h_H2O_sat=Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hv_p(p);
    SI.SpecificEnthalpy h_H2O=max(h_H2O_sat, Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p,T))
      "to make sure it is gaseous";
  //  SI.SpecificEnthalpy h_CO2=Modelica.Media.IdealGases.SingleGases.CO2.h_T(Modelica.Media.IdealGases.SingleGases.CO2.data,T);
    SI.SpecificEnthalpy h_CO2=SingleGasNasa.h_T(data=SingleGases.CO2.data,T=T);
    SI.SpecificEnthalpy h_N2=SingleGasNasa.h_T(data=SingleGases.N2.data,T=T);
    SI.SpecificEnthalpy h_CH4=SingleGasNasa.h_T(data=SingleGases.CH4.data,T=T);
    SI.SpecificEnthalpy[:] h_vec={h_CO2,h_N2,h_CH4,h_H2O};

  algorithm
    if debugmode then
      print("Running specificEnthalpy_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;

    if not min(X)>0 then
      print("No gas composition, assuming water vapour.(BrineProp.BrineGas_3Gas.specificEnthalpy_pTX)");
    end if;

  //  h := h_vec*X "mass weighted average";
    h := h_vec * cat(1,X[1:end-1],{if min(X)>0 then X[end] else 1});

  /*  print("h_CO2: "+String(h_CO2)+" J/kg");
  print("h_N2: "+String(h_N2)+" J/kg");
  print("h_CH4: "+String(h_CH4)+" J/kg");
  print("h_H2O: "+String(h_H2O)+" J/kg");
  */
  end specificEnthalpy_pTX;

  redeclare function extends specificEnthalpy
    "water-saturated specific enthalpy of gas phase"
  algorithm
       h := specificEnthalpy_pTX(
          p=state.p,
          T=state.T,
          X= if waterSaturated then
       waterSaturatedComposition_pTX(state.p,state.T,state.X)
    else state.X);
  //  else state.X[end - nX + 1:end]);

  end specificEnthalpy;
end BrineGas_3Gas;
