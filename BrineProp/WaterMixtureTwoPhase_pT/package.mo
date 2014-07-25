within BrineProp;
package WaterMixtureTwoPhase_pT "(incomplete) Water model from Modelica.Media compatible to PartialMixtureTwoPhaseMedium (Example use)"

constant Integer nX_salt = 0;
constant Integer nX_gas = 0;


 extends PartialMixtureTwoPhaseMedium(
    final mediumName="TwoPhaseMixtureWater",
    final substanceNames={"water"},
    final reducedX=true,
    final singleState=false,
    reference_X=cat(
        1,
        fill(0, nX - 1),
        {1}),
    fluidConstants=BrineConstants);
//   final extraPropertiesNames={"gas enthalpy","liquid enthalpy"},

  constant Modelica.SIunits.MolarMass M_H2O = 0.018015 "[kg/mol] TODO";


 redeclare model extends BaseProperties "Base properties of medium"

   Real GVF=x*d/d_g "gas void fraction";
   Modelica.SIunits.Density d_l = Modelica.Media.Water.IF97_Utilities.rhol_T(T);
   Modelica.SIunits.Density d_g = Modelica.Media.Water.IF97_Utilities.rhov_T(T);
 /*  Modelica.SIunits.Density d_l = Modelica.Media.Water.IF97_Utilities.rhol_p(p);
  Modelica.SIunits.Density d_g = Modelica.Media.Water.IF97_Utilities.rhov_p(p);*/
 /*  Modelica.SIunits.Density d_l = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.rhol_p(p);
  Modelica.SIunits.Density d_g = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.rhov_p(p);*/
   Modelica.SIunits.SpecificEnthalpy h_l = bubbleEnthalpy(sat);
   Modelica.SIunits.SpecificEnthalpy h_g = dewEnthalpy(sat);
  // Real x = Medium.vapourQuality(props.state);
   Real x = min(max((h - h_l)/(h_g - h_l+ 1e-18), 0), 1)
    "(min=0,max=1) gas phase mass fraction";
 //  Integer phase_out "calculated phase";
   //END no gas case
 equation
   u = h - p/d;
   MM = M_H2O;
   R  = Modelica.Constants.R/MM;

 //End GVF

 //DENSITY
 //  q = vapourQuality(state);
     d = Modelica.Media.Water.WaterIF97_base.density_ph(p,h);
 //  d = d_l/(1-q*(1-d_l/d_g));
 //End DENSITY

 //ENTHALPY
   h = specificEnthalpy_pTX(p,T,X);
 /*
      if (p_H2O>p) then
    h_H2O_g = Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p,T,1);
  else
    h_H2O_g = Modelica.Media.Water.WaterIF97_base.dewEnthalpy(Modelica.Media.Water.WaterIF97_base.setSat_p(p));
  end if;
  h_gas_dissolved = 0;
  Delta_h_solution = solutionEnthalpy(T) ""only valid for saturated solution";
*/
 //assert(abs(((1-q)*h_l + q*h_g-h)/h) < 1e-3,"Enthalpie stimmt nicht! h_calc="+String((1-q)*h_l + q*h_g)+"<>h="+String(h));
 //End ENTHALPY

   s=0 "TODO";

   state = ThermodynamicState(
     p=p,
     T=T,
     X=X,
     h=h,
     x=x,
     s=0,
     d=d,
     X_l=X,
     d_g=d_g,
     d_l=d_l,
     phase=0) "phase_out";
 //    GVF=GVF,

 /*  if dT_explicit then
    p = pressure_dT(d, T, phase);
    h = specificEnthalpy_dT(d, T, phase);
    sat.Tsat = T;
    sat.psat = saturationPressure(T);
  elseif ph_explicit then
    d = density_ph(p, h, phase);
    T = temperature_ph(p, h, phase);
    sat.Tsat = saturationTemperature(p);
    sat.psat = p;
  else*/
     sat.psat = p;
     sat.Tsat = saturationTemperature(p);
     sat.X = X;
 //  end if;
   annotation (Documentation(info="<html></html>"),
               Documentation(revisions="<html>

</html>"));
 end BaseProperties;


  redeclare function specificEnthalpy_pTX
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:]=fill(0,0) "mass fraction m_NaCl/m_Sol";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Modelica.SIunits.SpecificEnthalpy h=Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p,T);
  algorithm
  //  Modelica.Utilities.Streams.print("specificEnthalpy_pTX("+String(p)+","+String(T)+")");
    annotation(LateInline=true,inverse(T = temperature_phX(p=p,h=h,X=X,phase=phase)));
  end specificEnthalpy_pTX;


  redeclare function temperature_phX
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.SpecificEnthalpy h;
    input MassFraction X[:]=fill(0,0) "mass fraction m_XCl/m_Sol";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Modelica.SIunits.Temp_K T=Modelica.Media.Water.WaterIF97_base.temperature_ph(p,h);
  algorithm
  //  Modelica.Utilities.Streams.print("temperature_phX("+String(p)+","+String(h)+")");

     annotation(LateInline=true,inverse(h = specificEnthalpy_pTX(p=p,T=T,phase=phase,X=X)));
  end temperature_phX;


redeclare record extends ThermodynamicState
  "a selection of variables that uniquely defines the thermodynamic state"
/*  AbsolutePressure p "Absolute pressure of medium";
  Temperature T(unit="K") "Temperature of medium";
  MassFraction X[nX] "Mass fraction of NaCl in kg/kg";*/
  SpecificEnthalpy h "Specific enthalpy";
  SpecificEntropy s "Specific entropy";
  Density d(start=300) "density";
//  Real GVF "Gas Void Fraction";
//  Density d_l(start=300) "density liquid phase";
//  Density d_g(start=300) "density gas phase";
  Real x "vapor quality on a mass basis [mass vapor/total mass]";

   annotation (Documentation(info="<html>

</html>"));
end ThermodynamicState;


  redeclare function extends dewEnthalpy "dew curve specific enthalpy of water"
  algorithm
    hv :=  Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hv_p(sat.psat);
  end dewEnthalpy;


  redeclare function extends bubbleEnthalpy
  "boiling curve specific enthalpy of water"
  algorithm
    hl := Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hl_p(sat.psat);
  end bubbleEnthalpy;


  redeclare function extends saturationTemperature
  algorithm
     //T := Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
     T := Modelica.Media.Water.WaterIF97_base.saturationTemperature(p);
  end saturationTemperature;


  redeclare function extends dynamicViscosity
  algorithm
    eta := Modelica.Media.Water.WaterIF97_base.dynamicViscosity(state);
  end dynamicViscosity;


redeclare function extends specificEntropy "specific entropy of water"
algorithm
    s := Modelica.Media.Water.IF97_Utilities.s_ph(state.p, state.h, state.phase);
end specificEntropy;


redeclare function specificEnthalpy_ps
  "Computes specific enthalpy as a function of pressure and temperature"
    extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEnthalpy h "specific enthalpy";
algorithm
  h := Modelica.Media.Water.IF97_Utilities.h_ps(p, s, phase);
end specificEnthalpy_ps;


redeclare function extends setState_psX
  "Return thermodynamic state of water as function of p and s"
algorithm
  state := ThermodynamicState(
    d=density_ps(p,s),
    T=temperature_ps(p,s),
    phase=0,
    h=specificEnthalpy_ps(p,s),
    p=p,
    X=X,
    s=s,
    q=-1,
    GVF=-1,
    d_l=-1,
    d_g=-1);
end setState_psX;


  redeclare function extends temperature "return temperature of ideal gas"
  algorithm
    T := state.T;
  end temperature;


redeclare function extends setState_pTX
  "Return thermodynamic state of water as function of p and T"

protected
    constant SpecificEnthalpy eps = 1e-8;
  SpecificEnthalpy h=specificEnthalpy_pTX(p,T);
  SpecificEnthalpy hl=bubbleEnthalpy(setSat_pX(p,X));
  Real x =  min(max((h - hl)
  /(dewEnthalpy(setSat_pX(p,X)) - hl
   + eps), 0), 1);

algorithm
  state := ThermodynamicState(
    d=density_pTX(p,T),
    T=T,
    phase=0,
    h=h,
    p=p,
    X=X,
    X_l=X,
    s=specificEntropy_pTX(p,T),
    x=x,
    d_l=  Modelica.Media.Water.IF97_Utilities.rhol_T(T),
    d_g=  Modelica.Media.Water.IF97_Utilities.rhov_T(T));
//    GVF=-1,
end setState_pTX;


redeclare function specificEntropy_pTX
  "Computes specific entropy as a function of pressure and temperature"
    extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Temperature T "Specific entropy";
  input MassFraction X[:]=fill(0,0) "mass fraction m_XCl/m_Sol";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output SpecificEntropy s "specific enthalpy";
algorithm
  s := Modelica.Media.Water.IF97_Utilities.s_pT(p, T, phase);
end specificEntropy_pTX;


  redeclare function extends thermalConductivity
  "Thermal conductivity of water"
  algorithm
    lambda := Modelica.Media.Water.IF97_Utilities.thermalConductivity(
        state.d,
        state.T,
        state.p,
        state.phase);
  end thermalConductivity;


  redeclare function extends specificHeatCapacityCp
  "specific heat capacity at constant pressure of water"

  algorithm
    if Modelica.Media.Water.WaterIF97_base.dT_explicit then
      cp := Modelica.Media.Water.IF97_Utilities.cp_dT(
          state.d,
          state.T,
          state.phase);
    elseif Modelica.Media.Water.WaterIF97_base.pT_explicit then
      cp := Modelica.Media.Water.IF97_Utilities.cp_pT(state.p, state.T);
    else
      cp := Modelica.Media.Water.IF97_Utilities.cp_ph(
          state.p,
          state.h,
          state.phase);
    end if;
      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
  end specificHeatCapacityCp;


  redeclare function extends saturationPressure
  algorithm
     p := Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
  end saturationPressure;


  redeclare function extends specificHeatCapacityCv
  "specific heat capacity at constant pressure of water"

  algorithm
    if Modelica.Media.Water.WaterIF97_base.dT_explicit then
      cv := Modelica.Media.Water.IF97_Utilities.cv_dT(
          state.d,
          state.T,
          state.phase);
    elseif Modelica.Media.Water.WaterIF97_base.pT_explicit then
      cv := Modelica.Media.Water.IF97_Utilities.cv_pT(state.p, state.T);
    else
      cv := Modelica.Media.Water.IF97_Utilities.cv_ph(
          state.p,
          state.h,
          state.phase);
    end if;
      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
  end specificHeatCapacityCv;


 redeclare function extends dynamicViscosity_liq
 algorithm
 //  eta := Modelica.Media.Water.WaterIF97_base.dynamicViscosity(state);
   eta := Modelica.Media.Water.IF97_Utilities.dynamicViscosity(state.d, saturationTemperature(state.p), state.p+1, 1);
 end dynamicViscosity_liq;


 redeclare function extends dynamicViscosity_gas
 algorithm
 //  eta := Modelica.Media.Water.WaterIF97_base.dynamicViscosity(state);
   eta := Modelica.Media.Water.IF97_Utilities.dynamicViscosity(state.d, saturationTemperature(state.p), state.p-1, 1);
 end dynamicViscosity_gas;


  redeclare function extends surfaceTension
  "Surface tension in two phase region of water"
  algorithm
  sigma := Modelica.Media.Water.IF97_Utilities.surfaceTension(sat.Tsat);
  end surfaceTension;


redeclare replaceable partial function extends setState_phX
  "Calculates medium properties from p,h,X"
//      input String fluidnames;
protected
  constant SpecificEnthalpy eps = 1e-8;
  SpecificEnthalpy hl=bubbleEnthalpy(setSat_pX(p,X));
  Real x =  min(max((h - hl)
  /(dewEnthalpy(setSat_pX(p,X)) - hl
   + eps), 0), 1);
algorithm

  if BrineProp.debugmode then
    Modelica.Utilities.Streams.print("Running setState_phX(" + String(p/1e5) + " bar,"
       + String(h) + " J/kg,X)...");
  end if;
/*  state := setState_pTX(p,temperature_phX(p,h,X,phase),
    X,
    phase) ",fluidnames)";*/
  state := ThermodynamicState(
    d=density_phX(p,h),
    T=temperature_phX(p,h),
    phase=0,
    h=h,
    p=p,
    X=X,
    X_l=X,
    s=0,
    x=x,
    d_l=  Modelica.Media.Water.IF97_Utilities.rhol_p(p),
    d_g=  Modelica.Media.Water.IF97_Utilities.rhov_p(p));
//    s=specificEntropy_pTX(p,h)

end setState_phX;


  redeclare function density_pTX
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:]=fill(0,0) "mass fraction m_NaCl/m_Sol";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Modelica.SIunits.Density d=Modelica.Media.Water.WaterIF97_base.density_pT(p,T);
  algorithm
  //  Modelica.Utilities.Streams.print("density_pTX("+String(p)+","+String(T)+")");
  //  annotation(LateInline=true,inverse(T = temperature_phX(p=p,h=h,X=X,phase=phase)));
  end density_pTX;


  redeclare function density_phX
    input Modelica.SIunits.Pressure p;
    input SpecificEnthalpy h;
    input MassFraction X[:]=fill(0,0) "mass fraction m_NaCl/m_Sol";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    output Modelica.SIunits.Density d=Modelica.Media.Water.WaterIF97_base.density_ph(p,h);
  algorithm
  //  Modelica.Utilities.Streams.print("density_phX("+String(p)+","+String(h)+")");
  //  annotation(LateInline=true,inverse(T = temperature_phX(p=p,h=h,X=X,phase=phase)));
  end density_phX;


  redeclare function vapourQuality "Return vapour quality"
    input ThermodynamicState state "Thermodynamic state record";
    output MassFraction x= state.x "Vapour quality";
  algorithm
  //  x := state.x;
    annotation(Documentation(info="<html></html>"));
  end vapourQuality;


 annotation (Documentation(info="<html>
  <h1>Water_MixtureTwoPhase_pT</h1>
  This is a an example use of PartialMixtureTwoPhaseMedium.
  It is a (incomplete) water model using the template PartialMixtureTwoPhaseMedium. It uses the property functions from Modelica.Media.Water.<br/>
  
<h3> Created by</h3>
Henning Francke<br/>
Helmholtz Centre Potsdam<br/>
GFZ German Research Centre for Geosciences<br/>
Telegrafenberg, D-14473 Potsdam<br/>
Germany
<p>
<a href=mailto:info@xrg-simulation.de>francke@gfz-potsdam.de</a>
  </html>
"));
end WaterMixtureTwoPhase_pT;
