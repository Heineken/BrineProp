within BrineProp;
partial package PartialBrine_ngas_Newton "Template package for aqueous solutions of m Salts and n Gases, VLE solved by Newton's method"
  //definition of molar masses
 //constant Integer phase=0;

constant String explicitVars = "ph"
  "set of variables the model is explicit for, may be set to all combinations of ph or pT, setting pT should speed up the model in pT cases";


 replaceable package Salt_data = BrineProp.SaltData;

  import Partial_Units;


 extends BrineProp.PartialGasData;

 constant Modelica.SIunits.MolarMass[:] MM_gas;
 constant Integer[:] nM_gas "number of ions per molecule";

 constant Modelica.SIunits.MolarMass[:] MM_salt;
 constant Integer[:] nM_salt "number of ions per molecule";

 constant Modelica.SIunits.MolarMass[:] MM_vec = cat(1,MM_salt, MM_gas, {M_H2O});
 constant Modelica.SIunits.MolarMass[:] nM_vec = cat(1,nM_salt, nM_gas, {1});

//TWO-PHASE-STUFF
constant String saltNames[:]={""};
constant String gasNames[:]={""};

constant Integer nX_salt = size(saltNames, 1) "Number of salt components"   annotation(Evaluate=true);
constant Integer nX_gas = size(gasNames, 1) "Number of gas components" annotation(Evaluate=true);
//TWO-PHASE-STUFF

constant FluidConstants[nS] BrineConstants(
     each chemicalFormula = "H2O+NaCl+KCl+CaCl2+MgCl2+SrCl2+CO2+N2+CH4",
     each structureFormula="H2O+NaCl+KCl+CaCl2+MgCl2+SrCl2+CO2+N2+CH4",
     each casRegistryNumber="007",
     each iupacName="Geothermal Brine",
     each molarMass=0.1,
     each criticalTemperature = 600,
     each criticalPressure = 300e5,
     each criticalMolarVolume = 1,
     each acentricFactor = 1,
     each triplePointTemperature = 273.15,
     each triplePointPressure = 1e5,
     each meltingPoint = 1,
     each normalBoilingPoint = 1,
     each dipoleMoment = 1);


 extends MediaTwoPhaseMixture.PartialMixtureTwoPhaseMedium(
   final mediumName="TwoPhaseMixtureMedium",
   final substanceNames=cat(1,saltNames,gasNames,{"water"}),
   final reducedX =  true,
   final singleState=false,
   reference_X=cat(1,fill(0,nX-1),{1}),
   fluidConstants = BrineConstants);
//   final extraPropertiesNames={"gas enthalpy","liquid enthalpy"},

/*
,
   AbsolutePressure(
     min=1,
     max=500),
   Temperature(
     min=283.15,
     max=623.15,
     start=323)
   final fixedX = false
     */

//  type Pressure_bar = Real(final quantity="Pressure", final unit="bar") "pressure in bar";


 redeclare model extends BaseProperties "Base properties of medium"
 //  MassFraction[nX] X_l(start=cat(1,fill(0,nXi),{1}))  "= X cat(1,X_salt/(X[end]*m_l),X[nX_salt+1:end])";

   Modelica.SIunits.Density d_l;
   Modelica.SIunits.Density d_g;

   Real GVF "=x*d/d_ggas void fraction";
   Real x "(start=0) (min=0,max=1) gas phase mass fraction";
   Modelica.SIunits.Pressure p_H2O;
   Modelica.SIunits.Pressure[nX_gas] p_gas;
   Modelica.SIunits.Pressure p_degas;
 //  Modelica.SIunits.Pressure p_check=sum(p_gas)+p_H2O;
 //  Real k[nX_gas];
   parameter Real[nX_gas+1] n_g_norm_start = fill(0.5,
                                                     nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
 //  Real[nX_gas+1] n_g_norm;
protected
   MassFraction[nX_salt] X_salt = X[1:nX_salt];
   MassFraction[nX_gas] X_gas = X[nX_salt+1:end-1];
 //  Modelica.SIunits.Temperature T_corr = max(273.16,min(400,T)) "TODO";
 //  Modelica.SIunits.Pressure p_corr = max(1e5,min(455e5,p)) "TODO";
   Real y_vec[:]=massFractionsToMoleFractions(X,MM_vec);
   Integer pp(start=0)=state.phase
    "just to get rid of initialization problem warning";
 equation
    //   assert(nX_gas==2,"Wrong number of gas mass fractions specified (2 needed - CO2,N2)");
 //  assert(max(X)<=1 and min(X)>=0, "X out of range [0...1] = "+PowerPlant.vector2string(X)+" (saturationPressure_H2O())");
   u = h - p/d;
 //  MM = (X_salt*MM_salt + X_gas*MM_gas + X[end]*M_H2O);
   MM = y_vec*MM_vec;
   R  = Modelica.Constants.R/MM;

 //  (h,x,d,d_g,d_l) = specificEnthalpy_pTX(p,T,X) Leider nicht invertierbar;

   if explicitVars=="pT" or explicitVars=="pT" then
     h=state.h;
   else
     h = specificEnthalpy_pTX(p,T,X,phase,n_g_norm_start);
   end if;
 //  (x,d,d_g,d_l,p_H2O,p_gas,X_l,p_degas,k)= quality_pTX(p_corr,T_corr,X);

   state = setState_pTX(p,T,X,phase,n_g_norm_start);
 //  X_l=state.X_l;
   GVF=state.GVF;
   x=state.x;
   s=state.s;
   d_g=state.d_g;
   d_l=state.d_l;
   d=state.d;
   p_H2O=state.p_H2O;
   p_gas=state.p_gas;
   p_degas=sum(state.p_degas);
  /* 
  (x,d,d_g,d_l,p_H2O,p_gas,X_l,p_degas)= quality_pTX(p,T,X,n_g_start);
  s = 0 "specificEntropy_phX(p,h,X) TODO";
 
state =  ThermodynamicState( p=p,
                              T=T,
                              X=X,
                              X_l=X_l,
                              h=h,
                              GVF=x*d/d_g,
                              x=x,
                              s=0,
                              d_g=d_g,
                              d_l=d_l,
                              d=d,
                              phase=0) "phase_out";*/
 /*                              GVF=GVF,
                              x=GVF*d_g/d,*/

   sat.psat = sum(state.p_degas);
   sat.Tsat = T;
   sat.X = X;
  // sat.p_degas=p_degas;

   annotation (Documentation(info="<html></html>"),
               Documentation(revisions="<html>

</html>"));
 end BaseProperties;

 /* Provide implementations of the following optional properties.
     If not available, delete the corresponding function.
     The record "ThermodynamicState" contains the input arguments
     of all the function and is defined together with the used
     type definitions in PartialMedium. The record most often contains two of the
     variables "p, T, d, h" (e.g. medium.T)
  */

// redeclare replaceable record ThermodynamicState


redeclare record extends ThermodynamicState
  "a selection of variables that uniquely defines the thermodynamic state"
/*  AbsolutePressure p "Absolute pressure of medium";
  Temperature T(unit="K") "Temperature of medium";*/
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
  AbsolutePressure p_gas[nX_gas];
  AbsolutePressure[nX_gas + 1] p_degas
    "should be in SatProp, but is calculated in setState which returns a state";
   annotation (Documentation(info="<html>

</html>"));
end ThermodynamicState;



  redeclare function extends dewEnthalpy "dew curve specific enthalpy of water"
  algorithm
    hv := 1000;
  end dewEnthalpy;


  redeclare function extends bubbleEnthalpy
  "boiling curve specific enthalpy of water"
  algorithm
    hl := 2000;
  end bubbleEnthalpy;


  redeclare function extends saturationTemperature "saturation temperature"
  algorithm
    T := 373.15;
  end saturationTemperature;


  replaceable partial function solutionEnthalpy
    input Modelica.SIunits.Temp_K T;
    output Modelica.SIunits.SpecificEnthalpy Delta_h_solution;
  end solutionEnthalpy;


  replaceable partial function solubilities_pTX
  "solubility calculation of gas in m_gas/m_H2O"
    //Stoffdaten auslagern
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Modelica.SIunits.MassFraction X[nX] "mass fractions m_x/m_Sol";
    input Modelica.SIunits.MassFraction X_l[nX] "mass fractions m_x/m_Sol";
    input Modelica.SIunits.Pressure[nX_gas] p_gas;
  //  input String gasname;
  //  input Modelica.SIunits.MolarMass MM[:] "=fill(0,nX)molar masses of components";
  //  output Molality[nX_gas] solu;
    output MassFraction solu[nX_gas] "gas concentration in kg_gas/kg_fluid";
  end solubilities_pTX;


  replaceable function fugacity_pTX
  "Calculation of nitrogen fugacity coefficient extracted from EES"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T_K;
    input MassFraction X[:]=reference_X "Mass fractions";
    input String substancename;
    output Real phi;
protected
    Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
  end fugacity_pTX;


 replaceable function density_liquid_pTX "Dichte der flüssigen Phase"
   input Modelica.SIunits.Pressure p "TODO: Rename to density_liq_pTX";
   input Modelica.SIunits.Temp_K T;
   input MassFraction X[nX] "mass fraction m_NaCl/m_Sol";
   input Modelica.SIunits.MolarMass MM[:]
    "=MM_vec =fill(0,nX) molar masses of components";
   output Modelica.SIunits.Density d;
 end density_liquid_pTX;


redeclare function vapourQuality
  "Returns vapour quality, needs to be defined to overload function defined in PartialMixtureTwoPhaseMedium"
  input ThermodynamicState state "Thermodynamic state record";
  output MassFraction x "Vapour quality";
algorithm
x := state.x;
end vapourQuality;


  redeclare function specificEnthalpy_pTX
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    input Real[nX_gas+1] n_g_norm_start=fill(0.5,
                                                nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
    output Modelica.SIunits.SpecificEnthalpy h;
  /*  output Real x "gas mass fraction";
  output Modelica.SIunits.Density d;
  output Modelica.SIunits.Density d_g;
  output Modelica.SIunits.Density d_l;
  output Modelica.SIunits.Pressure p_H2O;
  output Modelica.SIunits.Pressure[nX_gas] p_gas;*/
  //  Real M_CO2 = MM_gas[1];
  //  Real M_N2 = MM_gas[2];
  /*  Real x "gas mass fraction";
  Modelica.SIunits.Density d;
  Modelica.SIunits.Density d_g;
  Modelica.SIunits.Density d_l;
  Modelica.SIunits.Pressure p_H2O;
  Modelica.SIunits.Pressure[nX_gas] p_gas;

  Modelica.SIunits.MassFraction c_gas[nX_gas] = X[end-nX_gas:end-1] 
    "mass ratio m_CO2/m_ges";
  Modelica.SIunits.MassFraction c_H2O = X[end] "mass ratio m_CO2/m_ges";
  Modelica.SIunits.Density d_H2O_g = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.rhov_p(p) 
    "density of water vapor";
*/
  //  Real solu[nX_gas]=solubilities_pTX(p=p,T=T,X=X,p_gas) "gas solubilities";
  /*  Modelica.SIunits.Pressure[nX_gas+1] p_sat=saturationPressures(p,T,X,MM_vec) "= c_CO2/k_H_CO2";
  Modelica.SIunits.Pressure p_H2O_sat= p_sat[end];*/

  //  Real M_CO2 = Modelica.Media.IdealGases.SingleGases.CO2.data.MM;
  //  Modelica.SIunits.Density d_gas_g = p/(Modelica.Constants.R/M_CO2*T) "pure gas density in gas phase (ideal gas)";
  //  Real x "gas mass fraction";
  //  Integer case;
  algorithm
  //  assert(T>273.15,"T too low in PartialBrine_ngas_Newton.specificEnthalpy_pTX()");
    if debugmode then
  //     Modelica.Utilities.Streams.print("Running specificEnthalpy_pTX("+String(p)+","+String(T)+" K)");
        Modelica.Utilities.Streams.print("Running specificEnthalpy_pTX("+String(p/1e5)+","+String(T-273.15)+"°C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;
  /* x:=quality_pTX(p,T,X,n_g_start);
 h := x*specificEnthalpy_gas_pTX(p,T,X) + (1-x)*specificEnthalpy_liq_pTX(p,T,X);
*/
   h:=specificEnthalpy(setState_pTX(
      p,
      T,
      X,
      phase,n_g_norm_start));

  //Modelica.Utilities.Streams.print(String(p)+","+String(T)+" K->"+String(h)+" J/kg & (PartialBrine_Multi_TwoPhase_ngas.specificEnthalpy_pTX)");
   //,p=pressure_ThX(T,h,X);

   annotation(LateInline=true,inverse(T=temperature_phX(p,h,X,phase,n_g_norm_start)));
  end specificEnthalpy_pTX;


  redeclare function temperature_phX
  "iterative inversion of specificEnthalpy_pTX by regula falsi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[nX] "Mass fractions";
    input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
    input Real[nX_gas + 1] n_g_start=fill(0.5,
                                             nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
    output Temperature T "Temperature";
protected
    Modelica.SIunits.SpecificHeatCapacity c_p;
    Modelica.SIunits.Temperature T_a=273.16;
  //  Modelica.SIunits.Temperature T0_a=273.16;
    Modelica.SIunits.Temperature T_b=400;
  //  Modelica.SIunits.Temperature T0_b=400 "limit of N2 solubility";
  //  Modelica.SIunits.Temperature T_neu;
    Modelica.SIunits.SpecificEnthalpy h_a;
    Modelica.SIunits.SpecificEnthalpy h_b;/**/
    Modelica.SIunits.SpecificEnthalpy h_T;
    Integer z=0 "Loop counter";
  algorithm
    if debugmode then
       Modelica.Utilities.Streams.print("\ntemperature_phX("+String(p)+","+String(h)+")");
    end if;
    //Find temperature with h above given h ->T_b
    assert(h>specificEnthalpy_pTX(p,T_a,X),"h="+String(h/1e3)+" kJ/kg -> Enthalpy too low (< 0°C) (Brine.PartialBrine_ngas_Newton.temperature_phX)");
    while true loop
      h_T:=specificEnthalpy_pTX(p,T_b,X);
  //    Modelica.Utilities.Streams.print(String(p)+","+String(T_b)+" K->"+String(h_T)+" J/kg (PartialBrine_ngas_Newton.temperature_phX)");
      if h>h_T then
        T_a := T_b;
        T_b := T_b + 50;
      else
        break;
      end if;
    end while;

  //BISECTION - is schneller, braucht 13 Iterationen
    while (T_b-T_a)>1e-2 and abs(h-h_T/h)>1e-5 loop   //stop when temperatures or enthalpy are close
  //  while abs(h-h_T/h)>1e-5 loop
  //    Modelica.Utilities.Streams.print("T_b-T_a="+String(T_b-T_a)+", abs(h-h_T)/h="+String(abs(h-h_T)/h));
      T:=(T_a+T_b)/2 "Halbieren";
  //    Modelica.Utilities.Streams.print("T_neu="+String(T)+"K");
      h_T:=specificEnthalpy_pTX(p,T,X);
      if h_T > h then
        T_b:=T;
  //      Modelica.Utilities.Streams.print("T_b="+String(T)+"K -> dh="+String(h_T-h));
      else
        T_a:=T;
  //      Modelica.Utilities.Streams.print("T_a="+String(T)+"K -> dh="+String(h_T-h));
      end if;
      z:=z+1;
  //    Modelica.Utilities.Streams.print(String(z)+": "+String(T_a)+" K & "+String(T_b)+" K -> "+String((h-h_T)/h)+"(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)\n");
  //    Modelica.Utilities.Streams.print("h("+String(T_a)+")="+String(h_a-h)+" J/kg & h("+String(T_b)+")="+String(h_b-h)+" J/kg");
      assert(z<100,"Maximum number of iteration reached for temperature calculation. Something's wrong here. Cancelling...(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)");
    end while;
  // Modelica.Utilities.Streams.print("BISECTION " + String(z)+": "+String(T));

  /*
//REGULA FALSI - is langsamer, braucht 19 Iterationen
  z:=0;
  T_a:=T0_a;
  T_b:=T0_b "limit of N2 solubility";
  h_a := specificEnthalpy_pTX(p,T_a,X);
  h_b := specificEnthalpy_pTX(p,T_b,X);
  while abs(T_b-T_a)>1e-2 and abs(h_T-h)/h>1e-5 loop
//  while abs(T_b-T_a)/T_l>1e-4 loop
    Modelica.Utilities.Streams.print("h_a("+String(T_a)+")="+String(h_a)+" / h_b("+String(T_b)+")="+String(h_b));
    T:=max(T0_a,min(T0_b,T_a-(T_b-T_a)/(h_b-h_a)*(h_a-h))) "Regula falsi";
    h_T:=specificEnthalpy_pTX(p,T,X);
    Modelica.Utilities.Streams.print("T_neu="+String(T)+"K");
    if h_T > h then
      T_b:=T;
      h_b:=h_T;
    else
      T_a:=T;
      h_a:=h_T;
//      Modelica.Utilities.Streams.print("T_a="+String(T)+"K -> h="+String(h_T-h));
    end if;
    z:=z+1;
//    Modelica.Utilities.Streams.print(String(z)+": "+String(T_a)+" K & "+String(T_b)+" K -> "+String((h-h_T)/h)+"(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)\n");
//    Modelica.Utilities.Streams.print("h("+String(T_a)+")="+String(h_a-h)+" J/kg & h("+String(T_b)+")="+String(h_b-h)+" J/kg");
    assert(z<100,"Maximum number of iteration reached for temperature calculation. Something's wrong here. Cancelling...(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)");
  end while;
 Modelica.Utilities.Streams.print("REGULA FALSI " + String(z)+": "+String(T));
*/

  end temperature_phX;


 replaceable function specificEnthalpy_liq_pTX
  "Specific enthalpy of liquid phase"
   input Modelica.SIunits.Pressure p;
   input Modelica.SIunits.Temp_K T;
   input MassFraction X[nX] "mass fraction m_NaCl/m_Sol";
   input Modelica.SIunits.MolarMass MM[:]=fill(0,nX)
    "molar masses of components";
   output Modelica.SIunits.SpecificEnthalpy h;
protected
   SI.SpecificEnthalpy[nX_salt] h_vec;
 end specificEnthalpy_liq_pTX;


 replaceable function specificEnthalpy_gas_pTX
  "Specific enthalpy of gas in gas phase"
   input Modelica.SIunits.Pressure p;
   input Modelica.SIunits.Temp_K T;
 //  input Modelica.SIunits.MolarMass MM[:]=fill(0,nX)     "molar masses of components";
   input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
   output Modelica.SIunits.SpecificEnthalpy h;
 end specificEnthalpy_gas_pTX;


  replaceable partial function saturationPressures
  "Return saturationPressures for gases and water"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
    input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
    output Modelica.SIunits.Pressure[nX_gas] p_sat;
  end saturationPressures;



  redeclare replaceable partial function extends setState_pTX
  "finds the VLE iteratively by varying the normalized quantity of gas in the gasphase, calculates the densities"
  input Real[nX_gas + 1] n_g_norm_start "=fill(.1,nX_gas+1) 
    start value, all gas in gas phase, all water liquid, set in BaseProps";
  /*
//output Modelica.SIunits.Density d_g= if x>0 then (n_CO2_g*d_g_CO2 + n_N2_g*d_g_N2)/(n_CO2_g + n_H2O_g) else -1;
//output Real[nX_gas + 1] n_g_norm;
//output Real k[nX_gas];
// Modelica.SIunits.Density d_g_H2O = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.rhov_p(p) "density of water vapor";
*/
protected
    Modelica.SIunits.Density d;
    Modelica.SIunits.Density d_g;
    Modelica.SIunits.Density d_l;
    Modelica.SIunits.SpecificEnthalpy h_g;
    Modelica.SIunits.SpecificEnthalpy h_l;
    Modelica.SIunits.MassFraction[nX] X_l=X "start value";
    Modelica.SIunits.MassFraction[nX_gas+1] X_g;
    Modelica.SIunits.Pressure p_H2O;
    Modelica.SIunits.MassFraction x;
    Modelica.SIunits.Pressure[nX_gas + 1] p_degas;
    Modelica.SIunits.Pressure p_sat_H2O
    "= saturationPressure_H2O(p,T2,X,MM_vec,nM_vec)";
    Modelica.SIunits.Pressure p_H2O_0;
    Modelica.SIunits.Pressure[nX_gas + 1] f;
    Modelica.SIunits.Pressure[nX_gas + 1] p_sat;
    Modelica.SIunits.Pressure[nX_gas + 1] p_sat_test;
    Modelica.SIunits.Pressure[nX_gas + 1] p_gas "=fill(0,nX_gas)";
    Modelica.SIunits.MassFraction[nX_gas + 1] Delta_n_g_norm = fill(1e3,nX_gas+1);
  //  Modelica.SIunits.MassFraction[nX_gas + 1] c = {3.16407e-5,0,3.6e-8,.746547} "cat(1,fill(1e-4, nX_gas), {X[end]})fill(0, nX_gas+1)X[nX_salt+1:end]";
    Real k_H2O "Henry coefficient";
    Real k[nX_gas];
    Real[nX_gas + 1] n "Total mol numbers";
    Real[nX_gas + 1] n_l "mols in liquid phase per kg fluid";
    Real[nX_gas + 1] n_g "mols in   gas  phase per kg fluid";
    Real[nX_gas + 1] n_g_norm_test;
  //  Modelica.SIunits.MassFraction[nX] X;
    Real[nX_gas + 1] n_g_norm
    "= X[end-nX_gas:end-1]./MM_gas fill(0,nX_gas) - start value: all degassed";
    Real dp_gas_dng_norm;
    Real dcdng_norm;
    Real dp_degas_dng_norm;
    Real[nX_gas + 1] dfdn_g_norm;
    Integer z=0;
    Real sum_n_ion;
    constant Integer zmax=1000 "maximum number of iteration";
  //  Integer ju = nX_gas+1;
    Real[nX_gas + 1,nX_gas + 1] Grad_f;
    Real DeltaC=0.001;
    Modelica.SIunits.Temperature T2;
    SpecificHeatCapacity R_gas;
  algorithm
    if debugmode then
  //    Modelica.Utilities.Streams.print("Running setState_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+"°C,X)...");
        Modelica.Utilities.Streams.print("Running setState_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;
  //  assert(T-273.15<300,"Moooot!");

   assert(p>0,"p="+String(p/1e5)+"bar - Negative pressure is not yet supported ;-) (PartialBrine_ngas_Newton.quality_pTX())");
  /*  Modelica.Utilities.Streams.print("quality_pTX("+String(p)+","+String(T2)+","+PowerPlant.vector2string(X_l[1:end],false)+")");
  X[1:nX_salt] = X_[1:nX_salt];
  for i in nX_salt+1:nX-1 loop
    X[i]:=max(0,min(1e-3,X_[i]));
  end for;
  X[end]=1-sum(X[1:end-1]);
  X_l:=X;*/
  //  Modelica.Utilities.Streams.print("quality_pTX("+String(p)+","+String(T)+","+PowerPlant.vector2string(X[1:end],false)+")");

    assert(max(X)-1<=1e-8 and min(X)>=-1e-8, "X out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X]))+" (quality_pTX())");
  //  assert(T>273.15,"T too low in PartialBrine_ngas_Newton.()");
  //    Modelica.Utilities.Streams.print("\nn_g_start=" + PowerPlant.vector2string(n_g_start));

      if T<273.15 then
      Modelica.Utilities.Streams.print("T="+String(T)+" too low (<0°C), setting to 0°C in PartialBrine_ngas_Newton.quality_pTX()");
    end if;
    T2:= max(273.16,T);

    p_H2O := saturationPressure_H2O(p,T2,X,MM_vec,nM_vec);

  /*
    p_sat_H2O := saturationPressure_H2O(p,T2,X,MM_vec,nM_vec);
   p_degas := if phase==1 then 0 else sum(saturationPressures(p,T2,X,MM_vec)) + p_sat_H2O;
   if  p_degas < p then*/
     p_degas := cat(1,saturationPressures(p,T2,X,MM_vec), {p_H2O});

     if phase==1 or sum(p_degas) < p then
     if debugmode then
      Modelica.Utilities.Streams.print("1Phase-Liquid (PartialBrine_Multi_TwoPhase_ngas.quality_pTX("+String(p)+","+String(T2)+"))");
     end if;
      x:=0;
      p_H2O := p_sat_H2O;
    else
      assert(max(X[end-nX_gas:end-1])>0,"Phase equilibrium cannot be calculated without dissolved gas at "+String(p/1e5)+" bar, "+String(T2-273.15)+"°C with p_degas="+String(sum(p_degas)/1e5)+" bar.");
  //    Modelica.Utilities.Streams.print("2Phase (PartialBrine_Multi_TwoPhase_ngas.quality_pTX)");
  //    Modelica.Utilities.Streams.print("p="+String(p/1e5)+" bar");

      n:=X[nX_salt + 1:end] ./ MM_vec[nX_salt + 1:nX]
      "total mole numbers per kg brine";
  //    n_g_norm:=cat(1, fill(1,nX_gas), {0.01})
  //    n_g:={0.0175164,0.00204528,0.00435, 2.51075};
  //    n_g_norm :=n_g ./ n;
      n_g_norm:=n_g_norm_start .* sign(X[nX_salt + 1:nX])
      "switch off unused salts";

      while z<1 or max(abs(Delta_n_g_norm))>5e-5 loop
      //abbrechen wenn Druck-GG gefunden oder sehr geringer Gasanteil
        z:=z + 1;
        assert(z<=zmax,"Reached maximum number of iterations ("+String(z)+"/"+String(zmax)+") for solution equilibrium calculation. (quality_pTX("+String(p/1e5)+"bar,"+String(T2-273.16)+"°C))\nDeltaP="+String(max(abs(p_sat-p_gas))));

  //     Modelica.Utilities.Streams.print("\nn_g_norm=" + PowerPlant.vector2string(n_g_norm));
        n_g :=n_g_norm .* n;
  /*      if abs(Delta_n_g_norm[ju])<1e-3 then
         ju:=if ju == nX_gas+1 then 1 else ju + 1;
         Modelica.Utilities.Streams.print("Gas "+String(ju)+"!");
      end if;
*/
        n_l := n-n_g;
        x := n_g*MM_vec[nX_salt+1:nX];
  //      Modelica.Utilities.Streams.print("\n"+String(z)+": x="+String(x)+" Delta_n_g="+String(max(abs(Delta_n_g_norm))));
        X_l:=cat(1, X[1:nX_salt], n_l.*MM_vec[nX_salt+1:nX])/(1-x);
   /*      Modelica.Utilities.Streams.print("n_l=" + PowerPlant.vector2string(n_l));
      Modelica.Utilities.Streams.print("n_g=" + PowerPlant.vector2string(n_g));
*/
    //PARTIAL PRESSURE
          p_gas := p * n_g/sum(n_g);

    //Concentration for that partial pressure

    //DEGASSING PRESSURE
          (p_H2O,p_H2O_0):=saturationPressure_H2O(p,T2,X_l,MM_vec,nM_vec)
        "X_l ändert sich";
      if (p_H2O>p) then
          Modelica.Utilities.Streams.print("p_H2O(" + String(p/1e5) + "bar," +
            String(T2 - 273.15) + "°C, " + Modelica.Math.Matrices.toString(transpose([X])) + ") = "
             + String(p_H2O/1e5) + "bar>p ! (PartialBrine_ngas_Newton.quality_pTX)");
        x:=1;
        break;
      end if;

  //       Modelica.Utilities.Streams.print("p_H2O_0=" + String(p_H2O_0));

  //        k:=solubilities_pTX(p=p, T=T2, X_l=X_l, X=X, p_gas=fill(p/3,3)) ./ fill(p/3,3);
  //  Modelica.Utilities.Streams.print("X_l="+PowerPlant.vector2string(X_l[nX_salt+1:end]));
          k:=solubilities_pTX(p=p, T=T2, X_l=X_l, X=X, p_gas=p_gas[1:nX_gas]) ./ p_gas[1:nX_gas];
  //    Modelica.Utilities.Streams.print("k="+PowerPlant.vector2string(k)+" (PartialBrine_ngas_Newton.quality_pTX)");p,T,X_l,MM_vec,p_gas[1])
          for i in 1:nX_gas loop
            p_sat[i] := X_l[nX_salt+i]/ (if k[i]>0 then k[i] else 1e10)
          "Entlösedruck";
          end for;
          p_sat[nX_gas+1] := p_H2O;

          f :=  p_gas-p_sat;
  //       Modelica.Utilities.Streams.print("p_gas=" + PowerPlant.vector2string(p_gas) + "=>" + String(sum(p_gas)));
  //       Modelica.Utilities.Streams.print("p_sat=" + PowerPlant.vector2string(p_sat));

         sum_n_ion :=cat(1,X[1:nX_salt] ./ MM_vec[1:nX_salt],n_l)*nM_vec;

    //GRADIENT analytisch df[gamma]/dc[gamma]

          for gamma in 1:nX_gas+1 loop
              dp_gas_dng_norm:=p*n[gamma]* (sum(n_g)-n_g[gamma])/(sum(n_g))^2
          "partial pressure";
              if gamma == nX_gas+1 then
                dp_degas_dng_norm := p_H2O_0*n[end]*( (if gamma == nX_gas+1 then -sum_n_ion else 0)+(1-n_g_norm[end])*n[gamma]) / sum_n_ion^2;
              else
                  dcdng_norm := n[gamma]*MM_vec[nX_salt+gamma]*( (x-1) +(1 - n_g_norm[gamma])*n[gamma]*MM_vec[nX_salt+gamma])/(1 - x)^2;
                  dp_degas_dng_norm := dcdng_norm / (if k[gamma] > 0 then k[gamma] else 1e-10)
            "degassing pressure";
              end if;
              dfdn_g_norm[gamma] := dp_gas_dng_norm-dp_degas_dng_norm;
          end for;

  /*        
  //GRADIENT analytisch df[alpha]/dc[gamma]
       for gamma in 1:nX_gas+1 loop
          for alpha in 1:nX_gas+1 loop
            dp_gas_dng_norm:=p*n[gamma]*((if alpha == gamma then sum(n_g) else 0)-n_g[alpha])/(sum(n_g))^2 
            "partial pressure";
            if alpha == nX_gas+1 then
              dp_degas_dng_norm := p_H2O_0*n[end]*( (if gamma == nX_gas+1 then -sum_n_ion else 0)+(1-n_g_norm[end])*n[gamma]) / sum_n_ion^2;
            else
               if alpha == gamma then
                dcdng_norm := n[alpha]*MM_vec[nX_salt+alpha]*( (x-1) +(1 - n_g_norm[alpha])*n[gamma]*MM_vec[nX_salt+gamma])/(1 - x)^2;
                dp_degas_dng_norm := dcdng_norm /k[alpha] "degassing pressure";
              else
                dp_degas_dng_norm := 0 "degassing pressure";
              end if;
//            Modelica.Utilities.Streams.print("dcdng_norm("+String(alpha)+","+String(gamma)+")=" + String(dcdng_norm));
            end if;
            Grad_f[gamma,alpha] := dp_gas_dng_norm-dp_degas_dng_norm;

/*           Modelica.Utilities.Streams.print("dp_gas_dng_norm("+String(gamma)+","+String(alpha)+")=" + String(dp_gas_dng_norm));
           Modelica.Utilities.Streams.print("dp_degas_dng_norm("+String(gamma)+","+String(alpha)+")=" + String(dp_degas_dng_norm));
           * /

          end for;
//         Modelica.Utilities.Streams.print("Grad_f["+String(gamma)+",:] =" + PowerPlant.vector2string(Grad_f[gamma,:]));
        end for;
*/

  //       Modelica.Utilities.Streams.print("k=" + PowerPlant.vector2string(k));/**/
  //       Modelica.Utilities.Streams.print("dp_gas=" + PowerPlant.vector2string(p_sat - p_gas));

    //SOLVE NEWTON STEP
  //        Delta_n_g_norm := Modelica.Math.Matrices.solve(Grad_f, -f)         "solve Grad_f*Delta_n_g_norm=-f";
  //        n_g_norm := n_g_norm + Delta_n_g_norm;

  //        Modelica.Utilities.Streams.print("n_g_norm="+Modelica.Math.Matrices.toString({n_g_norm}));
          for alpha in 1 :nX_gas+1 loop
  //        for alpha in ju:ju loop
  //          Delta_n_g_norm[alpha] := -f[alpha]/Grad_f[alpha,alpha];
            Delta_n_g_norm[alpha] := if X[nX_salt+alpha]>0 then -f[alpha]/dfdn_g_norm[alpha] else 0;
  //          if alpha==ju then
  //            n_g_norm[alpha] := max(0,min(1,n_g_norm[alpha] + b[alpha]*Delta_n_g_norm[alpha]))
              n_g_norm[alpha] := max(1e-9,min(1,n_g_norm[alpha] + Delta_n_g_norm[alpha]))
          "new concentration limited by all dissolved/none dissolved, 1e-9 to avoid k=NaN";
  //          end if;
          end for;
  //       Modelica.Utilities.Streams.print("p_sat="+String(p_sat[1])+", solu="+String(solubility_CO2_pTX_Duan2006(p,T2,X_l,MM_vec,p_gas[1]))+", p_gas="+String(p_gas[1]));
  //         Modelica.Utilities.Streams.print("p="+String(p)+",T2="+String(T2)+",p_gas[1]="+String(p_gas[1]));
  /*        Modelica.Utilities.Streams.print("X_l="+Modelica.Math.Matrices.toString({X_l}));
        Modelica.Utilities.Streams.print("MM_vec="+Modelica.Math.Matrices.toString({MM_vec}));
*/
      end while;

    end if "p_degas< p";

  //DENSITY
   X_g:=if x>0 then (X[end-nX_gas:end]-X_l[end-nX_gas:end]*(1-x))/x else fill(0,nX_gas+1);
  /*Calculation here  R_gas :=if x > 0 then sum(Modelica.Constants.R*X_g ./ cat(1,MM_gas,{M_H2O})) else -1;
  d_g :=if x > 0 then p/(T2*R_gas) else -1;*/
  //  d_g:= if x>0 then p/(Modelica.Constants.R*T2)*(n_g*cat(1,MM_gas,{M_H2O}))/sum(n_g) else -1;
    d_g :=if x > 0 then BrineGas_3Gas.density_pTX(p,T, X_g[end - nX_gas:end]) else -1
    "calculation in MoistAirModel";

    d_l:=if not x<1 then -1 else density_liquid_pTX(p,T2,X_l,MM_vec)
    "gases are ignored anyway";
    d:=1/(x/d_g + (1 - x)/d_l);
  //  Modelica.Utilities.Streams.print(String(z)+" (p="+String(p_gas[1])+" bar)");

  // X_g:=if x>0 then (X-X_l*(1-x))/x else fill(0,nX);
   h_l:=specificEnthalpy_liq_pTX(p,T,X_l);
   h_g:=specificEnthalpy_gas_pTX(p,T,X_g);
   state :=ThermodynamicState(
      p=p,
      T=T,
      X=X,
      X_l=X_l,
      h_g=h_g,
      h_l=h_l,
      h=x*h_g + (1-x)*h_l,
      GVF=x*d/d_g,
      x=x,
      s=0,
      d_g=d_g,
      d_l=d_l,
      d=d,
      phase=if x>0 and x<1 then 2 else 1,
      p_H2O=p_H2O,
      p_gas=p_gas[1:nX_gas],
      p_degas=p_degas) "phase_out";
    annotation (Diagram(graphics={Text(
            extent={{-96,16},{98,-16}},
            lineColor={0,0,255},
            textStyle={TextStyle.Bold},
            textString="find static VLE")}));
  end setState_pTX;


redeclare replaceable partial function extends setState_phX
  "Calculates medium properties from p,h,X"
//      input String fluidnames;
algorithm

  if debugmode then
    Modelica.Utilities.Streams.print("Running setState_phX("+String(p/1e5)+" bar,"+String(h)+" J/kg,X)...");
  end if;
  state := setState_pTX(p,
    temperature_phX(p,h,X,phase),
    X,
    phase) ",fluidnames)";
end setState_phX;


  replaceable function dynamicViscosity_pTX_unused "viscosity calculation"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
    output Modelica.SIunits.DynamicViscosity eta;
  //  constant Real M_NaCl=0.058443 "molar mass in [kg/mol]";
  end dynamicViscosity_pTX_unused;


  redeclare replaceable function extends specificHeatCapacityCp
  "numeric calculation of specific heat capacity at constant pressure"
protected
    Modelica.SIunits.SpecificHeatCapacity cp_liq=specificHeatCapacityCp_liq(state);
    Modelica.SIunits.SpecificHeatCapacity cp_gas=specificHeatCapacityCp_gas(state);
  algorithm
    cp:=state.x*cp_gas + (1-state.x)*cp_liq;

  //  assert(cp>0 and cp<5000,"T="+String(state.T-273.15)+"K, p="+String(state.p/1e5)+"bar, x="+String(state.x)+", cp_liq="+String(cp_liq)+"J(kgK), cp_gas="+String(cp_gas)+"J(kgK)");

  //  Modelica.Utilities.Streams.print("c_p_liq("+String(state.T)+"°C)="+String(p)+" J/(kg·K)");
      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
  end specificHeatCapacityCp;


  replaceable function specificHeatCapacityCp_liq
  //extends specificHeatCapacityCp;SHOULD WORK WITH THIS!
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificHeatCapacity cp
    "Specific heat capacity at constant pressure";

   /*protected 
  constant Modelica.SIunits.TemperatureDifference dT=.1;
algorithm 
//    cp := Modelica.Media.Water.IF97_Utilities.cp_pT(state.p, state.T) "TODO";
    cp:=(specificEnthalpy_pTX(state.p,state.T+dT,state.X)-state.h)/dT;
    */

  end specificHeatCapacityCp_liq;


  replaceable function specificHeatCapacityCp_gas
  //extends specificHeatCapacityCp;SHOULD WORK WITH THIS!
    extends Modelica.Icons.Function;
    input ThermodynamicState state "thermodynamic state record";
    output SpecificHeatCapacity cp
    "Specific heat capacity at constant pressure";
  /*protected 
  constant Modelica.SIunits.TemperatureDifference dT=.1;
algorithm 
//    cp := Modelica.Media.Water.IF97_Utilities.cp_pT(state.p, state.T) "TODO";
    cp:=(specificEnthalpy_pTX(state.p,state.T+dT,state.X)-state.h)/dT;
    */

protected
    Modelica.SIunits.MassFraction[nX] X_g=if state.x>0 then (state.X-state.X_l*(1-state.x))/state.x
   else
       fill(-1,nX);
    Modelica.SIunits.SpecificHeatCapacity cp_vec[nX_gas+1];
  end specificHeatCapacityCp_gas;


  annotation (Documentation(info="<html>
<p><b>PartialBrine_ngas_Newton</b> is template package for aqueous solution of m Salts and n Gases. The vapour-liquid-equilibrium (VLE) is defined by the water vapour pressure and the gas solubilites.</p><p>The VLE is solved by Newton&apos;s method.</p><p>Explicit functions for density and enthalpy of the phases are not specified in this package, as it is just a template.</p>
<p>The package has to be in one file, because it extends a MSL package (DYMOLA limitiation???).</p>
<p><b><font style=\"font-size: 12pt; \">Fluid model assumptions</b></p>
<p><ul>
<li>The fluid consists of water, Ns salts and Nggases.</li>
<li>Its total composition is given by vector of mass fractions X.</li>
<li>There are one or two phases: liquid and, if absolute pressure is low enough, gas.</li>
<li>The salts are completely dissolved in and limited to the liquid phase.</li>
<li>The gas phase is formed by water vapour and gases. </li>
<li>Water and gases are exchanged between the liquid and the gas phase by means of dissolution or evaporation/condensation.</li>
<li>Mass and energy conservation are fulfilled.</li>
<li>Gases dissolve in liquid depending on their respective solubility, which depends on temperature and salt content, but not on the content of other gases.</li>
<li>The saturation pressure of water is reduced by the salt content.</li>
<li>Boundary surface enthalpies are neglected. </li>
</ul></p>
<p><b>Specific Enthalpy</b></p>
<pre>h:=x&middot;h_G + (1-x)&middot;h_L</pre>
<p><b>Density</b></p>
<p>The total density <i>d</i> of the fluid is calculated by combining the densities of both phases (<i>dg</i> and <i>dl</i>) according to their volume fractions. The gas phase is assumed to be an Density of the gas phase is assumed to be an ideal mixture of ideal gases. </p>
<pre>d:=1/(x/d_g + (1 - x)/d_l)</pre>
<p>All files in this library, including the C source files are released under the Modelica License 2. </p>
<p><b>TODO:</b></p>
<p><b></font><font style=\"font-size: 10pt; \">Created by</b></p>
<p>Henning Francke</p><p>Helmholtz Centre Potsdam</p><p>GFZ German Research Centre for Geosciences</p><p>Telegrafenberg, D-14473 Potsdam</p><p>Germany </p>
<p><a href=\"mailto:francke@gfz-potsdam.de\">francke@gfz-potsdam.de</a> </p>
</html>",
 revisions="<html>

</html>"));
end PartialBrine_ngas_Newton;
