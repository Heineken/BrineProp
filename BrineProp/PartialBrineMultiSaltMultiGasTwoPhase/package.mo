within BrineProp;
partial package PartialBrineMultiSaltMultiGasTwoPhase "Template medium for aqueous solutions of m Salts and n Gases, VLE solved by Newton's method"

  //constant String explicitVars = "ph" "set of variables the model is explicit for, may be set to all combinations of ph or pT, setting pT should speed up the model in pT cases";
    constant Boolean input_ph=true
  "activates inversion for state definition by ph, slows calculation down";
    constant Boolean input_dT=false
  "activates inversion for state definition by dT, slows calculation down";

//   extends BrineProp.PartialFlags_(final nX_salt=size(saltNames, 1));
  constant String saltNames[:]={""} "TODO replace by saltorder enumeration";
   constant Integer nX_salt = size(saltNames, 1) "Number of salt components"   annotation(Evaluate=true);


   replaceable package Salt_data = BrineProp.SaltData;
   constant Salt_data.SaltConstants[:] saltConstants;


   extends Densities(final nX_salt_=nX_salt);


   extends SpecificEnthalpies(final nX_salt_=nX_salt);


   extends Viscosities(final nX_salt_=nX_salt);
    //definition of molar masses

   constant SI.MolarMass[:] MM_salt;
   constant Integer[:] nM_salt "number of ions per molecule";

  //TWO-PHASE-STUFF
   constant SI.MolarMass[:] MM_vec = cat(1,MM_salt, MM_gas, {M_H2O});
   constant Integer[:] nM_vec = cat(1,nM_salt, nM_gas, {1});

   constant String gasNames[:]={""};

   extends BrineProp.GasData(final nX_salt_=nX_salt);

   constant Integer nX_gas = size(gasNames, 1) "Number of gas components" annotation(Evaluate=true);
   constant SI.MolarMass[:] MM_gas;
   constant Integer[:] nM_gas "number of ions per molecule";


   extends PartialMixtureTwoPhaseMedium(
     final mediumName="TwoPhaseBrine",
     final substanceNames=cat(1,saltNames,gasNames,{"water"}),
     final reducedX =  true,
     final singleState=false,
     reference_X=cat(1,fill(0,nX-1),{1}));


   redeclare model extends BaseProperties(d(start=1000))
  "Base properties of medium"
     import BrineProp;
   //  MassFraction[nX] X_l(start=cat(1,fill(0,nXi),{1}))  "= X cat(1,X_salt/(X[end]*m_l),X[nX_salt+1:end])";
     SI.Density d_l "density of liquid phase";
     SI.Density d_g "density of gas phase";

     Real GVF "gas void fraction";
     Real x "gas phase mass fraction";
     SI.Pressure p_H2O;
     SI.Pressure[nX_gas + 1] p_gas "partial gas pressures in gas phase";
     SI.Pressure p_degas[nX_gas + 1] "degassing pressure in liquid phase";
     //  MassFraction[nX_gas+1] X_g "mass fractions in gas phase";
     MassFraction[nX] X_l "mass fractions in liquid phase";

     SI.SpecificEnthalpy h_l "specific enthalpy of liquid phase";
     SI.SpecificEnthalpy h_g "specific enthalpy of gas phase";

     parameter Real[nX_gas + 1] n_g_norm_start = fill(0.5,nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
   //  Real[nX_gas+1] n_g_norm_start "uncomment for Examples.InitializationTest & comment line above";

     SI.MoleFraction y_vec[:]=Utilities.massToMoleFractions(X,MM_vec);
     Real y_g[:]= p_gas/p "mole fractions in gas phase";

     Real[nX_gas + 1] x_vec = { (if X[nX_salt+i]>0 then state.X_g[i]*x/X[nX_salt+i] else 0) for i in 1:nX_gas+1}
    "Fractions of gas mass in gas phase";

   //  Real[nX_gas+1] n_g_norm;
   //  Integer z "Number of iterations in VLE algorithm";
   equation
   //  assert(max(X)<=1 and min(X)>=0, "X out of range [0...1] = "+PowerPlant.vector2string(X)+" (saturationPressure_H2O())");

     u = h - p/d;
     MM = y_vec*MM_vec;
     R  = Modelica.Constants.R/MM;

     if input_ph then
       h = specificEnthalpy_pTX(p,T,X,phase,n_g_norm_start);
     else
       h=state.h;
     end if;

     if input_dT then
       d = density_pTX(p,T,X,phase,n_g_norm_start);
     else
       d=state.d;
     end if;

   //  (state,GVF,h_l,h_g,p_gas,p_H2O,p_degas,z,n_g_norm) = setState_pTX(p,T,X,phase,n_g_norm_start);
     (state,GVF,h_l,h_g,p_gas,p_H2O,p_degas) = setState_pTX(p,T,X,phase,n_g_norm_start);
     X_l=state.X_l;
     x=state.x;
     s=state.s;
     d_g=state.d_g;
     d_l=state.d_l;

     sat.psat = sum(p_degas);
     sat.Tsat = T;
     sat.X = X_l;
   end BaseProperties;


  redeclare record extends ThermodynamicState(phase(start=0),X(start=cat(1,fill(0,nXi),{1})))
  "a selection of variables that uniquely defines the thermodynamic state"
    SpecificEnthalpy h(start=4e5) "Specific enthalpy";
    SpecificEntropy s(start=0) "Specific entropy";
    Density d(start=300) "density";
    Real x(start=0) "vapor quality on a mass basis [mass vapor/total mass]";
    MassFraction[nX_gas+1] X_g(start=cat(1,fill(0,nX_gas),{1}));
     annotation (Documentation(info="<html>

</html>"));
  end ThermodynamicState;


    redeclare function extends dewEnthalpy
  "dew curve specific enthalpy of water"
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
      input SI.Temp_K T;
      output SI.SpecificEnthalpy Delta_h_solution;
    end solutionEnthalpy;


    replaceable partial function solubilities_pTX
  "solubility calculation of gas in m_gas/m_H2O"
      input SI.Pressure p;
      input SI.Temp_K T;
      input SI.MassFraction X_l[nX] "mass fractions m_x/m_Sol";
      input SI.MassFraction X[nX] "mass fractions m_x/m_Sol";
      input SI.Pressure[nX_gas] p_gas;
      input Boolean ignoreTlimit=false "activated by temperature_phX";
    //  input SI.MolarMass MM[:] "=fill(0,nX)molar masses of components";
    //  output Molality[nX_gas] solu;
      output MassFraction solu[nX_gas] "gas concentration in kg_gas/kg_fluid";
    end solubilities_pTX;


    replaceable function fugacity_pTX
  "Calculation of nitrogen fugacity coefficient extracted from EES"
      input SI.Pressure p;
      input SI.Temp_K T;
      input MassFraction X[:]=reference_X "Mass fractions";
      input String substancename;
      output Real phi;
protected
    Types.Pressure_bar p_bar=SI.Conversions.to_bar(p);
    end fugacity_pTX;


   replaceable function density_liq_pTX "Density of the liquid phase"
     input SI.Pressure p "TODO: Rename to density_liq_pTX";
     input SI.Temp_K T;
     input MassFraction X[nX] "mass fraction m_NaCl/m_Sol";
     input SI.MolarMass MM[:] "=MM_vec =fill(0,nX) molar masses of components";
     output SI.Density d;
   end density_liq_pTX;


  redeclare function vapourQuality
  "Returns vapour quality, needs to be defined to overload function defined in PartialMixtureTwoPhaseMedium"
    input ThermodynamicState state "Thermodynamic state record";
    output MassFraction x "Vapour quality";
  algorithm
  x := state.x;
  end vapourQuality;


    redeclare function specificEnthalpy_pTX "wrapper to extract h from state"
      //necessary for declaration of inverse function T(p,h)
      input SI.Pressure p;
      input SI.Temp_K T;
      input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
      input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
      input Real[nX_gas+1] n_g_norm_start=fill(0.5,nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
      input Boolean ignoreTlimit=false;
      output SI.SpecificEnthalpy h;

    algorithm
      if debugmode then
          print("Running specificEnthalpy_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" C, ignoreTlimit="+String(ignoreTlimit)+", X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
      end if;
     h:=specificEnthalpy(setState_pTX(
        p,
        T,
        X,
        phase,
        n_g_norm_start,
        ignoreTlimit));

    //print(String(p)+","+String(T)+" K->"+String(h)+" J/kg & (PartialBrine_Multi_TwoPhase_ngas.specificEnthalpy_pTX)");
     //,p=pressure_ThX(T,h,X);

     annotation(LateInline=true,inverse(T=temperature_phX(p,h,X,phase,n_g_norm_start,ignoreTlimit)));
    end specificEnthalpy_pTX;


    redeclare function density_pTX "wrapper to extract d from state"
      //necessary for declaration of inverse function p(T,d)
      input SI.Pressure p;
      input SI.Temp_K T;
      input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
      input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
      input Real[nX_gas+1] n_g_norm_start=fill(0.5,nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
      output SI.Density d;

    algorithm
      if debugmode then
          print("Running density_pTX("+String(p/1e5)+","+String(T-273.15)+"degC, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
      end if;
     d:=density(setState_pTX(p,T,X,phase,n_g_norm_start));

     annotation(LateInline=true,inverse(p=pressure_dTX(d,T,X,phase,n_g_norm_start)));
    end density_pTX;


    redeclare function temperature_phX
  "iterative inversion of specificEnthalpy_pTX by regula falsi"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEnthalpy h "Specific enthalpy";
      input MassFraction X[nX] "Mass fractions";
      input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
      input Real[nX_gas + 1] n_g_start=fill(0.5,nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
      input Boolean ignoreTlimit=true
    "To avoid warnings due to low temperature";
      output Temperature T "Temperature";
protected
      SI.SpecificHeatCapacity c_p;
      SI.Temperature T_a=273.16 "lower window temperature";
      SI.Temperature T_b=min(400,Modelica.Media.Water.WaterIF97_pT.saturationTemperature(p)-1)
    "upper window temperature";

      SI.SpecificEnthalpy h_a "h at lower limit";
      SI.SpecificEnthalpy h_b "h at upper limit";
      SI.SpecificEnthalpy h_T "calculated h";
      Integer z=0 "Loop counter";
    algorithm
      if debugmode then
         print("temperature_phX("+String(p)+","+String(h)+")");
      end if;

      if not max(X[1:end-1])>0 then
        if debugmode then
          print("2Phase water (temperature_phX("+String(p)+","+String(h)+"))");
        end if;
          T :=Modelica.Media.Water.WaterIF97_pT.temperature_ph(p, h);
          //T = Modelica.Media.Water.IF97_Utilities.T_ph(p,h);
      else

        //Find temperature with h above given h ->T_b
        assert(h>specificEnthalpy_pTX(p,T_a,X,0,n_g_start,ignoreTlimit=true),"h="+String(h/1e3)+" kJ/kg -> Enthalpy below lower bisection limit ("+String(T_a-273.15)+" C)");
        while true loop
          h_T:=specificEnthalpy_pTX(p,T_b,X,0,n_g_start,ignoreTlimit=true);
          //print(String(p)+","+String(T_b)+" K->"+String(h_T)+" J/kg (PartialBrine_ngas_Newton.temperature_phX)");
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
      //    print("T_b-T_a="+String(T_b-T_a)+", abs(h-h_T)/h="+String(abs(h-h_T)/h));
          T:=(T_a+T_b)/2 "Halbieren";
      //    print("T_neu="+String(T)+"K");
          h_T:=specificEnthalpy_pTX(p,T,X,0,n_g_start,ignoreTlimit=true);
          if h_T > h then
            T_b:=T;
      //      print("T_b="+String(T)+"K -> dh="+String(h_T-h));
          else
            T_a:=T;
      //      print("T_a="+String(T)+"K -> dh="+String(h_T-h));
          end if;
          z:=z+1;
      //    print(String(z)+": "+String(T_a)+" K & "+String(T_b)+" K -> "+String((h-h_T)/h)+"(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)\n");
      //    print("h("+String(T_a)+")="+String(h_a-h)+" J/kg & h("+String(T_b)+")="+String(h_b-h)+" J/kg");
          assert(z<100,"Maximum number of iteration reached for temperature calculation. Something's wrong here. Cancelling...(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)");
        end while;
      // print("BISECTION " + String(z)+": "+String(T));

          /*
      //REGULA FALSI - is langsamer, braucht 19 Iterationen
        z:=0;
        T_a:=T0_a;
        T_b:=T0_b "limit of N2 solubility";
        h_a := specificEnthalpy_pTX(p,T_a,X);
        h_b := specificEnthalpy_pTX(p,T_b,X);
        while abs(T_b-T_a)>1e-2 and abs(h_T-h)/h>1e-5 loop
      //  while abs(T_b-T_a)/T_l>1e-4 loop
          print("h_a("+String(T_a)+")="+String(h_a)+" / h_b("+String(T_b)+")="+String(h_b));
          T:=max(T0_a,min(T0_b,T_a-(T_b-T_a)/(h_b-h_a)*(h_a-h))) "Regula falsi";
          h_T:=specificEnthalpy_pTX(p,T,X);
          print("T_neu="+String(T)+"K");
          if h_T > h then
            T_b:=T;
            h_b:=h_T;
          else
            T_a:=T;
            h_a:=h_T;
      //      print("T_a="+String(T)+"K -> h="+String(h_T-h));
          end if;
          z:=z+1;
      //    print(String(z)+": "+String(T_a)+" K & "+String(T_b)+" K -> "+String((h-h_T)/h)+"(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)\n");
      //    print("h("+String(T_a)+")="+String(h_a-h)+" J/kg & h("+String(T_b)+")="+String(h_b-h)+" J/kg");
          assert(z<100,"Maximum number of iteration reached for temperature calculation. Something's wrong here. Cancelling...(PartialBrine_Multi_TwoPhase_ngas.temperature_phX)");
        end while;
        //   print("REGULA FALSI " + String(z)+": "+String(T));
      */
      end if;
    end temperature_phX;


    redeclare function pressure_dTX
  "iterative inversion of density_pTX by regula falsi"
      extends Modelica.Icons.Function;
      input SI.Density d;
      input SI.Temperature T;
      input MassFraction X[nX] "Mass fractions";
      input FixedPhase phase=0
    "2 for two-phase, 1 for one-phase, 0 if not known";
      input Real[nX_gas + 1] n_g_start=fill(0.5, nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
      output SI.Pressure p;
protected
      SI.SpecificHeatCapacity c_p;
      SI.Pressure p_a=1e5 "lower limit";
      SI.Pressure p_b=999e5 "lower limit";

      SI.Density d_a "density at lower limit";
      SI.Density d_b "density at upper limit";
      SI.Density d_p "density calculated for tentative pressure";
      Integer z=0 "Loop counter";

    algorithm
      if debugmode then
         print("pressure_TdX("+String(T)+","+String(d)+")");
      end if;
      //Find temperature with h above given h ->T_b
      assert(d>density_pTX(p_a,T,X),"d="+String(d)+" kg/m^3 -> density too low (< d(1 bar))");
      while true loop
        d_p:=density_pTX(p_b,T,X);
        //print(String(p_b/1e5)+" bar,"+String(T)+" K->"+String(d_p)+" kg/m^3 (PartialBrine_ngas_Newton.pressure_TdX)");
        if d>d_p then
          p_a := p_b;
          p_b := p_b + 100;
        else
          break;
        end if;
      end while;

    //BISECTION - is schneller, braucht 13 Iterationen
      while (p_b-p_a)>1e-2 and abs(d-d_p/d)>1e-5 loop   //stop when temperatures or enthalpy are close
    //  while abs(h-h_T/h)>1e-5 loop
    //    print("T_b-T_a="+String(T_b-T_a)+", abs(h-h_T)/h="+String(abs(h-h_T)/h));
        p:=(p_a+p_b)/2 "Halbieren";
    //    print("T_neu="+String(T)+"K");
        d_p:=density_pTX(p,T,X);
        if d_p > d then
          p_b:=p;
    //      print("T_b="+String(T)+"K -> dh="+String(h_T-h));
        else
          p_a:=p;
    //      print("T_a="+String(T)+"K -> dh="+String(h_T-h));
        end if;
        z:=z+1;
    //    print(String(z)+": "+String(p_a)+" K & "+String(p_b)+" K -> "+String((h-h_T)/h)+"(PartialBrine_Multi_TwoPhase_ngas.pressure_TdX)\n");
    //    print("h("+String(p_a)+")="+String(d_a-h)+" J/kg & h("+String(p_b)+")="+String(d_b-d)+" J/kg");
        assert(z<100,"Maximum number of iteration reached for pressure calculation. Something's wrong here. Cancelling...");
      end while;
    // print("BISECTION " + String(z)+": "+String(T));

    end pressure_dTX;


   replaceable function specificEnthalpy_liq_pTX
  "Specific enthalpy of liquid phase"
     input SI.Pressure p;
     input SI.Temp_K T;
     input MassFraction X[nX] "mass fraction m_NaCl/m_Sol";
   //  input SI.MolarMass MM[:]=fill(0,nX) "molar masses of components";
     input Boolean ignoreTlimit=false "activated by temperature_phX";
     output SI.SpecificEnthalpy h;
   end specificEnthalpy_liq_pTX;


   replaceable function specificEnthalpy_gas_pTX
  "Specific enthalpy of gas in gas phase"
     input SI.Pressure p;
     input SI.Temp_K T;
   //  input SI.MolarMass MM[:]=fill(0,nX)     "molar masses of components";
     input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
     output SI.SpecificEnthalpy h;
   end specificEnthalpy_gas_pTX;


    replaceable partial function saturationPressures
  "Return saturationPressures for gases and water"
      extends Modelica.Icons.Function;
      input SI.Pressure p;
      input SI.Temp_K T;
      input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
      input SI.MolarMass MM[:] "molar masses of components";
      output SI.Pressure[nX_gas] p_sat;
    end saturationPressures;



    redeclare replaceable partial function extends setState_pTX
  "finds the VLE iteratively by varying the normalized quantity of gas in the gasphase, calculates the densities"
    input Real[nX_gas + 1] n_g_norm_start= fill(0.1,nX_gas+1)
    "start value, all gas in gas phase, all water liquid, set in BaseProps";
    /*
//output SI.Density d_g= if x>0 then (n_CO2_g*d_g_CO2 + n_N2_g*d_g_N2)/(n_CO2_g + n_H2O_g) else -1;
//output Real[nX_gas + 1] n_g_norm;
//output Real k[nX_gas];
// SI.Density d_g_H2O = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.rhov_p(p) "density of water vapor";
*/
      input Boolean ignoreTlimit=false "activated by temperature_phX";
      output Real GVF;
      output SI.SpecificEnthalpy h_l;
      output SI.SpecificEnthalpy h_g;
      output SI.Pressure[nX_gas + 1] p_gas "partial pressures of gases";
      output SI.Pressure p_H2O "water vapour pressure TODO is in p_gas drin";
      output SI.Pressure[nX_gas + 1] p_degas;
    //  output Integer z "number of iterations";
protected
      Integer z "number of iterations";
      SI.MassFraction[nX_gas+1] X_g;
      SI.MassFraction[nX] X_l=X "start value";
      SI.Density d;
      SI.Density d_l;
      SI.Density d_g;
      SI.MassFraction x;
      //SI.Pressure p_sat_H2O "water vapour pressure considering salinity";
      SI.Pressure p_H2O_0 "pure water vapour pressure";
      SI.Pressure[nX_gas + 1] f
    "componentwise pressure disbalance (to become zero)";
      SI.Pressure[nX_gas + 1] p_sat=fill(0,nX_gas + 1);
      SI.MassFraction[nX_gas + 1] Delta_n_g_norm = fill(1e3,nX_gas+1)
    "large initial value to enter while loop";
    //  SI.MassFraction[nX_gas + 1] c = {3.16407e-5,0,3.6e-8,.746547} "cat(1,fill(1e-4, nX_gas), {X[end]})fill(0, nX_gas+1)X[nX_salt+1:end]";
      Real solu[nX_gas];
      Real k[nX_gas];
      Real[nX_gas + 1] n "Total mol numbers";
      Real[nX_gas + 1] n_l "mols in liquid phase per kg fluid";
      Real[nX_gas + 1] n_g "mols in   gas  phase per kg fluid";
      Real[nX_gas + 1] n_g_norm_test;
    //  SI.MassFraction[nX] X;
      Real[nX_gas + 1] n_g_norm
    "= X[end-nX_gas:end-1]./MM_gas fill(0,nX_gas) - start value: all degassed";
      Real dp_gas_dng_norm;
      Real dcdng_norm;
      Real dp_degas_dng_norm;
      Real[nX_gas + 1] dfdn_g_norm;
      Real sum_n_ion;
      constant Integer zmax=1000 "maximum number of iterations";
    //  Integer ju = nX_gas+1;
      Real[nX_gas + 1,nX_gas + 1] Grad_f;
      Real DeltaC=0.001;
      SI.Temperature T2;
      SpecificHeatCapacity R_gas;
      Boolean isTwoPhaseWater=false;
    algorithm
      if debugmode then
          print("\nRunning setState_pTX("+String(p/1e5)+" bar,"+String(min(1000,T)-273.15)+" degC, ignoreTlimit="+String(ignoreTlimit)+", X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
      end if;

     assert(p>0,"p="+String(p/1e5)+" bar - Negative pressure is not yet supported.");
     assert(max(X)-1<=1e-8 and min(X)>=-1e-8, "X out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X])));

      if String(T)== "-1.#IND" then
      assert(false, "T is NaN, probably division by zero."+(if AssertLevel==1 then "Setting T=300 K." else ""),aLevel);
        T2:=350;
      else
        assert(T>273.15,"T="+String(T-273.15)+" C too low (<0 C)."+(if AssertLevel==1 then " setting to 0 C" else ""),aLevel);
        T2:= max(273.16,T);
      end if;

      p_H2O := saturationPressure_H2O(p,T2,X,MM_vec,nM_vec);

    /*     p_degas := cat(1,saturationPressures(p,T2,X,MM_vec), {p_H2O}); 60% slower*/
        p_gas :=fill(p/(nX_gas + 1), nX_gas + 1);
        solu :=solubilities_pTX(p,T,X_l,X,p_gas[1:nX_gas],ignoreTlimit=ignoreTlimit);

        k :=solu ./ p_gas[1:nX_gas];
        for i in 1:nX_gas loop
            p_degas[i] :=X_l[nX_salt + i]/(if k[i] > 0 then k[i] else 1^10)
      "Degassing pressure";
        end for;
        p_degas[nX_gas + 1] :=p_H2O;

       if phase==1 or sum(p_degas) < p then
        if debugmode then
          print("1Phase-Liquid (setState_pTX("+String(p)+","+String(T2)+"))");
        end if;
        x:=0;
        //p_H2O := p_sat_H2O;
       elseif not max(X[1:end-1])>0 then
        isTwoPhaseWater:=true;
          if debugmode then
            print("2Phase water (PartialBrine_Multi_TwoPhase_ngas.setState_pTX("+String(p)+","+String(T2)+"))");
          end if;
          X_l:=X_l;
       // h:=Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p, T);
        h_l :=Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hl_p(p);
        h_g :=Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hv_p(p);
        x :=min(max((Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p, T) - h_l)/(h_g - h_l + 1e-18), 0), 1);
        d_l :=Modelica.Media.Water.IF97_Utilities.rhol_T(T);
        d_g :=Modelica.Media.Water.IF97_Utilities.rhov_T(T);

       else
        assert(max(X[end-nX_gas:end-1])>0,"Phase equilibrium cannot be calculated without dissolved gas at "+String(p/1e5)+" bar, "+String(T2-273.15)+"degC with p_degas="+String(sum(p_degas)/1e5)+" bar.");

        n:=X[nX_salt + 1:end] ./ MM_vec[nX_salt + 1:nX]
      "total mole numbers per kg brine";
        n_g_norm:=n_g_norm_start .* sign(X[nX_salt + 1:nX])
      "switch off unused salts";

        z:=0;
        while z<1 or max(abs(Delta_n_g_norm))>1e-3 loop
        //stop iteration when p-equlibrium is found or gas fraction is very low
          z:=z + 1 "count iterations";
          assert(z<=zmax,"Reached maximum number of iterations ("+String(z)+"/"+String(zmax)+") for solution equilibrium calculation. (setState_pTX("+String(p/1e5)+"bar,"+String(T2-273.16)+"degC))\nDeltaP="+String(max(abs(p_sat-p_gas))));

          n_g :=n_g_norm .* n;
          n_l := n-n_g;
          x := n_g*MM_vec[nX_salt+1:nX];
    //      print("\n"+String(z)+": x="+String(x)+" Delta_n_g="+String(max(abs(Delta_n_g_norm))));
          X_l:=cat(1, X[1:nX_salt], n_l.*MM_vec[nX_salt+1:nX])/(1-x);
       /*      print("n_l=" + PowerPlant.vector2string(n_l));
        print("n_g=" + PowerPlant.vector2string(n_g));
        */
      //PARTIAL PRESSURE
            p_gas := p * n_g/sum(n_g);

      //DEGASSING PRESSURE
            (p_H2O,p_H2O_0):=saturationPressure_H2O(p,T2,X_l,MM_vec,nM_vec)
        "X_l changes";
        if (p_H2O>p) then
            print("p_H2O(" + String(p/1e5) + " bar," +
              String(T2 - 273.15) + "degC, " + Modelica.Math.Matrices.toString(transpose([X])) + ") = "
               + String(p_H2O/1e5) + " bar>p ! (PartialBrineMultiSaltMultiGasTwoPhase.setState_pTX())");
          x:=1;
          break;
        end if;

     //      print("p_gas[1]=" + String(p_gas[1]));

            solu:=solubilities_pTX(p=p, T=T2, X_l=X_l, X=X, p_gas=p_gas[1:nX_gas],ignoreTlimit=ignoreTlimit);

            for i in 1:nX_gas loop
              k[i]:=if p_gas[i] > 0 then solu[i]/p_gas[i] else 1e10;
            end for;
            p_sat[1:nX_gas] := X_l[nX_salt+1:nX-1]./ k;
            p_sat[nX_gas+1] := p_H2O;

            f :=  p_gas-p_sat;
    //       print("p_gas=" + PowerPlant.vector2string(p_gas) + "=>" + String(sum(p_gas)));

           sum_n_ion :=cat(1,X[1:nX_salt] ./ MM_vec[1:nX_salt],n_l)*nM_vec;

      //GRADIENT analytisch df[gamma]/dc[gamma]

            for gamma in 1:nX_gas+1 loop
                dp_gas_dng_norm:=p*n[gamma]* (sum(n_g)-n_g[gamma])/(sum(n_g))^2
          "partial pressure";
                if gamma == nX_gas+1 then
                  dp_degas_dng_norm := p_H2O_0*n[end]*( -sum_n_ion +(1-n_g_norm[end])*n[gamma]) / sum_n_ion^2;
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
//            print("dcdng_norm("+String(alpha)+","+String(gamma)+")=" + String(dcdng_norm));
            end if;
            Grad_f[gamma,alpha] := dp_gas_dng_norm-dp_degas_dng_norm;

/*           print("dp_gas_dng_norm("+String(gamma)+","+String(alpha)+")=" + String(dp_gas_dng_norm));
           print("dp_degas_dng_norm("+String(gamma)+","+String(alpha)+")=" + String(dp_degas_dng_norm));
           * /

          end for;
//         print("Grad_f["+String(gamma)+",:] =" + PowerPlant.vector2string(Grad_f[gamma,:]));
        end for;
*/

    //       print("k=" + PowerPlant.vector2string(k));/**/
    //       print("dp_gas=" + PowerPlant.vector2string(p_sat - p_gas));

      //SOLVE NEWTON STEP
    //        Delta_n_g_norm := Modelica.Math.Matrices.solve(Grad_f, -f)         "solve Grad_f*Delta_n_g_norm=-f";
    //        n_g_norm := n_g_norm + Delta_n_g_norm;

    //        print("n_g_norm="+Modelica.Math.Matrices.toString({n_g_norm}));
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
    //       print("p_sat="+String(p_sat[1])+", solu="+String(solubility_CO2_pTX_Duan2006(p,T2,X_l,MM_vec,p_gas[1]))+", p_gas="+String(p_gas[1]));
    //         print("p="+String(p)+",T2="+String(T2)+",p_gas[1]="+String(p_gas[1]));
    /*        print("X_l="+Modelica.Math.Matrices.toString({X_l}));
        print("MM_vec="+Modelica.Math.Matrices.toString({MM_vec}));
*/
        end while;

      end if "p_degas< p";

      if not isTwoPhaseWater then
    //DENSITY
       X_g:=if x>0 then (X[end-nX_gas:end]-X_l[end-nX_gas:end]*(1-x))/x else fill(0,nX_gas+1);
      /*Calculation here  R_gas :=if x > 0 then sum(Modelica.Constants.R*X_g ./ cat(1,MM_gas,{M_H2O})) else -1;
    d_g :=if x > 0 then p/(T2*R_gas) else -1;*/
      //  d_g:= if x>0 then p/(Modelica.Constants.R*T2)*(n_g*cat(1,MM_gas,{M_H2O}))/sum(n_g) else -1;
        if x > 0 then
          d_g :=BrineGas3Gas.density_pTX(p,T,X_g);
          h_g := specificEnthalpy_gas_pTX(p,T,X_g);
        else
          d_g := -1;
          h_g := -1;
        end if;
        d_l :=if not x < 1 then -1 else density_liq_pTX(
              p,
              T2,
              X_l,
              MM_vec) "no 1-phase gas";
        h_l := specificEnthalpy_liq_pTX(p,T,X_l,ignoreTlimit);
      end if "TwoPhaseWater";

      d:=1/(x/d_g + (1 - x)/d_l);
    //  d:=1/(x/max(1e-8,d_g) + (1 - x)/d_l);
    //  print(String(z)+" (p="+String(p_gas[1])+" bar)");

    // X_g:=if x>0 then (X-X_l*(1-x))/x else fill(0,nX);
     GVF:=x*d/d_g;
     state :=ThermodynamicState(
        p=p,
        T=T,
        X=X,
        X_l=X_l,
        X_g=X_g,
        h=x*h_g + (1-x)*h_l,
        x=x,
        s=0,
        d=d,
        d_l=d_l,
        d_g=max(0,d_g),
        phase=if x>0 and x<1 then 2 else 1) "phase_out";
    /*    h_g=h_g,
    h_l=h_l,
    p_H2O=p_H2O,
    p_gas=p_gas[1:nX_gas],
    p_degas=p_degas
    X_l={max(0,X_l[i]) for i in 1:nX}, didn't help!
    X_g={max(0,X_g[i]) for i in 1:nX_gas+1}, didn't help"
*/

      annotation (Diagram(graphics={Text(
              extent={{-96,16},{98,-16}},
              lineColor={0,0,255},
              textStyle={TextStyle.Bold},
              textString="find static VLE")}));
    end setState_pTX;


  redeclare replaceable partial function extends setState_phX
  "Calculates medium properties from p,h,X"
  //      input String fluidnames;
protected
    SI.SpecificEnthalpy h_l;
    SI.SpecificEnthalpy h_v;
    SI.Density d_l;
    SI.Density d_g;
    Real x;
    constant SI.SpecificEnthalpy eps = 1e-8;
  algorithm

    if debugmode then
      print("Running setState_phX("+String(p/1e5)+" bar,"+String(h)+" J/kg,X)...");
    end if;

    if not max(X[1:end-1])>0 then
      if debugmode then
        print("2Phase water (setState_phX("+String(p)+","+String(h)+"))");
      end if;
    h_l :=Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hl_p(p);
    h_v:=Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hv_p(p);
    d_l:=Modelica.Media.Water.IF97_Utilities.rhol_p(p);
    d_g:=Modelica.Media.Water.IF97_Utilities.rhov_p(p);
    x:=min(max((h - h_l)/(h_v - h_l + 1e-18), 0), 1);
    state := ThermodynamicState(
      d=1/(x/d_g + (1 - x)/d_l),
      T=temperature_phX(p,h,X),
      phase=0,
      h=h,
      p=p,
      X=X,
      X_l=X,
      X_g=fill(0,nX_gas+1),
      s=0,
      x=x,
      d_l=  d_l,
      d_g=  d_g);
    else
      state := setState_pTX(p,
        temperature_phX(p,h,X,phase),
        X,
        phase) ",fluidnames)";
    end if;
  end setState_phX;


    redeclare replaceable function extends specificHeatCapacityCp
  "numeric calculation of specific heat capacity at constant pressure"
protected
      SI.SpecificHeatCapacity cp_liq=specificHeatCapacityCp_liq(state);
      SI.SpecificHeatCapacity cp_gas=specificHeatCapacityCp_gas(state);
    algorithm
      cp:=state.x*cp_gas + (1-state.x)*cp_liq;

    //  assert(cp>0 and cp<5000,"T="+String(state.T-273.15)+"K, p="+String(state.p/1e5)+"bar, x="+String(state.x)+", cp_liq="+String(cp_liq)+"J(kgK), cp_gas="+String(cp_gas)+"J(kgK)");
    //  print("c_p_liq("+String(state.T)+"degC)="+String(p)+" J/(kg.K)");

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
  constant SI.TemperatureDifference dT=.1;
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
  constant SI.TemperatureDifference dT=.1;
algorithm 
//    cp := Modelica.Media.Water.IF97_Utilities.cp_pT(state.p, state.T) "TODO";
    cp:=(specificEnthalpy_pTX(state.p,state.T+dT,state.X)-state.h)/dT;
    */

protected
      SI.MassFraction[nX] X_g=if state.x>0 then (state.X-state.X_l*(1-state.x))/state.x
     else
         fill(-1,nX);
      SI.SpecificHeatCapacity cp_vec[nX_gas+1];
    end specificHeatCapacityCp_gas;


    replaceable function isobaricExpansionCoefficient_liq
    //  extends isobaricExpansionCoefficient;
      input ThermodynamicState state;
      input SI.Density d_l;
      output SI.LinearTemperatureCoefficient beta;
protected
      constant SI.Temperature Delta_T= 1;
    algorithm
    //  beta :=d_l*(1/d_l - 1/(density_liquid_pTX(state.p,state.T - Delta_T,state.X,MM_vec)))/Delta_T;
      beta :=(1 - d_l/(density_liq_pTX(
            state.p,
            state.T - Delta_T,
            state.X,
            MM_vec)))/Delta_T;
    end isobaricExpansionCoefficient_liq;


  redeclare replaceable function specificEnthalpy
    input ThermodynamicState state "Thermodynamic state record";
    output SI.SpecificEnthalpy h;
  algorithm
    h := state.h;
  end specificEnthalpy;


  redeclare function extends density "return density of ideal gas"
  algorithm
    d := state.d;
  end density;


    annotation (Documentation(info="<html>
<ul>
<li><b>PartialBrine_ngas_Newton</b> is based on <code>PartialMixtureTwoPhaseMedium</code>, an extension to the <code>Modelica.Media</code> library. This extension was necessary because <code>Modelica.Media</code> supports mixtures and two-phase media, but not both combined.</li>
<li>The vapour-liquid-equilibrium (VLE) is defined by the water vapour pressure and the gas solubilites. It is determined using Newton&apos;s method.</li>
<li>Explicit material functions are not specified in this package, as it is just a template. They need to be provided in the instantiating package (e.g.<code> BrineProp.Brine_5salts_TwoPhase_3gas</code>).</li>
<li>The model calculates properties for a thermodynamic state specified by <i>p</i> and <i>T</i>, <i>p</i> and <i>h</i>, <i>T</i> and <i>d</i> or <i>p</i> and <i>d</i>.</li>
</ul>
<h4>Fluid model assumptions</h4>
<ul>
<li>The fluid consists of water, <i>N<sub>s</sub></i> salts and <i>N<sub>g</sub></i> gases.</li>
<li>Its total composition is given by a vector of mass fractions X.</li>
<li>There are one or two phases: liquid and, if absolute pressure is low enough, gas (no solid phase).</li>
<li>The gas phase is an ideal mixture of water vapour and gases.</li>
<li>The salts are completely dissolved in and limited to the liquid phase (no precipitation/evaporation). </li>
<li>Water and gases are exchanged between both phases by degassing/dissolution or evaporation/condensation, taking into account mass and energy conservation.</li>
<li>Gases dissolve in liquid depending on their respective solubility, which depends on temperature and salt content, but not on the content of other gases.</li>
<li>The saturation pressure of water is reduced by the salt content.</li>
<li>Water evaporation and condensation depend on its saturation pressure, which depends on temperature and the salt content according to Raoult&apos;s law.</li>
<li>In two-phase state degassing pressures equal the respective partial pressures.</li>
<li>Both phases are assumed to be in thermodynamic equilibrium, i.e. they have the same pressure and temperature. The vapour-liquid equilibrium is instantly reached.</li>
<li>Boundary surface enthalpies are neglected.</li>
</ul>
<p>See <a href=\"http://nbn-resolving.de/urn:nbn:de:kobv:83-opus4-47126\">PhD thesis</a> for more details.</p>
<h4>Details</h4>
<p>Brine is a mixture of components and, if gases are involved, potentially has two phases. As this is not possible with the components provided in Modelica.Media a new Medium template had to be created by merging Modelica.Media.Interfaces.PartialMixtureMedium and Modelica.Media.Interfaces.PartialTwoPhaseMedium of the Modelica Standard Library 3.1. </p>
<p>The model is explicit for p and T, but for h(p,T) the inverse function T(p,h) is defined. T(p,h) is inverts h(p,T) numerically by bisection, stopping at a given tolerance.</p>
<p>In order to calculate h(p,T), the vapour-liquid-equilibrium (VLE) is determined, i.e. the gas mass fraction q and the compositions of the liquid phases X_l. Only h is returned, due to the limitation of DYMOLA/Modelica not allowing inverse functions of functions that are returning an array. As x (gas mass fraction) and X_l (composition of liquid phase) are of interest themselves and required to calculate density and viscosity, the VLE calculation is conducted one more time, this time with T known. This additiona unnecessary calculation doubles the workload when p,h are given. When p,T are given, however, it adds only one more calculation to the multiple iterations of the bisection algorithm. </p>
<h4><span style=\"color:#008000\">Usage</span></h4>
<p>This is a partial medium and cannot be used as is.</p>
<p>See <code>BrineProp.Examples.BrineProps2phase</code> for usage example.</p>
<p>See <code><a href=\"Modelica://BrineProp.Examples.BrineProps2phase\">BrineProp.Examples.BrineProps2phase</a></code> or info of <code><a href=\"Modelica://BrineProp.Brine_5salts_TwoPhase_3gas\">BrineProp.Brine_5salts_TwoPhase_3gas</a></code> for more usage examples.</p>

<h5>TODO:</h5>
<h5>Known Issues:</h5>
<ul>
<li>The package is in one file, because it extends a MSL package (DYMOLA limitiation?).</li>
</ul>
<h5>Created by</h5>
<div>Henning Francke<br/>
Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences<br/>
Telegrafenberg, D-14473 Potsdam<br/>
Germany</div>
<p><a href=\"mailto:info@xrg-simulation.de\">francke@gfz-potsdam.de</a></p>
</html>",
   revisions="<html>

</html>"));
end PartialBrineMultiSaltMultiGasTwoPhase;
