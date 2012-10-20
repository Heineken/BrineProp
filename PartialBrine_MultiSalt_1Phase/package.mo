within BrineProp;
partial package PartialBrine_MultiSalt_1Phase "Template for one-phase (liquid) brine based on PartialMediaMixtureMedium"
  extends Modelica.Media.Interfaces.PartialMixtureMedium(
   final mediumName="TwoPhaseMixtureMedium",
   final substanceNames=cat(1,saltNames,{"water"}),
   final reducedX =  true,
   final singleState=false,
   reference_X=cat(1,fill(0,nX-1),{1}),
   fluidConstants = BrineConstants);

  constant FluidConstants[nS] BrineConstants(
     each chemicalFormula = "H2O+NaCl+KCl+CaCl2+MgCl2+SrCl2",
     each structureFormula="H2O+NaCl+KCl+CaCl2+MgCl2+SrCl2",
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

  constant String explicitVars = "ph"
  "set of variables the model is explicit for, may be set to all combinations of ph or pT, setting pT should speed up the model in pT cases";


 replaceable package Salt_data = BrineProp.SaltData;

  import Partial_Units;

 constant Real[:] MM_salt;
 constant Integer[:] nM_salt "number of ions per molecule";

 constant Modelica.SIunits.MolarMass MM_vec = cat(1,MM_salt, {M_H2O});
 constant Modelica.SIunits.MolarMass nM_vec = cat(1,nM_salt, {1});

 constant String saltNames[:]={""};

  constant Integer nX_salt = size(saltNames, 1) "Number of salt components"   annotation(Evaluate=true);


  redeclare function extends dynamicViscosity "viscosity calculation"
  algorithm
    eta:=dynamicViscosity_pTX(
      state.p,
      state.T,
      state.X);
  end dynamicViscosity;

  replaceable function dynamicViscosity_pTX "viscosity calculation"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
    output Modelica.SIunits.DynamicViscosity eta;
  //  constant Real M_NaCl=0.058443 "molar mass in [kg/mol]";
  end dynamicViscosity_pTX;


 redeclare model extends BaseProperties "Base properties of medium"

    //    PowerPlant.Media.TableLookup Table;
    //  protected
    /*     constant Modelica.SIunits.MolarMass M_H2O = PartialBrine.M_H2O "[kg/mol]";
     constant Modelica.SIunits.MolarMass M_NaCl = PartialBrine.M_NaCl 
        "[kg/mol]";*/
   Real y_vec[:]=massFractionsToMoleFractions(X,MM_vec);
 equation
   d = density_pTX(p,T,X);
   h = specificEnthalpy_pTX(p,T,X);
 //  T = temperature_phX(p,h,X);
   u = 1 "h - p/d";
   MM = y_vec*MM_vec;
   R  = 8.3144/MM;

   state.p = p;
   state.T = T;
   state.X = X;

 //  state.s = 0 "specificEntropy_phX(p,h,X)";
 //  state.h = h;
 //  state.d = d;

   annotation (Documentation(revisions="<html>

</html>"));
 end BaseProperties;


  redeclare replaceable partial function density_pTX
  "Return density from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:] "Mass fractions";
    input Modelica.SIunits.MolarMass MM[:]={1} "molar masses of components";
    output Density d "Density";
    annotation(Documentation(info="<html></html>"));
  end density_pTX;


  redeclare replaceable function specificEnthalpy_pTX
     input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
    output Modelica.SIunits.SpecificEnthalpy h;

  /*algorithm 
  h := 4*T;
*/
  end specificEnthalpy_pTX;


  redeclare function temperature_phX
  "iterative inversion of specificEnthalpy_pTX by regula falsi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[nX] "Mass fractions";
    input Real[nX_gas + 1] n_g_start=fill(.5,nX_gas+1)
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
  // Modelica.Utilities.Streams.print(String(p)+","+String(T_b)+" K->"+String(h_T)+" J/kg (PartialBrine_ngas_Newton.temperature_phX)");
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
  //      Modelica.Utilities.Streams.print("T_b="+String(T)+"K -> h="+String(h_T-h));
      else
        T_a:=T;
  //      Modelica.Utilities.Streams.print("T_a="+String(T)+"K -> h="+String(h_T-h));
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


  redeclare replaceable partial function extends setState_pTX
  "finds the VLE iteratively by varying the normalized quantity of gas in the gasphase, calculates the densities"
  input Real[PartialBrine_ngas_Newton.nX_gas + 1]
                       n_g_start=fill(.5,nX_gas+1)
    "start value, all gas in gas phase, all water liquid";
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
    Modelica.SIunits.MassFraction[nX] X_l=X "start value";
    Modelica.SIunits.Pressure p_H2O;
    Modelica.SIunits.MassFraction x;
    Modelica.SIunits.Pressure p_degas;
    Modelica.SIunits.Pressure p_sat_H2O
    "= saturationPressure_H2O(p,T2,X,MM_vec,nM_vec)";
    Modelica.SIunits.Pressure p_H2O_0;
    Modelica.SIunits.Pressure[PartialBrine_ngas_Newton.nX_gas + 1]
                                        f;
    Modelica.SIunits.Pressure[PartialBrine_ngas_Newton.nX_gas + 1]
                                          p_sat;
    Modelica.SIunits.Pressure[PartialBrine_ngas_Newton.nX_gas + 1]
                                          p_sat_test;
    Modelica.SIunits.Pressure[PartialBrine_ngas_Newton.nX_gas + 1]
                                        p_gas "=fill(0,nX_gas)";
    Modelica.SIunits.MassFraction[PartialBrine_ngas_Newton.nX_gas + 1]
                                              Delta_n_g_norm = fill(1e3,nX_gas+1);
  //  Modelica.SIunits.MassFraction[nX_gas + 1] c = {3.16407e-5,0,3.6e-8,.746547} "cat(1,fill(1e-4, nX_gas), {X[end]})fill(0, nX_gas+1)X[nX_salt+1:end]";
    Real k_H2O "Henry coefficient";
    Real k[PartialBrine_ngas_Newton.nX_gas];
    Real[PartialBrine_ngas_Newton.nX_gas + 1]
                     n "Total mol numbers";
    Real[PartialBrine_ngas_Newton.nX_gas + 1]
                     n_l "mols in liquid phase per kg fluid";
    Real[PartialBrine_ngas_Newton.nX_gas + 1]
                     n_g "mols in   gas  phase per kg fluid";
    Real[PartialBrine_ngas_Newton.nX_gas + 1]
                     n_g_norm_test;
  //  Modelica.SIunits.MassFraction[nX] X;
    Real[PartialBrine_ngas_Newton.nX_gas + 1]
                     n_g_norm
    "= X[end-nX_gas:end-1]./MM_gas fill(0,nX_gas) - start value: all degassed";
    Real dp_gas_dng_norm;
    Real dcdng_norm;
    Real dp_degas_dng_norm;
    Real[PartialBrine_ngas_Newton.nX_gas + 1]
                   dfdn_g_norm;
    Integer z=0;
    Real sum_n_ion;
    constant Integer zmax=1000 "maximum number of iteration";
  //  Integer ju = nX_gas+1;
    Real[PartialBrine_ngas_Newton.nX_gas + 1,PartialBrine_ngas_Newton.nX_gas + 1]
                            Grad_f;
    Real DeltaC=.001;
    Modelica.SIunits.Temperature T2;
  algorithm
    if debugmode then
      Modelica.Utilities.Streams.print("Running setState_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+"°C,X)...");
    end if;

  assert(p>0,"p="+String(p/1e5)+"bar - Negative pressure is not yet supported ;-) (PartialBrine_ngas_Newton.quality_pTX())");
  /*  Modelica.Utilities.Streams.print("quality_pTX("+String(p)+","+String(T2)+","+PowerPlant.vector2string(X_l[1:end],false)+")");
  X[1:nX_salt] = X_[1:nX_salt];
  for i in nX_salt+1:nX-1 loop
    X[i]:=max(0,min(1e-3,X_[i]));
  end for;
  X[end]=1-sum(X[1:end-1]);
  X_l:=X;*/
  //  Modelica.Utilities.Streams.print("quality_pTX("+String(p)+","+String(T)+","+PowerPlant.vector2string(X[1:end],false)+")");

    assert(max(X)<=1 and min(X)>=0, "X out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X]))+" (quality_pTX())");
  //  assert(T>273.15,"T too low in PartialBrine_ngas_Newton.()");
  //    Modelica.Utilities.Streams.print("\nn_g_start=" + PowerPlant.vector2string(n_g_start));

      if T<273.15 then
      Modelica.Utilities.Streams.print("T="+String(T)+" too low (<0°C), setting to 0°C in PartialBrine_ngas_Newton.quality_pTX()");
    end if;
    T2:= max(273.16,T);

    p_sat_H2O := saturationPressure_H2O(p,T2,X,MM_vec,nM_vec);
    p_degas := if phase==1 then 0 else sum(saturationPressures(p,T2,X,MM_vec)) + p_sat_H2O;

     if  p_degas< p then
  //    Modelica.Utilities.Streams.print("1Phase-Liquid (PartialBrine_Multi_TwoPhase_ngas.quality_pTX("+String(p)+","+String(T2)+"))");
      x:=0;
      p_H2O := p_sat_H2O;
    else
      assert(max(X[end-nX_gas:end-1])>0,"Phase equilibrium cannot be calculated without dissolved gas ("+String(p/1e5)+" bar, "+String(T2-273.15)+"°C).");
  //    Modelica.Utilities.Streams.print("2Phase (PartialBrine_Multi_TwoPhase_ngas.quality_pTX)");
  //    Modelica.Utilities.Streams.print("p="+String(p/1e5)+" bar");

      n:=X[nX_salt + 1:end] ./ MM_vec[nX_salt + 1:nX] "total mole numbers";
  //    n_g_norm:=cat(1, fill(1,nX_gas), {0.01})
  //    n_g:={0.0175164,0.00204528,0.00435, 2.51075};
  //    n_g_norm :=n_g ./ n;
      n_g_norm:=n_g_start .* sign(X[nX_salt + 1:nX]);

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
  //    Modelica.Utilities.Streams.print("k="+PowerPlant.vector2string(k)+" (PartialBrine_ngas_Newton.quality_pTX)");
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
  //        Modelica.Utilities.Streams.print("Delta_n_g_norm="+PowerPlant.vector2string(Delta_n_g_norm));

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

      end while;

    end if "p_degas< p";

  //  assert(x>=0 and (sum(n_gas_g) + n_H2O_g)>=0 or (not x>=0 and not (sum(n_gas_g) + n_H2O_g)>=0),"WEIRD!");
    d_g:= if x>0 then p/(Modelica.Constants.R*T2)*(n_g*cat(1,MM_gas,{M_H2O}))/sum(n_g) else -1;

    d_l:=if not x<1 then -1 else density_liquid_pTX(p,T2,X_l,MM_vec)
    "gases are ignored anyway";
    d:=1/(x/d_g + (1 - x)/d_l);
  //  Modelica.Utilities.Streams.print(String(z)+" (p="+String(p)+" bar)");

   state :=ThermodynamicState(
      p=p,
      T=T,
      X=X,
      X_l=X_l,
      h=x*specificEnthalpy_gas_pTX(p,T,X) + (1-x)*specificEnthalpy_liq_pTX(p,T,X),
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
  state := setState_pTX(p,temperature_phX(p,h,X,phase),X,phase) ",fluidnames)";
end setState_phX;
end PartialBrine_MultiSalt_1Phase;
