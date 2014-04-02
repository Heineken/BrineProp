within BrineProp.PartialBrine_ngas_Newton;
function setState_pTX_cEQ
  "finds the VLE iteratively by varying the normalized quantity of gas in the gasphase, calculates the densities"
//redeclare replaceable partial function extends setState_pTX
input Real[nX_gas + 1] n_g_norm_start "=fill(.1,nX_gas+1) 
    start value, all gas in gas phase, all water liquid, set in BaseProps";
/*
//output SI.Density d_g= if x>0 then (n_CO2_g*d_g_CO2 + n_N2_g*d_g_N2)/(n_CO2_g + n_H2O_g) else -1;
//output Real[nX_gas + 1] n_g_norm;
//output Real k[nX_gas];
// SI.Density d_g_H2O = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.rhov_p(p) "density of water vapor";
*/
  output Real GVF;
  output SI.SpecificEnthalpy h_l;
  output SI.SpecificEnthalpy h_g;
  output SI.Pressure[nX_gas + 1] p_gas "=fill(0,nX_gas)";
  output SI.Pressure p_H2O "water vapour pressure TODO is in p_gas drin";
  output SI.Pressure[nX_gas + 1] p_degas;
  output Integer z "number of iterations";
//  output
//  Integer z=0;
protected
  SI.MassFraction[nX_gas+1] X_g;
  SI.MassFraction[nX] X_l=X "start value";
  SI.Density d;
  SI.Density d_l;
  SI.Density d_g;
  SI.MassFraction x;
  SI.Pressure p_sat_H2O "water vapour pressure considering salinity";
  SI.Pressure p_H2O_0 "pure water vapour pressure";
  SI.Pressure[nX_gas + 1] f
    "componentwise pressure disbalance (to become zero)";
  SI.MassFraction[nX_gas] X_l_sat "start value";
  SI.Pressure[nX_gas + 1] p_sat;
  SI.MassFraction[nX_gas + 1] Delta_n_g_norm = fill(1e3,nX_gas+1)
    "large initial value to enter while loop";
//  SI.MassFraction[nX_gas + 1] c = {3.16407e-5,0,3.6e-8,.746547} "cat(1,fill(1e-4, nX_gas), {X[end]})fill(0, nX_gas+1)X[nX_salt+1:end]";
  Real k_H2O "Henry coefficient";
  Real k[nX_gas];
  Real[nX_gas + 1] n "Total mol numbers";
  Real[nX_gas + 1] n_l "mols in liquid phase per kg fluid";
  Real[nX_gas + 1] n_g "mols in   gas  phase per kg fluid";
  Real[nX_gas + 1] n_g_norm_test;
//  SI.MassFraction[nX] X;
  Real[nX_gas + 1] n_g_norm
    "= X[end-nX_gas:end-1]./MM_gas fill(0,nX_gas) - start value: all degassed";
  Real dXl_dng_norm;
  Real dXl_sat_dng_norm;
  Real[nX_gas + 1] dfdn_g_norm;
  Real sum_n_ion;
  constant Integer zmax=1000 "maximum number of iteration";
//  Integer ju = nX_gas+1;
  Real[nX_gas + 1,nX_gas + 1] Grad_f;
  Real DeltaC=0.001;
  SI.Temperature T2;
  SpecificHeatCapacity R_gas;
  Real y_l_H2O;
  Real y_l_H2O_sat;
  Boolean debugmode = true;
algorithm
  if debugmode then
      print("Running setState_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
  end if;

 assert(p>0,"p="+String(p/1e5)+"bar - Negative pressure is not yet supported ;-) (PartialBrine_ngas_Newton.setState_pTX())");
 assert(max(X)-1<=1e-8 and min(X)>=-1e-8, "X out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X]))+" (setState_pTX())");
//  assert(T>273.15,"T too low in PartialBrine_ngas_Newton.()");

  if T<273.15 then
    print("T="+String(T)+" too low (<0°C), setting to 0°C in PartialBrine_ngas_Newton.setState_pTX()");
  end if;
  T2:= max(273.16,T);

  p_H2O := saturationPressure_H2O(p,T2,X,MM_vec,nM_vec);
  p_degas := cat(1,saturationPressures(p,T2,X,MM_vec), {p_H2O});

   if phase==1 or sum(p_degas) < p then
    if debugmode then
      print("1Phase-Liquid (PartialBrine_Multi_TwoPhase_ngas.setState_pTX("+String(p)+","+String(T2)+"))");
    end if;
    x:=0;
    p_H2O := p_sat_H2O;
  else
    assert(max(X[end-nX_gas:end-1])>0,"Phase equilibrium cannot be calculated without dissolved gas at "+String(p/1e5)+" bar, "+String(T2-273.15)+"°C with p_degas="+String(sum(p_degas)/1e5)+" bar.");

    n:=X[nX_salt + 1:end] ./ MM_vec[nX_salt + 1:nX]
      "total mole numbers per kg brine";
    n_g_norm:=n_g_norm_start .* sign(X[nX_salt + 1:nX])
      "switch off unused salts";

    while z<1 or max(abs(Delta_n_g_norm))>1e-3 loop
    //stop iteration when p-equlibrium is found or gas fraction is very low
      z:=z + 1 "count iterations";
      assert(z<=zmax,"Reached maximum number of iterations ("+String(z)+"/"+String(zmax)+") for solution equilibrium calculation. (setState_pTX("+String(p/1e5)+"bar,"+String(T2-273.16)+"°C))\nDeltaP="+String(max(abs(p_sat-p_gas))));

//     print("\nn_g_norm=" + PowerPlant.vector2string(n_g_norm));
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
        "X_l ändert sich";
    if (p_H2O>p) then
        print("p_H2O(" + String(p/1e5) + "bar," +
          String(T2 - 273.15) + "°C, " + Modelica.Math.Matrices.toString(transpose([X])) + ") = "
           + String(p_H2O/1e5) + "bar>p ! (PartialBrine_ngas_Newton.setState_pTX)");
      x:=1;
      break;
    end if;

//       print("p_H2O_0=" + String(p_H2O_0));

//        k:=solubilities_pTX(p=p, T=T2, X_l=X_l, X=X, p_gas=fill(p/3,3)) ./ fill(p/3,3);
//  print("X_l="+PowerPlant.vector2string(X_l[nX_salt+1:end]));
        X_l_sat:=solubilities_pTX(
        p=p,
        T=T2,
        X_l=X_l,
        X=X,
        p_gas=p_gas[1:nX_gas]);
        k:=X_l_sat ./ p_gas[1:nX_gas];
//    print("k="+PowerPlant.vector2string(k)+" (PartialBrine_ngas_Newton.setState_pTX)");p,T,X_l,MM_vec,p_gas[1])
/*        for i in 1:nX_gas loop
          p_sat[i] := X_l[nX_salt+i]/ (if k[i]>0 then k[i] else 1e10) 
          "Degassing pressure";
        end for;
        p_sat[nX_gas+1] := p_H2O;

        f :=  p_gas-p_sat;*/

        y_l_H2O:=n_l[end]/sum(n_l);
        y_l_H2O_sat:=p_gas[end]/p_H2O;
        f[1:nX_gas] :=  X_l[nX_salt+1:end-1]-X_l_sat;
        f[end] :=y_l_H2O-y_l_H2O_sat;

//       print("p_gas=" + PowerPlant.vector2string(p_gas) + "=>" + String(sum(p_gas)));
//       print("p_sat=" + PowerPlant.vector2string(p_sat));

       sum_n_ion :=cat(1,X[1:nX_salt] ./ MM_vec[1:nX_salt],n_l)*nM_vec;

  //GRADIENT analytisch df[gamma]/dc[gamma]

        for gamma in 1:nX_gas+1 loop
            dXl_dng_norm:=-n[gamma]*MM_vec[nX_salt+gamma] / (1- n_g*MM_vec[nX_salt+1:nX])
          "partial pressure";
            if gamma == nX_gas+1 then
              dXl_sat_dng_norm := n[end]/ sum_n_ion;
            else
                dXl_sat_dng_norm := p * n[gamma]/sum(n_g) *k[gamma];
            end if;
            dfdn_g_norm[gamma] := dXl_dng_norm-dXl_sat_dng_norm;
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
            n_g_norm[alpha] := max(1e-9,min(1,n_g_norm[alpha] + 0.01*Delta_n_g_norm[alpha]))
          "new concentration limited by all dissolved/none dissolved, 1e-9 to avoid k=NaN";
//          end if;
        end for;
//       print("p_sat="+String(p_sat[1])+", solu="+String(solubility_CO2_pTX_Duan2006(p,T2,X_l,MM_vec,p_gas[1]))+", p_gas="+String(p_gas[1]));
         print("p="+String(p)+",T2="+String(T2)+",p_gas="+Modelica.Math.Matrices.toString(transpose([p_gas])));
/*        print("X_l="+Modelica.Math.Matrices.toString({X_l}));
        print("MM_vec="+Modelica.Math.Matrices.toString({MM_vec}));
*/
    end while;

  end if "p_degas< p";

//DENSITY
 X_g:=if x>0 then (X[end-nX_gas:end]-X_l[end-nX_gas:end]*(1-x))/x else fill(0,nX_gas+1);
/*Calculation here  R_gas :=if x > 0 then sum(Modelica.Constants.R*X_g ./ cat(1,MM_gas,{M_H2O})) else -1;
  d_g :=if x > 0 then p/(T2*R_gas) else -1;*/
//  d_g:= if x>0 then p/(Modelica.Constants.R*T2)*(n_g*cat(1,MM_gas,{M_H2O}))/sum(n_g) else -1;
  d_g :=if x > 0 then BrineGas_3Gas.density_pTX(p,T, X_g) else -1
    "calculation in MoistAirModel";

  d_l:=if not x<1 then -1 else density_liquid_pTX(p,T2,X_l,MM_vec)
    "no 1-phase gas";
  d:=1/(x/d_g + (1 - x)/d_l);
//  print(String(z)+" (p="+String(p_gas[1])+" bar)");

// X_g:=if x>0 then (X-X_l*(1-x))/x else fill(0,nX);
 h_l:=specificEnthalpy_liq_pTX(p,T,X_l);
 h_g:=specificEnthalpy_gas_pTX(p,T,X_g);
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
    d_g=d_g,
    phase=if x>0 and x<1 then 2 else 1) "phase_out";
/*    h_g=h_g,
    h_l=h_l,
    p_H2O=p_H2O,
    p_gas=p_gas[1:nX_gas],
    p_degas=p_degas

*/

  annotation (Diagram(graphics={Text(
          extent={{-96,16},{98,-16}},
          lineColor={0,0,255},
          textStyle={TextStyle.Bold},
          textString="find static VLE")}));
end setState_pTX_cEQ;
