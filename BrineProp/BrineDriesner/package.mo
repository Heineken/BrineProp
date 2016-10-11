within BrineProp;
package BrineDriesner "NaCl solution using Driesner density and enthalpy function"
  extends BrineProp.PartialBrineMultiSaltOnePhase(
      redeclare package Salt_data = BrineProp.SaltData,
      final saltNames = {"sodium chloride"},
      saltConstants = {Salt_data.saltConstants_NaCl},
      final MM_salt = {Salt_data.M_NaCl},
      final nM_salt = {Salt_data.nM_NaCl},
      final iNaCl=1);

  redeclare function extends density_pTX
  algorithm
      d:=Densities.density_Driesner2007_pTX(p,T,X);
  end density_pTX;


  redeclare function extends specificEnthalpy_pTX
  "enthalpy calculation according to Driesner 2007 et al: 0-1000degC; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  //Pressure limited to 100 MPa by Modelica Water property function
  /*  input SI.Pressure p;
  input SI.Temp_K T;
  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output SI.SpecificEnthalpy h;*/

  algorithm
    h := specificEnthalpy_pTX_Driesner(p,T,X[1]);

  //  print("Brine_Driesner.specificEnthalpy_pTX: "+String(p*1e-5)+"bar."+String(T_Scale_h)+"degC->"+String(h)+" J/kg");
  end specificEnthalpy_pTX;
end BrineDriesner;
