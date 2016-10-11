within BrineProp;
package BrineDuan "NaCl solution using Duan density"
  extends BrineProp.PartialBrineMultiSaltOnePhase(
      redeclare package Salt_data = BrineProp.SaltDataDuan,
      final saltNames = {"sodium chloride"},
      saltConstants = {Salt_data.saltConstants_NaCl},
      final MM_salt = {Salt_data.M_NaCl},
      final nM_salt = {Salt_data.nM_NaCl},
      final iNaCl=1);



  redeclare function extends specificEnthalpy_pTX
  "enthalpy calculation according to Driesner et al: 0-1000degC; 0.1-500MPa"
  algorithm
    h :=BrineDriesner.specificEnthalpy_pTX(p,T,X);
  end specificEnthalpy_pTX;


  redeclare function extends density_pTX
  algorithm
    d:=density_Duan2008_pTX(p,T,X,MM_vec,saltConstants);
  end density_pTX;


  redeclare function extends dynamicViscosity_pTX
  algorithm
  //  eta :=dynamicViscosity_Duan2009NaCl_pTX(p,T,X);
    eta :=dynamicViscosity_Duan2009NaClKCl_pTX(p,T,X,MM_vec,saltConstants);
  end dynamicViscosity_pTX;

  annotation (Documentation(info="<html>
<p>Implementation of property functions (h,rho,eta) for NaCl solution by Duan.</p>
<p>Based on multi-salt template.</p>
</html>"));
end BrineDuan;
