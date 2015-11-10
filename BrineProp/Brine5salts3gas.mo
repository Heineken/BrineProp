within BrineProp;
package Brine5salts3gas
  "Two-phase aqueous solution of NaCl, KCl, CaCl2, N2, CO2, CH4"
  extends Brine3salts3gas(
    redeclare package Salt_data = BrineProp.SaltDataDuan,
    final saltConstants = {
      saltConstants_NaCl,
      saltConstants_KCl,
      saltConstants_CaCl2,
      saltConstants_MgCl2,
      saltConstants_SrCl2},
    final iNaCl=1,
    final iKCl=2,
    final iCaCl2=3,
    final iMgCl2=4,
    final iSrCl2=5,
    final saltNames = {"sodium chloride","potassium chloride","calcium chloride","magnesium chloride","strontium chloride"},
    final MM_salt = {M_NaCl,M_KCl,M_CaCl2,M_MgCl2,M_SrCl2},
    final nM_salt = {nM_NaCl,nM_KCl,nM_CaCl2,nM_MgCl2,nM_SrCl2});
end Brine5salts3gas;
