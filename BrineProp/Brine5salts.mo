within BrineProp;
package Brine5salts
  extends Brine3salts(
    final saltNames = {"sodium chloride","potassium chloride","calcium chloride","magnesium chloride","strontium chloride"},
    final saltConstants = {
      saltConstants_NaCl,
      saltConstants_KCl,
      saltConstants_CaCl2,
      saltConstants_MgCl2,
      saltConstants_SrCl2},
    final MM_salt = {M_NaCl,M_KCl,M_CaCl2,M_MgCl2,M_SrCl2},
    final nM_salt = {nM_NaCl,nM_KCl,nM_CaCl2,nM_MgCl2,nM_SrCl2},
      final iMgCl2=4,
      final iSrCl2=5);
end Brine5salts;
