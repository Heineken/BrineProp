within BrineProp;
package SaltData "Molar masses and mole numbers of the contained salts"

   constant Modelica.SIunits.MolarMass M_NaCl = 0.058443 "[kg/mol]";
   constant Integer nM_NaCl = 2 "[ion moles/mol]";
   constant Modelica.SIunits.MolarMass M_KCl = 0.074551 "[kg/mol]";
   constant Integer nM_KCl = 2 "[ion moles/mol]";
   constant Modelica.SIunits.MolarMass M_CaCl2 = 0.1109840 "[kg/mol]";
   constant Integer nM_CaCl2 = 3 "[ion moles/mol]";
   constant Modelica.SIunits.MolarMass M_MgCl2 = 0.095236 "[kg/mol]";
   constant Integer nM_MgCl2 = 3 "[ion moles/mol]";
   constant Modelica.SIunits.MolarMass M_SrCl2 = 0.158536 "[kg/mol]";
   constant Integer nM_SrCl2 = 3 "[ion moles/mol]";
//   constant Modelica.SIunits.MolarMass M_H2O = 0.018015 "[kg/mol]";

  constant Real[:] MM_salt = {
    M_NaCl,
    M_KCl,
    M_CaCl2,
    M_MgCl2,
    M_SrCl2};
//    ,    M_H2O};

  constant Real[:] nM_salt = {
    nM_NaCl,
    nM_KCl,
    nM_CaCl2,
    nM_MgCl2,
    nM_SrCl2};


  replaceable record SaltConstants
    extends Modelica.Icons.Record;
    Modelica.SIunits.MolarMass M_salt "Molar Mass in kg/mol";
    //    Integer nM_salt "Molar Mass in kg/mol";
    String name;
    annotation (Documentation(info="<html></html>"));
  end SaltConstants;
//    ,    M_H2O};

end SaltData;
