within Brine;
partial package Partial_Units "Additional units"
  type Pressure_bar = Real(final quantity="Pressure", final unit="bar")
    "pressure in bar";
  type Pressure_kPa = Real(final quantity="Pressure", final unit="kPa")
    "pressure in kPa";
  type Pressure_MPa = Real(final quantity="Pressure", final unit="MPa")
    "pressure in MPa";
  type Pressure_GPa = Real(final quantity="Pressure", final unit="GPa")
    "pressure in GPa";
  type Molality = Real (final quantity="molality", final unit=
     "mol/kg");
  type PartialMolarEnthalpy = Real (final quantity="PartialMolarEnthalpy", final unit=
     "J/mol");
  type PartialMolarHeatCapacity = Real (final quantity="PartialMolarHeatCapacity", final unit=
     "J/(mol.K)");
end Partial_Units;
