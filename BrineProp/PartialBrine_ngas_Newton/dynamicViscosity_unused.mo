within BrineProp.PartialBrine_ngas_Newton;
function dynamicViscosity_unused
//redeclare function extends dynamicViscosity_unused
algorithm
  eta :=dynamicViscosity_pTX_unused(
      state.p,
      state.T,
      state.X);
end dynamicViscosity_unused;
