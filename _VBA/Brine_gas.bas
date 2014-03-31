Attribute VB_Name = "Brine_gas"
' Gas phase properties
' (density, viscosity, specific heat capacity, specific enthalpy)
Option Explicit
Option Base 1

Public Const nX_gas = 3
Public Const nX = nX_gas + 1
Public Const ignoreLimitN2_T = False

Const r = 8.314472
'Const M_CO2 = 0.0440095 '[kg/mol]
Const nM_CO2 = 1 'number of ions per molecule
'Const M_N2 = 0.0280134 '[kg/mol]
Const nM_N2 = 1 'number of ions per molecule
'Const M_CH4 = 0.01604246 '[kg/mol]
Const nM_CH4 = 1 'number of ions per molecule
'Function gas_vec()
'    gas_vec = Array(M_CO2, M_N2, M_CH4, M_H2O)
'End Function

'Function MM_vec()
'' generates variant vector of molar masses
'    MM_vec = Array(M_CO2, M_N2, M_CH4, M_H2O) '[first index=0 !]
'    Dim v(1 To 4) As Variant
'    v(1) = M_CO2
'    v(2) = M_N2
'    v(3) = M_CH4
'    v(4) = M_H2O
'End Function

Function MM_vec()
 'generates double vector of molar masses
    Dim V(1 To nX_gas + 1) As Double
    V(1) = M_CO2
    V(2) = M_N2
    V(3) = M_CH4
    V(4) = M_H2O
    MM_vec = V
End Function

Function nM_vec() As Double()
 'generates double vector of molar masses
    Dim V(1 To nX_gas + 1) As Double
    V(1) = CO2.nM
    V(2) = N2.nM
    V(3) = CH4.nM
    V(4) = H2O.nM
    nM_vec = V
End Function

Function R_gas(Xi)
    Dim X
    X = CheckMassVector(Xi, nX)
    If VarType(X) = vbString Then
        R_gas = X
        Exit Function
   End If
    R_gas = ScalProd(X, Array(CO2.r, N2.r, CH4.r, H2O.r))
    ' R_gas = Application.SumProduct(X, Array(CO2.R, N2.R, CH4.R, H2O.R)) 'fails for arrays from workbook
End Function

'BELOW GENERIC PART

Private Function MM_gas(X) As Double
    MM_gas = CheckMassVector(X)
    If VarType(MM_gas) <> vbBoolean Then
        Exit Function
    End If
    MM_gas = ScalProd(X, MM_vec)
End Function

Function density(p As Double, T As Double, Xin)
    'Density of an ideal mixture of ideal gases
    If DebugMode Then
        Debug.Print "Running BrineGas.Density(" & p / 10 ^ 5 & " bar," & T - 273.15 & " °C"
    End If
    
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        density = X & " (Brine_gas.density)"
        Exit Function
    End If
    
    density = p / (T * R_gas(X))
End Function

Private Function SingleGasNasa_h_T(data As DataRecord, T As Double, Optional exclEnthForm As Boolean = True, Optional ZeroAt0K As Boolean = True, Optional h_off As Double = 0) As Double
    SingleGasNasa_h_T = _
    IIf(T < data.Tlimit, data.r * ((-data.alow(1) + T * (data.blow(1) + data.alow(2) * Math.Log(T) + T * (1# * data.alow(3) + T * (0.5 * data.alow(4) + T * (1 / 3 * data.alow(5) + T * (0.25 * data.alow(6) + 0.2 * data.alow(7) * T)))))) / T) _
        , data.r * ((-data.ahigh(1) + T * (data.bhigh(1) + data.ahigh(2) * Math.Log(T) + T * (1# * data.ahigh(3) + T * (0.5 * data.ahigh(4) + T * (1 / 3 * data.ahigh(5) + T * (0.25 * data.ahigh(6) + 0.2 * data.ahigh(7) * T)))))) / T)) _
    + IIf(exclEnthForm, -data.Hf, 0#) + IIf(ZeroAt0K, data.h0, 0#) + h_off
End Function

Function specificEnthalpy(p As Double, T As Double, Xin)
    'calculation of specific enthalpy of gas mixture
    'specificEnthalpy = checkMassVector(X)
    'If VarType(specificEnthalpy) <> vbBoolean Then
    '    Exit Function
    'End If
    
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        specificEnthalpy = X & " (Brine_gas.specificEnthalpy)"
        Exit Function
    End If
    
    Dim h_H2O_sat As Double, h_H2O As Double, h_CO2 As Double, h_N2 As Double, h_CH4 As Double
    h_H2O_sat = Waterh_satv_p(p) 'Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.hv_p(p)
    Dim h_vec()
    h_vec() = Array( _
    SingleGasNasa_h_T(CO2, T), _
    SingleGasNasa_h_T(N2, T), _
    SingleGasNasa_h_T(CH4, T), _
    Application.Max(h_H2O_sat, SpecificEnthalpy_pT(p, T)) _
    ) 'to make sure it is gaseous TODO:Take regions directly

    If DebugMode Then
      Debug.Print "Running BrineGas.SpecificEnthalpy(" & p / 10 ^ 5 & " bar," & T - 273.15 & " °C, X=" & Vector2String(X) & ")"
    End If

    'If Not Application.Min(X) > 0 Then
    '  Debug.Print "No gas composition, assuming water vapour.(Brine_gas.SpecificHeatCapacity_pTX)"
    'End If

    specificEnthalpy = ScalProd(h_vec, X) 'mass weighted average
    
End Function
    
Private Function SingleGasNasa_cp_T(data As DataRecord, T As Double, Optional exclEnthForm As Boolean = True, Optional ZeroAt0K As Boolean = True, Optional h_off As Double = 0) As Double
    SingleGasNasa_cp_T = _
    IIf(T < data.Tlimit, _
        data.r * (1 / (T * T) * (data.alow(1) + T * (data.alow(2) + T * (1# * data.alow(3) + T * (data.alow(4) + T * (data.alow(5) + T * (data.alow(6) + data.alow(7) * T))))))) _
        , data.r * (1 / (T * T) * (data.ahigh(1) + T * (data.ahigh(2) + T * (1# * data.ahigh(3) + T * (data.ahigh(4) + T * (data.ahigh(5) + T * (data.ahigh(6) + data.ahigh(7) * T))))))) _
    )
End Function

Function specificHeatCapacityCp(p As Double, T As Double, Xin) 'calculation of specific enthalpy of gas mixture
    'Argument X is either X or XI (mass vector with or without water)
    If DebugMode Then
      Debug.Print "Running specificHeatCapacityCp_gas(" & p / 10 ^ 5 & " bar," & T - 273.15 & " °C)"
    End If
    
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        specificHeatCapacityCp = X & " (Brine_gas.specificHeatCapacityCp)"
        Exit Function
    End If
    
    Dim cp_H2O_sat As Double, cp_H2O As Double, cp_CO2 As Double, cp_N2 As Double, cp_CH4 As Double
    cp_CO2 = SingleGasNasa_cp_T(CO2, T)
    cp_N2 = SingleGasNasa_cp_T(N2, T)
    cp_CH4 = SingleGasNasa_cp_T(CH4, T)
    'cp_H2O = SingleGasNasa_cp_T(H2O, T)
    cp_H2O = IAPWS.SpecificHeatCapacityCp_pT(Application.Min(p, IAPWS.Waterpsat_T(T) - 1), T)
    Dim cp_vec()
    cp_vec() = Array(cp_CO2, cp_N2, cp_CH4, cp_H2O)

    If DebugMode Then
        Debug.Print "cp_CO2: " & cp_CO2
        Debug.Print "cp_N2: " & cp_N2
        Debug.Print "cp_CH4: " & cp_CH4
        Debug.Print "cp_H2O: " & cp_H2O
    End If

    specificHeatCapacityCp = ScalProd(cp_vec, X) 'mass weighted average
End Function

Function dynamicViscosity(p As Double, T As Double, Optional X)
    dynamicViscosity = GasData.MoistAirDynamicViscosity(T)
End Function
