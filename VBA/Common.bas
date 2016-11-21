Attribute VB_Name = "Common"
' Functions used by both phases
' Functions to provide vector calculation functionality
' and more

' by Henning Francke francke@gfz-potsdam.de
' 2014 GFZ Potsdam

Option Explicit
Option Base 1

Public Function isArrayEmpty(vec As Variant) As Boolean
'Returns false if not an array or dynamic array that has not been initialised (ReDim) or has been erased (Erase)

    If IsArray(vec) = False Then isArrayEmpty = True
    
    On Error Resume Next
    If vec.Count > 0 Then isArrayEmpty = False: Exit Function 'Range
    If UBound(vec) < LBound(vec) Then isArrayEmpty = True: Exit Function Else: isArrayEmpty = False 'Array
End Function

Function Vector2String(vec)
'    Vector2String = "{" & Join(vec, ";") & "}" Works only for string arrays
    If isArrayEmpty(vec) Then
        Vector2String = "{}"
        Exit Function
    End If
        
    Dim val
    Vector2String = "{"
    For Each val In vec
        Vector2String = Vector2String & CStr(val) & ";"
    Next val
    Vector2String = Left(Vector2String, Len(Vector2String) - 1) & "}" ' crop last semicolon
End Function

Function String2Vector(Xi_string, Optional ByRef n As Integer) 'As Double() Converts composition string to vector of doubles'
    ' input as {1.1;2.2;3.3} or 1.1;2.2;3.3 or {1,1;2.2;3.3} or 1.1;2.2;3.3
    If Len(Xi_string) = 0 Or Xi_string = "{}" Then
        String2Vector = "#Stringvector is empty"
        Exit Function
    End If
        
    If Left(Xi_string, 1) = "#" Then 'ErrorMsg
        String2Vector = Xi_string
        Exit Function
    End If
    
    
    Dim Xi_str As String: Xi_str = Trim(Xi_string)
    If Left(Xi_str, 1) = "{" And Right(Xi_string, 1) = "}" Then
        Xi_str = Mid(Xi_str, 2, Len(Xi_str) - 2) ' remove curly braces
    Else
        String2Vector = "String format incorrect"
        Exit Function
    End If
    
    'Change decimal separator
    If Mid(CCur(1.1), 2, 1) = "," Then ' If Excel is set to dot as decimal separator, but system is not
        Xi_str = Replace(Xi_str, ".", ",")
    Else
        Xi_str = Replace(Xi_str, ",", ".")
    End If
     
    Dim Xi_vec() As String:
    Xi_vec = Split(Xi_str, ";")
       
    n = UBound(Xi_vec) + 1 'Split returns vector starting at 0
    If Len(Xi_vec(n - 1)) = 0 Then 'if last field is empty then crop
        n = n - 1
    End If
    Dim Xi() As Double
    ReDim Xi(1 To n)
    
    Dim i As Integer
    For i = 1 To n
        If Not IsNumeric(Xi_vec(i - 1)) Then
            String2Vector = "#Error in input string"
            Exit Function
        End If
        Xi(i) = CDbl(Xi_vec(i - 1))
    Next i
    String2Vector = Xi
End Function


Function GetValueFromJSON(JSON As String, PropertyName As String) 'Parse JSON String
    Dim equations() As String
    Dim equation
    Dim keyval() As String
    equations = Split(Mid(JSON, 2, Len(JSON) - 2), ",") ' remove curly braces
    For Each equation In equations
        keyval = Split(equation, ":")
        If Trim(keyval(0)) = Trim(PropertyName) Then
            GetValueFromJSON = keyval(1)
            GoTo Done
        End If
    Next equation
    GetValueFromJSON = "#Not found: '" & PropertyName & "'"
Done:
End Function

Function SubArray(sourceArray, Indexfrom As Integer, IndexTo As Integer)
    Dim b() As Double, i As Integer
    Dim n As Integer
    
    n = Length(sourceArray)
    If Indexfrom < 1 Or IndexTo < 1 Or Indexfrom > n Or IndexTo > n Then
        SubArray = "#Index out of valid range for input array in SubArray"
        Exit Function
    End If
    ReDim b(1 To IndexTo - Indexfrom + 1)
    For i = 0 To IndexTo - Indexfrom
        b(1 + i) = sourceArray(Indexfrom + i)
    Next i
    SubArray = b
End Function


Private Sub TestMulVecElwise()
    Dim a(3) As Double
    a(1) = 1
    a(2) = 2
    a(3) = 3
    Dim b() As Double: b = MulVecElwise(a, a)
End Sub

Function Length(vec, Optional ByRef offset As Integer) As Integer
    Dim vt As Integer
    offset = 0
    vt = VarType(vec) 'http://www.java2s.com/Code/VBA-Excel-Access-Word/Data-Type/ValuesreturnedbytheVarTypefunction.htm
    If vt < 2 Then
        Length = 0
    ElseIf vt < 12 Then
        Length = 1
    ElseIf IsObject(vec) Then
        Length = vec.Count
    Else
        On Error Resume Next 'return length=0 for empty error
        offset = LBound(vec) - 1
        Length = UBound(vec) - LBound(vec) + 1 ' gives error for empty array
    End If
End Function


Function ToDouble(vec, Optional ByRef n As Integer, Optional reduce = False) 'Typecast scalar/array to double scalar/array with index starting at 1
    If VarType(vec) = vbString Then
        ToDouble = "Type error (ToDouble)"
        Exit Function
    End If
        
    Dim offset As Integer
    n = Length(vec, offset)
    Dim vt As Integer
    vt = VarType(vec)
    If vt < 12 Then ' if scalar
        ToDouble = vec
    ElseIf n = 1 And reduce Then
        ToDouble = vec(1) ' reduce 1-element-array to scalar
    Else
        Dim dbl() As Double
        If n = 0 Then
            Exit Function
        End If
        ReDim dbl(1 To n)
        Dim i As Integer
 '       Dim offset As Integer: offset = LBound(vec) - 1

'Is passed array 1D or 2D (when given anot as String or Range, but as {1,2,3})
On Error GoTo TwoD: 'coz there is no function to query the number of dimensions in VBA (https://stackoverflow.com/questions/6901991/how-to-return-the-number-of-dimensions-of-a-variant-variable-passed-to-it-in-v)
        For i = 1 To n
            dbl(i) = vec(i + offset) ' Try 1D Array
        Next i
        ToDouble = dbl
    End If
    Exit Function
    
TwoD:
    For i = 1 To n
        dbl(i) = vec(i + offset, 1) ' 2D Array
    Next i
    ToDouble = dbl
End Function

Function VecAbs(vec) 'As Double()
    Dim c() As Double, n As Integer
    n = Length(vec)
    ReDim c(1 To n)
    
    Dim i As Integer
    For i = 1 To n
        c(i) = Abs(vec(i))
    Next i
    VecAbs = c
End Function

Function VecSgn(vec) 'As Double()
    Dim c() As Double, n As Integer
    n = Length(vec)
    ReDim c(1 To n)
    
    Dim i As Integer
    For i = 1 To n
        c(i) = Math.Sgn(vec(i))
    Next i
    VecSgn = c
End Function

Function VecSum(a, b) 'As Double()
    VecSum = VecOp(a, b, "add")
End Function

Function VecProd(a, b) 'As Double()
    VecProd = VecOp(a, b, "multiply")
End Function

Function VecDiv(a, b) 'As Double()
    VecDiv = VecOp(a, b, "divide")
End Function

Function VecDiff(a, b) 'As Double()
    VecDiff = VecOp(a, b, "substract")
End Function

Function ScalProd(a, b) As Double 'scalar product of two vectors
    ScalProd = VecOp(a, b, "scalarProduct")
End Function

Function VecOp(A_, B_, what) 'As Double()
    Dim i As Integer, n_a As Integer, n_b As Integer
    'n_a = Length(a)
    'n_b = Length(b)
    
    Dim a, b
    a = ToDouble(A_, n_a, True)
    b = ToDouble(B_, n_b, True)
    
    If VarType(a) = vbString Then
        VecOp = a
        Exit Function
    End If
    If VarType(b) = vbString Then
        VecOp = b
        Exit Function
    End If

    If what = "substract" Then
        If n_a <> n_b And n_b <> 1 Then
                VecOp = "#2nd Factor must be scalar or both factors must be vectors of same length for division."
                Exit Function
        End If
    Else
        If n_a <> n_b And n_a <> 1 And n_b <> 1 Then
                VecOp = "#Factors must be scalar or vectors of same length for multiplication."
                Exit Function
        End If
    End If
    
    Dim n As Integer
    n = Application.Max(n_a, n_b)
    
    If what = "scalarProduct" Then
        For i = 1 To n
            VecOp = VecOp + a(i) * b(i)
        Next i
    Else
        Dim c() As Double
        ReDim c(1 To n)
        
        Select Case what
        Case "multiply"
            If n_a = 1 Then
                For i = 1 To n
                    c(i) = a * b(i)
                Next i
            ElseIf n_b = 1 Then
                For i = 1 To n
                    c(i) = a(i) * b
                Next i
            Else
                For i = 1 To n
                    c(i) = a(i) * b(i)
                Next i
            End If
        Case "divide"
            If n_b = 1 Then
                For i = 1 To n_a ' divide all elements by scalar
                    c(i) = a(i) / b
                Next i
            Else
                For i = 1 To n_a
                    If b(i) > 0 Then
                        c(i) = a(i) / b(i) ' divide elementwise
                    Else
                        VecOp = "Division by zero in VecOp"
                        Exit Function
                    End If
                Next i
            End If
        Case "add"
            If n_a = 1 Then
                For i = 1 To n
                    c(i) = a + b(i)
                Next i
            ElseIf n_b = 1 Then
                For i = 1 To n
                    c(i) = a(i) + b
                Next i
            Else
                For i = 1 To n
                    c(i) = a(i) + b(i)
                Next i
            End If
        Case "substract"
            If n_a = 1 Then
                For i = 1 To n
                    c(i) = a - b(i)
                Next i
            ElseIf n_b = 1 Then
                For i = 1 To n
                    c(i) = a(i) - b
                Next i
            Else
                For i = 1 To n
                    c(i) = a(i) - b(i)
                Next i
            End If
        Case Else
                VecOp = "#Don't know what to do (VecOp)"
        End Select
        VecOp = c
    End If
End Function


Function cat(a, b)
    If VarType(a) = vbString Then
        cat = a
        Exit Function
    End If
    If VarType(b) = vbString Then
        cat = b
        Exit Function
    End If

    Dim n_a As Integer, n_b As Integer
    n_a = Length(a)
    n_b = Length(b)
    'b = ToDouble(b)
    Dim c() As Double, i As Integer
    c = ToDouble(a)
    ReDim Preserve c(1 To n_a + n_b)
    
    For i = 1 To n_b
        c(n_a + i) = b(i)
    Next i
    cat = c
End Function

Function fill(val, n) As Double()
    Dim vec() As Double
    ReDim vec(1 To n)
    Dim i As Integer
    For i = 1 To n
        vec(i) = val
    Next i
    fill = vec
End Function

Function FullMassVector(Xi, Optional ByRef nX As Integer) 'As Double()
    Dim nXi As Integer
    Dim x '() As Double
    x = ToDouble(Xi, nXi)
    If VarType(x) = vbString Or VarType(x) = vbError Or IsEmpty(x) Then
        FullMassVector = x
        Exit Function
    End If
    nX = nXi + 1
    ReDim Preserve x(1 To nX)
  
    x(nX) = 1 - Application.Sum(Xi)
    If x(nX) > 1 Or x(nX) < 0 Then 'removed X(nX) <= 0 to allow for pure gases
        'X(1) = -1
        FullMassVector = "Mass vector is wrong"
        Exit Function
    End If
    FullMassVector = x
End Function

Function massFractionsToMolalities(x, MM) 'Calculate molalities (mole_i per kg H2O) from mass fractions X
  Dim molalities, nX As Integer, nM As Integer
  nX = Length(x)
  nM = Length(MM)
  ReDim molalities(1 To nX) 'Molalities moles/m_H2O
    If nX <> nM Then
        massFractionsToMolalities = "#Inconsistent vectors for mass fraction(" & nX & ") and molar masses(" & nM & ")"
    End If

    Dim i As Integer
    For i = 1 To nX
        If x(nX) > 0 Then
            If x(i) > 10 ^ -6 Then 'to prevent division by zero
                molalities(i) = x(i) / (MM(i) * x(nX)) 'numerical errors may create X[i]>0 for non-present salts, this prevents it
           'Else
           '    molalities(i) = 0
            End If
        Else
           molalities(i) = -1 ' High number that will exceed any molality limit
        End If
    Next i
  massFractionsToMolalities = molalities
End Function

Function massFractionToMolality(x, X_H2O, MM) 'Calculate molalities (mole_i per kg H2O) from mass fractions X
'used in worksheet
  Dim nX As Integer: nX = Length(x)

    If X_H2O > 0 Then
        If x > 10 ^ -6 Then
            massFractionToMolality = x / (MM * X_H2O) 'numerical errors my create X[i]>0, this prevents it
       'Else
       '    molalities(i) = 0
        End If
    Else
       massFractionToMolality = -1
    End If
End Function


Function CheckMassVector(x, nX_must) As Variant
    Dim nX As Integer, msg As String
    Dim Xout, s2v As Boolean
If VarType(x) = vbString Then
        Xout = String2Vector(x, nX) 'make sure first index is 1
        If VarType(Xout) = vbString Then
            CheckMassVector = Xout
            Exit Function
        End If
        s2v = True 'stupid flag to avoid having to recheck or copy Xout=X
    Else
        nX = Length(x)
        ' Xout = X Doesn't work
        s2v = False
    End If
    
    If nX = nX_must - 1 Then 'without water
        Xout = FullMassVector(IIf(s2v, Xout, x), nX) 'make sure first index is 1
        'If VarType(Xout) = vbString Then
        '    CheckMassVector = Xout
        'Else
        If VarType(Xout) = vbError Then ' Or IIf(s2v, Xout(1), X(1)) = -1 'when does that happen?
            CheckMassVector = "#Mass vector is wrong"
        Else 'also if error string was returned
            CheckMassVector = Xout
        End If
'    ElseIf nX = nX_salt + 1 Then 'Full mass vector with water
    ElseIf nX = nX_must Then 'Full mass vector with water
        If Abs(Application.Sum(IIf(s2v, Xout, x)) - 1) > 10 ^ -6 Then
            CheckMassVector = "#Mass vector does not add up to 1"
        Else
            CheckMassVector = ToDouble(IIf(s2v, Xout, x)) 'to prevent adding a dimension
        End If
    Else
        CheckMassVector = "#Mass vector has wrong number of elements (" & nX & " instead of " & nX_must - 1 & " or " & nX_must & " )"
    End If
    
    'If Len(msg) <> 0 Then
    '    CheckMassVector = msg
    '    Exit Function
    'End If
End Function
