Attribute VB_Name = "Common"
' Functions used by both phases
' Functions to provide vector calculation functionality
' and more

' by Henning Francke francke@gfz-potsdam.de
' 2014 GFZ Potsdam

Option Explicit
Option Base 1

Function Vector2String(vec)
    Dim val
    For Each val In vec
        Vector2String = Vector2String & CStr(val) & ";"
    Next val
End Function

Function String2Vector(Xi_string, Optional ByRef n As Integer) 'As Double() Converts composition string to vector of doubles'
    If Len(Xi_string) = 0 Then
        String2Vector = "#Xi is empty"
        Exit Function
    End If
        
    Dim Xi_vec() As String:
    If Mid(CCur(1.1), 2, 1) = "," Then ' If Excel is set to dot as decimal separator, but system is not
        Xi_vec = Split(Replace(Xi_string, ".", ","), ";")
    Else
        Xi_vec = Split(Replace(Xi_string, ",", "."), ";")
     End If
        
    n = UBound(Xi_vec) + 1 'Split returns vector starting at 0
    If Len(Xi_vec(n - 1)) = 0 Then 'if last field is empty then crop
        n = n - 1
    End If
    Dim Xi() As Double
    ReDim Xi(1 To n)
    
    Dim i As Integer
    For i = 1 To n
        If Not IsNumeric(Xi_vec(i - 1)) Then
            String2Vector = "#Error in composition string"
            Exit Function
        End If
        Xi(i) = CDbl(Xi_vec(i - 1))
    Next i
    String2Vector = Xi
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

Function Length(vec) As Integer
    Dim vt As Integer
    vt = VarType(vec) 'http://www.java2s.com/Code/VBA-Excel-Access-Word/Data-Type/ValuesreturnedbytheVarTypefunction.htm
    If vt < 2 Then
        Length = 0
    ElseIf vt < 12 Then
        Length = 1
    ElseIf IsObject(vec) Then
        Length = vec.Count
    Else
        On Error Resume Next 'return length=0 for empty error
        Length = UBound(vec) ' gives error for empty array
    End If
End Function


Function ToDouble(vec, Optional ByRef n As Integer, Optional reduce = False) 'Typecast scalar/array to double scalar/array with index starting at 1
    If VarType(vec) = vbString Then
        ToDouble = "Type error (ToDouble)"
        Exit Function
    End If
        
    n = Length(vec)
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
        For i = 1 To n
            'If vt = 8204 Then
            '    dbl(i) = vec(i, 1)
            'Else
                dbl(i) = vec(i)
            'End If
        Next i
        ToDouble = dbl
    End If
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
    Dim X '() As Double
    X = ToDouble(Xi, nXi)
    If VarType(X) = vbString Or VarType(X) = vbError Or IsEmpty(X) Then
        FullMassVector = X
        Exit Function
    End If
    nX = nXi + 1
    ReDim Preserve X(1 To nX)
  
    X(nX) = 1 - Application.Sum(Xi)
    If X(nX) > 1 Or X(nX) <= 0 Then
        X(1) = -1
        ' MsgBox "Mass vector is wrong"
    End If
    FullMassVector = X
End Function

Function massFractionsToMolalities(X, MM) 'Calculate molalities (mole_i per kg H2O) from mass fractions X
  Dim molalities, nX As Integer, nM As Integer
  nX = Length(X)
  nM = Length(MM)
  ReDim molalities(1 To nX) 'Molalities moles/m_H2O
    If nX <> nM Then
        massFractionsToMolalities = "#Inconsistent vectors for mass fraction(" & nX & ") and molar masses(" & nM & ")"
    End If

    Dim i As Integer
    For i = 1 To nX
        If X(nX) > 0 Then
            If X(i) > 10 ^ -6 Then
                molalities(i) = X(i) / (MM(i) * X(nX)) 'numerical errors my create X[i]>0, this prevents it
           'Else
           '    molalities(i) = 0
            End If
        Else
           molalities = -1
        End If
    Next i
  massFractionsToMolalities = molalities
End Function

Function massFractionToMolality(X, X_H2O, MM) 'Calculate molalities (mole_i per kg H2O) from mass fractions X
  Dim nX As Integer: nX = Length(X)

    If X_H2O > 0 Then
        If X > 10 ^ -6 Then
            massFractionToMolality = X / (MM * X_H2O) 'numerical errors my create X[i]>0, this prevents it
       'Else
       '    molalities(i) = 0
        End If
    Else
       massFractionToMolality = -1
    End If
End Function


Function CheckMassVector(X, nX_must) As Variant
    Dim nX As Integer, msg As String
    Dim Xout, s2v As Boolean
If VarType(X) = vbString Then
        Xout = String2Vector(X, nX) 'make sure first index is 1
        s2v = True 'stupid flag to avoid having to recheck or copy Xout=X
    Else
        nX = Length(X)
        Xout = X
        s2v = False
    End If
    
    If nX = nX_must - 1 Then 'without water
        Xout = FullMassVector(IIf(s2v, Xout, X), nX) 'make sure first index is 1
        If VarType(Xout) = vbString Then
            CheckMassVector = Xout
        ElseIf VarType(Xout) = vbError Or Xout(1) = -1 Then
            CheckMassVector = "#Mass vector is wrong"
        Else
            CheckMassVector = Xout
        End If
'    ElseIf nX = nX_salt + 1 Then 'Full mass vector with water
    ElseIf nX = nX_must Then 'Full mass vector with water
        If Abs(Application.Sum(IIf(s2v, Xout, X)) - 1) > 10 ^ -8 Then
            CheckMassVector = "#Mass vector does not add up to 1"
        Else
            CheckMassVector = ToDouble(IIf(s2v, Xout, X)) 'to prevent adding a dimension
        End If
    Else
        CheckMassVector = "#Mass vector has wrong number of elements (" & nX & " instead of " & nX_must - 1 & " or " & nX_must & " )"
    End If
    
    'If Len(msg) <> 0 Then
    '    CheckMassVector = msg
    '    Exit Function
    'End If
End Function
