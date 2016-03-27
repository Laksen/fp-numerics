unit ops;

interface

{$mode objfpc}
{$modeswitch advancedrecords}

uses
  sysutils, base, math, matrix;

type
  Decompose = record
    type
      QR = record
        class procedure GramSchmidt(const AMatrix: TMatrix; out AQ, AR: TMatrix); static;
        class procedure Householder(const AMatrix: TMatrix; out AQ, AR: TMatrix); static;
      end;

    class procedure LU(const AMatrix: TMatrix; out AL, AU: TMatrix); static;
    class function Cholesky(const AMatrix: TMatrix): TMatrix; static;
  end;

  Solve = record
    class function LU(const AL, AU, AB: TMatrix): TMatrix; static;
  end;

  Invert = record
    class function GaussJordan(const AMatrix: TMatrix): TMatrix; static;
    class function LU(const AL, AU: TMatrix): TMatrix; static;
  end;

  Determinant = record
    class function Gen(const AMatrix: TMatrix): double; static;
    class function Tridiagonal(const AMatrix: TMatrix): double; static;
    class function Triangular(const AMatrix: TMatrix): double; static;
    class function Unitary(const AMatrix: TMatrix): double; static;
  end;

  Helpers = record helper for TMatrix
  strict private
    function GetDeterminant: double;
    function GetInverse: TMatrix;
  public
    property Inverse: TMatrix read GetInverse;
    property Determinant: double read GetDeterminant;
  end;

operator +(const a,b: TVector): TVector;
operator -(const a,b: TVector): TVector;
operator *(const a,b: TVector): TVector;
operator /(const a,b: TVector): TVector;

operator +(const a: TVector; const b: double): TVector;
operator +(const a: double; const b: TVector): TVector;
operator -(const a: TVector; const b: double): TVector;
operator -(const a: double; const b: TVector): TVector;
operator *(const a: TVector; const b: double): TVector;
operator *(const a: double; const b: TVector): TVector;
operator /(const a: TVector; const b: double): TVector;
operator /(const a: double; const b: TVector): TVector;

operator **(const a, b: double): double;
operator **(const a: TVector; const b: double): TVector;

operator :=(const a: TVector): ansistring;

implementation

operator +(const a,b: TVector): TVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]+b[i];
end;

operator -(const a,b: TVector): TVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]-b[i];
end;

operator *(const a,b: TVector): TVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]*b[i];
end;

operator /(const a,b: TVector): TVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]/b[i];
end;

operator +(const a: TVector; const b: double): TVector;
var i: longint;
begin
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]+b;
end;

operator +(const a: double; const b: TVector): TVector;
var i: longint;
begin
   setlength(result, length(b));
   
   for i := 0 to high(b) do
      result[i] := a+b[i];
end;

operator -(const a: TVector; const b: double): TVector;
var i: longint;
begin
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]-b;
end;

operator -(const a: double; const b: TVector): TVector;
var i: longint;
begin
   setlength(result, length(b));
   
   for i := 0 to high(b) do
      result[i] := a-b[i];
end;

operator *(const a: TVector; const b: double): TVector;
var i: longint;
begin
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]*b;
end;

operator *(const a: double; const b: TVector): TVector;
var i: longint;
begin
   setlength(result, length(b));
   
   for i := 0 to high(b) do
      result[i] := a*b[i];
end;

operator /(const a: TVector; const b: double): TVector;
var i: longint;
begin
   assert(b<>0, 'Cannot divide by zero');
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := a[i]/b;
end;

operator /(const a: double; const b: TVector): TVector;
var i: longint;
begin
   setlength(result, length(b));

   for i := 0 to high(b) do
   begin
      if b[i] = 0 then
         result[i] := 0
      else
         result[i] := a/b[i];
   end;
end;

operator **(const a, b: double): double;
begin
   result := power(a,b);
end;

operator **(const a: TVector; const b: double): TVector;
var i: longint;
begin
   setlength(result, length(a));
   
   for i := 0 to high(a) do
      result[i] := power(a[i],b);
end;

operator :=(const a: TVector): ansistring;
var i: longint;
begin
   result := '[';

   for i := 0 to high(a) do
   begin
      if i > 0 then
         result := result + format(', %3.4F',[a[i]])
      else
         result := result + format('%3.4F',[a[i]]);
   end;

   result := result + ']';
end;

function Helpers.GetDeterminant: double;
begin
  result:=Ops.Determinant.Gen(self);
end;

function Helpers.GetInverse: TMatrix;
begin
  result:=Invert.GaussJordan(self);
end;

class function Determinant.Gen(const AMatrix: TMatrix): double;
begin
  result:=Matrix.Determinant(AMatrix);
end;

class function Determinant.Tridiagonal(const AMatrix: TMatrix): double;
begin
  result:=DeterminantTriangular(AMatrix);
end;

class function Determinant.Triangular(const AMatrix: TMatrix): double;
begin
  result:=DeterminantTriangular(AMatrix);
end;

class function Determinant.Unitary(const AMatrix: TMatrix): double;
begin
  result:=DeterminantUnitary(AMatrix);
end;

class procedure Decompose.QR.GramSchmidt(const AMatrix: TMatrix; out AQ, AR: TMatrix);
begin
  QRDecomposeGS(AMatrix, AQ,AR);
end;

class procedure Decompose.QR.Householder(const AMatrix: TMatrix; out AQ, AR: TMatrix);
begin
  QRDecomposeHR(AMatrix, AQ,AR);
end;

class procedure Decompose.LU(const AMatrix: TMatrix; out AL, AU: TMatrix);
begin
  LUDecompose(AMatrix, al, au);
end;

class function Decompose.Cholesky(const AMatrix: TMatrix): TMatrix;
begin
  result:=CholeskyDecompose(AMatrix);
end;

class function Invert.GaussJordan(const AMatrix: TMatrix): TMatrix;
begin
  result:=InvertGJ(AMatrix);
end;

class function Invert.LU(const AL, AU: TMatrix): TMatrix;
begin
  result:=InvertLU(AL, AU);
end;

class function Solve.LU(const AL, AU, AB: TMatrix): TMatrix;
begin
  result:=SolveLU(al,au,ab);
end;

end.
