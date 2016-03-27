unit gmatrix;

{$mode objfpc}
{$modeswitch advancedrecords}
{$Macro on}
{$h+}

interface

uses
  SysUtils,
  gbase,
  UComplex;

{$define IsInterface}
type
  {$undef IsComplex}
  {$define TFloatType:=double}
  {$define TGMatrix:=TGMatrixD}
  {$define TGVector:=TGVectorD}
  TGMatrixD = {$i gmatrixh.inc}

type
  {$define TFloatType:=single}
  {$define TGMatrix:=TGMatrixS}
  {$define TGVector:=TGVectorS}
  TGMatrixS = {$i gmatrixh.inc}

type
  {$define IsComplex}
  {$define TFloatType:=complex}
  {$define TGMatrix:=TGMatrixC}
  {$define TGVector:=TGVectorC}
  TGMatrixC = {$i gmatrixh.inc}

{$undef IsInterface}

type
  TMatrix = TGMatrixD;

function RandomUniformD(ARows, AColumns: longint; ALow, AHigh: double): TGMatrixD;

implementation

uses
  Math, gops;

{$define IsImplementation}

  {$undef IsComplex}
  {$define vzeros:=ZerosD}
  {$define vsqrt:=system.sqrt}
  {$define vsqr:=system.sqr}
  {$define TFloatType:=double}
  {$define TGMatrix:=TGMatrixD}
  {$define TGVector:=TGVectorD}
  {$i gmatrixh.inc}

  {$define vzeros:=ZerosS}
  {$define TFloatType:=single}
  {$define TGMatrix:=TGMatrixS}
  {$define TGVector:=TGVectorS}
  {$i gmatrixh.inc}

  {$define IsComplex}
  {$define vzeros:=ZerosC}
  {$define vsqrt:=cmod}
  {$define vsqr:=csqr}
  {$define TFloatType:=complex}
  {$define TGMatrix:=TGMatrixC}
  {$define TGVector:=TGVectorC}
  {$i gmatrixh.inc}

{$undef IsImplementation}

function RandomUniformD(ARows, AColumns: longint; ALow, AHigh: double): TGMatrixD;
var
  i, i2: longint;
begin
  result:=TGMatrixD.Create(arows, AColumns);

  for i:=0 to ARows-1 do
    for i2:=0 to AColumns-1 do
      result[i,i2]:=RandUniform(ALow, AHigh);
end;

end.

