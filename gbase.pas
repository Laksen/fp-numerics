unit gbase;

interface

{$mode delphi}
{$macro on}

uses
  sysutils, math, UComplex;

type
 TIndices = array of longint;

function Indices(const v: array of longint): TIndices;

{$define IsInterface}

type
  {$define comparable}
  {$define TDataType:=longint}
  {$define TFloatType:=double}
  {$define TGVector:=TGVectorI}
  TGVectorI = type {$i gbaseh.inc}

type
  {$define comparable}
  {$define isreal}
  {$define TDataType:=double}
  {$define TFloatType:=double}
  {$define TGVector:=TGVectorD}
  TGVectorD = type {$i gbaseh.inc}

type
  {$define comparable}
  {$define TDataType:=single}
  {$define TFloatType:=single}
  {$define TGVector:=TGVectorS}
  TGVectorS = type {$i gbaseh.inc}

type
  {$undef comparable}
  {$define TDataType:=complex}
  {$define TFloatType:=complex}
  {$define TGVector:=TGVectorC}
  TGVectorC = type {$i gbaseh.inc}

type
  TGVectorsI = array of TGVectorI;
  TGVectorsD = array of TGVectorD;
  TGVectorsS = array of TGVectorS;
  TGVectorsC = array of TGVectorC;

{$undef IsInterface}

function ZerosI(len: longint): TGVectorI; overload;
function ZerosS(len: longint): TGVectorS; overload;
function ZerosD(len: longint): TGVectorD; overload;
function ZerosC(len: longint): TGVectorC; overload;

function RandomUniform: double;
function RandUniform(ALow, AHigh: Double): Double;

implementation

uses
  ops;

const
  LCGMul = 1013904223;
  LCGInc = 1664525;

var
  LCGNext: longword = 0;

function RandomLCG: longword;
begin
  result:=LCGNext;
  LCGNext:=longword((LCGNext*LCGMul)+LCGInc);
end;

function RandomUniform: double;
begin
  result:=RandomLCG/high(LongWord);
end;

function RandUniform(ALow, AHigh: Double): Double;
begin
  result:=RandomUniform()*(AHigh-ALow)+ALow;
end;

function ZerosI(len: longint): TGVectorI; overload;
var i: longint;
begin
   if len < 0 then exit;

   SetLength(result, len);
   For i := 0 to len-1 do
      result[i] := 0;
end;

function ZerosS(len: longint): TGVectorS; overload;
var i: longint;
begin
   if len < 0 then exit;

   SetLength(result, len);
   For i := 0 to len-1 do
      result[i] := 0;
end;

function ZerosD(len: longint): TGVectorD; overload;
var i: longint;
begin
   if len < 0 then exit;

   SetLength(result, len);
   For i := 0 to len-1 do
      result[i] := 0;
end;

function ZerosC(len: longint): TGVectorC; overload;
var i: longint;
begin
   if len < 0 then exit;

   SetLength(result, len);
   For i := 0 to len-1 do
      result[i] := 0;
end;

function Indices(const v: array of longint): TIndices;
var i: longint;
begin
   setlength(result, length(v));
   for i := 0 to high(v) do
      result[i] := v[i];
end;

{$define IsImplementation}

  {$undef IsComplex}
  {$undef isreal}
  {$define VSqrt:=system.sqrt}
  {$define VSqr:=system.sqr}
  {$define VZeros:=zerosI}
  {$define VSin:=sin}
  {$define VCos:=cos}
  {$define VAbs:=system.abs}
  {$define comparable}
  {$define TDataType:=longint}
  {$define TFloatType:=double}
  {$define TGVector:=TGVectorI}
  {$i gbaseh.inc}

  {$undef IsComplex}
  {$define isreal}
  {$define VSqrt:=system.sqrt}
  {$define VSqr:=system.sqr}
  {$define VZeros:=zerosd}
  {$define VSin:=sin}
  {$define VCos:=cos}
  {$define VAbs:=system.abs}
  {$define comparable}
  {$define TDataType:=double}
  {$define TFloatType:=double}
  {$define TGVector:=TGVectorD}
  {$i gbaseh.inc}

  {$define comparable}
  {$define VZeros:=zeross}
  {$define TDataType:=single}
  {$define TFloatType:=single}
  {$define TGVector:=TGVectorS}
  {$i gbaseh.inc}

  {$define IsComplex}
  {$define VSqrt:=csqrt}
  {$define VSqr:=csqr}
  {$define VZeros:=zerosc}
  {$define VSin:=csin}
  {$define VCos:=ccos}
  {$define VAbs:=cmod}
  {$undef comparable}
  {$define TDataType:=complex}
  {$define TFloatType:=complex}
  {$define TGVector:=TGVectorC}
  {$i gbaseh.inc}

{$undef IsImplementation}

initialization
   DecimalSeparator := '.';
   randomize;
   LCGNext:=random(high(longint));

end.
