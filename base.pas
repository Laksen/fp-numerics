unit base;

interface

{$mode objfpc}

uses sysutils, math;

type
 TVector = array of double;
 TVectors = array of TVector;

 TIndices = array of longint;

function Vector(const v: array of double): TVector;
function Indices(const v: array of longint): TIndices;

function Sequence(AStart, AStop: longint): TVector;

function Zeros(len: longint): TVector;
function Concat(const a,b: TVector): TVector;
function Pick(const a: TVector; const I: TIndices): TVector;
function Sort(const a: TVector): TVector;

function Sum(const t: TVector): double;
function Dot(const a,b: TVector): double;
function Magnitude(const a: TVector): double;

function Diff(const y: TVector): TVector;

function Mean(const t: TVector): double;
function Median(const t: TVector): double;
function Variance(const t: TVector): double;
//function Covariance(const t: TMatrix): TMatrix;
function AutoCorrelation(const t: TVector): TVector;
function Hist(const t: TVector; l,h: double): TVector;

function Largest(const T: TVector; N: longint): TIndices;

function Max(const T: TVector): double; overload;
function Max(const T: TVector; out Pos: longint): double; overload;
function Min(const T: TVector): double; overload;
function Min(const T: TVector; out Pos: longint): double; overload;

function MaxPos(const T: TVector): longint;
function MinPos(const T: TVector): longint;

function Max(const T: array of TVector): double; overload;
function Min(const T: array of TVector): double; overload;

function RoundVector(const t: TVector): TVector;
procedure Normalize(var T: TVector);
procedure Normalize(var T: array of TVector);
function Clamp(const T: TVector; l,h: double): TVector;
function Threshold(const T: TVector; Thres, l, h: double): TVector;

function Sine(const T: TVector): TVector;
function Cosine(const T: TVector): TVector;

function EvalPoly(const T: TVector; X: double): double;

function Rand(const Len: longint): TVector;

function Sqr(const T: TVector): TVector;
function Sqrt(const T: TVector): TVector;
function Abs(const T: TVector): TVector;
function Ln(const t: TVector): TVector;
function Log(const t: TVector): TVector;
function Exp(const T: TVector): TVector;

function SubChunk(const T: TVector; Index, Len: longint): TVector;
function Reverse(const T: TVector): TVector;

function Pow2(p: double): longint;

function DumpVector(const v: TVector): ansistring;
//function DumpMatrix(const v: TMatrix): ansistring;

implementation

uses ops;

function Vector(const v: array of double): TVector;
begin
   setlength(result, length(v));
   
   move(v[0], result[0], length(v)*sizeof(double));
end;

function Indices(const v: array of longint): TIndices;
var i: longint;
begin
   setlength(result, length(v));
   for i := 0 to high(v) do
      result[i] := v[i];
end;

function Sequence(AStart, AStop: longint): TVector;
var
  i: longint;
begin
  setlength(result,AStop-AStart+1);
  for i:=0 to high(result) do
     result[i]:=i+AStart;
end;

function Zeros(len: longint): TVector;
var i: longint;
begin
   if len < 0 then exit;

   SetLength(result, len);
   For i := 0 to len-1 do
      result[i] := 0;
end;

function Concat(const a, b: TVector): TVector;
var i: longint;
begin
   setlength(result, length(a)+length(b));

   for i := 0 to high(a) do
      result[i] := a[i];
   for i := 0 to high(b) do
      result[i+length(a)] := b[i];
end;

function Pick(const a: TVector; const I: TIndices): TVector;
var i2: longint;
begin
   setlength(result, length(i));
   for i2 := 0 to high(i) do
      result[i2] := a[i[i2]];
end;

procedure QuickSort(var A: TVector; iLo, iHi: Integer) ;
var Lo, Hi: Integer;
    pivot,t: double;
begin
   if (iHi-iLo) <= 0 then exit;

   Lo := iLo;
   Hi := iHi;
   Pivot := A[(Lo + Hi) div 2];

   repeat
      while A[Lo] < Pivot do Inc(Lo);
      while A[Hi] > Pivot do Dec(Hi);

      if Lo <= Hi then
      begin
         T := A[Lo];
         A[Lo] := A[Hi];
         A[Hi] := T;
         Inc(Lo) ;
         Dec(Hi) ;
      end;
   until Lo > Hi;
   if Hi > iLo then QuickSort(A, iLo, Hi) ;
   if Lo < iHi then QuickSort(A, Lo, iHi) ;
end;

function Sort(const a: TVector): TVector;
begin
   result := a;
   QuickSort(result, 0, high(result));
end;

function Sum(const t: TVector): double;
var tmp: double;
    i: longint;
begin
   tmp := 0;
   for i := 0 to high(t) do
      tmp := tmp + t[i];
   result := tmp;
end;

function Dot(const a, b: TVector): double;
var i: longint;
    t: double;
begin
   t := 0;
   for i := 0 to high(a) do
      t := t + a[i]*b[i];
   result := t;
end;

function Magnitude(const a: TVector): double;
begin
  result:=system.sqrt(dot(a,a));
end;

function Diff(const y: TVector): TVector;
var i: longint;
begin
   if length(y) <= 1 then
      exit(y);

   setlength(result, length(y)-1);

   for i := 0 to high(y)-1 do
      result[i] := y[i+1]-y[i];
end;

function Mean(const t: TVector): double;
var i,n: longint;
    tmp: double;
begin
   n := length(t);
   if n <= 0 then exit(0);
   
   tmp := 0;
   for i := 0 to n-1 do
      tmp := tmp + t[i];
   result := tmp / n;
end;

function Median(const t: TVector): double;
var h: TVector;
begin
   if length(t) <= 0 then exit(0);

   h := sort(t);

   if odd(length(t)) then
      result := h[((length(t)+1) div 2)-1]
   else
      result := (h[length(t) div 2-1]+h[length(t) div 2])/2;
end;

function Variance(const t: TVector): double;
var i,n: longint;
    tmp,m,h: double;
begin
   n := length(t);
   if n <= 1 then exit(0);
   
   m := Mean(t);
   tmp := 0;
   
   for i := 0 to n-1 do
   begin
      h := t[i]-m;
      tmp := tmp + h*h;
   end;
   
   result := tmp/(n-1);
end;

{function Covariance(const t: TMatrix): TMatrix;
var n,i,j: longint;
begin
   n := length(t);
   setlength(result, n);

   for i := 0 to n-1 do
   begin
      setlength(result[i], n);
      for j := 0 to n-1 do
         result[i, j] := 1/(n-1)*sum((t[i]-mean(t[i]))*(t[j]-mean(t[j])));
   end;
end;}

function AutoCorrelation(const t: TVector): TVector;
var n,i,i2: longint;
    tmp: double;
begin
   n := length(t)-1;

   result := zeros(n+1);

   for i := 0 to n do
   begin
      tmp := 0;
      for i2 := 0 to n-i do
         tmp := tmp + t[i2]*t[i2+i];
      result[i] := tmp;
   end;
end;

function Hist(const t: TVector; l,h: double): TVector;
var i,x: longint;
    fs: double;
begin
   result := zeros(length(t));

   fs := (h-l)/high(t);

   for i := 0 to high(t) do
   begin
      x := system.round((t[i]-l)/fs);
      result[x] := result[x] + 1;
   end;
end;

function Largest(const T: TVector; N: longint): TIndices;
var h: TVector;
    i: longint;
begin
   setlength(result, n);

   h := copy(t);

   for i := 0 to n-1 do
   begin
      max(h, result[i]);
      h[result[i]] := min(h);
   end;
end;

function Max(const T: TVector): double;
var i: longint;
begin
   result := Max(t, i);
end;

function Max(const T: TVector; out Pos: longint): double;
var tmp: Double;
    i, bi: longint;
begin
   assert(length(t) > 0, 'Cannot find minimum of zero-length sequence');

   tmp := t[0];
   bi := 0;

   for i := 0 to high(t) do
      if t[i] > tmp then
      begin
         bi := i;
         tmp := t[i];
      end;

   result := tmp;
   pos := bi;
end;

function Min(const T: TVector): double;
var i: longint;
begin
   result := Min(t, i);
end;

function Min(const T: TVector; out Pos: longint): double;
var tmp: Double;
    i, bi: longint;
begin
   assert(length(t) > 0, 'Cannot find minimum of zero-length sequence');

   tmp := t[0];
   bi := 0;

   for i := 1 to high(t) do
      if t[i] < tmp then
      begin
         bi := i;
         tmp := t[i];
      end;

   result := tmp;
   Pos := bi;
end;

function MaxPos(const T: TVector): longint;
begin
   Max(t, result);
end;

function MinPos(const T: TVector): longint;
begin
   Min(t, result);
end;

function Max(const T: array of TVector): double;
var i: longint;
begin
   assert(length(t) > 0);

   result := Max(t[0]);

   for i := 1 to high(t) do
      result := max(result, max(t[i]));
end;

function Min(const T: array of TVector): double;
var i: longint;
begin
   assert(length(t) > 0);

   result := Min(t[0]);

   for i := 1 to high(t) do
      result := Min(result, Min(t[i]));
end;

function RoundVector(const t: TVector): TVector;
var i: longint;
begin
   setlength(result, length(t));

   for i := 0 to high(t) do
      result[i] := round(t[i]);
end;

procedure Normalize(var T: TVector);
var i: longint;
    fl,fh: double;
begin
   fl := min(t);
   fh := max(t);

   for i := 0 to high(t) do
      t[i] := (t[i]-fl)/(fh-fl);
end;

procedure Normalize(var T: array of TVector);
var i,i2: longint;
    fl,fh: double;
begin
   fl := min(t);
   fh := max(t);

   if (fh-fl) = 0 then exit;

   for i := 0 to high(t) do
      for i2 := 0 to high(t[i]) do
         t[i,i2] := (t[i,i2]-fl)/(fh-fl);
end;

function Clamp(const T: TVector; l, h: double): TVector;
var i: longint;
    x: Double;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
   begin
      x := t[i];
      if x > h then
         result[i] := h
      else if x < l then
         result[i] := l
      else
         result[i] := x;
   end;
end;

function Threshold(const T: TVector; Thres, l, h: double): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
   begin
      if t[i] >= thres then
         result[i] := h
      else
         result[i] := l;
   end;
end;

function Sine(const T: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := sin(t[i]);
end;

function Cosine(const T: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := cos(t[i]);
end;

function EvalPoly(const T: TVector; X: double): double;
var i: longint;
begin
   result := t[0];

   for i := 1 to high(t) do
   begin
      result := result + t[i]*x;
      x := x*x;
   end;
end;

function Rand(const Len: longint): TVector;
var i: longint;
begin
   setlength(result, len);

   for i := 0 to len-1 do
      result[i] := random(1000)/500-1;
end;

function Sqr(const T: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := t[i]*t[i];
end;

function Sqrt(const T: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := System.sqrt(t[i]);
end;

function Abs(const T: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := system.abs(t[i]);
end;

function Ln(const t: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
   begin
      if t[i] <= 0 then
         result[i] := 0
      else
         result[i] := system.ln(t[i]);
   end;
end;

function Log(const t: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
   begin
      if t[i] <= 0 then
         result[i] := 0
      else
         result[i] := log10(t[i]);
   end;
end;

function Exp(const T: TVector): TVector;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := system.exp(t[i]);
end;

function SubChunk(const T: TVector; Index, Len: longint): TVector;
var i: longint;
begin
   if index < 0 then
   begin
      len := len+index;
      index := 0;
   end;

   result := Zeros(len);

   if (index+Len) > length(t) then
   begin
      for i := 0 to (high(t)-index)-1 do
         result[i] := t[index+i];
   end
   else
   begin
      for i := 0 to len-1 do
         result[i] := t[index+i];
   end;
end;

function Reverse(const T: TVector): TVector;
var i: longint;
begin
   setlength(result, length(t));

   for i := 0 to high(result) do
      result[i] := t[high(t)-i];
end;

function Pow2(p: double): longint;
begin
   result := 1 shl ceil(log2(p));
end;

function DumpVector(const v: TVector): ansistring;
var m,i: longint;
  x: Double;
begin
   x:=max(v);
   if x<0 then
     m := floor(log10(-x))+6
   else if x>0 then
     m := floor(log10(x))+5
   else
     m := 5;

   result := format(format('%%%D.2F',[m]),[v[0]]);
   for i := 1 to high(v) do
      result := result + format(format(', %%%D.2F',[m]),[v[i]]);
end;

{function DumpMatrix(const v: TMatrix): ansistring;
var m,i,i2: longint;
begin
   m := trunc(log10(max(v)))+5;

   result := '';
   for i2 := 0 to high(v) do
   begin
      result := result + format(format('%%%D.2F',[m]),[v[i2,0]]);
      for i := 1 to high(v[i2]) do
         result := result + format(format(', %%%D.2F',[m]),[v[i2,i]]);
      if i2 <> high(v) then result := result + lineending;
  end;
end;}

initialization
   DecimalSeparator := '.';

end.
