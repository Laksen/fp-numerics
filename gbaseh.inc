{$ifdef IsInterface}
array of TDataType;

function Vector(const v: array of TDataType): TGVector; overload;

function Sequence(AStart, AStop: longint): TGVector; overload;

function Concat(const a,b: TGVector): TGVector; overload;
function Pick(const a: TGVector; const I: TIndices): TGVector; overload;
{$ifdef comparable}
function Sort(const a: TGVector): TGVector; overload;
{$endif comparable}

{$ifdef isreal}
function NormEuclid(const t: TGVector): TFloatType; overload;
function Magnitude(const a: TGVector): TFloatType; overload;
{$endif isreal}

function Sum(const t: TGVector): TDataType; overload;
function Dot(const a,b: TGVector): TDataType; overload;
function DiffSqr(const a,b: TGVector): TDataType; overload;

function Diff(const y: TGVector): TGVector; overload;

{$ifdef isreal}
function Mean(const t: TGVector): TFloatType; overload;
{$ifdef comparable}
function Median(const t: TGVector): TFloatType; overload;
{$endif comparable}
function Variance(const t: TGVector): TFloatType; overload;
//function Covariance(const t: TMatrix): TMatrix;
function AutoCorrelation(const t: TGVector): TGVector; overload;
{$ifdef comparable}
function Hist(const t: TGVector; l,h: TFloatType): TGVector; overload;
{$endif comparable}
{$endif isreal}

{$ifdef comparable}
function Largest(const T: TGVector; N: longint): TIndices; overload;

function Max(const T: TGVector): TDataType; overload;
function Max(const T: TGVector; out Pos: longint): TDataType; overload;
function Min(const T: TGVector): TDataType; overload;
function Min(const T: TGVector; out Pos: longint): TDataType; overload;

function MaxPos(const T: TGVector): longint; overload;
function MinPos(const T: TGVector): longint; overload;

function Max(const T: array of TGVector): TDataType; overload;
function Min(const T: array of TGVector): TDataType; overload;

function RoundVector(const t: TGVector): TGVector; overload;
{$ifdef isreal}
procedure Normalize(var T: TGVector); overload;
procedure Normalize(var T: array of TGVector); overload;
{$endif isreal}
function Clamp(const T: TGVector; l,h: TDataType): TGVector; overload;
function Threshold(const T: TGVector; Thres, l, h: TDataType): TGVector; overload;
{$endif comparable}

{$ifdef isreal}
function Sine(const T: TGVector): TGVector; overload;
function Cosine(const T: TGVector): TGVector; overload;
{$endif isreal}

function EvalPoly(const T: TGVector; X: TDataType): TDataType; overload;

function Rand(const Len: longint): TGVector; overload;

function Conjugate(const t: TGVector): TGVector; overload;

function Sqr(const T: TGVector): TGVector; overload;
{$ifdef isreal}
function Sqrt(const T: TGVector): TGVector; overload;
{$endif isreal}
function Abs(const T: TGVector): TGVector; overload;
{$ifdef comparable}
{$ifdef isreal}
function Ln(const t: TGVector): TGVector; overload;
function Log(const t: TGVector): TGVector; overload;
function Exp(const T: TGVector): TGVector; overload;
{$endif isreal}
{$endif comparable}

function SubChunk(const T: TGVector; Index, Len: longint): TGVector; overload;
function Reverse(const T: TGVector): TGVector; overload;

{$ifdef comparable}
function Pow2(p: TDataType): longint; overload;

function DumpVector(const v: TGVector): ansistring; overload;
{$endif comparable}
//function DumpMatrix(const v: TMatrix): ansistring;
{$else}
function Vector(const v: array of TDataType): TGVector; overload;
begin
   setlength(result, length(v));

   move(v[0], result[0], length(v)*sizeof(TDataType));
end;

function Sequence(AStart, AStop: longint): TGVector; overload;
var
  i: longint;
begin
  setlength(result,AStop-AStart+1);
  for i:=0 to high(result) do
     result[i]:=i+AStart;
end;

function Concat(const a, b: TGVector): TGVector;
var i: longint;
begin
   setlength(result, length(a)+length(b));

   for i := 0 to high(a) do
      result[i] := a[i];
   for i := 0 to high(b) do
      result[i+length(a)] := b[i];
end;

function Pick(const a: TGVector; const I: TIndices): TGVector;
var i2: longint;
begin
   setlength(result, length(i));
   for i2 := 0 to high(i) do
      result[i2] := a[i[i2]];
end;

{$ifdef comparable}
procedure QuickSort(var A: TGVector; iLo, iHi: Integer); overload;
var Lo, Hi: Integer;
    pivot,t: TDataType;
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

function Sort(const a: TGVector): TGVector;
begin
   result := a;
   QuickSort(result, 0, high(result));
end;
{$endif}

function Conjugate(const t: TGVector): TGVector; overload;
var i: longint;
begin
{$ifdef IsComplex}
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := cong(t[i]);
{$else}
   result:=system.copy(t);
{$endif}
end;

function NormEuclid(const t: TGVector): TFloatType; overload;
var tmp: TDataType;
    i: longint;
begin
   tmp := 0;
   for i := 0 to high(t) do
      tmp := tmp + vsqr(t[i]);
{$ifdef comparable}
  if tmp=0 then exit(0);
{$endif}
   result := vsqrt(tmp);
end;

function Sum(const t: TGVector): TDataType;
var tmp: TDataType;
    i: longint;
begin
   tmp := 0;
   for i := 0 to high(t) do
      tmp := tmp + t[i];
   result := tmp;
end;

function Dot(const a, b: TGVector): TDataType;
var i: longint;
    t: TDataType;
begin
   t := 0;
   for i := 0 to high(a) do
      t := t + a[i]*b[i];
   result := t;
end;

function DiffSqr(const a, b: TGVector): TDataType;
var i: longint;
    t: TDataType;
begin
   t := 0;
   for i := 0 to high(a) do
      t := t + vsqr(a[i]-b[i]);
   result := t;
end;

function Magnitude(const a: TGVector): TFloatType; overload;
begin
  result:=NormEuclid(a);
end;

function Diff(const y: TGVector): TGVector;
var i: longint;
begin
   if length(y) <= 1 then
      exit(y);

   setlength(result, length(y)-1);

   for i := 0 to high(y)-1 do
      result[i] := y[i+1]-y[i];
end;

function Mean(const t: TGVector): TFloatType; overload;
var i,n: longint;
    tmp: TDataType;
begin
   n := length(t);
   if n <= 0 then exit(0);

   tmp := 0;
   for i := 0 to n-1 do
      tmp := tmp + t[i];
   result := tmp / n;
end;

{$ifdef comparable}
function Median(const t: TGVector): TFloatType; overload;
var h: TGVector;
begin
   if length(t) <= 0 then exit(0);

   h := sort(t);

   if odd(length(t)) then
      result := h[((length(t)+1) div 2)-1]
   else
      result := (h[length(t) div 2-1]+h[length(t) div 2])/2;
end;
{$endif comparable}

function Variance(const t: TGVector): TFloatType; overload;
var i,n: longint;
    tmp,m,h: TFloatType;
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

function AutoCorrelation(const t: TGVector): TGVector; overload;
var n,i,i2: longint;
    tmp: TDataType;
begin
   n := length(t)-1;

   result := (vzeros(n+1));

   for i := 0 to n do
   begin
      tmp := 0;
      for i2 := 0 to n-i do
         tmp := tmp + t[i2]*t[i2+i];
      result[i] := tmp;
   end;
end;

{$ifdef comparable}
function Hist(const t: TGVector; l,h: TFloatType): TGVector; overload;
var i,x: longint;
    fs: TFloatType;
begin
   result := (vzeros(length(t)));

   fs := (h-l)/high(t);

   for i := 0 to high(t) do
   begin
      x := system.round((t[i]-l)/fs);
      result[x] := result[x] + 1;
   end;
end;

function Largest(const T: TGVector; N: longint): TIndices;
var h: TGVector;
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

function Max(const T: TGVector): TDataType;
var i: longint;
begin
   result := Max(t, i);
end;

function Max(const T: TGVector; out Pos: longint): TDataType;
var tmp: TDataType;
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

function Min(const T: TGVector): TDataType;
var i: longint;
begin
   result := Min(t, i);
end;

function Min(const T: TGVector; out Pos: longint): TDataType;
var tmp: TDataType;
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

function MaxPos(const T: TGVector): longint;
begin
   Max(t, result);
end;

function MinPos(const T: TGVector): longint; overload;
begin
   Min(t, result);
end;

function Max(const T: array of TGVector): TDataType; overload;
var i: longint;
begin
   assert(length(t) > 0);

   result := Max(t[0]);

   for i := 1 to high(t) do
      result := max(result, max(t[i]));
end;

function Min(const T: array of TGVector): TDataType; overload;
var i: longint;
begin
   assert(length(t) > 0);

   result := Min(t[0]);

   for i := 1 to high(t) do
      result := Min(result, Min(t[i]));
end;

function RoundVector(const t: TGVector): TGVector; overload;
var i: longint;
begin
   setlength(result, length(t));

   for i := 0 to high(t) do
      result[i] := round(t[i]);
end;

{$ifdef isreal}
procedure Normalize(var T: TGVector); overload;
var i: longint;
    fl,fh: TDataType;
begin
   fl := min(t);
   fh := max(t);

   for i := 0 to high(t) do
      t[i] := (t[i]-fl)/(fh-fl);
end;

procedure Normalize(var T: array of TGVector); overload;
var i,i2: longint;
    fl,fh: TDataType;
begin
   fl := min(t);
   fh := max(t);

   if (fh-fl) = 0 then exit;

   for i := 0 to high(t) do
      for i2 := 0 to high(t[i]) do
         t[i,i2] := (t[i,i2]-fl)/(fh-fl);
end;
{$endif isreal}

function Clamp(const T: TGVector; l, h: TDataType): TGVector; overload;
var i: longint;
    x: TDataType;
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

function Threshold(const T: TGVector; Thres, l, h: TDataType): TGVector; overload;
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
{$endif comparable}

{$ifdef IsReal}
function Sine(const T: TGVector): TGVector; overload;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := vsin(t[i]);
end;

function Cosine(const T: TGVector): TGVector; overload;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := vcos(t[i]);
end;
{$endif IsReal}

function EvalPoly(const T: TGVector; X: TDataType): TDataType; overload;
var i: longint;
begin
   result := t[0];

   for i := 1 to high(t) do
   begin
      result := result + t[i]*x;
      x := x*x;
   end;
end;

function Rand(const Len: longint): TGVector; overload;
var i: longint;
begin
   setlength(result, len);

{$ifdef IsReal}
   for i := 0 to len-1 do
      result[i] := random(1000)/500-1;
{$else}
   for i := 0 to len-1 do
      result[i] := random(high(tdatatype));
{$endif IsReal}
end;

function Sqr(const T: TGVector): TGVector; overload;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := t[i]*t[i];
end;

{$ifdef IsReal}
function Sqrt(const T: TGVector): TGVector; overload;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := VSqrt(t[i]);
end;
{$endif IsReal}

function Abs(const T: TGVector): TGVector; overload;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := vabs(t[i]);
end;

{$ifdef comparable}
{$ifdef IsReal}
function Ln(const t: TGVector): TGVector; overload;
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

function Log(const t: TGVector): TGVector; overload;
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

function Exp(const T: TGVector): TGVector; overload;
var i: longint;
begin
   SetLength(result, length(t));

   for i := 0 to high(t) do
      result[i] := system.exp(t[i]);
end;
{$endif isreal}
{$endif comparable}

function SubChunk(const T: TGVector; Index, Len: longint): TGVector; overload;
var i: longint;
begin
   if index < 0 then
   begin
      len := len+index;
      index := 0;
   end;

   result := (VZeros(len));

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

function Reverse(const T: TGVector): TGVector; overload;
var i: longint;
begin
   setlength(result, length(t));

   for i := 0 to high(result) do
      result[i] := t[high(t)-i];
end;

{$ifdef comparable}
function Pow2(p: TDataType): longint; overload;
begin
   result := 1 shl ceil(log2(p));
end;

function DumpVector(const v: TGVector): ansistring; overload;
var m,i: longint;
  x: TDataType;
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
{$endif comparable}
{$endif}
