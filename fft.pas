unit fft;

interface

{$mode objfpc}

uses base;

type
 TFFT = class
 private
  fPoints: PtrInt;
  fIndex: array of PtrInt;
  fTwiddleR, fTwiddleI: array of single;
  procedure CalcFFT(var OutR, OutI: TVector);
 public
  procedure CalcFFTRealAbs(const InpR: TVector; var OutR: TVector);
  procedure CalcFFTReal(const InpR: TVector; var OutR, OutI: TVector);
  procedure CalcFFTComplex(const InpR, InpI: TVector; var OutR, OutI: TVector);
  
  constructor Create(Points: PtrInt; Inverse: boolean = false);
  destructor Destroy; override;
  
  property Points: PtrInt read fPoints;
 end;

function FFTRealReal(const T: TVector): TVector; // Transform a vector of only real components, and return the real components of the FFT
function FFTRealRealSqr(const T: TVector): TVector; // Transform a vector of only real components, and return the square of the FFT result
function FFTRealRealAbs(const T: TVector): TVector; // Transform a vector of only real components, and return the magnitude of the FFT result

function IFFTRealReal(const T: TVector): TVector;
function IFFTRealRealSqr(const T: TVector): TVector;
function IFFTRealRealAbs(const T: TVector): TVector;

implementation

function BitReverse(X, N: PtrInt): PtrInt;
begin
   result := 0;
   while n > 1 do
   begin
      result := result shl 1;
      if (x and 1) = 1 then
         result := result or 1;
      x := x shr 1;
      n := n shr 1;
   end;
end;

function FFTRealReal(const T: TVector): TVector;
var f: TFFT;
    tmp: TVector;
begin
   setlength(result, length(t));
   setlength(tmp, length(t));

   f := TFFT.Create(length(t));
   f.CalcFFTReal(t, result, tmp);
   f.Free;
end;

function FFTRealRealSqr(const T: TVector): TVector;
var f: TFFT;
begin
   setlength(result, length(t));

   f := TFFT.Create(length(t));
   f.CalcFFTRealAbs(t, result);
   f.Free;

   result := sqr(result);
end;

function FFTRealRealAbs(const T: TVector): TVector;
var f: TFFT;
begin
   setlength(result, length(t));

   f := TFFT.Create(length(t));
   f.CalcFFTRealAbs(t, result);
   f.Free;
end;

function IFFTRealReal(const T: TVector): TVector;
var f: TFFT;
    tmp: TVector;
begin
   setlength(result, length(t));
   setlength(tmp, length(t));

   f := TFFT.Create(length(t), true);
   f.CalcFFTReal(t, result, tmp);
   f.Free;
end;

function IFFTRealRealSqr(const T: TVector): TVector;
var f: TFFT;
begin
   setlength(result, length(t));

   f := TFFT.Create(length(t), true);
   f.CalcFFTRealAbs(t, result);
   f.Free;

   result := sqr(result);
end;

function IFFTRealRealAbs(const T: TVector): TVector;
var f: TFFT;
begin
   setlength(result, length(t));

   f := TFFT.Create(length(t), true);
   f.CalcFFTRealAbs(t, result);
   f.Free;
end;

procedure TFFT.CalcFFT(var OutR, OutI: TVector);
var n,nb, pti,i: PtrInt;
    xr,xi,sr,si,yr,yi,tr,ti: single;
begin
   n := 2;
   pti := 0;
   nb := 1;
   repeat
      pti := 0;
      
      repeat
         
         for i := 0 to nb-1 do
         begin
            xr := outr[pti+i];
            xi := outi[pti+i];
            yr := outr[pti+i+nb];
            yi := outi[pti+i+nb];
            
            tr := fTwiddleR[fPoints*i div n];
            ti := fTwiddleI[fPoints*i div n];
            
            sr := yr*tr-yi*ti;
            si := yr*ti+yi*tr;
            
            outr[pti+i] := (xr+sr)/2;
            outi[pti+i] := (xi+si)/2;
            outr[pti+i+nb] := (xr-sr)/2;
            outi[pti+i+nb] := (xi-si)/2;
         end;
         
         inc(pti, n);
      until pti >= fPoints;
      
      nb := n;
      n := n shl 1;
   until n > fPoints;
end;

procedure TFFT.CalcFFTRealAbs(const InpR: TVector; var OutR: TVector);
var i: ptrint;
    tmpi: TVector;
begin
   //Bit reverse input
   setlength(tmpi, length(inpr));
   
   for i := 0 to fPoints-1 do
   begin
      OutR[i] := Inpr[fIndex[i]];
      tmpi[i] := 0;
   end;
   
   CalcFFT(OutR, tmpi);
   
   for i := 0 to fPoints-1 do
      OutR[i] := system.sqrt(outR[i]*outR[i]+tmpi[i]*tmpi[i]);
end;

procedure TFFT.CalcFFTReal(const InpR: TVector; var OutR, OutI: TVector);
var i: ptrint;
begin
   //Bit reverse input
   for i := 0 to fPoints-1 do
   begin
      OutR[i] := InpR[fIndex[i]];
      OutI[i] := 0;
   end;
   
   CalcFFT(OutR, OutI);
end;

procedure TFFT.CalcFFTComplex(const InpR, InpI: TVector; var OutR, OutI: TVector);
var i: ptrint;
begin
   //Bit reverse input
   for i := 0 to fPoints-1 do
   begin
      OutR[i] := InpR[fIndex[i]];
      OutI[i] := InpI[fIndex[i]];
   end;
   
   CalcFFT(OutR, OutI);
end;

constructor TFFT.Create(Points: PtrInt; Inverse: boolean);
var i, sign: ptrint;
begin
   inherited Create;
   assert(system.sqrt(points*points)=points, 'FFT must be of a power of 2');

   fPoints := Points;
   
   setlength(fIndex, points);
   for i := 0 to Points-1 do fIndex[i] := BitReverse(i, Points);
   
   setlength(fTwiddleR, points div 2);
   setlength(fTwiddleI, points div 2);
   
   if inverse then
      sign := 1
   else
      sign := -1;
   
   for i := 0 to (points div 2)-1 do
   begin
      fTwiddleR[i] := cos(sign*i/points*2*pi);
      fTwiddleI[i] := sin(sign*i/points*2*pi);
   end;
end;

destructor TFFT.Destroy;
begin
   setlength(fIndex, 0);
   setlength(fTwiddleR, 0);
   setlength(fTwiddleI, 0);
   inherited Destroy;
end;

end.
