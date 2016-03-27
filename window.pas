unit window;

{$mode objfpc}{$H+}

interface

uses base;

function Hamming(Width: longint): TVector;
function Hann(Width: longint): TVector;
function Cosine(Width: longint): TVector;
function Lanczos(Width: longint): TVector;
function Tukey(Width: longint; alpha: double): TVector;
function Blackman(Width: longint; alpha: double = 0.16): TVector;

implementation

function Hamming(Width: longint): TVector;
var i: longint;
begin
   setlength(result, width);

   for i := 0 to width-1 do
      result[i] := 0.54-0.46*cos((2*pi*i)/(width-1));
end;

function Hann(Width: longint): TVector;
var i: longint;
begin
   setlength(result, width);
   
   if width <= 0 then
      exit;
   
   result[0] := 0.5*(1+cos((2*pi*width)/(width-1)));
   for i := 1 to width-1 do
      result[i] := 0.5*(1-cos((2*pi*width)/(width-1)));
end;

function Cosine(Width: longint): TVector;
var i: longint;
begin
   setlength(result, width);

   for i := 0 to width-1 do
      result[i] := sin((pi*i)/(width-1));
end;

function sinc(x: double): double;
begin
   result := sin(pi*x)/(pi*x);
end;

function Lanczos(Width: longint): TVector;
var i: longint;
begin
   setlength(result, width);

   for i := 0 to width-1 do
      result[i] := sinc((2*i)/(width-1)-1);
end;

function Tukey(Width: longint; alpha: double): TVector;
var i, an2: longint;
begin
   setlength(result, width);

   an2 := system.round(width*alpha/2);

   for i := 0 to width-1 do
   begin
      if i <= an2 then
         result[i] := 0.5*(1+cos(pi*(i/an2-1)))
      else if i >= (width-an2) then
         result[i] := 0.5*(1+cos(pi*(i/an2-2/alpha+1)))
      else
         result[i] := 1;
   end;
end;

function Blackman(Width: longint; alpha: double): TVector;
var i: longint;
    a0,a1,a2: double;
begin
   a0 := (1-alpha)/2;
   a1 := 1/2;
   a2 := alpha/2;

   setlength(result, width);

   for i := 0 to width-1 do
      result[i] := a0-a1*cos(2*pi*i/(width-1))+a2*cos(4*pi*i/(width-1));
end;

end.

