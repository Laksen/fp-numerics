{$ifdef IsInterface}
operator +(const a,b: TGVector): TGVector;
operator -(const a,b: TGVector): TGVector;
operator *(const a,b: TGVector): TGVector;
operator /(const a,b: TGVector): TGVector;

operator +(const a: TGVector; const b: TFloatType): TGVector;
operator +(const a: TFloatType; const b: TGVector): TGVector;
operator -(const a: TGVector; const b: TFloatType): TGVector;
operator -(const a: TFloatType; const b: TGVector): TGVector;
operator *(const a: TGVector; const b: TFloatType): TGVector;
operator *(const a: TFloatType; const b: TGVector): TGVector;
operator /(const a: TGVector; const b: TFloatType): TGVector;
operator /(const a: TFloatType; const b: TGVector): TGVector;

operator **(const a, b: TFloatType): TFloatType ;
operator **(const a: TGVector; const b: TFloatType ): TGVector;
{$else}
operator +(const a,b: TGVector): TGVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]+b[i];
end;

operator -(const a,b: TGVector): TGVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]-b[i];
end;

operator *(const a,b: TGVector): TGVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]*b[i];
end;

operator /(const a,b: TGVector): TGVector;
var i: longint;
begin
   assert(length(a)=length(b));
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]/b[i];
end;

operator +(const a: TGVector; const b: TFloatType ): TGVector;
var i: longint;
begin
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]+b;
end;

operator +(const a: TFloatType ; const b: TGVector): TGVector;
var i: longint;
begin
   setlength(result, length(b));

   for i := 0 to high(b) do
      result[i] := a+b[i];
end;

operator -(const a: TGVector; const b: TFloatType ): TGVector;
var i: longint;
begin
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]-b;
end;

operator -(const a: TFloatType ; const b: TGVector): TGVector;
var i: longint;
begin
   setlength(result, length(b));

   for i := 0 to high(b) do
      result[i] := a-b[i];
end;

operator *(const a: TGVector; const b: TFloatType ): TGVector;
var i: longint;
begin
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]*b;
end;

operator *(const a: TFloatType ; const b: TGVector): TGVector;
var i: longint;
begin
   setlength(result, length(b));

   for i := 0 to high(b) do
      result[i] := a*b[i];
end;

operator /(const a: TGVector; const b: TFloatType ): TGVector;
var i: longint;
begin
   assert(b<>0, 'Cannot divide by zero');
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := a[i]/b;
end;

operator /(const a: TFloatType ; const b: TGVector): TGVector;
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

operator **(const a, b: TFloatType ): TFloatType ;
begin
   result := vpower(a,b);
end;

operator **(const a: TGVector; const b: TFloatType ): TGVector;
var i: longint;
begin
   setlength(result, length(a));

   for i := 0 to high(a) do
      result[i] := vpower(a[i],b);
end;
{$endif}
