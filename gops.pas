unit gops;

{$mode objfpc}
{$macro on}

interface

uses
  gbase, ucomplex;

{$define IsInterface}

  {$define TFloatType:=double}
  {$define TGVector:=TGVectorD}
  {$i gopsh.inc}

  {$define TFloatType:=single}
  {$define TGVector:=TGVectorS}
  {$i gopsh.inc}

  {$define TFloatType:=complex}
  {$define TGVector:=TGVectorC}
  {$i gopsh.inc}

{$undef IsInterface}

implementation

uses
  sysutils, math;

function cpower(z1, z2 : complex): complex;
{ exp : z := z1 ** z2 }
begin
  result:= cexp(z2*cln(z1));
end;

{$define IsImplementation}

  {$define VPower:=power}
  {$define TFloatType:=double}
  {$define TGVector:=TGVectorD}
  {$i gopsh.inc}

  {$define TFloatType:=single}
  {$define TGVector:=TGVectorS}
  {$i gopsh.inc}

  {$define VPower:=cpower}
  {$define TFloatType:=complex}
  {$define TGVector:=TGVectorC}
  {$i gopsh.inc}

{$undef IsImplementation}

end.

