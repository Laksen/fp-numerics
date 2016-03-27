unit graphs;

{$mode objfpc}
{$modeswitch advancedrecords}
{$H+}

interface

uses
  gbase, gops,
  gmatrix;

type
  TGraphBuilder = record
  strict private
    Data: TGMatrixD;
    Nodes: longint;

    Edges: longint;
  public
    procedure GetEdges(out ASuccessors, APredecessors: TGMatrixD);
    function GetNodes: TGVectorI;

    function AddNode: LongInt;
    procedure AddEdge(ASource, ASink, ACapacity: longint);
  end;

function LayoutSpring(const ANodes: TGVectorI; AEdges: TGMatrixD): TGMatrixD;

implementation

procedure QuickSort(var A: TMatrix; iLo, iHi: Integer; AColumn: longint);
var Lo, Hi: Integer;
    pivot: double;
    T: TGVectorD;
begin
   if (iHi-iLo) <= 0 then exit;

   Lo := iLo;
   Hi := iHi;
   Pivot := A[(Lo + Hi) div 2, AColumn];

   repeat
      while A[Lo, AColumn] < Pivot do Inc(Lo);
      while A[Hi, AColumn] > Pivot do Dec(Hi);

      if Lo <= Hi then
      begin
         T := A.Row[Lo];
         A.Row[Lo] := A.Row[Hi];
         A.Row[Hi] := T;

         Inc(Lo) ;
         Dec(Hi) ;
      end;
   until Lo > Hi;
   if Hi > iLo then QuickSort(A, iLo, Hi, AColumn);
   if Lo < iHi then QuickSort(A, Lo, iHi, AColumn);
end;

function LayoutSpring(const ANodes: TGVectorI; AEdges: TGMatrixD): TGMatrixD;
const
  c1 = 2;
  c2 = 1;
  c3 = 1;
  c4 = 0.1;

  function IsEdge(ASrc, ADst: longint): boolean;
  var
    i: longint;
  begin
    for i:=0 to AEdges.rows-1 do
      begin
        if (AEdges[i,0]=asrc) and (AEdges[i, 1]=ADst) then
          exit(true);
      end;
    exit(false);
  end;

var
  s,t,rnd, i,i2: longint;
  Force: TGMatrixD;
  d: double;
  ps,pt,diff: TGVectorD;
begin
  result:=RandomUniformD(length(ANodes), 2, 0, 1);

  Force:=TMatrix.Create(result.Rows, 2);

  for rnd:=0 to 9 do
    begin
      FillChar(force.Data[0], length(force.data)*sizeof(double), 0);

      for i := 0 to high(ANodes) do
        for i2 := 0 to high(ANodes) do
          if i<>i2 then
            begin
              diff:=Result.Row[i]-Result.row[i2];
              d:=Magnitude(diff);

              if d=0 then
                continue;

              if IsEdge(i, i2) then
                force.row[i2]:=force.row[i2]+(diff/d*c1*ln(d/c2))
              else
                force.row[i2]:=force.row[i2]+(diff/d*(c3/sqr(d)))
            end;

      for i := 0 to force.rows-1 do
        result.row[i]:=result.row[i]+c4*force.Row[i];
    end;
end;

procedure TGraphBuilder.GetEdges(out ASuccessors, APredecessors: TMatrix);
begin
  ASuccessors:=Data.Submatrix(0,0,Edges,3);

  QuickSort(ASuccessors, 0, ASuccessors.Rows-1,0);
  APredecessors:=ASuccessors.Copy;
  QuickSort(APredecessors, 0, ASuccessors.Rows-1,1);
end;

function TGraphBuilder.GetNodes: TGVectorI;
begin
  result:=Sequence(0, nodes-1);
end;

function TGraphBuilder.AddNode: LongInt;
begin
  result:=Nodes;
  inc(nodes);
end;

procedure TGraphBuilder.AddEdge(ASource, ASink, ACapacity: longint);
begin
  if data.rows=0 then
    data:=TMatrix.Create(4,3)
  else if edges>=Data.Rows then
    data:=Data.Extend(0,0,data.rows*4 div 3, 0, false);

  data[edges,0]:=ASource;
  data[edges,1]:=ASink;
  data[edges,2]:=ACapacity;
  inc(edges);
end;

end.

