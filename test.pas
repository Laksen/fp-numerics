program test;

{$mode delphi}

uses
  gmatrix, gbase,
  optimization,
  graphs,
  plot,
  sysutils;

procedure TTTT;
var
  c: Pointer;
  gb: TGraphBuilder;
  s, t, a, b: LongInt;
  nodes: TGVectorI;
  pre, edge, loc: TGMatrixD;
begin
  c:=CreateWindow(500,500);

  FillChar(gb, sizeof(gb), 0);
  s:=gb.AddNode;
  t:=gb.AddNode;

  a:=gb.AddNode;
  b:=gb.AddNode;

  gb.AddEdge(s,a,1);
  gb.AddEdge(s,b,1);
  gb.AddEdge(a,b,1);
  gb.AddEdge(a,t,1);
  gb.AddEdge(b,t,1);

  nodes:=gb.GetNodes;
  gb.GetEdges(edge, pre);

  //loc:=RandomUniformD(length(nodes), 2, 0, 499);
  loc:=LayoutSpring(nodes, edge);

  PlotGraph(c, loc, edge, nodes, 500, 500, $FF0000FF, 4);

  Present(c);

  WaitWindow(c);
end;

const
  clusters = 9;

var
  x, ct, cl, ct2: TGMatrixD;
  pts: TGVectorI;
  cnt, i, cli, i5: Integer;
  c: Pointer;
begin
  tttt;
  exit;

  x:=RandomUniformD(20000, 2, 0,499);

  c:=CreateWindow(500,500);

  ct2:=RandomUniformD(Clusters, x.Columns,0,499);

  for i5:=1 to 1000 do
    begin
      ct:=ct2.Copy;

      Clear(c, 0);

      pts:=ClusterKMeansInitialized(x, clusters, i5, ct);

      for cli:=0 to clusters-1 do
        begin
          cnt:=0;
          for i:=0 to x.Rows-1 do
            if pts[i]=cli then
              inc(cnt);

          if cnt<=0 then continue;

          cl:=TGMatrixD.Create(cnt, 2);

          cnt:=0;
          for i:=0 to x.Rows-1 do
            if pts[i]=cli then
              begin
                cl.row[cnt]:=x.row[i];
                inc(cnt);
              end;

          PlotData(c, cl, PlotColors[cli], 0, 2);
        end;

      PlotData(c, ct, $FFFFFFFF, 0, 2);
      present(c);

      if HasExit(c) then
        exit;

      sleep(50);
    end;

  WaitWindow(c);
end.
