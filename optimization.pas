unit optimization;

{$mode objfpc}{$H+}

interface

uses
  gbase, gops,
  gmatrix;

type
  TOptVector = TGVectorD;

  TObjectiveFunc = function(AState: TOptVector; AData: pointer): double;

function OptimizeSA(AFunc: TObjectiveFunc; AData: pointer; const AStart, ALow, AHigh: TOptVector; AAcceptPropability: double; AMaxIter: longint): TOptVector;

function ClusterKMeansInitialized(const APoints: TGMatrixD; AClusters, AMaxIterations: longint; var ACentroids: TGMatrixD): TGVectorI;
function ClusterKMeans(const APoints: TGMatrixD; AClusters, AMaxIterations: longint; out ACentroids: TGMatrixD): TGVectorI;

implementation

uses
  math;

function Clamp(AValue, ALow, AHigh: double): double;
begin
  if AValue<ALow then
    exit(ALow)
  else if AValue>AHigh then
    exit(AHigh)
  else
    exit(AValue);
end;

function OptimizeSA(AFunc: TObjectiveFunc; AData: pointer; const AStart, ALow, AHigh: TOptVector; AAcceptPropability: double; AMaxIter: longint): TOptVector;
var
  temp,cooling, newbest, best: Double;
  new: TOptVector;
  i, i2, cnt: longint;
begin
  result:=copy(AStart);
  new:=copy(AStart);

  temp:=1;
  cooling:=10/AMaxIter;

  best:=AFunc(AStart, AData);

  for i:=0 to AMaxIter-1 do
    begin
      cnt:=0;
      for i2:=0 to high(new) do
        if random()<temp then
          begin
            new[i2]:=Clamp(new[i2]+0.2*(random()*2-1)*(ahigh[i2]-alow[i2]), alow[i2], ahigh[i2]);
            cnt:=1;
          end;

      if cnt=0 then
        begin
          i2:=random(length(new));
          new[i2]:=Clamp(new[i2]+0.2*(random()*2-1)*(ahigh[i2]-alow[i2]), alow[i2], ahigh[i2]);
        end;

      newbest:=AFunc(new, AData);

      if (best>newbest) or
         (random()<AAcceptPropability) then
        begin
          move(new[0], result[0], length(new)*sizeof(new[0]));
          best:=newbest;
        end
      else
        begin
          move(result[0], new[0], length(new)*sizeof(new[0]));
        end;

      temp:=temp-temp*cooling;
    end;
end;

function ClusterKMeansInitialized(const APoints: TGMatrixD; AClusters, AMaxIterations: longint; var ACentroids: TGMatrixD): TGVectorI;
var
  pc, pt: TGVectorD;
  i, i2, i3, cnt, bi: longint;
  bestd, nd, t: Double;
begin
  // Find nearest centroids
  result:=ZerosI(APoints.Rows);
  for i:=0 to high(result) do
    begin
      pc:=APoints.Row[i];
      pt:=ACentroids.Row[0]-pc;
      bestd:=dot(pt,pt);
      bi:=0;

      for i2:=1 to ACentroids.Rows-1 do
        begin
          pt:=ACentroids.Row[i2]-pc;
          nd:=dot(pt,pt);
          if nd<bestd then
            begin
              bestd:=nd;
              bi:=i2;
            end;
        end;

      result[i]:=bi;
    end;

  while AMaxIterations>0 do
    begin
      // Update centroids
      for i:=0 to ACentroids.Rows-1 do
        for i2:=0 to ACentroids.Columns-1 do
          begin
            t:=0;
            cnt:=0;
            for i3:=0 to high(result) do
              if result[i3]=i then
                begin
                  t:=t+APoints[i3,i2];
                  inc(cnt);
                end;

            if cnt<>0 then
              t:=t/cnt;

            ACentroids[i,i2]:=t;
          end;

      // Update clusters
      for i:=0 to high(result) do
        begin
          pc:=APoints.Row[i];
          pt:=ACentroids.Row[0]-pc;
          bestd:=dot(pt,pt);
          bi:=0;

          for i2:=1 to ACentroids.Rows-1 do
            begin
              pt:=ACentroids.Row[i2]-pc;
              nd:=dot(pt,pt);
              if nd<bestd then
                begin
                  bestd:=nd;
                  bi:=i2;
                end;
            end;

          result[i]:=bi;
        end;

      dec(AMaxIterations);
    end;
end;

function ClusterKMeans(const APoints: TGMatrixD; AClusters, AMaxIterations: longint; out ACentroids: TGMatrixD): TGVectorI;
var
  v: TGVectorD;
  i, i2: longint;
  vmin, vmax: Double;
begin
  ACentroids:=TGMatrixD.Create(AClusters, APoints.Columns);

  for i:=0 to APoints.Columns-1 do
    begin
      v:=APoints.Column[i];

      vmin:=min(v);
      vmax:=max(v);

      for i2:=0 to AClusters-1 do
        ACentroids[i2,i]:=RandomUniform()*(vmax-vmin)+vmin;
    end;

  result:=ClusterKMeansInitialized(APoints, AClusters, AMaxIterations, ACentroids);
end;

end.

