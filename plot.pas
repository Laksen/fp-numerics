unit plot;

{$mode objfpc}{$H+}

interface

uses
  sdl,
  sdl_gfx,
  gbase, gmatrix;

const
  PlotColors: array[0..16] of longword = (
    $FF0000FF,
    $FFFF00FF,
    $FF00FFFF,
    $00FF00FF,
    $00FFFFFF,
    $0000FFFF,
    $7F0000FF,
    $7F7F00FF,
    $7F007FFF,
    $007F00FF,
    $007F7FFF,
    $FF007FFF,
    $7F00FFFF,
    $00FF7FFF,
    $007FFFFF,
    $FF700FFF,
    $7FF00FFF
  );

function CreateWindow(AWidth, AHeight: longint): pointer;
function HasExit(AWindow: Pointer): boolean;
procedure WaitWindow(AWindow: Pointer);

procedure Present(AWindow: pointer);
procedure Clear(AWindow: pointer; AColor: longword);
procedure PlotData(AWindow: Pointer; AData: TGMatrixD; AColor, ALines, APoints: longword);
procedure PlotData(AWindow: Pointer; AX,AY: TGVectorD; AColor, ALines, APoints: longword);

procedure PlotGraph(AWindow: Pointer; ALocations, AEdges: TGMatrixD; ANodes: TGVectorI; AWidth, AHeight, AColor, APoints: longword);

implementation

uses
  math;

function CreateWindow(AWidth, AHeight: longint): pointer;
begin
  SDL_VideoInit(nil,0);
  result:=SDL_SetVideoMode(AWidth,AHeight,32,SDL_HWSURFACE or SDL_HWPALETTE or SDL_DOUBLEBUF);
end;

function HasExit(AWindow: Pointer): boolean;
var
  event : TSDL_Event;
begin
  result:=false;
  while SDL_PollEvent(@event)<>0 do
    begin
      case event.type_ of
        SDL_KEYDOWN:;
        SDL_QUITEV:
          begin
            result:=true;
            break;
          end;
      end;
    end;
end;

procedure WaitWindow(AWindow: Pointer);
var
  event : TSDL_Event;
begin
  while SDL_WaitEvent(@event)<>0 do
    begin
      case event.type_ of
        SDL_KEYDOWN:;
        SDL_QUITEV:break;
      end;
    end;
end;

procedure Present(AWindow: pointer);
begin
  SDL_Flip(AWindow);
end;

procedure Clear(AWindow: pointer; AColor: longword);
begin
  SDL_FillRect(AWindow, nil, AColor);
end;

procedure PlotData(AWindow: Pointer; AData: TGMatrixD; AColor, ALines, APoints: longword);
begin
  PlotData(AWindow, AData.Column[0], AData.Column[1], AColor, ALines, APoints);
end;

procedure PlotData(AWindow: Pointer; AX, AY: TGVectorD; AColor, ALines, APoints: longword);
var
  m,i,
  xl,xn,yl,yn: longint;
begin
  m:=min(length(ax), length(ay));

  xl:=round(ax[0]);
  yl:=round(ay[0]);

  for i:=0 to m-1 do
    begin
      xn:=round(ax[i]);
      yn:=round(ay[i]);

      if (ALines>0) and (i>0) then
        lineColor(AWindow, xl,yl,xn,yn, AColor);

      if APoints>0 then
        filledCircleColor(AWindow, xn,yn, APoints, AColor);

      xl:=xn;
      yl:=yn;
    end;
end;

procedure PlotGraph(AWindow: Pointer; ALocations, AEdges: TGMatrixD; ANodes: TGVectorI; AWidth, AHeight, AColor, APoints: longword);
var
  i,s,t, xl,xn,yl,yn: longint;
  vx, vy: TGVectorD;
  mxi, mxa, myi, mya, dx, dy: Double;
begin
  vx:=ALocations.Column[0];
  vy:=ALocations.Column[1];

  mxi:=(min(vx));
  mxa:=(max(vx));
  myi:=(min(vy));
  mya:=(max(vy));

  dx:=mxa-mxi;
  dy:=mya-myi;

  for i:=0 to AEdges.Rows-1 do
    begin
      s:=round(AEdges[i, 0]);
      t:=round(AEdges[i, 1]);

      xl:=round((ALocations[s,0]-mxi)/dx*AWidth);
      yl:=round((ALocations[s,1]-myi)/dy*AHeight);

      xn:=round((ALocations[t,0]-mxi)/dx*AWidth);
      yn:=round((ALocations[t,1]-myi)/dy*AHeight);

      lineColor(AWindow, xl,yl,xn,yn, AColor);

      if APoints>0 then
        filledCircleColor(AWindow, xn,yn, APoints, AColor);
    end;

  for s:=0 to ALocations.Rows-1 do
    begin
      xn:=round((ALocations[s,0]-mxi)/dx*AWidth);
      yn:=round((ALocations[s,1]-myi)/dy*AHeight);

      filledCircleColor(AWindow, xn,yn, APoints, AColor);
    end;
end;

initialization
  SDL_Init(SDL_INIT_VIDEO);
finalization
  SDL_Quit;

end.

