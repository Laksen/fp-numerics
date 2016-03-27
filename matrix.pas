unit matrix;

{$mode objfpc}
{$modeswitch advancedrecords}
{$h+}

interface

uses
  SysUtils,
  base;

type
  TMatrix = record
  strict private
    function GetColumn(AColumn: longint): TVector;
    function GetCopy: TMatrix;
    function GetHermitian: TMatrix;
    function GetMinorColumn(AStartRow, AStartColumn, AColumn: longint): TVector;
    function GetRow(ARow: longint): TVector;
    function GetTranspose: TMatrix;
    function GetValue(ARow, AColumn: longint): double;
    procedure SetColumn(AColumn: longint; AValue: TVector);
    procedure SetRow(ARow: longint; AValue: TVector);
    procedure SetValue(ARow, AColumn: longint; AValue: double);
  public
    Data: array of double;
    Rows, Columns: longint;

    function Submatrix(AStartRow,AStartColumn: longint): TMatrix;
    function Submatrix(AStartRow,AStartColumn,ARows,AColumns: longint): TMatrix;
    function Extend(ARowS,AColumnS, ARowE,AColumnE: longint; AIdent: boolean): TMatrix;

    constructor Identity(ARows, AColumns: longint);
    constructor Identity(N: longint);
    constructor Create(ARows, AColumns: longint);
    constructor Create(ARows, AColumns: longint; const AValues: array of double);

    procedure SubScaled(ASrc, ADst: longint; AValue: double);
    procedure ScaleRow(ARow: longint; AValue: double);
    procedure Scale(AValue: double);

    function Dump: string;

    property Transpose: TMatrix read GetTranspose;
    property Hermitian: TMatrix read GetHermitian;
    property Copy: TMatrix read GetCopy;

    property Value[ARow, AColumn: longint]: double read GetValue write SetValue; default;

    property Row[ARow: longint]: TVector read GetRow write SetRow;
    property Column[AColumn: longint]: TVector read GetColumn write SetColumn;
    property MinorColumn[AStartRow,AStartColumn,AColumn: longint]: TVector read GetMinorColumn;
  end;

operator explicit(const A: TVector): TMatrix;
operator explicit(const A: TMatrix): TVector;

operator -(const A: TMatrix): TMatrix;

operator +(const A, B: TMatrix): TMatrix;
operator -(const A, B: TMatrix): TMatrix;
operator * (const A, B: TMatrix): TMatrix;
operator * (const A: TMatrix; const B: double): TMatrix;
operator * (const B: double; const A: TMatrix): TMatrix;

procedure LUDecompose(const AMatrix: TMatrix; out AL, AU: TMatrix);

procedure QRDecomposeGS(const AMatrix: TMatrix; out AQ, AR: TMatrix); // Gram-Schmidt
procedure QRDecomposeHR(const AMatrix: TMatrix; out AQ, AR: TMatrix); // Householder reflection

function ReduceHR(const AMatrix: TMatrix): TMatrix; // Reduction to tridiagonal by Householder reflection

function CholeskyDecompose(const AMatrix: TMatrix): TMatrix; // Returns L where AMatrix=L*L.Transpose

procedure ArnoldiIteration(const A: TMatrix; out H, Q: TMatrix);

function SolveLU(const AL, AU, Ab: TMatrix): TMatrix;

function InvertGJ(const A: TMatrix): TMatrix;
function InvertLU(const AL, AU: TMatrix): TMatrix;

function Determinant(const AMatrix: TMatrix): double;
function DeterminantUnitary(const AMatrix: TMatrix): double;
function DeterminantTriangular(const AMatrix: TMatrix): double;

implementation

uses
  Math, ops;

operator explicit(const A: TVector): TMatrix;
begin
  result:=TMatrix.Create(length(a), 1, a);
end;

operator explicit(const A: TMatrix): TVector;
begin
  Assert((a.Rows = 1) or (a.Columns = 1));

  Result := copy(a.Data);
end;

operator-(const A: TMatrix): TMatrix;
var
  pa, pr: PDouble;
  i: longint;
begin
  pa := @a.Data[0];

  Result := TMatrix.Create(a.Rows, a.Columns);

  pr := @Result.Data[0];

  for i := 0 to high(a.Data) do
  begin
    pr^ := -pa^;
    Inc(pa);
    Inc(pr);
  end;
end;

operator +(const A, B: TMatrix): TMatrix;
var
  pa, pb, pr: PDouble;
  i: longint;
begin
  Assert(a.Rows = b.Rows);
  Assert(a.Columns = b.Columns);

  pa := @a.Data[0];
  pb := @b.Data[0];

  Result := TMatrix.Create(a.Rows, a.Columns);

  pr := @Result.Data[0];

  for i := 0 to high(a.Data) do
  begin
    pr^ := pa^ + pb^;
    Inc(pa);
    Inc(pb);
    Inc(pr);
  end;
end;

operator -(const A, B: TMatrix): TMatrix;
var
  pa, pb, pr: PDouble;
  i: longint;
begin
  Assert(a.Rows = b.Rows);
  Assert(a.Columns = b.Columns);

  pa := @a.Data[0];
  pb := @b.Data[0];

  Result := TMatrix.Create(a.Rows, a.Columns);

  pr := @Result.Data[0];

  for i := 0 to high(a.Data) do
  begin
    pr^ := pa^ - pb^;
    Inc(pa);
    Inc(pb);
    Inc(pr);
  end;
end;

operator * (const A, B: TMatrix): TMatrix;
var
  pa, pb, pr: PDouble;
  i, i2, i3:  longint;
  sum:        double;
begin
  Assert(a.Columns = b.Rows);

  Result := TMatrix.Create(a.Rows, b.Columns);

  pr := @Result.Data[0];

  for i := 0 to a.Rows - 1 do
    for i2 := 0 to b.Columns - 1 do
    begin
      sum := 0;

      pa := @a.Data[i * a.Columns];
      pb := @b.Data[i2];

      for i3 := 0 to a.Columns - 1 do
      begin
        sum := sum + pa^ * pb^;
        Inc(pa);
        Inc(pb, b.Columns);
      end;

      pr^ := sum;
      Inc(pr);
    end;
end;

operator*(const A: TMatrix; const B: double): TMatrix;
begin
  result:=A.Copy;
  result.Scale(b);
end;

operator*(const B: double; const A: TMatrix): TMatrix;
begin
  result:=A.Copy;
  result.Scale(b);
end;

procedure LUDecompose(const AMatrix: TMatrix; out AL, AU: TMatrix);
var
  i, i2: longint;
  t: Double;
  Ln: TMatrix;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in LUDecompose');

  AL := TMatrix.Identity(AMatrix.Rows);
  AU:=AMatrix;

  for i:=0 to AMatrix.Rows-2 do
    begin
      Ln:=TMatrix.Identity(AMatrix.Rows);

      t:=AU[i,i];
      assert(t<>0, 'Zero entry in diagonal at iteration '+inttostr(i));
      t:=1/t;

      for i2:=i+1 to AMatrix.Rows-1 do
        Ln[i2,i]:=-AU[i2,i]*t;

      AL:=Ln*AL;
      AU:=Ln*AU;
    end;

  for i:=0 to AL.Rows-2 do
    begin
      for i2:=i+1 to AL.Rows-1 do
          AL[i2,i]:=-AL[i2,i];
    end;
end;

procedure QRDecomposeGS(const AMatrix: TMatrix; out AQ, AR: TMatrix);

  function Proj(e,a: TVector): TVector;
  begin
    result:=e*(Dot(e,a)/Dot(e,e));
  end;

var
  um: TMatrix;
  u, e, a: TVector;
  i, i2: longint;
  t: double;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in QRDecomposeGS');

  um:=TMatrix.Create(AMatrix.Rows, AMatrix.Columns);

  e:=Zeros(AMatrix.Columns);
  for i:=0 to AMatrix.Columns-1 do
    begin
      a:=AMatrix.Column[i];
      u:=a;
      for i2:=0 to i-1 do
        u:=u-proj(um.Column[i2],a);
      um.Column[i]:=u;
      t:=system.sqrt(dot(u,u));
      assert(t<>0, 'Matrix is singular');
      e[i]:=1/t;
    end;

  AQ:=TMatrix.Create(AMatrix.Rows, AMatrix.Columns);
  for i:=0 to AMatrix.Rows-1 do
    for i2:=0 to AMatrix.Columns-1 do
      AQ[i,i2]:=um[i,i2]*e[i2];

  AR:=AQ.Transpose*AMatrix;
end;

function Sgn(d: double): double;
begin
  if d<0 then
    exit(-1)
  else
    exit(1);
end;

procedure QRDecomposeHR(const AMatrix: TMatrix; out AQ, AR: TMatrix);
var
  x, u, v: TVector;
  t, alpha: Double;
  Qn, A: TMatrix;
  i: longint;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in QRDecomposeHR');

  AQ:=TMatrix.Identity(AMatrix.Rows);
  A:=AMatrix;

  for i:=0 to AMatrix.Rows-2 do
    begin
      x:=A.MinorColumn[i,i,0];
      alpha:=-Sgn(x[1])*system.sqrt(dot(x,x));
      u:=Zeros(length(x));
      u[0]:=1;
      u:=x-alpha*u;

      t:=system.sqrt(dot(u,u));
      if t=0 then continue;
      //if dot(x,x)=0 then continue;
      //assert(t<>0, 'Matrix is singular');
      v:=u/t;

      Qn:=TMatrix.Identity(A.Rows-I) - (TMatrix.Create(length(v),1, 2*v)*TMatrix.Create(1,length(v),v));
      Qn:=Qn.Extend(i,i,0,0, true);

      AQ:=AQ*Qn.Transpose;
      A:=Qn*A;
    end;

  AR:=AQ.Transpose*AMatrix;
end;

function ReduceHR(const AMatrix: TMatrix): TMatrix;
var
  x, v: TVector;
  alpha: Double;
  Qn, A: TMatrix;
  i, i2, n: longint;
  r, rr: double;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in ReduceHR');

  n:=AMatrix.Columns;

  A:=AMatrix;

  for i:=0 to n-2 do
    begin
      x:=A.MinorColumn[i,i,0];
      x[0]:=0;
      alpha:=-Sgn(x[1])*system.sqrt(dot(x,x));

      r:=system.sqrt(0.5*(system.sqr(alpha)-x[1]*alpha));
      if r=0 then continue;
      rr:=1/(2*r);

      v:=zeros(length(x));
      v[0]:=0;
      v[1]:=(x[1]-alpha)*rr;

      for i2:=2 to high(x) do
        v[i2]:=x[i2]*rr;

      Qn:=TMatrix.Identity(A.Rows-I) - (TMatrix.Create(length(v),1, 2*v)*TMatrix.Create(1,length(v),v));
      Qn:=Qn.Extend(i,i,0,0, true);

      A:=Qn*A*Qn;
    end;

  result:=A;
end;

function CholeskyDecompose(const AMatrix: TMatrix): TMatrix;
var
  Am: TMatrix;
  aii,t: Double;
  i, i2: longint;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in CholeskyDecompose');

  am:=AMatrix.Copy;
  result:=TMatrix.Identity(AMatrix.Rows);

  for i:=0 to result.Rows-1 do
    begin
      t:=am[i,i];

      assert(t<>0, 'Matrix is not invertible');
      assert(t>0, 'Resulting matrix is not real');

      aii:=system.sqrt(t);
      result[i,i]:=aii;
      aii:=1/aii;

      am.ScaleRow(i,system.sqr(aii));

      for i2:=i+1 to result.rows-1 do
        begin
          t:=am[i2,i];
          am.SubScaled(i,i2,t);
          result[i2,i]:=t*aii;
        end;
    end;
end;

procedure ArnoldiIteration(const A: TMatrix; out H, Q: TMatrix);
var
  n, j: longint;
  b, v: TVector;
begin
  b:=Rand(a.Rows);

  h:=TMatrix.Identity(a.Rows);
  q:=TMatrix.Identity(a.Rows);
  q.Column[0]:=b/Magnitude(b);

  for n:=1 to a.Rows-1 do
    begin
      v:=(A*TMatrix(q.Column[n-1])).data;

      for j:=0 to n-1 do
        begin
          h[j,n-1]:=Dot(q.Column[j],q.Column[n]);
          v:=v-h[j,n]*q.Column[j];
        end;

      h[n,n-1]:=Magnitude(v);
      q.Column[n]:=v/h[n,n-1];
    end;
end;

function InvertGJ(const A: TMatrix): TMatrix;
var
  t: TMatrix;
  n, i, i2: LongInt;
  x: Double;
begin
  assert(a.Columns = a.Rows, 'Only square matrices supported in InvertGJ');

  n:=A.rows;

  t:=a.Copy;
  result:=TMatrix.Identity(n);

  for i:=0 to n-1 do
    begin
      x:=t[i,i];
      result.ScaleRow(i,x);
      t.ScaleRow(i,x);

      for i2:=i+1 to n-1 do
        begin
          x:=t[i2,i];
          result.SubScaled(i,i2,x);
          t.SubScaled(i,i2,x);
        end;
    end;

  for i:=n-1 downto 1 do
    begin
      for i2:=i-1 downto 0 do
        begin
          x:=t[i2,i];
          result.SubScaled(i,i2,x);
          t.SubScaled(i,i2,x);
        end;
    end;
end;

function SolveLU(const AL, AU, Ab: TMatrix): TMatrix;
var
  y: TMatrix;
  n, c, i, i2: LongInt;
  t: double;
begin
  n:=AL.Rows;

  result:=TMatrix.Create(n,Ab.Columns);

  y:=TMatrix.Create(n,Ab.Columns);
  for c:=0 to ab.Columns-1 do
    begin
      for i:=0 to n-1 do
        begin
          t:=0;

          for i2:=0 to i-1 do
            t:=t+y[i2,0]*al[i,i2];

          y[i,0]:=(ab[i,c]-t);
        end;

      for i:=n-1 downto 0 do
        begin
          t:=0;

          for i2:=n-1 downto i+1 do
            t:=t+result[i2,c]*au[i,i2];

          result[i,c]:=(y[i,0]-t)/au[i,i];
        end;
    end;
end;

function InvertLU(const AL, AU: TMatrix): TMatrix;
var
  i,i2,n: LongInt;
  t: double;
  Y, b: TMatrix;
begin
  assert(al.Columns = al.Rows, 'Only square matrices supported in InvertLU');

  n:=al.rows;

  b:=TMatrix.Identity(n);

  result:=SolveLU(al, au, b);
end;

function Determinant(const AMatrix: TMatrix): double;
var
  i,i2: longint;
  t: double;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in Determinant');

  result:=0;
  for i:=0 to AMatrix.Columns-1 do
    begin
      t:=1;
      for i2:=0 to AMatrix.Rows-1 do
        t:=t*AMatrix[i2,(i+i2) mod AMatrix.Columns];
      result:=result+t;

      t:=1;
      for i2:=0 to AMatrix.Rows-1 do
        t:=t*AMatrix[AMatrix.Rows-i2-1, (i+i2) mod AMatrix.Columns];
      result:=result-t;
    end;
end;

function DeterminantUnitary(const AMatrix: TMatrix): double;
begin
  result:=1;
end;

function DeterminantTriangular(const AMatrix: TMatrix): double;
var
  i: longint;
begin
  result:=1;
  for i:=0 to min(AMatrix.Rows,AMatrix.Columns)-1 do
    result:=result*AMatrix[i,i];
end;

function TMatrix.GetColumn(AColumn: longint): TVector;
var
  pd: PDouble;
  i: longint;
begin
  Assert((AColumn<Columns) and (AColumn>=0), 'Invalid column index');

  setlength(Result,Rows);

  pd:=@Data[AColumn];

  for i:=0 to rows-1 do
    begin
      result[i]:=pd^;
      inc(pd,Columns);
    end;
end;

function TMatrix.GetCopy: TMatrix;
begin
  result:=TMatrix.Create(Rows,Columns,system.Copy(Data));
end;

function TMatrix.GetHermitian: TMatrix;
var
  i, i2: longint;
begin
{$ifdef IsComplex}
  result:=TMatrix.Create(Columns,Rows);
  for i:=0 to rows-1 do
    for i2:=0 to Columns-1 do
      result[i2,i]:=cong(self[i,i2]);
{$else}
  result:=Transpose;
{$endif}
end;

function TMatrix.GetMinorColumn(AStartRow, AStartColumn, AColumn: longint): TVector;
var
  pd: PDouble;
  i: longint;
begin
  Assert((AStartRow<Rows) and (AStartRow>=0), 'Invalid start row index');
  Assert((AStartColumn<Columns) and (AStartColumn>=0), 'Invalid start column index');

  Assert((AColumn<(Columns-AStartColumn)) and (AColumn>=0), 'Invalid column index');

  setlength(Result,Rows-AStartRow);

  pd:=@Data[AStartRow*Columns+AStartColumn+AColumn];

  for i:=0 to rows-AStartRow-1 do
    begin
      result[i]:=pd^;
      inc(pd,Columns);
    end;
end;

function TMatrix.GetRow(ARow: longint): TVector;
begin
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  setlength(Result,Columns);
  move(Data[Columns*ARow],result[0], Columns*sizeof(double));
end;

function TMatrix.GetTranspose: TMatrix;
var
  i, i2: longint;
begin
  result:=TMatrix.Create(Columns,Rows);
  for i:=0 to rows-1 do
    for i2:=0 to Columns-1 do
      result[i2,i]:=self[i,i2];
end;

function TMatrix.GetValue(ARow, AColumn: longint): double;
begin
  Assert((AColumn<Columns) and (AColumn>=0), 'Invalid column index');
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  result:=Data[AColumn+ARow*Columns];
end;

procedure TMatrix.SetColumn(AColumn: longint; AValue: TVector);
var
  pd: PDouble;
  i: longint;
begin
  Assert((AColumn<Columns) and (AColumn>=0), 'Invalid column index');

  pd:=@Data[AColumn];

  for i:=0 to rows-1 do
    begin
      pd^:=AValue[i];
      inc(pd,Columns);
    end;
end;

procedure TMatrix.SetRow(ARow: longint; AValue: TVector);
begin
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  move(AValue[0], Data[Columns*ARow], Columns*sizeof(double));
end;

procedure TMatrix.SetValue(ARow, AColumn: longint; AValue: double);
begin
  Assert((AColumn<Columns) and (AColumn>=0), 'Invalid column index');
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  Data[AColumn+ARow*Columns]:=AValue;
end;

function TMatrix.Submatrix(AStartRow, AStartColumn: longint): TMatrix;
var
  i, i2: longint;
begin
  Assert((AStartColumn<Columns) and (AStartColumn>=0), 'Invalid column index');
  Assert((AStartRow<Rows) and (AStartRow>=0), 'Invalid row index');

  result:=TMatrix.Create(Rows-AStartRow, Columns-AStartColumn);
  for i:=0 to result.Rows-1 do
    for i2:=0 to result.Columns-1 do
      result[i,i2]:=self[i+AStartRow,i2+AStartColumn];
end;

function TMatrix.Submatrix(AStartRow, AStartColumn, ARows, AColumns: longint): TMatrix;
var
  i, i2: longint;
begin
  Assert((AStartColumn<Columns) and (AStartColumn>=0), 'Invalid column index');
  Assert((AStartRow<Rows) and (AStartRow>=0), 'Invalid row index');

  Assert((AStartColumn+AColumns)<=Columns, 'Invalid number of columns');
  Assert((AStartRow+ARows)<=Rows, 'Invalid number of rows');

  result:=TMatrix.Create(ARows, AColumns);
  for i:=0 to result.Rows-1 do
    for i2:=0 to result.Columns-1 do
      result[i,i2]:=self[i+AStartRow,i2+AStartColumn];
end;

function TMatrix.Extend(ARowS, AColumnS, ARowE, AColumnE: longint; AIdent: boolean): TMatrix;
var
  i, i2: longint;
begin
  if AIdent then
    begin
      result:=TMatrix.Identity(ARowS+rows+ARowE);
    end
  else
    begin
      result:=TMatrix.Create(Rows+ARowS+ARowE, Columns+AColumnS+AColumnE);
      FillChar(result.Data[0], sizeof(double)*result.Rows*result.Columns, 0);
    end;

  for i:=0 to Rows-1 do
    for i2:=0 to Columns-1 do
      result[i+ARowS,i2+AColumnS]:=self[i,i2];
end;

constructor TMatrix.Identity(ARows, AColumns: longint);
var
  pd:    PDouble;
  i, i2: longint;
begin
  setlength(Data, ARows*AColumns);
  Rows := ARows;
  Columns := AColumns;

  pd := @Data[0];

  for i := 0 to ARows - 1 do
    for i2 := 0 to AColumns - 1 do
    begin
      if i = i2 then
        pd^ := 1
      else
        pd^ := 0;
      Inc(pd);
    end;
end;

constructor TMatrix.Identity(N: longint);
var
  pd:    PDouble;
  i, i2: longint;
begin
  setlength(Data, n * n);
  Rows := n;
  Columns := n;

  pd := @Data[0];

  for i := 0 to n - 1 do
    for i2 := 0 to n - 1 do
    begin
      if i = i2 then
        pd^ := 1
      else
        pd^ := 0;
      Inc(pd);
    end;
end;

constructor TMatrix.Create(ARows, AColumns: longint);
begin
  Data:=Zeros(ARows*AColumns);
  Rows := ARows;
  Columns := AColumns;
end;

constructor TMatrix.Create(ARows, AColumns: longint; const AValues: array of double);
begin
  Rows := ARows;
  Columns := AColumns;

  setlength(Data, ARows * AColumns);
  move(AValues[0], Data[0], length(AValues) * sizeof(double));
end;

procedure TMatrix.SubScaled(ASrc, ADst: longint; AValue: double);
var
  pd, ps: PDouble;
  i:  longint;
begin
  ps := @Data[ASrc*Columns];
  pd := @Data[ADst*Columns];
  for i := 0 to Columns-1 do
  begin
    pd^ := pd^-ps^ * AValue;
    Inc(pd);
    Inc(ps);
  end;
end;

procedure TMatrix.ScaleRow(ARow: longint; AValue: double);
var
  pd: PDouble;
  i:  longint;
begin
  pd := @Data[ARow*Columns];
  for i := 0 to Columns-1 do
  begin
    pd^ := pd^ * AValue;
    Inc(pd);
  end;
end;

procedure TMatrix.Scale(AValue: double);
var
  pd: PDouble;
  i:  longint;
begin
  pd := @Data[0];
  for i := 0 to high(Data) do
  begin
    pd^ := pd^ * AValue;
    Inc(pd);
  end;
end;

function TMatrix.Dump: string;
var
  m, i, i2: longint;
  x:        double;
begin
  x := max(Data);
  if x<0 then
    m := trunc(log10(-x)) + 6
  else if x<>0 then
    m := trunc(log10(x)) + 5
  else
    m := 1;

  Result := '';
  for i2 := 0 to Rows - 1 do
  begin
    Result := Result + format(format('%%%D.2F', [m]), [Data[i2 * Columns + 0]]);
    for i := 1 to Columns - 1 do
      Result := Result + format(format(', %%%D.2F', [m]), [Data[i2 * Columns + i]]);
    if i2 <> (rows - 1) then
      Result := Result + lineending;
  end;
end;

end.
