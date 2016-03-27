{$ifdef IsInterface}
record
  strict private
    function GetColumn(AColumn: longint): TGVector;
    function GetCopy: TGMatrix;
    function GetHermitian: TGMatrix;
    function GetMinorColumn(AStartRow, AStartColumn, AColumn: longint): TGVector;
    function GetRow(ARow: longint): TGVector;
    function GetTranspose: TGMatrix;
    function GetValue(ARow, AColumn: longint): TFloatType;
    procedure SetColumn(AColumn: longint; AValue: TGVector);
    procedure SetRow(ARow: longint; AValue: TGVector);
    procedure SetValue(ARow, AColumn: longint; AValue: TFloatType);
  public
    Data: array of TFloatType;
    Rows, Columns: longint;

    function Submatrix(AStartRow,AStartColumn: longint): TGMatrix;
    function Submatrix(AStartRow,AStartColumn,ARows,AColumns: longint): TGMatrix;
    function Extend(ARowS,AColumnS, ARowE,AColumnE: longint; AIdent: boolean): TGMatrix;

    constructor Identity(ARows, AColumns: longint);
    constructor Identity(N: longint);
    constructor Create(ARows, AColumns: longint);
    constructor Create(ARows, AColumns: longint; const AValues: array of TFloatType);

    procedure SubScaled(ASrc, ADst: longint; AValue: TFloatType);
    procedure ScaleRow(ARow: longint; AValue: TFloatType);
    procedure Scale(AValue: TFloatType);

    function Dump: string;

    property Transpose: TGMatrix read GetTranspose;
    property Hermitian: TGMatrix read GetHermitian;

    property Copy: TGMatrix read GetCopy;

    property Value[ARow, AColumn: longint]: TFloatType read GetValue write SetValue; default;

    property Row[ARow: longint]: TGVector read GetRow write SetRow;
    property Column[AColumn: longint]: TGVector read GetColumn write SetColumn;
    property MinorColumn[AStartRow,AStartColumn,AColumn: longint]: TGVector read GetMinorColumn;
  end;

operator explicit(const A: TGVector): TGMatrix;
operator explicit(const A: TGMatrix): TGVector;

operator -(const A: TGMatrix): TGMatrix;

operator +(const A, B: TGMatrix): TGMatrix;
operator -(const A, B: TGMatrix): TGMatrix;
operator * (const A, B: TGMatrix): TGMatrix;
operator * (const A: TGMatrix; const B: TFloatType): TGMatrix;
operator * (const B: TFloatType; const A: TGMatrix): TGMatrix;

procedure LUDecompose(const AMatrix: TGMatrix; out AL, AU: TGMatrix);

procedure QRDecomposeGS(const AMatrix: TGMatrix; out AQ, AR: TGMatrix); // Gram-Schmidt
procedure QRDecomposeHR(const AMatrix: TGMatrix; out AQ, AR: TGMatrix); // Householder reflection

function ReduceHR(const AMatrix: TGMatrix): TGMatrix; // Reduction to tridiagonal by Householder reflection

function CholeskyDecompose(const AMatrix: TGMatrix): TGMatrix; // Returns L where AMatrix=L*L.Transpose

procedure ArnoldiIteration(const A: TGMatrix; out H, Q: TGMatrix);

function SolveLU(const AL, AU, Ab: TGMatrix): TGMatrix;

function InvertGJ(const A: TGMatrix): TGMatrix;
function InvertLU(const AL, AU: TGMatrix): TGMatrix;
function InvertLT(const AL: TGMatrix): TGMatrix; // Invert lower triangular matrix
function InvertUT(const AU: TGMatrix): TGMatrix; // Invert upper triangular matrix

function Determinant(const AMatrix: TGMatrix): TFloatType;
function DeterminantUnitary(const AMatrix: TGMatrix): TFloatType;
function DeterminantTriangular(const AMatrix: TGMatrix): TFloatType;
{$else}
operator explicit(const A: TGVector): TGMatrix;
begin
  result:=TGMatrix.Create(length(a), 1, a);
end;

operator explicit(const A: TGMatrix): TGVector;
begin
  Assert((a.Rows = 1) or (a.Columns = 1));

  Result := copy(a.Data);
end;

operator-(const A: TGMatrix): TGMatrix;
var
  pa, pr: ^TFloatType;
  i: longint;
begin
  pa := @a.Data[0];

  Result := TGMatrix.Create(a.Rows, a.Columns);

  pr := @Result.Data[0];

  for i := 0 to high(a.Data) do
  begin
    pr^ := -pa^;
    Inc(pa);
    Inc(pr);
  end;
end;

operator +(const A, B: TGMatrix): TGMatrix;
var
  pa, pb, pr: ^TFloatType;
  i: longint;
begin
  Assert(a.Rows = b.Rows);
  Assert(a.Columns = b.Columns);

  pa := @a.Data[0];
  pb := @b.Data[0];

  Result := TGMatrix.Create(a.Rows, a.Columns);

  pr := @Result.Data[0];

  for i := 0 to high(a.Data) do
  begin
    pr^ := pa^ + pb^;
    Inc(pa);
    Inc(pb);
    Inc(pr);
  end;
end;

operator -(const A, B: TGMatrix): TGMatrix;
var
  pa, pb, pr: ^TFloatType;
  i: longint;
begin
  Assert(a.Rows = b.Rows);
  Assert(a.Columns = b.Columns);

  pa := @a.Data[0];
  pb := @b.Data[0];

  Result := TGMatrix.Create(a.Rows, a.Columns);

  pr := @Result.Data[0];

  for i := 0 to high(a.Data) do
  begin
    pr^ := pa^ - pb^;
    Inc(pa);
    Inc(pb);
    Inc(pr);
  end;
end;

operator * (const A, B: TGMatrix): TGMatrix;
var
  pa, pb, pr: ^TFloatType;
  i, i2, i3:  longint;
  sum:        TFloatType;
begin
  Assert(a.Columns = b.Rows);

  Result := TGMatrix.Create(a.Rows, b.Columns);

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

operator*(const A: TGMatrix; const B: TFloatType): TGMatrix;
begin
  result:=A.Copy;
  result.Scale(b);
end;

operator*(const B: TFloatType; const A: TGMatrix): TGMatrix;
begin
  result:=A.Copy;
  result.Scale(b);
end;

procedure LUDecompose(const AMatrix: TGMatrix; out AL, AU: TGMatrix);
var
  i, i2: longint;
  t: TFloatType;
  Ln: TGMatrix;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in LUDecompose');

  AL := TGMatrix.Identity(AMatrix.Rows);
  AU:=AMatrix;

  for i:=0 to AMatrix.Rows-2 do
    begin
      Ln:=TGMatrix.Identity(AMatrix.Rows);

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

procedure QRDecomposeGS(const AMatrix: TGMatrix; out AQ, AR: TGMatrix);

  function Proj(e,a: TGVector): TGVector;
  var
    ec: tgvector;
  begin
{$ifdef IsComplex}
    ec:=Conjugate(e);
    result:=e*(dot(ec,a)/Dot(ec,e));
{$else}
    result:=e*(Dot(e,a)/Dot(e,e));
{$endif IsComplex}
  end;

var
  um: TGMatrix;
  u, e, a: TGVector;
  i, i2: longint;
  t: TFloatType;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in QRDecomposeGS');

  um:=TGMatrix.Create(AMatrix.Rows, AMatrix.Columns);

  e:=(VZeros(AMatrix.Columns));
  for i:=0 to AMatrix.Columns-1 do
    begin
      a:=AMatrix.Column[i];
      u:=a;
      for i2:=0 to i-1 do
        u:=u-proj(um.Column[i2],a);
      um.Column[i]:=u;
      t:=NormEuclid(u);
      assert(t<>0, 'Matrix is singular');
      e[i]:=1/t;
    end;

  AQ:=TGMatrix.Create(AMatrix.Rows, AMatrix.Columns);
  for i:=0 to AMatrix.Rows-1 do
    for i2:=0 to AMatrix.Columns-1 do
      AQ[i,i2]:=um[i,i2]*e[i2];

  AR:=AQ.Transpose*AMatrix;
end;

{$ifndef IsComplex}
function Sgn(d: TFloatType): TFloatType;
begin
  if d<0 then
    exit(-1)
  else
    exit(1);
end;
{$endif IsComplex}

procedure QRDecomposeHR(const AMatrix: TGMatrix; out AQ, AR: TGMatrix);
var
  x, u, v: TGVector;
  t, alpha: TFloatType;
  Qn, A: TGMatrix;
  i: longint;
  w: TFloatType;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in QRDecomposeHR');

  AQ:=TGMatrix.Identity(AMatrix.Rows);
  A:=AMatrix;

  for i:=0 to AMatrix.Rows-2 do
    begin
      x:=A.MinorColumn[i,i,0];
{$ifdef IsComplex}
      alpha:=-cexp(ucomplex.i*carg(x[1]))*NormEuclid(x);
{$else}
      alpha:=-Sgn(x[1])*NormEuclid(x);
{$endif IsComplex}
      u:=(VZeros(length(x)));
      u[0]:=1;
      u:=x-alpha*u;

      t:=NormEuclid(u);
      if t=0 then continue;
      v:=u/t;

{$ifdef IsComplex}
      w:=Dot(Conjugate(x), v)/Dot(Conjugate(v), x);
      Qn:=TGMatrix.Identity(A.Rows-I) - (1+w)*(TGMatrix.Create(length(v),1, v)*TGMatrix.Create(length(v),1,v).hermitian);
{$else}
      Qn:=TGMatrix.Identity(A.Rows-I) - (TGMatrix.Create(length(v),1, 2*v)*TGMatrix.Create(1,length(v),v));
{$endif IsComplex}
      Qn:=Qn.Extend(i,i,0,0, true);

      AQ:=AQ*Qn.Transpose;
      A:=Qn*A;
    end;

  AR:=AQ.Transpose*AMatrix;
end;

function ReduceHR(const AMatrix: TGMatrix): TGMatrix;
var
  x, v: TGVector;
  alpha: TFloatType;
  Qn, A: TGMatrix;
  i, k, i2, n: longint;
  r, rr: TFloatType;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in ReduceHR');

  n:=AMatrix.Columns;

  A:=AMatrix;

  for i:=0 to n-2 do
    begin
      x:=A.MinorColumn[i,i,0];
      x[0]:=0;
{$ifdef IsComplex}
      k:=1;
      alpha:=-cexp(ucomplex.i*carg(x[k]))*NormEuclid(x);
{$else}
      alpha:=-Sgn(x[1])*NormEuclid(x);
{$endif IsComplex}

      r:=vsqrt(0.5*(vsqr(alpha)-x[1]*alpha));
      if r=0 then continue;
      rr:=1/(2*r);

      v:=(Vzeros(length(x)));
      v[0]:=0;
      v[1]:=(x[1]-alpha)*rr;

      for i2:=2 to high(x) do
        v[i2]:=x[i2]*rr;

      Qn:=TGMatrix.Identity(A.Rows-I) - (TGMatrix.Create(length(v),1, 2*v)*TGMatrix.Create(1,length(v),v));
      Qn:=Qn.Extend(i,i,0,0, true);

      A:=Qn*A*Qn;
    end;

  result:=A;
end;

function CholeskyDecompose(const AMatrix: TGMatrix): TGMatrix;
var
  Am: TGMatrix;
  aii,t: TFloatType;
  i, i2: longint;
begin
  assert(AMatrix.Columns = AMatrix.Rows, 'Only square matrices supported in CholeskyDecompose');

  am:=AMatrix.Copy;
  result:=TGMatrix.Identity(AMatrix.Rows);

  for i:=0 to result.Rows-1 do
    begin
      t:=am[i,i];

      assert(t<>0, 'Matrix is not invertible');
{$ifndef IsComplex}
      assert(t>0, 'Resulting matrix is not real');
{$endif IsComplex}

      aii:=vsqrt(t);
      result[i,i]:=aii;
      aii:=1/aii;

      am.ScaleRow(i,vsqr(aii));

      for i2:=i+1 to result.rows-1 do
        begin
          t:=am[i2,i];
          am.SubScaled(i,i2,t);
          result[i2,i]:=t*aii;
        end;
    end;
end;

procedure ArnoldiIteration(const A: TGMatrix; out H, Q: TGMatrix);
var
  n, j: longint;
  b, v: TGVector;
begin
  b:=TGVector(Rand(a.Rows));

  h:=TGMatrix.Identity(a.Rows);
  q:=TGMatrix.Identity(a.Rows);
  q.Column[0]:=b/Magnitude(b);

  for n:=1 to a.Rows-1 do
    begin
      v:=(A*TGMatrix(q.Column[n-1])).data;

      for j:=0 to n-1 do
        begin
          h[j,n-1]:=Dot(q.Column[j],q.Column[n]);
          v:=v-h[j,n]*q.Column[j];
        end;

      h[n,n-1]:=Magnitude(v);
      q.Column[n]:=v/h[n,n-1];
    end;
end;

function InvertGJ(const A: TGMatrix): TGMatrix;
var
  t: TGMatrix;
  n, i, i2: LongInt;
  x: TFloatType;
begin
  assert(a.Columns = a.Rows, 'Only square matrices supported in InvertGJ');

  n:=A.rows;

  t:=a.Copy;
  result:=TGMatrix.Identity(n);

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

function SolveLU(const AL, AU, Ab: TGMatrix): TGMatrix;
var
  y: TGMatrix;
  n, c, i, i2: LongInt;
  t: TFloatType;
begin
  n:=AL.Rows;

  result:=TGMatrix.Create(n,Ab.Columns);

  y:=TGMatrix.Create(n,Ab.Columns);
  for c:=0 to ab.Columns-1 do
    begin
      for i:=0 to n-1 do
        begin
          t:=0;

          for i2:=0 to i-1 do
            t:=t+y[i2,0]*al[i,i2];

          y[i,0]:=(ab[i,c]-t)/al[i,i];
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

function InvertLU(const AL, AU: TGMatrix): TGMatrix;
var
  i,i2,n: LongInt;
  t: TFloatType;
  Y, b: TGMatrix;
begin
  assert(al.Columns = al.Rows, 'Only square matrices supported in InvertLU');

  n:=al.rows;

  b:=TGMatrix.Identity(n);

  result:=SolveLU(al, au, b);
end;

function InvertLT(const AL: TGMatrix): TGMatrix;
var
  i,i2,c,n: LongInt;
  t: TFloatType;
  Y, b: TGMatrix;
begin
  assert(al.Columns = al.Rows, 'Only square matrices supported in InvertLU');

  n:=al.rows;

  b:=TGMatrix.Identity(n);

  n:=AL.Rows;

  result:=TGMatrix.Create(n,b.Columns);

  y:=TGMatrix.Create(n,b.Columns);
  for c:=0 to b.Columns-1 do
    begin
      for i:=0 to n-1 do
        begin
          t:=0;

          for i2:=0 to i-1 do
            t:=t+y[i2,0]*al[i,i2];

          y[i,0]:=(b[i,c]-t)/al[i,i];
        end;
    end;
end;

function InvertUT(const AU: TGMatrix): TGMatrix;
var
  i,i2,c,n: LongInt;
  t: TFloatType;
  Y, b: TGMatrix;
begin
  assert(AU.Columns = AU.Rows, 'Only square matrices supported in InvertLU');

  n:=AU.rows;

  b:=TGMatrix.Identity(n);

  n:=AU.Rows;

  result:=TGMatrix.Create(n,b.Columns);

  y:=TGMatrix.Create(n,b.Columns);
  for c:=0 to b.Columns-1 do
    begin
      for i:=n-1 downto 0 do
        begin
          t:=0;

          for i2:=n-1 downto i+1 do
            t:=t+result[i2,c]*au[i,i2];

          result[i,c]:=(y[i,0]-t)/au[i,i];
        end;
    end;
end;

function Determinant(const AMatrix: TGMatrix): TFloatType;
var
  i,i2: longint;
  t: TFloatType;
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

function DeterminantUnitary(const AMatrix: TGMatrix): TFloatType;
begin
  result:=1;
end;

function DeterminantTriangular(const AMatrix: TGMatrix): TFloatType;
var
  i: longint;
begin
  result:=1;
  for i:=0 to min(AMatrix.Rows,AMatrix.Columns)-1 do
    result:=result*AMatrix[i,i];
end;

function TGMatrix.GetColumn(AColumn: longint): TGVector;
var
  pd: ^TFloatType;
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

function TGMatrix.GetCopy: TGMatrix;
begin
  result:=TGMatrix.Create(Rows,Columns,system.Copy(Data));
end;

function TGMatrix.GetHermitian: TGMatrix;
var
  i, i2: longint;
  pd: ^TFloatType;
begin
{$ifdef IsComplex}
  result:=TGMatrix.Create(Columns,Rows);
  pd:=@data[0];
  for i:=0 to rows-1 do
    for i2:=0 to Columns-1 do
      begin
        result[i2,i]:=cong(pd^);
        inc(pd);
      end;
{$else}
  result:=Transpose;
{$endif}
end;

function TGMatrix.GetMinorColumn(AStartRow, AStartColumn, AColumn: longint): TGVector;
var
  pd: ^TFloatType;
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

function TGMatrix.GetRow(ARow: longint): TGVector;
begin
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  setlength(Result,Columns);
  move(Data[Columns*ARow],result[0], Columns*sizeof(TFloatType));
end;

function TGMatrix.GetTranspose: TGMatrix;
var
  i, i2: longint;
begin
  result:=TGMatrix.Create(Columns,Rows);
  for i:=0 to rows-1 do
    for i2:=0 to Columns-1 do
      result[i2,i]:=self[i,i2];
end;

function TGMatrix.GetValue(ARow, AColumn: longint): TFloatType;
begin
  Assert((AColumn<Columns) and (AColumn>=0), 'Invalid column index');
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  result:=Data[AColumn+ARow*Columns];
end;

procedure TGMatrix.SetColumn(AColumn: longint; AValue: TGVector);
var
  pd: ^TFloatType;
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

procedure TGMatrix.SetRow(ARow: longint; AValue: TGVector);
begin
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  move(AValue[0], Data[Columns*ARow], Columns*sizeof(TFloatType));
end;

procedure TGMatrix.SetValue(ARow, AColumn: longint; AValue: TFloatType);
begin
  Assert((AColumn<Columns) and (AColumn>=0), 'Invalid column index');
  Assert((ARow<Rows) and (ARow>=0), 'Invalid row index');

  Data[AColumn+ARow*Columns]:=AValue;
end;

function TGMatrix.Submatrix(AStartRow, AStartColumn: longint): TGMatrix;
var
  i, i2: longint;
begin
  Assert((AStartColumn<Columns) and (AStartColumn>=0), 'Invalid column index');
  Assert((AStartRow<Rows) and (AStartRow>=0), 'Invalid row index');

  result:=TGMatrix.Create(Rows-AStartRow, Columns-AStartColumn);
  for i:=0 to result.Rows-1 do
    for i2:=0 to result.Columns-1 do
      result[i,i2]:=self[i+AStartRow,i2+AStartColumn];
end;

function TGMatrix.Submatrix(AStartRow, AStartColumn, ARows, AColumns: longint): TGMatrix;
var
  i, i2: longint;
begin
  Assert((AStartColumn<Columns) and (AStartColumn>=0), 'Invalid column index');
  Assert((AStartRow<Rows) and (AStartRow>=0), 'Invalid row index');

  Assert((AStartColumn+AColumns)<=Columns, 'Invalid number of columns');
  Assert((AStartRow+ARows)<=Rows, 'Invalid number of rows');

  result:=TGMatrix.Create(ARows, AColumns);
  for i:=0 to result.Rows-1 do
    for i2:=0 to result.Columns-1 do
      result[i,i2]:=self[i+AStartRow,i2+AStartColumn];
end;

function TGMatrix.Extend(ARowS, AColumnS, ARowE, AColumnE: longint; AIdent: boolean): TGMatrix;
var
  i, i2: longint;
begin
  if AIdent then
    begin
      result:=TGMatrix.Identity(ARowS+rows+ARowE);
    end
  else
    begin
      result:=TGMatrix.Create(Rows+ARowS+ARowE, Columns+AColumnS+AColumnE);
      FillChar(result.Data[0], sizeof(TFloatType)*result.Rows*result.Columns, 0);
    end;

  for i:=0 to Rows-1 do
    for i2:=0 to Columns-1 do
      result[i+ARowS,i2+AColumnS]:=self[i,i2];
end;

constructor TGMatrix.Identity(ARows, AColumns: longint);
var
  pd:    ^TFloatType;
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

constructor TGMatrix.Identity(N: longint);
var
  pd:    ^TFloatType;
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

constructor TGMatrix.Create(ARows, AColumns: longint);
begin
  Data:=(VZeros(ARows*AColumns));
  Rows := ARows;
  Columns := AColumns;
end;

constructor TGMatrix.Create(ARows, AColumns: longint; const AValues: array of TFloatType);
begin
  Rows := ARows;
  Columns := AColumns;

  setlength(Data, ARows * AColumns);
  move(AValues[0], Data[0], length(AValues) * sizeof(TFloatType));
end;

procedure TGMatrix.SubScaled(ASrc, ADst: longint; AValue: TFloatType);
var
  pd, ps: ^TFloatType;
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

procedure TGMatrix.ScaleRow(ARow: longint; AValue: TFloatType);
var
  pd: ^TFloatType;
  i:  longint;
begin
  pd := @Data[ARow*Columns];
  for i := 0 to Columns-1 do
  begin
    pd^ := pd^ * AValue;
    Inc(pd);
  end;
end;

procedure TGMatrix.Scale(AValue: TFloatType);
var
  pd: ^TFloatType;
  i:  longint;
begin
  pd := @Data[0];
  for i := 0 to high(Data) do
  begin
    pd^ := pd^ * AValue;
    Inc(pd);
  end;
end;

function TGMatrix.Dump: string;
var
  m, i, i2: longint;
  x:        TFloatType;
begin
{$ifndef IsComplex}
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
{$else}
  Result := '';
  for i2 := 0 to Rows - 1 do
  begin
    Result := Result + cstr(Data[i2 * Columns + 0], 5,2);
    for i := 1 to Columns - 1 do
      Result := Result + ', '+ cstr(Data[i2 * Columns + i], 5,2);
    if i2 <> (rows - 1) then
      Result := Result + lineending;
  end;
{$endif IsComplex}
end;
{$endif}
