function Numeric()

    fig = uifigure('Name','NumericLab Solver','Position',[80 60 840 660],...
        'Color',[0.95 0.95 0.97]);

    uilabel(fig,'Text','NumericLab Solver',...
        'Position',[20 625 420 30],'FontSize',19,'FontWeight','bold',...
        'FontColor',[0.10 0.40 0.70]);
    uilabel(fig,'Text','Numerical Methods & Matrix Algebra  —  MATLAB GUI App',...
        'Position',[20 603 600 18],'FontSize',11,'FontColor',[0.45 0.45 0.45]);

    tg = uitabgroup(fig,'Position',[10 10 820 585]);
    t1 = uitab(tg,'Title','   Root Finding   ');
    t2 = uitab(tg,'Title','   Matrix Operations   ');

    buildRootTab(t1);
    buildMatrixTab(t2);
end

%%==========================================================================
%% ROOT FINDING TAB
%%==========================================================================
function buildRootTab(t)
    uilabel(t,'Text','Equation  f(x) = 0  — use x as variable',...
        'Position',[14 530 500 18],'FontSize',11,'FontColor',[0.3 0.3 0.3]);
    fxField = uieditfield(t,'text','Position',[14 505 440 28],'FontSize',12,...
        'Placeholder','e.g.   x^3 - x - 2     or     sin(x) - x/2');

    uilabel(t,'Text','Method:','Position',[465 513 58 18],'FontSize',11);
    methodDrop = uidropdown(t,...
        'Items',{'Incremental','Bisection','Regula-Falsi','Newton-Raphson','Secant'},...
        'Position',[525 505 175 28],'FontSize',11);

    % Parameter fields
    aLabel = uilabel(t,'Text','a (lower):','Position',[14 472 72 18],'FontSize',11);
    aField = uieditfield(t,'numeric','Position',[90 468 70 26],'Value',-3);
    bLabel = uilabel(t,'Text','b (upper):','Position',[170 472 72 18],'FontSize',11);
    bField = uieditfield(t,'numeric','Position',[246 468 70 26],'Value',3);
    x0Label= uilabel(t,'Text','x0 (init):','Position',[14 472 72 18],'FontSize',11,'Visible','off');
    x0Field= uieditfield(t,'numeric','Position',[90 468 70 26],'Value',1,'Visible','off');
    x1Label= uilabel(t,'Text','x1 (2nd):','Position',[170 472 72 18],'FontSize',11,'Visible','off');
    x1Field= uieditfield(t,'numeric','Position',[246 468 70 26],'Value',2,'Visible','off');

    uilabel(t,'Text','Tolerance:','Position',[330 472 72 18],'FontSize',11);
    tolField = uieditfield(t,'numeric','Position',[406 468 80 26],'Value',1e-4);
    uilabel(t,'Text','Max Iter:','Position',[496 472 66 18],'FontSize',11);
    maxField = uieditfield(t,'numeric','Position',[566 468 55 26],'Value',50);

    methodDrop.ValueChangedFcn = @(s,~) toggleParams(s.Value,...
        aLabel,aField,bLabel,bField,x0Label,x0Field,x1Label,x1Field);

    solveBtn = uibutton(t,'Text','SOLVE','Position',[636 465 150 32],...
        'FontSize',13,'FontWeight','bold',...
        'BackgroundColor',[0.10 0.42 0.72],'FontColor','white');

    resLabel = uilabel(t,'Text','Enter equation above then press SOLVE',...
        'Position',[14 438 790 20],'FontSize',11,'FontColor',[0.2 0.5 0.2]);

    tbl = uitable(t,'Position',[14 200 395 228],...
        'ColumnName',{'i','a','b','c','f(c)','|error|'},...
        'FontSize',10,'RowName',{});

    ax = uiaxes(t,'Position',[420 195 385 235]);
    title(ax,'Graph of f(x)'); xlabel(ax,'x'); ylabel(ax,'f(x)');
    grid(ax,'on'); ax.FontSize=10;

    uilabel(t,'Text',...
        'Syntax: x^n  sin(x)  cos(x)  tan(x)  exp(x)  log(x)  sqrt(x)  abs(x)  pi',...
        'Position',[14 178 760 18],'FontSize',10,'FontColor',[0.55 0.55 0.55]);

    solveBtn.ButtonPushedFcn = @(~,~) doSolve(fxField,methodDrop,...
        aField,bField,x0Field,x1Field,tolField,maxField,resLabel,tbl,ax);
end

function toggleParams(m,aL,aF,bL,bF,x0L,x0F,x1L,x1F)
    bracket = any(strcmp(m,{'Incremental','Bisection','Regula-Falsi'}));
    secant  = strcmp(m,'Secant');
    newton  = strcmp(m,'Newton-Raphson');
    aL.Visible  = onoff(bracket||secant);
    aF.Visible  = onoff(bracket||secant);
    bL.Visible  = onoff(bracket||secant);
    bF.Visible  = onoff(bracket||secant);
    x0L.Visible = onoff(newton);
    x0F.Visible = onoff(newton);
    x1L.Visible = 'off'; x1F.Visible = 'off';
    if secant, aL.Text='x0:'; bL.Text='x1:';
    else,      aL.Text='a (lower):'; bL.Text='b (upper):'; end
end

function doSolve(fxField,methodDrop,aField,bField,x0Field,~,tolField,maxField,resLabel,tbl,ax)
    raw = strtrim(fxField.Value);
    if isempty(raw)
        resLabel.Text='ERROR: Please enter an equation.'; resLabel.FontColor=[.8 .1 .1]; return
    end
    try
        f = str2func(['@(x)' prepExpr(raw)]); f(1.0);
    catch
        resLabel.Text='ERROR: Invalid equation syntax.'; resLabel.FontColor=[.8 .1 .1]; return
    end

    method = methodDrop.Value;
    tol    = max(tolField.Value, 1e-14);
    maxIt  = max(round(maxField.Value),1);
    root=NaN; conv=false; hdr={}; rows={};

    switch method
      case 'Incremental'
        a=aField.Value; b=bField.Value; dx=(b-a)/400;
        hdr={'i','x','f(x)','f_prev * f(x)'};
        prevF=f(a);
        for i=1:maxIt*8
            x=a+i*dx; if x>b, break; end
            fx=f(x);
            rows(end+1,:)={i,rnd(x),rnd(fx),rnd(prevF*fx)};
            if prevF*fx<0, root=x-dx/2; conv=true; break; end
            prevF=fx;
        end

      case 'Bisection'
        a=aField.Value; b=bField.Value;
        if f(a)*f(b)>=0
            resLabel.Text='ERROR: f(a)*f(b) must be negative (opposite signs).';
            resLabel.FontColor=[.8 .1 .1]; return
        end
        hdr={'i','a','b','c=(a+b)/2','f(c)','|b-a|'};
        for i=1:maxIt
            c=(a+b)/2; fc=f(c);
            rows(end+1,:)={i,rnd(a),rnd(b),rnd(c),rnd8(fc),rnd8(abs(b-a))};
            if abs(fc)<tol||abs(b-a)<tol, root=c; conv=true; break; end
            if f(a)*fc<0, b=c; else, a=c; end
        end

      case 'Regula-Falsi'
        a=aField.Value; b=bField.Value;
        if f(a)*f(b)>=0
            resLabel.Text='ERROR: f(a)*f(b) must be negative (opposite signs).';
            resLabel.FontColor=[.8 .1 .1]; return
        end
        hdr={'i','a','b','c (false pos)','f(c)','|f(c)|'};
        for i=1:maxIt
            fa=f(a); fb=f(b);
            c=a-fa*(b-a)/(fb-fa); fc=f(c);
            rows(end+1,:)={i,rnd(a),rnd(b),rnd(c),rnd8(fc),rnd8(abs(fc))};
            if abs(fc)<tol, root=c; conv=true; break; end
            if fa*fc<0, b=c; else, a=c; end
        end

      case 'Newton-Raphson'
        x=x0Field.Value;
        hdr={'i','x_n','f(x_n)','df(x_n)','x_{n+1}','|error|'};
        for i=1:maxIt
            fx=f(x); dfx=(f(x+1e-7)-f(x-1e-7))/(2e-7);
            if abs(dfx)<1e-14, resLabel.Text='Stopped: derivative near zero.'; break; end
            xn=x-fx/dfx;
            rows(end+1,:)={i,rnd(x),rnd8(fx),rnd8(dfx),rnd(xn),rnd8(abs(xn-x))};
            if abs(xn-x)<tol, root=xn; conv=true; break; end
            x=xn;
        end

      case 'Secant'
        x0=aField.Value; x1=bField.Value;
        hdr={'i','x_{n-1}','x_n','f(x_{n-1})','f(x_n)','x_{n+1}','|error|'};
        for i=1:maxIt
            f0=f(x0); f1=f(x1);
            if abs(f1-f0)<1e-14, break; end
            x2=x1-f1*(x1-x0)/(f1-f0);
            rows(end+1,:)={i,rnd(x0),rnd(x1),rnd8(f0),rnd8(f1),rnd(x2),rnd8(abs(x2-x1))};
            if abs(x2-x1)<tol, root=x2; conv=true; break; end
            x0=x1; x1=x2;
        end
    end

    tbl.ColumnName = hdr;
    if ~isempty(rows), tbl.Data=rows; end

    if conv && ~isnan(root)
        resLabel.Text = sprintf('[%s]  Root ≈ %.8f   |   f(root) = %.3e   |   %d iterations',...
            method, root, f(root), size(rows,1));
        resLabel.FontColor=[0.05 0.42 0.12];
    else
        resLabel.Text = sprintf('[%s]  Did not fully converge — %d iterations completed. Adjust range/tolerance.',...
            method, size(rows,1));
        resLabel.FontColor=[0.70 0.30 0.05];
    end

    % --- Plot ---
    if strcmp(method,'Newton-Raphson')
        xv=linspace(x0Field.Value-4, x0Field.Value+4, 400);
    else
        xv=linspace(aField.Value, bField.Value, 400);
    end
    yv=arrayfun(f,xv);
    cla(ax);
    plot(ax,xv,yv,'b-','LineWidth',2); hold(ax,'on');
    yline(ax,0,'k--','LineWidth',1);
    if ~isnan(root)
        plot(ax,root,f(root),'ro','MarkerSize',10,'MarkerFaceColor','r');
        legend(ax,{'f(x)','y = 0',sprintf('Root \\approx %.5f',root)},'Location','best');
    end
    title(ax,['f(x) = ' raw],'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'f(x)'); grid(ax,'on');
    yf=yv(isfinite(yv));
    if ~isempty(yf) && max(abs(yf))>0
        ylim(ax,[-max(abs(yf))*1.4, max(abs(yf))*1.4]);
    end
    hold(ax,'off');
end

function e=prepExpr(r)
    e=strrep(r,'^','.^');
    e=strrep(e,'*','.*');
    e=strrep(e,'/','./');
    e=regexprep(e,'(\d)(x)','$1.*x');
end
function v=rnd(x);  v=round(x,6);  end
function v=rnd8(x); v=round(x,8);  end

%%==========================================================================
%% MATRIX OPERATIONS TAB
%%==========================================================================
function buildMatrixTab(t)
    bg=[0.95 0.95 0.97];

    uilabel(t,'Text','Rows (max 6):','Position',[14 530 90 18],'FontSize',11);
    rF=uieditfield(t,'numeric','Position',[108 527 42 24],'Value',3,'Limits',[1 6]);
    uilabel(t,'Text','Cols (max 6):','Position',[160 530 90 18],'FontSize',11);
    cF=uieditfield(t,'numeric','Position',[254 527 42 24],'Value',3,'Limits',[1 6]);
    uilabel(t,'Text','Operation:','Position',[308 530 72 18],'FontSize',11);
    opD=uidropdown(t,'Items',{'Addition (A+B)','Multiplication (A×B)',...
        'Transpose of A','Determinant of A','Inverse of A','Adjoint of A',...
        'Power of A  (A^n)','Linear Equations (Ax=b)'},...
        'Position',[384 524 220 26],'FontSize',11);
    uilabel(t,'Text','n =','Position',[612 530 28 18],'FontSize',11);
    nF=uieditfield(t,'numeric','Position',[642 527 40 24],'Value',2,'Limits',[1 20]);

    % Build button
    buildBtn=uibutton(t,'Text','Build / Reset Grids','Position',[696 524 110 26],...
        'FontSize',10,'BackgroundColor',[0.82 0.82 0.90],'FontColor','Black');

    % Static matrix grid labels
    lblA=uilabel(t,'Text','Matrix A','Position',[14 510 80 18],'FontSize',11,'FontWeight','bold');  %#ok<NASGU>
    lblB=uilabel(t,'Text','Matrix B','Position',[310 510 80 18],'FontSize',11,'FontWeight','bold');
    lblVec=uilabel(t,'Text','Vector b','Position',[310 502 80 18],'FontSize',11,'FontWeight','bold','Visible','off');

    % Preallocate 6x6 grids A, B, Vec
    gA = makeGrid(t, 14, 490, 'gA', bg);
    gB = makeGrid(t, 310, 490, 'gB', bg);
    gV = makeGrid(t, 310, 490, 'gV', bg);
    showGrid(gA,3,3,true); showGrid(gB,3,3,false); showGrid(gV,3,1,false);

    resArea=uitextarea(t,'Position',[14 10 790 200],'Editable','off',...
        'FontSize',11,'FontName','Courier New');
    resArea.Value={'Matrix values will appear here after CALCULATE.'};

    calcBtn=uibutton(t,'Text','CALCULATE','Position',[696 214 110 28],...
        'FontSize',12,'FontWeight','bold',...
        'BackgroundColor',[0.10 0.42 0.72],'FontColor','white');

    % Callbacks
    buildBtn.ButtonPushedFcn=@(~,~) onBuild(rF,cF,opD,gA,gB,gV,lblB,lblVec);
    calcBtn.ButtonPushedFcn =@(~,~) onCalc(rF,cF,opD,nF,gA,gB,gV,resArea);
end

function g=makeGrid(parent,x,y,tag,~)
    cw=52; ch=24; gap=3;
    g=cell(6,6);
    for i=1:6
        for j=1:6
            xp=x+(j-1)*(cw+gap);
            yp=y-(i-1)*(ch+gap);
            v=double(i==j);
            g{i,j}=uieditfield(parent,'numeric','Position',[xp yp cw ch],...
                'Value',v,'FontSize',10,'Tag',sprintf('%s_%d_%d',tag,i,j),'Visible','off');
        end
    end
end

function showGrid(g,rows,cols,vis)
    for i=1:6
        for j=1:6
            g{i,j}.Visible=onoff(vis && i<=rows && j<=cols);
        end
    end
end

function onBuild(rF,cF,opD,gA,gB,gV,lblB,lblVec)
    r=round(rF.Value); c=round(cF.Value); op=opD.Value;
    needB  = any(strcmp(op,{'Addition (A+B)','Multiplication (A×B)'}));
    needV  = strcmp(op,'Linear Equations (Ax=b)');
    showGrid(gA,r,c,true);
    showGrid(gB,r,c,needB);
    showGrid(gV,r,1,needV);
    lblB.Visible  =onoff(needB);
    lblVec.Visible=onoff(needV);
end

function onCalc(rF,cF,opD,nF,gA,gB,gV,resArea)
    r=round(rF.Value); c=round(cF.Value); op=opD.Value;
    A=readG(gA,r,c);
    lines={};

    switch op
      case 'Addition (A+B)'
        B=readG(gB,r,c); R=A+B;
        lines=fmtMat('A + B',R);

      case 'Multiplication (A×B)'
        B=readG(gB,r,c);
        if size(A,2)~=size(B,1)
            resArea.Value={'ERROR: For A*B, cols(A) must equal rows(B).'}; return
        end
        lines=fmtMat('A × B', A*B);

      case 'Transpose of A'
        lines=fmtMat('Transpose  (A^T)', A');

      case 'Determinant of A'
        if r~=c, resArea.Value={'ERROR: Square matrix required.'}; return; end
        lines={sprintf('det(A)  =  %.8f', det(A))};

      case 'Inverse of A'
        if r~=c, resArea.Value={'ERROR: Square matrix required.'}; return; end
        if abs(det(A))<1e-12
            resArea.Value={'ERROR: Matrix is singular — inverse does not exist.'}; return
        end
        lines=fmtMat('Inverse  A^(-1)', inv(A));

      case 'Adjoint of A'
        if r~=c, resArea.Value={'ERROR: Square matrix required.'}; return; end
        lines=fmtMat('Adjoint of A', adjMat(A));

      case 'Power of A  (A^n)'
        if r~=c, resArea.Value={'ERROR: Square matrix required.'}; return; end
        n=round(nF.Value);
        lines=fmtMat(sprintf('A^%d',n), A^n);

      case 'Linear Equations (Ax=b)'
        if r~=c, resArea.Value={'ERROR: Square coefficient matrix required.'}; return; end
        b=readG(gV,r,1);
        if abs(det(A))<1e-12
            resArea.Value={'ERROR: Singular system — no unique solution.'}; return
        end
        x=A\b;
        lines={'Solution  x  of  Ax = b:', ' '};
        for i=1:length(x)
            lines{end+1}=sprintf('   x%d  =  %12.8f', i, x(i));
        end
        lines{end+1}=' ';
        lines{end+1}=sprintf('Verification  ||Ax - b||  =  %.4e', norm(A*x-b));
    end
    resArea.Value=lines;
end

function M=readG(g,rows,cols)
    M=zeros(rows,cols);
    for i=1:rows
        for j=1:cols
            if ~isempty(g{i,j}), M(i,j)=g{i,j}.Value; end
        end
    end
end

function lines=fmtMat(label,M)
    lines={[label ':'],' '};
    for i=1:size(M,1)
        row='  |';
        for j=1:size(M,2)
            row=[row sprintf('  %11.5f', M(i,j))];
        end
        lines{end+1}=[row '  |'];
    end
    lines{end+1}=' ';
end

function A=adjMat(M)
    n=size(M,1); C=zeros(n);
    for i=1:n
        for j=1:n
            C(i,j)=(-1)^(i+j)*det(M([1:i-1,i+1:end],[1:j-1,j+1:end]));
        end
    end
    A=C';
end

%% Shared utility — single definition
function s=onoff(b)
    if b, s='on'; else, s='off'; end
end