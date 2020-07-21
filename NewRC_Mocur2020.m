function varargout = NewRC_Mocur2020(varargin)
% NEWRC_MOCUR2020 MATLAB code for NewRC_Mocur2020.fig
%      NEWRC_MOCUR2020, by itself, creates a new NEWRC_MOCUR2020 or raises the existing
%      singleton*.
%
%      H = NEWRC_MOCUR2020 returns the handle to a new NEWRC_MOCUR2020 or the handle to
%      the existing singleton*.
%
%      NEWRC_MOCUR2020('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWRC_MOCUR2020.M with the given input arguments.
%
%      NEWRC_MOCUR2020('Property','Value',...) creates a new NEWRC_MOCUR2020 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NewRC_Mocur2020_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NewRC_Mocur2020_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NewRC_Mocur2020

% Last Modified by GUIDE v2.5 20-Apr-2020 11:10:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @NewRC_Mocur2020_OpeningFcn, ...
    'gui_OutputFcn',  @NewRC_Mocur2020_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NewRC_Mocur2020 is made visible.
function NewRC_Mocur2020_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NewRC_Mocur2020 (see VARARGIN)

% Choose default command line output for NewRC_Mocur2020
handles.output = hObject;
handles.legend_list = [];
% Update handles structure
guidata(hObject, handles);

set(handles.axes5,'visible','off');
axes(handles.axes5);
image=imread('school.jpg');
imshow(image);

axes(handles.axes1)
set(gca,'XColor',[0.94 0.94 0.94]);
set(gca,'YColor',[0.94 0.94 0.94]);
xlabel('b','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' h','color','k','FontSize',12,'FontName','Times');
hold on

axes(handles.axes2)
xlabel('{\epsilon }','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' fc (kgf/cm2)','color','k','FontSize',12,'FontName','Times');
hold on
title('Concrete','color','k','FontSize',13,'FontName','Times');

axes(handles.axes3)
xlabel('{\epsilon }','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' fs (kgf/cm2)','color','k','FontSize',12,'FontName','Times');
hold on
title('Steel','color','k','FontSize',13,'FontName','Times');

axes(handles.axes4)
xlabel('Curvature (1/cm)','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' Moment (kgf-cm)','color','k','FontSize',12,'FontName','Times');
hold on
title('Moment - Curvature','color','k','FontSize',13,'FontName','Times');



% UIWAIT makes NewRC_Mocur2020 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NewRC_Mocur2020_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	param = refactor_read_data(handles);
    
    % 讀取參數
    % unpack all the parameters (temparary)
    handles = param.handles;
    
    pathname = param.pathname;
    filename = param.filename;
    name = param.name;
    concretetype = param.concretetype;
        
    h = param.h;
    b = param.b;
    cover = param.cover;
    P = param.P;
    stnum = param.stnum;
    
    fy = param.fy;
    fu = param.fu;
    Es = param.Es;
    esh = param.esh;
    esu = param.esu;
    power = param.power;
    fcr = param.fcr;

    s = param.s;
    ss = param.ss;
    number = param.number;
    mat = param.mat;
    sstyp = param.sstyp;
    ssrow = param.ssrow;
    
    if concretetype == 0 || concretetype == 1
        mander = param.mander;
        fc = mander.fc;
        Ec = mander.Ec;
        K = mander.K;
        ecu = mander.ecu;
        ft = mander.ft;
    elseif concretetype == 2
        hung = param.hung;
        SIGtc = hung.SIGtc;
        SIGtp = hung.SIGtp;
        EPSILONtc = hung.EPSILONtc;
        EPSILONtp = hung.EPSILONtp;
        EPSILONtu = hung.EPSILONtu;
        SIGcp = hung.SIGcp;
        SIGcu = hung.SIGcu;
        EPSILONcp = hung.EPSILONcp;
        EPSILONcu = hung.EPSILONcu;
    end


p=1;
for i=1:1:mat(1)
    areaAs(i)=0;j=1;
    for j=1:1:(mat(2)-2)
        switch ss(i,j)
            case 0
                As=0;
            case 3
                As=0.71;
            case 4
                As=1.27;
            case 5
                As=1.99;
            case 6
                As=2.87;
            case 7
                As=3.87;
            case 8
                As=5.07;
            case 9
                As=6.47;
            case 10
                As=8.14;
            case 11
                As=10.07;
            case 12
                As=12.19;
            case 14
                As=14.52;
            case 16
                As=19.79;
            case 18
                As=25.79;
        end
        sAs(i,j)=As;
        if As>0
            ssAs(1,p)=As;
            p=p+1;
        end
        areaAs(i)=areaAs(i)+sAs(i,j);
        As=0;
    end
end
steelnum=(p-1);
totalAs=0;dtmax=1;
for n=1:1:mat(1)
    d(n)=s(n,1);
    totalAs=totalAs+areaAs(n);
    if d(n)>dtmax
        dtmax=d(n);
    end
end

steelratio=roundn(totalAs/(b*h),-4);
set(handles.text28,'String',name);
set(handles.text4,'String',cover);
set(handles.text6,'String',steelratio);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%畫斷面圖
x=0:0.1:b;
y=0:0.1:h;
axes(handles.axes1)
if b>=h
    axis([0 b 0 b]);
else
    axis([0 h 0 h]);
end
set(gca,'XColor','k');
set(gca,'YColor','k');
xlabel('b','color','k','FontSize',10,'FontWeight','bold');
hold on
ylabel('h','color','k','FontSize',10,'FontWeight','bold');
hold on
axis equal;
plot(x,0,'k.');
hold on;
plot(x,h,'k.');
hold on;
plot(0,y,'k.');
hold on;
plot(b,y,'k.');
hold on;
axis equal;


for i=1:1:mat(1)
    number(i)=number(i)-2;
    if number(i)>0
        sxfirst(i,1)=h-dtmax;
        sx(i,1)=sxfirst(i,1);
        syfirst(1,1)=sxfirst(1,1);
        if number(i)> 1
            sxend(i,(number(i)))=b-sxfirst(i,1);
            sx(i,(number(i)))=sxend(i,(number(i)));
        end
        sy(i)=s(i,1);
        axes(handles.axes1)
        plot(sxfirst(i,1),h-sy(i),'ko');hold on;%%%%%%
        if (number(i))> 1
            plot(sxend(i,(number(i))),h-sy(i),'ko');hold on;%%%%%%
        end
        if (number(i))>2
            for j=2:1:(number(i)-1)
                sx(i,j)=sxfirst(i,1)+(j-1)*(((sxend(i,(number(i))))-sxfirst(i,1))/(number(i)-1));
                axes(handles.axes1)
                plot(sx(i,j),h-sy(i),'ko');hold on;%%%%%
            end
        end
    end
end

p=1;
for i=1:1:mat(1)
    for f=1:1:(number(i))
        scoordinate(:,p)=[sx(i,f) sy(i)];
        p=p+1;
    end
end





betal=1.05-(0.05*(fc/70));
if betal<0.65
    betal=0.65;
end
if betal>0.85
    betal=0.85;
end


%%%%%%%%%%%%%%%% MANDER Model %%%%%%%%%%%%%%%%%
if concretetype==0||concretetype==1
    unfcc=fc;
    if unfcc<=204
        unecc=0.002;
    elseif unfcc>204||unfcc<1020
        unecc=0.002+0.001*(unfcc-204)/816;
    else
        unecc=0.003;
    end
    unEsec=unfcc/unecc;
    unr=Ec/(Ec-unEsec);
    fcc=K*fc;
    ecc=unecc*(1+5*(K-1));
    Esec=fcc/ecc;
    r=Ec/(Ec-Esec);
    
    v=1;u=1;
    for xx=0:0.0001:ecu
        confinedyy(u,v)= (r*(xx/ecc))/((r-1)+(xx/ecc).^r)*fcc;
        v=v+1;
    end
    v=1;u=1;
    for xx=0:0.0001:0.004
        unconfinedyy(u,v)= (unr*(xx/unecc))/((unr-1)+(xx/unecc).^unr)*unfcc;
        v=v+1;
    end
    unconfined=unconfinedyy(u,v-1);
    axes(handles.axes2)
    ectension=ft/Ec;
    xx=0:0.0001:ecu;
    plot(xx,confinedyy,'r');hold on;xtickformat('%.3f');
    xx=0:0.0001:0.004;
    plot(xx,unconfinedyy,'r');hold on;xtickformat('%.3f');
    line([0.004,0.006],[unconfined,0],'Color','red');hold on;
    line([0,-ectension],[0,-ft],'Color','blue');grid on;hold on;
end

%%%%%%%%%%%%%%%%%%%%% HUNG Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if concretetype==2
    v=1;u=1;
    for xx=0:0.0001:EPSILONcp
        HPFRCCyy(u,v)= SIGcp*(2*(xx/EPSILONcp)-(xx/EPSILONcp)^2);
        v=v+1;
    end
    axes(handles.axes2)
    line([0,EPSILONtc],[0,SIGtc],'Color','blue');grid on;hold on;
    line([EPSILONtc,EPSILONtp],[SIGtc,SIGtp],'Color','blue');hold on;
    line([EPSILONtp,EPSILONtu],[SIGtp,0],'Color','blue');hold on;
    xx=0:0.0001:EPSILONcp;
    plot(-xx,-HPFRCCyy,'r');hold on;xtickformat('%.3f');
    line([-EPSILONcp,-EPSILONcu],[-SIGcp,-SIGcu],'Color','red');hold on;
    line([-EPSILONcu,-0.1],[-SIGcu,-SIGcu],'Color','red');hold on;
end
%%%%%%%%%%%%%%%% HUNG Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for num=1:1:stnum
    steely=[];
    axes(handles.axes3)
    line([0,fy(num)/Es(num)],[0,fy(num)],'Color','blue');hold on;
    line([fy(num)/Es(num),esh(num)],[fy(num),fy(num)],'Color','blue');hold on;
    f=1;
    for x=esh(num):0.01:esu(num)
        steely(1,f)=fu(num)+(fy(num)-fu(num))*((esu(num)-x)/(esu(num)-esh(num))).^power(num);
        f=f+1;
    end
    axes(handles.axes3)
    x=esh(num):0.01:esu(num);
    plot(x,steely,'b');hold on;grid on;
end

for num=1:1:stnum
    steelxx=[];steelyy=[];
    axes(handles.axes3)
    line([0,-(fy(num)/Es(num))],[0,-fy(num)],'Color','red');hold on;
    %%line([-0.004,-(fy(num)/Es(num))],[-fy(num),-fy(num)],'Color','red');hold on;
    if fcr(num)>fy(num)
        line([-(fy(num)/Es(num)),-esh(num)],[-fy(num),-fy(num)],'Color','red');hold on;
    end
    f=1;
    if fcr(num)>fu(num)
        for x=esh(num):0.01:esu(num)
            steelyy(1,f)=-(fu(num)+(fy(num)-fu(num))*((esu(num)-x)/(esu(num)-esh(num))).^power(num));
            f=f+1;
        end
        axes(handles.axes3)
        x=esh(num):0.01:esu(num);
        plot(-x,steelyy,'r');hold on;grid on;
        
    elseif fcr(num)<=fu(num)&&fcr(num)>=fy(num)
        ex=esu(num)-(((fcr(num)-fu(num))/(fy(num)-fu(num)))^(1/power(num)))*(esu(num)-esh(num));
        for x=esh(num):0.01:ex
            steelxx(1,f)=-x;
            steelyy(1,f)=-(fu(num)+(fy(num)-fu(num))*((esu(num)-x)/(esu(num)-esh(num))).^power(num));
            f=f+1;
        end
        axes(handles.axes3)
        plot(steelxx,steelyy,'r');hold on;grid on;
    end
    AXES=handles.axes3;
    plot(AXES,0,0,'r');
    AXES.YTickLabel = string(AXES.YTick);
end
AXE=handles.axes2;
plot(AXE,0,0);
AXE.XTickLabel=string(AXE.XTick);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axes2)
plot(0,0,'r');

p=1;
for e=1:1:200%%%
    for g=1:1:200%%%
        cxy(:,p)=[(b/200)/2+(b/200)*(g-1) (h/200)/2+(h/200)*(e-1)]; %%%
        p=p+1;
    end
end

theta=0;
srotation=[cosd(theta) -sind(theta);sind(theta) cosd(theta)]*[scoordinate];
crotation=[cosd(theta) -sind(theta);sind(theta) cosd(theta)]*[cxy];
dtmax=1;
for K=1:1:steelnum
    if srotation(2,K)>dtmax
        dtmax=srotation(2,K);
        dtnum=K;
    end
end
handles.myData1= s;
guidata(hObject, handles);
handles.myData2= ss;
guidata(hObject, handles);
handles.myData3= ssrow;
guidata(hObject, handles);
handles.myData4= sstyp;
guidata(hObject, handles);


path=strcat(pathname,'Mocur_',filename);
fp=fopen(path,'w');
fprintf(fp,'Mn                        Curvature                        c                    MaxStrainofCoreConc                  fsmax                 fsmin\r\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Paxial=P*b*h*fc;
epsilonC=0.00001;Curvature=0;
Csum=-1;Cs=0;count=1;fsmax=-20000;fsmin=20000;
if Paxial>0
    while abs((-Csum/Paxial)-1)>0.03&&count<200
        Csum=0;Cs=0;Mpo=0;fs=0;
        for m=1:1:steelnum
            if epsilonC*Es(sstyp(m))>=fy(sstyp(m))
                fs(m)=-(fy(sstyp(m)));
                Cs=Cs+ssAs(1,m)*fs(m);
            else
                fs(m)=(-epsilonC)*Es(sstyp(m));
                Cs=Cs+ssAs(1,m)*fs(m);
            end
            Mpo=-(ssAs(1,m)*fs(m)*((h/2)-srotation(2,m)))+Mpo;
        end
        if concretetype==0||concretetype==1
            eachconcreteforce=-(unr*(epsilonC/unecc))/((unr-1)+(epsilonC/unecc).^unr)*unfcc;
        elseif concretetype==2
            eachconcreteforce=-SIGcp*(2*(-epsilonC/(-EPSILONcp))-(-epsilonC/(-EPSILONcp))^2);
        end
        
        Csum=((eachconcreteforce)*b*h+Cs+totalAs*(-eachconcreteforce));
        if abs((-Csum/Paxial)-1)>0.3&&(-Csum/Paxial)<1
            epsilonC=epsilonC+0.0001;
        elseif abs((-Csum/Paxial)-1)>0.1&&(-Csum/Paxial)<1
            epsilonC=epsilonC+0.00005;
        elseif abs((-Csum/Paxial)-1)>0.01&&(-Csum/Paxial)<1
            epsilonC=epsilonC+0.00002;
        elseif abs((-Csum/Paxial)-1)>0.01&&(-Csum/Paxial)>1
            epsilonC=epsilonC-0.00001;
        end
        count=count+1;
    end
    fsmax=fs(1);
    fsmin=fs(1);
    PepsilonC=epsilonC;
    c=Inf;
    MaxStrainofCoreConc=-epsilonC;
    xCurvature(1,1)=Curvature;
    yMn(1,1)=Mpo;
    fprintf(fp,'%e              %e                   %2.4f                          %+2.7f                       %+5.1f                %+5.1f\r\n',Mpo,Curvature,c,MaxStrainofCoreConc,fsmax,fsmin);
end

PepsilonC=epsilonC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=h/2;central=h/2;flag=0;flagn=0;MaxStrainofCoreConc=0;Mpeak=0;
v=2;epsilonTmax=0.001;Mnominal=0;Curvaturebefore=0;flagfcr=0;
while (-MaxStrainofCoreConc<=ecu)&(epsilonTmax<=esu(sstyp(dtnum)))
    Csum=-1;Tsum=0.1;count=1;Mxx=0;Myy=0;
    c=central;
    while abs((-Csum/Tsum)-1)>0.01&&count<80
        Mcx=0;Mcy=0;Cc=0;Tc=0;Cs=0;T=0;Csum=0;Tsum=0;
        if concretetype==0||concretetype==1
            for g=1:1:40000%%%
                epsilonc(g)=epsilonC*(crotation(2,g)-c)/c;   %斷面下壓上拉
                if epsilonc(g)<0 && epsilonc(g)> -ecu    %受壓混凝土
                    if crotation(1,g)<=cover || crotation(1,g)>= (b-cover) || crotation(2,g)<=cover || crotation(2,g)>=(h-cover) %保護層
                        if epsilonc(g)>-0.004
                            eachfc(g)=-(unr*(-epsilonc(g)/unecc))/((unr-1)+(-epsilonc(g)/unecc).^unr)*unfcc;
                        elseif epsilonc(g)<=-0.004 && epsilonc(g)>=-0.006
                            eachfc(g)=-(-unconfined/0.002*-epsilonc(g)+3*unconfined);
                        else
                            eachfc(g)=0;
                        end
                    else
                        eachfc(g)=-(r*(-epsilonc(g)/ecc))/((r-1)+(-epsilonc(g)/ecc).^r)*fcc;  %將epsilonc(g)加負號變回正值，算完後再將應力變回負值代表壓力
                    end
                elseif epsilonc(g)>0 && epsilonc(g)<=ectension       %受拉混凝土
                    eachfc(g)=epsilonc(g)*Ec;                        %受拉應力為正
                else
                    eachfc(g)=0;
                end
                if  eachfc(g)<0
                    Cc=Cc+eachfc(g)*b/200*h/200;%%%
                else
                    Tc=Tc+eachfc(g)*b/200*h/200;%%%
                end
                
                Mcx=-(eachfc(g)*b/200*h/200*((h/2)-cxy(2,g)))+Mcx;%%%
                Mcy=(eachfc(g)*b/200*h/200*((b/2)-cxy(1,g)))+Mcy;%%%
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% HUNG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if concretetype==2
            for g=1:1:40000%%%
                epsilonc(g)=epsilonC*(crotation(2,g)-c)/c;   %斷面下壓上拉
                if epsilonc(g)>= 0 && epsilonc(g)<= EPSILONtc                 %受拉混凝土
                    if  EPSILONtc==0||SIGtc==0
                        eachfc(g)=0;
                    else
                        eachfc(g)=epsilonc(g)/EPSILONtc*SIGtc;
                    end
                elseif epsilonc(g)>=EPSILONtc && epsilonc(g)<=EPSILONtp       %受拉混凝土
                    if  EPSILONtp==0||EPSILONtc==0||SIGtc==0||SIGtp==0
                        eachfc(g)=0;
                    else
                        eachfc(g)=SIGtc+(SIGtp-SIGtc)*(epsilonc(g)-EPSILONtc)/(EPSILONtp-EPSILONtc);
                    end
                elseif epsilonc(g)>=EPSILONtp && epsilonc(g)<= EPSILONtu       %受拉混凝土
                    if  EPSILONtu==0||EPSILONtp==0||SIGtp==0
                        eachfc(g)=0;
                    else
                        eachfc(g)=SIGtp*(1-(epsilonc(g)-EPSILONtp)/(EPSILONtu-EPSILONtp));
                    end
                    
                elseif epsilonc(g)>=EPSILONtu                                 %受拉混凝土
                    eachfc(g)=0;
                elseif epsilonc(g)>=-EPSILONcp && epsilonc(g)< 0              %受壓混凝土
                    eachfc(g)=-SIGcp*(2*(epsilonc(g)/(-EPSILONcp))-(epsilonc(g)/(-EPSILONcp))^2);
                elseif epsilonc(g)>=-EPSILONcu && epsilonc(g)< -EPSILONcp       %受壓混凝土
                    eachfc(g)=-SIGcp*(1-((epsilonc(g)-(-EPSILONcp))/((-EPSILONcu)-(-EPSILONcp)))*(1-SIGcu/SIGcp));
                else
                    eachfc(g)=-SIGcu;
                end
                if  eachfc(g)<0
                    Cc=Cc+eachfc(g)*b/200*h/200;%%%
                else
                    Tc=Tc+eachfc(g)*b/200*h/200;%%%
                end
                
                Mcx=-(eachfc(g)*b/200*h/200*((h/2)-cxy(2,g)))+Mcx;%%%
                Mcy=(eachfc(g)*b/200*h/200*((b/2)-cxy(1,g)))+Mcy;%%%
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%   Steel   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        epsilonTmax=epsilonC*(srotation(2,dtnum)-c)/c;
        
        
        Msx=0;Msy=0;Mn=0;fsmax=-20000;fsmin=20000;
        for m=1:1:steelnum
            epsilon(m)=epsilonC*(srotation(2,m)-c)/c;
            ety=fy(sstyp(m))/Es(sstyp(m));
            if epsilon(m)<ety&& epsilon(m)>=0
                fs(m)=epsilon(m)*Es(sstyp(m));
            elseif epsilon(m)<0&& epsilon(m)>-(ety)&&flagfcr==0
                fs(m)=epsilon(m)*Es(sstyp(m));
            elseif epsilon(m)<=esh(sstyp(m))&&epsilon(m)>=ety
                fs(m)=fy(sstyp(m));
            elseif epsilon(m)>=-esh(sstyp(m))&&epsilon(m)<=-(ety)
                if fcr(sstyp(m))<=fy(sstyp(m))%&&epsilon(m)>=0.004%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fs(m)=0;
                    flagfcr=flagfcr+1;
                elseif flagfcr==0
                    fs(m)=-fy(sstyp(m));
                end
            elseif epsilon(m)<-esh(sstyp(m))&&epsilon(m)>=-esu(sstyp(m))
                fs(m)=-(fu(sstyp(m))+(fy(sstyp(m))-fu(sstyp(m)))*((esu(sstyp(m))+epsilon(m))/(esu(sstyp(m))-esh(sstyp(m)))).^power(sstyp(m)));
                if fs(m)<=-fcr(sstyp(m))
                    fs(m)=0;
                    flagfcr=flagfcr+1;
                end
            elseif epsilon(m)>esh(sstyp(m))&&epsilon(m)<=esu(sstyp(m))
                fs(m)=fu(sstyp(m))+(fy(sstyp(m))-fu(sstyp(m)))*((esu(sstyp(m))-epsilon(m))/(esu(sstyp(m))-esh(sstyp(m)))).^power(sstyp(m));
            else
                fs(m)=0;
            end
            
            if fs(m)>fsmax
                fsmax=fs(m);
            end
            if fs(m)<fsmin
                fsmin=fs(m);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            epsilonc(m)=epsilonC*(srotation(2,m)-c)/c;
            if concretetype==0||concretetype==1
                if epsilonc(m)<0 && epsilonc(m)> -ecu
                    eachfc(m)=-(r*(-epsilonc(m)/ecc))/((r-1)+(-epsilonc(m)/ecc).^r)*fcc;
                end
            elseif concretetype==2
                if epsilonc(m)>=-EPSILONcp && epsilonc(m)< 0              %受壓混凝土
                    eachfc(m)=-SIGcp*(2*(epsilonc(m)/(-EPSILONcp))-(epsilonc(m)/(-EPSILONcp))^2);
                elseif epsilonc(m)>=-EPSILONcu && epsilonc(m)< -EPSILONcp       %受壓混凝土
                    eachfc(m)=-SIGcp*(1-((epsilonc(m)-(-EPSILONcp))/((-EPSILONcu)-(-EPSILONcp)))*(1-SIGcu/SIGcp));
                else
                    eachfc(m)=-SIGcu;
                end
            else
                eachfc(m)=0;
            end
            fs(m)=fs(m)-eachfc(m);        %eachfc(m)為負值 fs(m)為負值
            
            
            if  fs(m)<0
                Cs=Cs+fs(m)*ssAs(1,m);
            else
                T=T+fs(m)*ssAs(1,m);
            end
            Msx=-fs(m)*ssAs(1,m)*((h/2)-srotation(2,m))+Msx;
            Msy=fs(m)*ssAs(1,m)*((b/2)-srotation(1,m))+Msy;
        end
        Csum=Cc+Cs;
        Tsum=T+Tc+Paxial;       %%%%%%%%%%%%%%%%%Paxial
        Mxx=(Mcx+Msx);
        Myy=(Mcy+Msy);
        
        %%%%%%%%%fprintf(fp,'%2.4f           %e         %f\r\n',c,epsilonC,-Csum/Tsum);%%%%%%%%%%%
        if c>3*h
            if (abs((-Csum/Tsum)-1))>0.01
                c=c/(-Csum/Tsum);
            end
        else
            if (abs((-Csum/Tsum)-1))>0.01
                c=c/(((-Csum/Tsum)+5)/6);
            end
        end
        count=count+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    central=c;
    MaxStrainofCoreConc=-epsilonC/c*(c-syfirst(1,1));
    Curvature=epsilonC/c;
    MM=(Mxx*Mxx)+(Myy*Myy);
    Mn=sqrt(MM);
    epsilonTy=fy(sstyp(dtnum))/Es(sstyp(dtnum));
    
    
    if epsilonTmax<=esu(sstyp(dtnum)) & Curvature>=Curvaturebefore
        %%%%%%%%%fprintf(fp,'Mn                        Curvature                        c                    MaxStrainofCoreConc                  fsmax                 fsmin\r\n');
        fprintf(fp,'%e              %e                   %2.4f                     %+2.7f                       %+5.1f                %+5.1f\r\n',Mn,Curvature,c,MaxStrainofCoreConc,fsmax,fsmin);
        xCurvature(1,v)=Curvature;
        yMn(1,v)=Mn;
        if epsilonTmax>epsilonTy&&flag==0
            Myield=Mn;
            Curvatureyield=Curvature;
            cyield=c;
            flag=flag+1;
        end
        
        if epsilonC>=0.003&&flagn==0&&PepsilonC<0.003
            Mnominal=Mn;
            Curvaturenominal=Curvature;
            cnominal=c;
            flagn=flagn+1;
        end
        if Mn>Mpeak
            Mpeak=Mn;Curvaturepeak=Curvature;cpeak=c;
        end
        Mult=Mn;Curvatureult=Curvature;cult=c;
        v=v+1;Curvaturebefore=Curvature;
    end
    
    %%epsilonC=epsilonC+0.0005;
    if epsilonC<=0.003
        epsilonC=epsilonC+0.00003;
    elseif 0.003<epsilonC<=0.01
        epsilonC=epsilonC+0.0002;
    else
        epsilonC=epsilonC+0.0005;
    end
    
end



if flag==0
    set(handles.text16,'String','NA');
    set(handles.text17,'String','NA');
    set(handles.text18,'String','NA');
else
    set(handles.text16,'String',num2str(Curvatureyield,'%.3e'));
    set(handles.text17,'String',num2str(Myield,'%.3e'));
    set(handles.text18,'String',num2str(cyield,'%.2f'));
end


if flagn==0
    set(handles.text19,'String','NA');
    set(handles.text20,'String','NA');
    set(handles.text21,'String','NA');
else
    set(handles.text19,'String',num2str(Mnominal,'%.3e'));
    set(handles.text20,'String',num2str(Curvaturenominal,'%.3e'));
    set(handles.text21,'String',num2str(cnominal,'%.2f'));
end

set(handles.text31,'String',num2str(Mpeak,'%.3e'));
set(handles.text30,'String',num2str(Curvaturepeak,'%.3e'));
set(handles.text32,'String',num2str(cpeak,'%.2f'));

set(handles.text23,'String',num2str(Mult,'%.3e'));
set(handles.text24,'String',num2str(Curvatureult,'%.3e'));
set(handles.text25,'String',num2str(cult,'%.2f'));

axes(handles.axes4)
plot(xCurvature,yMn,'-o','MarkerSize',5);
legend_list = handles.legend_list;
legend(legend_list);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=handles.myData1;
ss=handles.myData2;
ssrow=handles.myData3;
sstyp=handles.myData4;

handles.legend_list = [];
guidata(hObject, handles);
axes(handles.axes1);
cla reset
axes(handles.axes2);
cla reset
axes(handles.axes3);
cla reset
axes(handles.axes4);

cla reset
s=0;ss=0;ssrow=0;sstyp=0;
handles.text4.String='';handles.text6.String='';
handles.text16.String='';handles.text17.String='';handles.text18.String='';
handles.text19.String='';handles.text20.String='';handles.text21.String='';
handles.text23.String='';handles.text24.String='';handles.text25.String='';
handles.text30.String='';handles.text31.String='';handles.text32.String='';
axes(handles.axes1)
set(gca,'XColor',[0.94 0.94 0.94]);
set(gca,'YColor',[0.94 0.94 0.94]);
xlabel('b','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' h','color','k','FontSize',12,'FontName','Times');
hold on

axes(handles.axes2)
xlabel('{\epsilon }','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' fc (kgf/cm2)','color','k','FontSize',12,'FontName','Times');
hold on
title('Concrete','color','k','FontSize',13,'FontName','Times');

axes(handles.axes3)
xlabel('{\epsilon }','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' fs (kgf/cm2)','color','k','FontSize',12,'FontName','Times');
hold on
title('Steel','color','k','FontSize',13,'FontName','Times');

axes(handles.axes4)
xlabel('Curvature (1/cm)','color','k','FontSize',12,'FontName','Times');
hold on
ylabel(' Moment (kgf-cm)','color','k','FontSize',12,'FontName','Times');
hold on
title('Moment-Curvature','color','k','FontSize',13,'FontName','Times');
