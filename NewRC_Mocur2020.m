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
end

function createAxes(handles)
    axes(handles.axes1);
    hold on;
%     set(gca,'XColor',[0.94 0.94 0.94]);
%     set(gca,'YColor',[0.94 0.94 0.94]);
    xlabel('b', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    ylabel(' h', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');

    axes(handles.axes2);
    hold on;
    xlabel('{\epsilon }', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    ylabel(' fc (kgf/cm2)', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    title('Concrete', 'color', 'k', 'FontSize', 13, 'FontName', 'Times');

    axes(handles.axes3);
    hold on;
    xlabel('{\epsilon }', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    ylabel(' fs (kgf/cm2)', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    title('Steel','color', 'k', 'FontSize', 13, 'FontName', 'Times');

    axes(handles.axes4);
    hold on;
    xlabel('Curvature (1/cm)', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    ylabel(' Moment (kgf-cm)', 'color', 'k', 'FontSize', 12, 'FontName', 'Times');
    title('Moment - Curvature', 'color', 'k', 'FontSize', 13, 'FontName', 'Times');
end

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

    createAxes(handles);
end


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
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	param = refactor_read_data(handles);
    
    % 弄把计
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
    
    steelParam = refactor_steel_conversion(param);
    
    % 弄 steel_conversion 把计
    % unpack parameters in steel_conversion
    handles = steelParam.handles;
    param.handles = handles;
    ssAs = steelParam.ssAs;
    steelnum = steelParam.steelnum;
    totalAs = steelParam.totalAs;
    dtmax = steelParam.dtmax;
    
    visualParam = refactor_visualization(param, dtmax);
    
    % 弄 visualization 把计
    % unpack parameters in visualization
    scoordinate = visualParam.scoordinate;
    handles = visualParam.handles;
    param.handles = handles;
    if concretetype == 0 || concretetype == 1
        unr = visualParam.mander.unr;
        unecc = visualParam.mander.unecc;
        unfcc = visualParam.mander.unfcc;
        r = visualParam.mander.r;
        ecc = visualParam.mander.ecc;
        fcc = visualParam.mander.fcc;
        ectension = visualParam.mander.ectension;
        unconfined = visualParam.mander.unconfined;
    end

%     p=1;
%     for e=1:1:200%%%
%         for g=1:1:200%%%
%             cxy(:,p)=[(b/200)/2+(b/200)*(g-1) (h/200)/2+(h/200)*(e-1)]; %%%
%             p=p+1;
%         end
%     end
%   DM: `x = 1:5; y = 10:10:50; cxy = [reshape(repmat(x, 1, 5), 1, []); reshape(repmat(y, 5, 1), 1, [])]`
    x = 1:200;
    cx = (b / 200) / 2 + (b / 200) * (x - 1);
    cy = (h / 200) / 2 + (h / 200) * (x - 1);
    cxy = [reshape(repmat(cx, 1, 200), 1, []); reshape(repmat(cy, 200, 1), 1, [])];

    theta = 0;
    srotation = [cosd(theta), -sind(theta); sind(theta), cosd(theta)] * scoordinate;
    crotation = [cosd(theta), -sind(theta); sind(theta), cosd(theta)] * cxy;
%     dtmax=1;
%     for K=1:1:steelnum
%         if srotation(2,K)>dtmax
%             dtmax=srotation(2,K);
%             dtnum=K;
%         end
%     end
    [dtmax, dtnum] = max(srotation(2, :));
    
    handles.myData1 = s;
    handles.myData2 = ss;
    handles.myData3 = ssrow;
    handles.myData4 = sstyp;
    guidata(hObject, handles);

    refactor_write_file(param, steelParam, visualParam, ...
        srotation, crotation, cxy, dtmax, dtnum);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%     s=handles.myData1;
%     ss=handles.myData2;
%     ssrow=handles.myData3;
%     sstyp=handles.myData4;

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
    
%     s=0;ss=0;ssrow=0;sstyp=0;
    handles.text4.String = '';
    handles.text6.String = '';
    handles.text16.String = '';
    handles.text17.String = '';
    handles.text18.String = '';
    handles.text19.String = '';
    handles.text20.String = '';
    handles.text21.String = '';
    handles.text23.String = '';
    handles.text24.String = '';
    handles.text25.String = '';
    handles.text30.String = '';
    handles.text31.String = '';
    handles.text32.String = '';

    createAxes(handles);
end