function varargout = DataHunter(varargin)
% DATAHUNTER MATLAB code for DataHunter.fig
%      DATAHUNTER, by itself, creates a new DATAHUNTER or raises the existing
%      singleton*.
%
%      H = DATAHUNTER returns the handle to a new DATAHUNTER or the handle to
%      the existing singleton*.
%
%      DATAHUNTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAHUNTER.M with the given input arguments.
%
%      DATAHUNTER('Property','Value',...) creates a new DATAHUNTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataHunter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataHunter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataHunter

% Last Modified by GUIDE v2.5 03-Mar-2019 00:26:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataHunter_OpeningFcn, ...
                   'gui_OutputFcn',  @DataHunter_OutputFcn, ...
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



% --- Executes just before DataHunter is made visible.
function DataHunter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataHunter (see VARARGIN)

% Choose default command line output for DataHunter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% Init function
% close all
% constant init
global mainDir filtDir
 mainDir = evalin('base','which(''DataHunter'');');
 mainDir = mainDir(1:end-12);
 addpath(genpath(mainDir));
 filtDir = [mainDir, 'dataFilter\'];
 filters = dir(filtDir);
 filterGroup = {};
 for i = 1:numel(filters)
     if ~filters(i).isdir
         filterGroupTemp = [filtDir filters(i).name];
         try
             % try load config mat for mdl, i.e. for mdl.mdl, load mdl.mat
             load(strrep(filters(i).name,'mdl','mat'));
             load(strrep(filters(i).name,'slx','mat'));
         end
         filterGroup = [filterGroup;{filterGroupTemp,boolean(1)}];
     end
 end
 set(handles.ModelTable,'Data',filterGroup);

% UIWAIT makes DataHunter wait for user response (see UIRESUME)
% uiwait(handles.DataHunterMain);


% --- Outputs from this function are returned to the command line.
function varargout = DataHunter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function PinOnTop_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to PinOnTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp=get(hObject,'CData');
set(hObject,'CData',get(hObject,'UserData'));
set(hObject,'UserData',tmp);
% Get JavaFrame of Figure.
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
fJFrame = get(handles.DataHunterMain,'JavaFrame');
% Set JavaFrame Always-On-Top-Setting.
% if strcmpi(get(hObject,'State'),'on')
%     fJFrame.fFigureClient.getWindow.setAlwaysOnTop(1);
% else
%     fJFrame.fFigureClient.getWindow.setAlwaysOnTop(0);
% end

% if verLessThan('matlab', '7.10')
%     figclient='fFigureClient';
% elseif verLessThan('matlab', '8.4')
%     figclient = 'fHG1Client';
% else
%     figclient = 'fHG2Client';
% end
% Set JavaFrame Always-On-Top-Setting.
fig = get(hObject,'CData');
% figInv = reshape([reshape(flipud(fig(:,:,1)'),1,256) ...
%                   reshape(flipud(fig(:,:,2)'),1,256)... 
%                   reshape(flipud(fig(:,:,3)'),1,256)],...
%                   16,16,3);
% assignin('base','a',figInv);
if strcmp(get(hObject,'State'),'on')
    if verLessThan('matlab', '7.10')
        figclient='fFigureClient';
    elseif verLessThan('matlab', '8.4')
        figclient = 'fHG1Client';
    else
        figclient = 'fHG2Client';
    end
    fJFrame.(figclient).getWindow.setAlwaysOnTop(1);
%     set(hObject,'CData',figInv);
else
    if verLessThan('matlab', '7.10')
        figclient='fFigureClient';
    elseif verLessThan('matlab', '8.4')
        figclient = 'fHG1Client';
    else
        figclient = 'fHG2Client';
    end
    fJFrame.(figclient).getWindow.setAlwaysOnTop(0);
%     set(hObject,'CData',figInv);
end

% --------------------------------------------------------------------
function OpenDataFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to OpenDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
defaultdir = evalin('base','cd;');
[dataName,dataDir]=uigetfile({'*.mdf;*.dat;*.ascii;*.xls;*.xlsx'},'Pick the files to read',defaultdir,'MultiSelect','on');
if isequal(dataName,0) && isequal(dataDir,0)
    ;
else
dataGroup ={};
% main data loop
if iscell(dataName)
    [~,cnt]=size(dataName);
else
    [cnt,~]=size(dataName);
end
for i = 1:cnt
    try 
        tgtDataName = [dataDir,dataName{i}];
        dataNameTemp = dataName{i};
    catch
        tgtDataName = [dataDir,dataName];
        dataNameTemp = dataName;
    end
    dataGroup = [dataGroup;{tgtDataName,boolean(1),''}];
end
set(handles.DataTable,'Data',dataGroup);
end

% --------------------------------------------------------------------
function DataHunt_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to DataHunt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save results to a var
persistent rsltList
rslt.crtTm = strrep(datestr(now,31),' ','_');
rslt.crtTm = strrep(rslt.crtTm,':','-');
rslt.author = getComputerName();
% close warning function
warning('off');

% init
sampleTm = 0.01;
initTm   = 0;
stopTm   = 9999;

mdlIdx = cell2mat(handles.ModelTable.Data(:,2));
mdlIdx = find(mdlIdx == 1);
mdlLst = handles.ModelTable.Data(mdlIdx,1);
dataGroup = handles.DataTable.Data;
% refresh outcome
[row,~] = size(dataGroup);
for i = 1:row
    dataGroup{i,3} = '';
end
set(handles.DataTable,'Data',dataGroup);
dataIdx = cell2mat(handles.DataTable.Data(:,2));
dataIdx = find(dataIdx == 1);
dataLst = handles.DataTable.Data(:,1);
    
for k = 1:numel(mdlIdx)
    fIdx = 1;
    targetSim = mdlLst{k,1}; % current solution
    
    load_system(targetSim);
    inports = find_system(bdroot(gcs),'FindAll','on','SearchDepth',1,'BlockType','Inport');
    varNames = {};
    for i = 1:numel(inports)
        set_param(inports(i),'OutDataTypeStr','double');
        inportNameTemp = get_param(inports(i),'name');
        varNames = [varNames;inportNameTemp];
    end
    
    % set simulink model parameters for sim
    cfg = getActiveConfigSet(gcs);
    cfg.Components(1).SolverType='Fixed-step';
    cfg.Components(1).Solver='FixedStepDiscrete';
    cfg.Components(1).FixedStep=num2str(sampleTm);
    cfg.Components(2).SaveOutput = 'off';
    cfg.Components(2).SaveTime = 'off';
    cfg.Components(2).SaveFormat = 'StructureWithTime';
    
    % close_system(targetSim);
    
    for i = 1:numel(dataIdx)
        tgtDataName = dataLst{dataIdx(i)};
        dataGroup = read_MDA_dat(tgtDataName);
        
        % prep single sim
        ExpData2Sim_time=[];
        ExpData2Sim_inport=[];
        ExpData2Sim_outport=[];
        
        VarDGs=[];
        
        for j = 1:numel(varNames)
            varname = varNames{j};
            [st,idx]=getVarSampleTime(varname,dataGroup);
            VarDGs=[VarDGs;dataGroup{idx}];
        end
        
        % get time
        [refcnt,refidx]=max([VarDGs.sampleCount]);
        full_time_rng=VarDGs(refidx).Time;
        ref_time=[0:sampleTm:full_time_rng(end)]';
        sampleCnt=numel(ref_time);
        ExpData2Sim_time=ref_time; % time check
        cfg.Components(1).StartTime=num2str(min(ref_time));
        cfg.Components(1).StopTime=num2str(max(ref_time));
        
        % get input signals;
        for j=1:numel(varNames)
            fD=VarDGs(j);
            varname=varNames{j};
            vardata=fD.ExpData{strmatch(varname,fD.VariableNames,'exact')};
            colData=interp1(fD.Time,vardata,ref_time,'linear','extrap');
            ExpData2Sim_inport=[ExpData2Sim_inport,colData];
        end
        [t,~,r] = sim(targetSim,ExpData2Sim_time,[],[ExpData2Sim_time ExpData2Sim_inport]);
        %     if norm(single(r)) > 999e-3
        figure(mdlIdx(k));
        subplot(numel(dataIdx),1,fIdx)
        plot(t,r,'LineWidth',3);
        title(tgtDataName)
        fIdx = fIdx + 1;
        %     end
        
        %   execuate if save result check
        %   find event time
        eventTime = getEventTime(t,r);
        rslt.runs.(getFileName(tgtDataName)).(getFileName(targetSim))= eventTime;
        
    end
    
    % save and quit
    save_system(targetSim);
    % close system
    close_system(targetSim)

end
% Restore warning
warning('on');
% Write results to excel use another function
saveRunResult(rslt)

% --- Executes when user attempts to close DataHunterMain.
function DataHunterMain_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to DataHunterMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function Refresh_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filtDir
filters = dir(filtDir);
 filterGroup = {};
 for i = 1:numel(filters)
     if ~filters(i).isdir
         filterGroupTemp = [filtDir filters(i).name];
         try
             % try load config mat for mdl, i.e. for mdl.mdl, load mdl.mat
             load(strrep(filters(i).name,'mdl','mat'));
             load(strrep(filters(i).name,'slx','mat'));
         end
         filterGroup = [filterGroup;{filterGroupTemp,boolean(1)}];
     end
 end
 set(handles.ModelTable,'Data',filterGroup);
 
function [st,idx]=getVarSampleTime(varname,DataGroup)
%Get the Sample Time of the variable
%find the first one in the data group list
st=[];idx=[];
for j=1:numel(DataGroup)
    if ismember(varname,DataGroup{j}.VariableNames)
        st=DataGroup{j}.SampleTime;
        idx=j;
        break;
    end
end

function ret = getEventTime(t,y)
temp = {};
y_old = 0;
for i = 1:numel(y)
    y_new = y(i);
    if abs(y_old) <= 1e-5 && y_new~=0
        temp = {temp{:},num2str(t(i))};
    end
    y_old = y_new;
end
ret = temp;

function ret = getFileName(rawName)
% get rid of filename and extension
temp = regexp(rawName,'\w*(?=\.)','match');
temp = strrep(temp,' ','_');
ret = temp{:};

function saveRunResult(rslt)
global mainDir
% starting excel
hExcel=actxserver('Excel.Application');
hExcel.visible = 1; % for debug
hWkbk = invoke(hExcel.Workbooks,'Add'); % create excel files;
hSht = hWkbk.Sheets.Item(1);
% writing results
% info
hSht.Range('A1').Value = 'DataHuntResults';
hSht.Range('A1').Font.Bold = 4;
hSht.Range('A1').Font.Size = 16;
hSht.Range('A3').Value = 'Author';
hSht.Range('A3').Font.Bold = 4;
hSht.Range('B3').Value = rslt.author;
hSht.Range('A4').Value = 'Date';
hSht.Range('A4').Font.Bold = 4;
hSht.Range('B4').Value = rslt.crtTm;

% data char(65) = 
rowIdx = 6; % starting for 6th row
datas = fields(rslt.runs);
for i = 1:numel(datas)
    hSht.Range(['A',num2str(rowIdx)]).Value = datas{i};
    hSht.Range(['A',num2str(rowIdx)]).Font.Bold = 4;
    hSht.Range(['A',num2str(rowIdx)]).Font.Size = 11;
    rowIdx = rowIdx + 1;
    mdls = fields(rslt.runs.(datas{i}));
    for j = 1:numel(mdls)
        hSht.Range(['A',num2str(rowIdx)]).Value = mdls{j};
        hSht.Range(['A',num2str(rowIdx)]).Font.Bold = 4;
        rowIdx = rowIdx + 1;
        rsltTmp = rslt.runs.(datas{i}).(mdls{j});
        hSht.Range(['A',num2str(rowIdx),':',char(numel(rsltTmp)+64),num2str(rowIdx)]).Value = rsltTmp;
        rowIdx = rowIdx + 2;
    end
end
% char 
% mkdir
try 
    mkdir([mainDir '\dataHuntResult']) 
end
% save and quit
hWkbk.SaveAs([mainDir 'dataHuntResult\DataHuntResult_',rslt.crtTm,'.xlsx']);
Quit(hExcel)
delete(hExcel)
