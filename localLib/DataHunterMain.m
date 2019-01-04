function DataHunterMain()
close all
% constant init
fIdx = 1;
sampleTm = 0.01;
initTm   = 0;
stopTm   = 9999;
% get main dir
 mainDir = evalin('base','which(''DataHunterMain'');');
 mainDir = mainDir(1:end-25);
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
         filterGroup = [filterGroup;filterGroupTemp];
     end
 end
% prep data filler
% use simulink model as data filter
targetSim = filterGroup{1};

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
% save and quit
save_system(targetSim);
% close_system(targetSim);

% prep data file
defaultdir = evalin('base','cd;');
[dataName,dataDir]=uigetfile({'*.mdf;*.dat;*.ascii;*.xls;*.xlsx'},'Pick the files to read',defaultdir,'MultiSelect','on');
dataGroup ={};
% main data loop
[cnt,~]=size(dataName);
for i = 1:cnt
    try 
        tgtDataName = [dataDir,dataName{i}];
        dataNameTemp = dataName{i};
    catch
        tgtDataName = [dataDir,dataName];
        dataNameTemp = dataName;
    end
    dataGroup = read_MDA_dat(tgtDataName);
    % target data acquired
    
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
    if norm(single(r)) > 999e-3
        figure(fIdx);
        plot(t,r)
        title(dataNameTemp)
        fIdx = fIdx + 1;
    end
end

% close system
close_system(targetSim)

end

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
end
