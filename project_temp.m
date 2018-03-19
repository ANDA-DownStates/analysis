pair = [3,5];

%%
rootFolder = 'C:\Users\superuser\Documents\ANDA\project\data';
dataFiles = ls([rootFolder, '\*.mat']);

% load file
fileID = pair(1);
load( fullfile(rootFolder, dataFiles(fileID+1,:)));

unitsNum = numel(block.segments{1}.spiketrains);
trialID = block.annotation_all_trial_ids;   %
trialNum = numel(trialID);

spikeTrains = cell(trialNum,1);   %
for i = 1:trialNum
    st_temp = cell(unitsNum,1);
    for j = 1:unitsNum
        st_temp{j,1} = block.segments{1,i}.spiketrains{1,j}.times;
    end
    spikeTrains{i,1} = st_temp;
end

%% get event name
evt2collect = [...
    'TS-ON      ';...
    'WS-ON      ';...
    'CUE-ON     ';...
    'CUE-OFF    ';...
    'GO-ON      ';...
    'SR         ';...
    'OT         ';...
    'DO         ';...
    %     'FSRplat-ON ';...
    %     'HEplat-ON  ';...
    %     'FSRplat-OFF';...
    %     'HEplat-OFF ';...
    %     'SR-REP     ';...
    'RW-ON      ';...
    'STOP       ';...
    %     'OBB        ';...
    %     'OR         ';...
    ];

tsON = eventTime(block, 'TS-ON');
wsON = eventTime(block, 'WS-ON');
cueON = eventTime(block, 'CUE-ON');
cueOFF = eventTime(block, 'CUE-OFF');
goON = eventTime(block, 'GO-ON');
SR = eventTime(block, 'SR');
OT = eventTime(block, 'OT');
DO = eventTime(block, 'SR');
rwON = eventTime(block, 'RW-ON');
stop = eventTime(block, 'STOP');



%% get trail type
trailType = char(zeros(trialNum,4));
for i = 1:trialNum
    trailType(i,:) = block.segments{1,i}.spiketrains{1,1}.annotation_belong_to_trialtype;
end

trailType_num = zeros(size(trailType,1),1);
typeName = {'PGHF'; 'SGHF'; 'PGLF'; 'SGLF'};
for t = 1:numel(trailType_num)
    trailType_num(t,1) = find(~cellfun(@isempty, strfind(typeName,trailType(t,:))));
end

% get information of single units: channels, unit # in the channel
unitNum = numel(spikeTrains{1,1});
unitInfo_temp = zeros( unitNum, 2, trialNum);   % format is [channel #, unit # in the channel, trial #]
for u1 = 1:unitNum
    for u3 = 1:trialNum
        unitInfo_temp(u1,1,u3) = block.segments{1,u3}.spiketrains{1,u1}.annotation_channel_id;
        unitInfo_temp(u1,2,u3) = block.segments{1,u3}.spiketrains{1,u1}.annotation_unit_id;
        unitInfo_temp(u1,1,u3) = u3;
    end
end


spikeTrains_chID_order = cell(trialNum,1);
for tri = 1:numel(spikeTrains)
    unitInfo_temp = zeros( unitNum, 2);   % format is [channel #, unit # in the channel]
    
    for u1 = 1:unitNum
        unitInfo_temp(u1,1) = block.segments{1,tri}.spiketrains{1,u1}.annotation_channel_id;
        unitInfo_temp(u1,2) = block.segments{1,tri}.spiketrains{1,u1}.annotation_unit_id;
    end
    [ui_sort,ui_ind] = sortrows(unitInfo_temp,[1,2]);
    spikeTrains_chID_order{tri,1} = spikeTrains{tri,1}(ui_ind);
end


figure; hold on;
trial = 50;
for i = 1:unitNum
    spk2plot = spikeTrains_chID_order{trial,1}{i,1};
    plot( spk2plot, i*ones(numel(spk2plot),1),'ko');
end



%% 
rootFolder = 'C:\Users\superuser\Documents\ANDA\project\data';
dataFiles = ls([rootFolder, '\*.mat']);

% load file
fileID = pair(2);
load( fullfile(rootFolder, dataFiles(fileID+1,:)));

unitsNum = numel(block.segments{1}.spiketrains);
trialID = block.annotation_all_trial_ids;   %
trialNum = numel(trialID);

spikeTrains = cell(trialNum,1);   %
for i = 1:trialNum
    st_temp = cell(unitsNum,1);
    for j = 1:unitsNum
        st_temp{j,1} = block.segments{1,i}.spiketrains{1,j}.times;
    end
    spikeTrains{i,1} = st_temp;
end

%% get event name
evt2collect = [...
    'TS-ON      ';...
    'WS-ON      ';...
    'CUE-ON     ';...
    'CUE-OFF    ';...
    'GO-ON      ';...
    'SR         ';...
    'OT         ';...
    'DO         ';...
    %     'FSRplat-ON ';...
    %     'HEplat-ON  ';...
    %     'FSRplat-OFF';...
    %     'HEplat-OFF ';...
    %     'SR-REP     ';...
    'RW-ON      ';...
    'STOP       ';...
    %     'OBB        ';...
    %     'OR         ';...
    ];

tsON = eventTime(block, 'TS-ON');
wsON = eventTime(block, 'WS-ON');
cueON = eventTime(block, 'CUE-ON');
cueOFF = eventTime(block, 'CUE-OFF');
goON = eventTime(block, 'GO-ON');
SR = eventTime(block, 'SR');
OT = eventTime(block, 'OT');
DO = eventTime(block, 'SR');
rwON = eventTime(block, 'RW-ON');
stop = eventTime(block, 'STOP');



%% get trail type
trailType = char(zeros(trialNum,4));
for i = 1:trialNum
    trailType(i,:) = block.segments{1,i}.spiketrains{1,1}.annotation_belong_to_trialtype;
end

trailType_num = zeros(size(trailType,1),1);
typeName = {'PGHF'; 'SGHF'; 'PGLF'; 'SGLF'};
for t = 1:numel(trailType_num)
    trailType_num(t,1) = find(~cellfun(@isempty, strfind(typeName,trailType(t,:))));
end

% get information of single units: channels, unit # in the channel
unitNum = numel(spikeTrains{1,1});
unitInfo_temp = zeros( unitNum, 2, trialNum);   % format is [channel #, unit # in the channel, trial #]
for u1 = 1:unitNum
    for u3 = 1:trialNum
        unitInfo_temp(u1,1,u3) = block.segments{1,u3}.spiketrains{1,u1}.annotation_channel_id;
        unitInfo_temp(u1,2,u3) = block.segments{1,u3}.spiketrains{1,u1}.annotation_unit_id;
        unitInfo_temp(u1,1,u3) = u3;
    end
end


spikeTrains_chID_order = cell(trialNum,1);
for tri = 1:numel(spikeTrains)
    unitInfo_temp = zeros( unitNum, 2);   % format is [channel #, unit # in the channel]
    
    for u1 = 1:unitNum
        unitInfo_temp(u1,1) = block.segments{1,tri}.spiketrains{1,u1}.annotation_channel_id;
        unitInfo_temp(u1,2) = block.segments{1,tri}.spiketrains{1,u1}.annotation_unit_id;
    end
    [ui_sort,ui_ind] = sortrows(unitInfo_temp,[1,2]);
    spikeTrains_chID_order{tri,1} = spikeTrains{tri,1}(ui_ind);
end


for i = 1:unitNum
    spk2plot = spikeTrains_chID_order{trial,1}{i,1};
    plot( spk2plot, i*ones(numel(spk2plot),1),'r.');
end

%%
title(['data' num2str(pair(1)), '(black circle) vs. data' num2str(pair(2)), ', trial ' num2str(trial) ', units are reordered by channel ID on Utah array'])
xlabel('time (s)')
ylabel('unit ID')

