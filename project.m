% clear; clc;

rootFolder = 'C:\Users\superuser\Documents\ANDA\project\data';
dataFiles = ls([rootFolder, '\*.mat']);

% load file
fileID = 4;
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

% % % plot all trial with critical time point
% % figure; hold on;
% % for e = 1:trialNum
% %     plot(tsON(e,1),e,'k.');
% %     plot(cueON(e,1),e,'r.');
% %     plot(cueOFF(e,1),e,'g.');
% %     plot(goON(e,1),e,'c.');
% %     plot(SR(e,1),e,'m.');
% %     plot(rwON(e,1),e,'b.');
% %     plot(stop(e,1),e,'k.');
% % end
% % xlabel('time (s)')
% % ylabel('trail number')
% % title('all trial with critical time points')
% % axis tight

%% rasterplot of one trial
% % trial = 15;
% % figure; hold on
% % for k = 1:numel(spikeTrains{trial,1})
% %     plot(spikeTrains{trial,1}{k,1},k*ones(numel(spikeTrains{trial,1}{k,1}),1), 'k.');
% % end
% % title(['trial ' num2str(trial) ', all units'])
% % xlabel('time (s)')
% % ylabel('unit ID')
% % 
% % yAxis = [1, k];
% % jetColor = jet(15);
% % p{1,1} = plot(tsON(trial,1)*[1,1], yAxis,'Color',jetColor(1,:),'LineWidth',2);
% % p{2,1} = plot(wsON(trial,1)*[1,1], yAxis,'Color',jetColor(2,:),'LineWidth',2);
% % p{3,1} = plot(cueON(trial,1)*[1,1], yAxis,'Color',jetColor(3,:),'LineWidth',2);
% % p{4,1} = plot(cueOFF(trial,1)*[1,1], yAxis,'Color',jetColor(4,:),'LineWidth',2);
% % p{5,1} = plot(goON(trial,1)*[1,1], yAxis,'Color',jetColor(5,:),'LineWidth',2);
% % p{6,1} = plot(SR(trial,1)*[1,1], yAxis,'Color',jetColor(6,:),'LineWidth',2);
% % p{7,1} = plot(OT(trial,1)*[1,1], yAxis,'Color',jetColor(7,:),'LineWidth',2);
% % p{8,1} = plot(DO(trial,1)*[1,1], yAxis,'Color',jetColor(8,:),'LineWidth',2);
% % p{9,1} = plot(rwON(trial,1)*[1,1], yAxis,'Color', jetColor(9,:),'LineWidth',2);
% % p{10,1} = plot(stop(trial,1)*[1,1], yAxis,'Color', jetColor(10,:),'LineWidth',2);
% % 
% % p_ind = find(~cellfun(@isempty, p));
% % legend([p{p_ind,1}], cellstr(evt2collect(p_ind,:)), 'Location', 'northeastoutside')
% % axis tight

%% re-organize the spike train/unit id

% get trail type
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

% % figure; hold on;
% % trial = 15;
% % for i = 1:unitNum
% %     spk2plot = spikeTrains_chID_order{trial,1}{i,1};
% %     plot( spk2plot, i*ones(numel(spk2plot),1),'k.');
% % end
% % title(['trial ' num2str(trial) ', units are reordered by channel ID on Utah array'])
% % xlabel('time (s)')
% % ylabel('unit ID')
% % 
% % yAxis = [1, unitNum];
% % jetColor = jet(10);
% % p{1,1} = plot(tsON(trial,1)*[1,1], yAxis,'Color',jetColor(1,:),'LineWidth',2);
% % p{2,1} = plot(wsON(trial,1)*[1,1], yAxis,'Color',jetColor(2,:),'LineWidth',2);
% % p{3,1} = plot(cueON(trial,1)*[1,1], yAxis,'Color',jetColor(3,:),'LineWidth',2);
% % p{4,1} = plot(cueOFF(trial,1)*[1,1], yAxis,'Color',jetColor(4,:),'LineWidth',2);
% % p{5,1} = plot(goON(trial,1)*[1,1], yAxis,'Color',jetColor(5,:),'LineWidth',2);
% % p{6,1} = plot(SR(trial,1)*[1,1], yAxis,'Color',jetColor(6,:),'LineWidth',2);
% % p{7,1} = plot(OT(trial,1)*[1,1], yAxis,'Color',jetColor(7,:),'LineWidth',2);
% % p{8,1} = plot(DO(trial,1)*[1,1], yAxis,'Color',jetColor(8,:),'LineWidth',2);
% % p{9,1} = plot(rwON(trial,1)*[1,1], yAxis,'Color', jetColor(9,:),'LineWidth',2);
% % p{10,1} = plot(stop(trial,1)*[1,1], yAxis,'Color', jetColor(10,:),'LineWidth',2);
% % 
% % p_ind = find(~cellfun(@isempty, p));
% % legend([p{p_ind,1}], cellstr(evt2collect(p_ind,:)), 'Location', 'northeastoutside')
% % axis tight


%% psth
trial = 1;
range = [cueOFF, goON];    %sec
binSize = 0.025; %sec

trailType_num_merge = trailType_num;
trailType_num_merge( trailType_num_merge==3) = 1;   % merge the condition 1 and 3 (both PG)
trailType_num_merge( trailType_num_merge==4) = 2;


spkTrain_inRange = cell(unitNum,numel(unique(trailType_num_merge)));
for i = 1:trialNum
    spkTra_sinlge_trial = spikeTrains_chID_order{i,1};  % all spk train in a single trial
    for j = 1:unitNum
        spkTrain_temp = spkTra_sinlge_trial{j,1};   %single unit spike train
        spkTrain_temp = spkTrain_temp(spkTrain_temp>=range(i,1) & spkTrain_temp<range(i,2));
        spkTrain_inRange{j,trailType_num_merge(i)} = [spkTrain_inRange{j,trailType_num_merge(i)}, spkTrain_temp];
    end
end

% use the fix bin size
psth_inRange = zeros(unitNum, 1);
psth_length = numel(range(trial,1):binSize:range(trial,2))-1;   % the -1 is because the psth function use the range as edge in histcount function, and the gap is 1 element less than edges


% the psth here is calculate across the same trial type in the required
% range. The format is this:
% [unit1-condition1 unit1-condition2 unit1-condition3 unit1-condition4;
%  unit2-condition1 unit2-condition2 unit2-condition3 unit2-condition4;
%  unit3-condition1 unit3-condition2 unit3-condition3 unit3-condition4;
%  .......
%                                                     unitN-condition4];
% 
% This part can only work if the range is a fix duration. If not, the psth
% will crash because the psth function reports different numbers (within a 
% fix binsize)

for i = 1:unitNum
    for j = 1:numel(unique(trailType_num_merge))  % loop through trial type (conditions)
        psth_inRange(i,(j-1)*psth_length+1:j*psth_length) = psth(spkTrain_inRange{i,j}, range(trial,:), binSize);
    end
end

% % % pca: manual eig calculation
psth_inRange_centered = psth_inRange - mean(psth_inRange,2);

C = cov(psth_inRange_centered');
[eigVec, eigVal] = eig(C);
nPC = 10;    % number of PC
PCZ   = eigVec(:,end-nPC+1:end)'*psth_inRange_centered; % extract first n PCs
PCX   = reshape( PCZ, nPC, numel(unique(trailType_num_merge)), psth_length);  % reshape back into three-d array: (# pcs) x (# conditions) x (# psth)
score = eigVec/psth_inRange_centered';  % ??

% create color map for plotting
mp = colormap;
mp = mp( round(linspace(1,64,6)), : );

% plot first K principal components (in order)
figure;
K = 4;
for i=1:K
    subplot(floor(sqrt(K)), ceil( K/floor(sqrt(K))),i)
    hold on;
    for k = 1:numel(unique(trailType_num_merge))     % loop through 4 conditions
        plot( 1:psth_length, squeeze( PCX((nPC-i+1), k, :)));
    end
    title(['PC' num2str(i)])
end
suptitle(['data' num2str(fileID)]);


%%
[~,mrFet,~] = pca(psth_inRange', 'Centered', 1, 'NumComponents', 3);

subGroupNum = 2;
group = ones(psth_length,2);
group = group.*[1:subGroupNum];
group = reshape(group, 1, numel(group));

colors = ['k', 'r', 'b', 'g', 'm', 'c'];
figure;
subplot(2,2,1)
hold on;
for i = 1:subGroupNum
    plot(mrFet(group==i,1), mrFet(group==i,2), '.', 'Color', colors(i))
end
xlabel('PC1')
ylabel('PC2')

subplot(2,2,2)
% plot(mrFet(:,1), mrFet(:,3), 'k.')
hold on;
for i = 1:subGroupNum
    plot(mrFet(group==i,1), mrFet(group==i,3), '.', 'Color', colors(i))
end
xlabel('PC1')
ylabel('PC3')

subplot(2,2,3)
% plot(mrFet(:,2), mrFet(:,3), 'k.')
hold on;
for i = 1:subGroupNum
    plot(mrFet(group==i,2), mrFet(group==i,3), '.', 'Color', colors(i))
end
xlabel('PC2')
ylabel('PC3')

subplot(2,2,4)
% plot3(mrFet(:,1), mrFet(:,2), mrFet(:,3),'k.')
hold on;
for i = 1:subGroupNum
    plot3(mrFet(group==i,1), mrFet(group==i,2), mrFet(group==i,3), '.', 'Color', colors(i))
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
view(3)

suptitle(['data' num2str(fileID)]);

% % % % %
% % figure;
% % plot3(mrFet(:,1), mrFet(:,2), mrFet(:,3),'k.')
% % xlabel('PC1')
% % ylabel('PC2')
% % zlabel('PC3')
% % title(['data' num2str(fileID)]);
% % 




%%
% % %% psth
% % trial = 15;
% % range = [cueON, cueOFF];    %sec
% % binSize = 0.030; %sec
% %
% % psth_inRange = cell(unitNum,1);
% % for p = 1:unitNum
% %     spkTrain_temp = spikeTrains_chID_order{trial,1}{p,1};
% %     psth_inRange{p,1} = psth(spkTrain_temp, range(trial,:), binSize);
% % end
% %
% % figure; hold on;
% % for i = 1:156
% %     plot(psth_inRange{i,1}, i*ones(numel(psth_inRange{i,1}),1), 'k.')
% % end


