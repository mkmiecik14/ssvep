% Epoching and Final Cleaning with MARA
% Matt Kmiecik
% Started 3 Feb 2020

vis_workspace_prep % Prepares workspace (see src/...)

% Initializations ---- (change these paths specific to your machine)
ica_data_path =     'C:\Users\pains\Desktop\matt-eeglab\data\3-icaweights\';
prepro_data_path =  'C:\Users\pains\Desktop\matt-eeglab\data\2-prepro\';
outpath =           'C:\Users\pains\Desktop\matt-eeglab\data\4-epochs\';

% Gathers subject numbers from subject info excel read ----
subjs = string({RAW{2:size(RAW,1),1}}); %  {'352'}

% For loop to process participants
for i = 1:length(subjs)

    % Creating variables ----
    this_subj = dir(strcat(prepro_data_path, "vis-", subjs(i), "*.set")); % grabs subject info
    this_subj.id = strcat('vis-',char(subjs(i))); % eeglab file name

    % Load in participant preprocessed and visually inspected file ----
    EEG = pop_loadset('filename',this_subj.name,'filepath',prepro_data_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    % Import ICA results (weights) and settings ----
    load(fullfile(ica_data_path, strcat(this_subj.id, '-ica-res.mat'))); %loads
    EEG.icachansind = res.icachansind; % assigns the channels used in ICA
    EEG = pop_editset(EEG, 'icaweights', res.icaweights); % assigns ICA weights
    eeg_checkset(EEG);

    % Dropping ICA components ----
    % res.comps2drop uses the threshold defined in ica_mara.m and reported
    % in res.thresh
    % use res.artcomps to use the default MARA threshold {.50}
    EEG = pop_subcomp(EEG, res.comps2drop, 0); 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on', 'gui', 'off');

    % If there were any bad channels identified during inspection, they are
    % imported here:
    interpchans = str2num(RAW{i+1,2}); % gathers bad channels from excel
    
     % Interpolating any bad channels ----
    if sum(size(interpchans)) > 0 % only interpolates if there are bad channels
        disp('Bad channel(s) detected...interpolating now...');
        EEG = pop_interp(EEG, interpchans, 'spherical');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on', 'gui', 'off');
        EEG.icachansind = res.icachansind; % assigns the channels used in ICA (needs to be set again)
    else
        disp("No channels to interpolate...");
    end

    % Extracting epochs surrounding each stimulation block [-7 22]
    % This extracts 7 seconds before each stimulation = baseline
    % 22 seconds after stimulation = stimulation
    % However, the events are only placed between -6 to 21 seconds; a
    % little buffer was needed because later epoching was dropping epochs
    % due to boundary events
    EEG = pop_rmdat( EEG, {'S  1' 'S  2' 'S  3' 'S  4' 'S  5'},[-7 22] ,0); % original was [-6 21]
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',strcat(this_subj.id, '-s-epochs'),'overwrite','on','gui','off');
    EEG = eeg_checkset( EEG );

    % Retrieving the latency (in seconds) of each block start
    s1_lat = (eeg_getepochevent(EEG, 'S  1', [], 'latency'))/1000;
    s2_lat = (eeg_getepochevent(EEG, 'S  2', [], 'latency'))/1000;
    s3_lat = (eeg_getepochevent(EEG, 'S  3', [], 'latency'))/1000;
    s4_lat = (eeg_getepochevent(EEG, 'S  4', [], 'latency'))/1000;
    s5_lat = (eeg_getepochevent(EEG, 'S  5', [], 'latency'))/1000;
    
    % Creates new events for baseline and stim
    % events are purposely spaced 1 second apart for overlap
    stim_lats = {s1_lat, s2_lat, s3_lat, s4_lat, s5_lat};
    for k = 1:length(stim_lats)
        
        this_lat = cell2mat(stim_lats(k)); % init loop latency
        this_baseline = strcat('S', num2str(k),'_bl'); % init baseline label
        this_stim = strcat('S', num2str(k), '_stim'); % init stim label
        
        nevents = length(EEG.event); % calculates the number of events
        
        for j = 1:5
            EEG.event(nevents+j).latency = (this_lat-7+j)*EEG.srate + 1;  % latency of marker
            EEG.event(nevents+j).type = this_baseline;                    % baseline
            EEG.event(nevents+j).duration = 1;                            % arbitrary length
        end
        
        % Check for consistency and reorder the events chronologically
        EEG = eeg_checkset(EEG, 'eventconsistency');
        
        nevents = length(EEG.event); % calculates the new number of events
        
        % Loop for stimulation
        % 20 2sec windows starting at stimulation : 1 sec overlap
        for j = 1:20
            EEG.event(nevents+j).latency = (this_lat-1+j)*EEG.srate + 1;  % latency of marker
            EEG.event(nevents+j).type = this_stim;                        % new marker name
            EEG.event(nevents+j).duration = 1;                            % arbitrary length
        end
        % Check for consistency and reorder the events chronologically
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end
    
    % Extracts 2s epochs for FFT ----
    trigs = {'S1_bl' 'S1_stim' 'S2_bl' 'S2_stim' 'S3_bl' 'S3_stim' 'S4_bl' 'S4_stim' 'S5_bl' 'S5_stim'};
    EEG = pop_epoch( EEG, trigs, [0 2], 'newname', strcat(this_subj.id, '-epochs'), 'epochinfo', 'yes'); %2
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

    % Epoch thresholding ----
    % threshold limits are -100 to +100 microvolts
    EEG = pop_eegthresh(EEG,1,[1:32],-100,100,0,2,0,0); % last argument toggles marking
    % More criteria can be placed here

    % Epoch rejection ----
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);   % All the 1's mean that all rej criteria are considered
    EEG = pop_rejepoch( EEG, EEG.reject.rejthresh, 0);      % EEG.reject.rejthresh is the trials exceeding threshold
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',strcat(this_subj.id, '-epochs'),'overwrite','on','gui','off');

    % Save out epoched data ----
    EEG = eeg_checkset( EEG );
    outname = strcat(this_subj.id,'-epochs.set');
    EEG = pop_saveset( EEG, 'filename', outname, 'filepath' ,outpath);
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    % Decided not to save out trials that were rejected
    % each subject should have 5 baseline and 20 stim epochs
    % you can calculate the % dropped in each one in the spectral analysis
    % script

end

eeglab('redraw'); % updates GUI
