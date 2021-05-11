% ICA Decompositions with MARA identifying and marking artifactual ICs
% Matt Kmiecik
% Started 24 Jan 2020

tic; % starts timer

vis_workspace_prep % Prepares workspace (see src/...)

% Initializations ---- (change these paths specific for your machine)
prepro_data_path = 'C:\Users\pains\Desktop\matt-eeglab\data\2-prepro\';
outpath = 'C:\Users\pains\Desktop\matt-eeglab\data\3-icaweights\';

% Gathers subject numbers from subject info excel read
subjs = string({RAW{2:size(RAW,1),1}}); % subjs = {'363'}; % Initializes subjects for batch processing (if applicable)

% ICA ----
for i = 1:length(subjs)
    
    % Creating variables ----
    this_subj = dir(strcat(prepro_data_path, "vis-", subjs(i), "*.set")); % grabs subject info
    this_subj.id = strcat('vis-',char(subjs(i))); % eeglab file name
    outname = strcat(this_subj.id,'-ica-res.mat'); % save out subject name

    % Loads in preprocessed data ----
    EEG = pop_loadset('filename',this_subj.name,'filepath',prepro_data_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

    % Runs ICA ----
    % Channel indices the ica will be performed on
    % Importantly skips the EOG channel, as it doesn't have the same reference
    % as EEG channels
    this_chan_index = str2num(RAW{i+1,3}); % ensures ICA does not perform on bad channels
    this_rank = length(this_chan_index)-1; % computes rank based on good channels (minus 1 bc of linked mastoid)
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on','pca',this_rank, 'chanind', this_chan_index);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    % If there were any bad channels identified during inspection, they are
    % imported here:
    interpchans = str2num(RAW{i+1,2}); % gathers bad channels from excel
    
    % Bad channels must be removed from the dataset in order for MARA to
    % work; MARA must have as many datachannels as components
    if sum(size(interpchans)) > 0 % checks to see if there are bad channels
        disp('Bad channel(s) detected...removing channels for MARA...');
        % origchanlocs = EEG.chanlocs; % for later interpolation (this line
        % not used, but kept because it's an interesting way to interpolate
        % chans)
        EEG = pop_select( EEG, 'nochannel', interpchans);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
    else
        disp('No bad channels detected...')
    end
    
    % The rawest form of MARA ----
    % artcomps = vector of artifactual ICs
    % info = posterior probs and normed features for all ICs
    [artcomps, info] = MARA(EEG); % Performs automatic IC rejection
    
    % MARA's threshold for artifactual ICs is quite liberal = .5 {default}
    % to make this more conservative, = .6, use the posterior priors:
    mara_thresh = .6; % increase this number to be more conservative
    comps2drop = 1:size(EEG.icaweights, 1); % gathers vector of IC comps
    comps2drop = comps2drop(info.posterior_artefactprob > mara_thresh); % use threshold to index
    
    icatime = toc; % time in seconds to finish ICA + MARA
    
    % Saving out results ----
    
    % ICA settings and results
    res.ss = subjs{i};                  % subject number
    res.icaweights = EEG.icaweights;    % ica weights matrix
    res.icachansind = this_chan_index;  % ica channel indices
    res.proctime = icatime;             % time in seconds to finish ICA computations
    
    % MARA settings and results
    res.artcomps = artcomps;            % MARA marked artifactual ICs
    res.info = info;                    % MARA info structure
    res.comps2drop = comps2drop;        % actual components droped using thresh
    res.thresh = mara_thresh;           % threshold used for posterior priors    

    % Saves out results to file ----
    save(fullfile(outpath, outname),'res'); % saves out as matlab struct

end

