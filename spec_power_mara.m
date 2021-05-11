% Power Spectrum Calculations with and without Surface Laplacian filter
% Matt Kmiecik
% Started 28 Jan 2020

vis_workspace_prep % Prepares workspace (see src/...)

% Initializations ---- (modify these paths specific to your machine)
epoch_data_path =   'C:\Users\pains\Desktop\matt-eeglab\data\4-epochs\';
outpath =           'C:\Users\pains\Desktop\matt-eeglab\data\5-spec-res\';
fig_path =          'C:\Users\pains\Desktop\matt-eeglab\figs\';

% Saving out figures {default = 0 does not save out figures}
% Change to = 1 if you want to save out figures of topo maps
save_figs = 0;

% Initializes subjects for batch processing (if applicable)
subjs = string({RAW{2:size(RAW,1),1}}); %  {'352'}

for i = 1:length(subjs)

    % Creating variables ----
    this_subj = dir(strcat(epoch_data_path, "vis-", subjs(i), "*epochs.set")); % grabs subject info
    this_subj.id = strcat('vis-',char(subjs(i))); % eeglab file name

    % Load in participant epoched and cleaned file ----
    EEG = pop_loadset('filename',this_subj.name,'filepath',this_subj.folder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % Compute surface laplacian spatial filter ----
    % Calculating Surface Laplacian via CSD toolbox functions
    % http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/tutorial.html#PrepareInput

    chan_mont = cell(32,1); % initializes cell array

    % Fills cell array with electrode labels
    for j = 1:size(chan_mont)
        chan_mont(j) = cellstr(EEG.chanlocs(j).labels);
    end

    % Derives spherical coordinates via CSD toolbox fx ExtractMontage()
    csd_mont = ExtractMontage('10-5-System_Mastoids_EGI129.csd', chan_mont);
    % To view: MapMontage(csd_mont)

    [G, H] = GetGH(csd_mont); % Calculates G and H matrices

    % Applies surface laplacian to all epochs (known also as current source
    % density)
    for k = 1:size(EEG.data, 3)
        EEG.data(:,:,k) = CSD(EEG.data(:,:,k), G, H, 1.0e-5, 10); % lambda left at default, 10 = cm head size, so units are microvolt/cm^2
    end
    
    % Spectral power calculations for baseline ----
    bl_labs = {'S1_bl', 'S2_bl', 'S3_bl', 'S4_bl', 'S5_bl'};
    EEG_bl = pop_selectevent(EEG,'type',bl_labs, 'deleteevents','off','deleteepochs','on','invertepochs','off');
    % Computes FFTs
    % Note: given sampling rate (256 Hz) and window length (2 seconds),
    % frequency step should be .5 Hz. Also, Hamming window is used to taper per
    % help spectopo
    % Note: spectra (mean power over epochs), in dB; freqs are in Hz
    [spectra, freqs] = spectopo(EEG_bl.data(:,:,:), 0, EEG_bl.srate, 'winsize', 2*EEG_bl.srate, 'plot', 'off'); % winsize is 2 seconds
    
    % Saving out results to structure ----
    spec_res.bl_spectra = spectra;
    spec_res.bl_freqs = freqs;
    spec_res.bl_trials = size(EEG_bl.data, 3); % Third dimension is number of epochs
    
    % Saving out a figure of this subject's spectral topo
    if save_figs > 0
        figure; pop_spectopo(EEG_bl,1,[],'EEG' ,'freq',[10 12 24 25 50],'freqrange',[0 75],'electrodes','on');
        saveas(gcf, fullfile(fig_path, strcat(this_subj.id,'-','baseline','-mara-csd.png')));
        close; % closes figure
    else
        disp('Baseline topo maps not saved...');
    end
    
    % Spectral power calculations for stimulation blocks ----
    
    % Initializations
    stim_labs = {'S1_stim', 'S2_stim', 'S3_stim', 'S4_stim', 'S5_stim'};
    stim_spectra = zeros(32, 257, length(stim_labs)); % 32 chans x 257 freq bins x 5 stim strengths
    stim_freqs = zeros(257, 1, length(stim_labs)); % 257 freq bins x 1 column of freqs x 5 stim strengths
    stim_trials = zeros(1, length(stim_labs)); % vector to store number of trials in each average (20 = max)
    
    for n = 1:length(stim_labs)
        % Selects epochs of a certain label
        this_EEG = pop_selectevent(EEG,'type',stim_labs(n),'deleteevents','off','deleteepochs','on','invertepochs','off');
        % Computes FFT and averages spectra across epochs
        [spectra, freqs] = spectopo(this_EEG.data(:,:,:), 0, this_EEG.srate, 'winsize', 2*this_EEG.srate,'plot','off'); % winsize is 2 seconds
        stim_spectra(:,:,n) = spectra;              % Stores spectra
        stim_freqs(:,:,n) = freqs;                  % Stores freq bands
        stim_trials(n) = size(this_EEG.data, 3);    % Stores accepted num of epochs
        % Saving out a figure of this subject's spectral topo
        if save_figs > 0
            figure; pop_spectopo(this_EEG,1,[],'EEG' ,'freq',[10 12 24 25 50],'freqrange',[0 75],'electrodes','on');
            saveas(gcf, fullfile(fig_path, strcat(this_subj.id,'-',stim_labs{n},'-mara-csd.png')));
            close; % closes figure
        else
            disp('Stimulation topo maps not saved...');
        end
        
    end
    
    % Storing stimulation results ----
    % stimulation results are stored differently with each index being the
    % stimulation strength in order: 1, 2, 3, 4, 5
    spec_res.stim_spectra   = stim_spectra; % 3D mat of spectra
    spec_res.stim_freqs     = stim_freqs;   % 3D mat of freqs bins
    spec_res.stim_trials    = stim_trials;  % vector of trials

    % Saving out all data ----
    spec_outname = strcat(this_subj.id, '-mara-csd-spectral-res.mat');
    save(fullfile(outpath, spec_outname),'spec_res'); % saves out as matlab struct

end
