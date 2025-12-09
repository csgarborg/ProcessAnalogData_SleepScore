clc; clear all; close all;

%% ------------------------------------------------------------------------
%  STEP 1 — PRE-PROCESS ANALOG TDMS DATA
%
%  This step generates:
%      • ProcData.mat files (per recording)
%      • Cleaned & downsampled ECoG, EMG, Force traces
%      • EMG power, EMG 10–100 Hz signal
%      • binForceSensor (thresholded movement)
%      • Spectrograms (SpecData.mat)
%
%  Usage:
%      Batch_Preprocess_TDMS_AnalogData;                  
%          → Normal run: uses existing spectrograms and force thresholds.
%
%      Batch_Preprocess_TDMS_AnalogData(true, false);
%          → Recompute ONLY the spectrograms.
%
     Batch_Preprocess_TDMS_AnalogData(true, true);
%          → Fully recompute spectrograms AND force sensor processing
%            (threshold GUI will pop up again).
%
%      Batch_Preprocess_TDMS_AnalogData(false, true);
%          → Recompute ONLY the force sensor + threshold.
%
%  Notes:
%      • This step must be run before baseline normalization.
%      • Processed files are saved as: <AnimalID>_<Timestamp>_ProcData.mat
% ------------------------------------------------------------------------
% Batch_Preprocess_TDMS_AnalogData;


%% ------------------------------------------------------------------------
%  STEP 2 — BASELINE NORMALIZATION (ANALOG)
%
%  This step:
%      • Prompts for a manual "baseline window" per file (once per file)
%      • Automatically detects resting epochs (binForce == 0, 2–20 s long)
%      • Computes ONE resting baseline per Animal × Day:
%            mean/std for:
%                 – forceSensor
%                 – EMG power
%                 – EMG 10–100 Hz signal (if present)
%                 – ECoG_DS
%                 – Spectrogram (per-frequency baseline)
%      • Writes normalized fields BACK into each ProcData file:
%            forceSensor_norm
%            EMG.emgPower_norm
%            EMG.emgSignal_norm
%            ECoG_DS_norm
%            SpecData.ECoG.normS and normS_dB
%
%  Default call:
%      BaselineNormalization_Analog;
%
%  Optional:
%      BaselineNormalization_Analog('ForceReEntry', true);
%         → Forces re-selection of manual baseline windows,
%           even if previously saved.
%
%  Notes:
%      • Must be run AFTER Step 1 (preprocessing).
%      • Does NOT create new files; it appends normalized data into
%        existing *_ProcData.mat and *_SpecData.mat files.
% ------------------------------------------------------------------------
clc; close all;
fprintf('\n=========================================\n');
fprintf(' Running Baseline Normalization (Analog) \n');
fprintf(' Folder: %s\n', pwd);
fprintf('=========================================\n\n');

BaselineNormalization_Analog;

fprintf('\n-----------------------------------------\n');
fprintf(' Baseline Normalization COMPLETED\n');
fprintf('-----------------------------------------\n\n');
