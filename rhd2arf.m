function rhd2arf(dirname, arf_filename, groupname, varargin)

%rhd2arf(DIRNAME,ARF_FILENAME,GROUPNAME) creates an .arf file with name ARF_FILENAME
%from the .rhd files in the directory DIRNAME.  A single dataset will be created
%for each channel, and all of them will be placed in the group GROUPNAME.
%
%rhd2arf(DIRNAME,ARF_FILENAME,GROUPNAME,PULSE_CHANNEL,PRESTIM,POSTSTIM) 
%creates an .arf file with a separate group for each stimulus presentation.
%This requires a channel that records a pulse at the onset of each
%stimulus.  The name of this channel should be specified by PULSE_CHANNEL.
%The base name of each group will be GROUPNAME.
%PRESTIM seconds before each stimulus onset that will be saved in the
%arf file and POSTSTIM seconds after.
%
%The names of the each channels should be used to specify information about
%the data.  If the channel name of the format NAME_DATATYPE.    

disp(length(varargin)==3 && ischar(varargin{1}))

p = inputParser;

addRequired(p,'dirname', @ischar);
addRequired(p,'arf_filename', @ischar);
addRequired(p,'groupname', @ischar);
addOptional(p,'pulse_channel', '', @(x) (length(varargin)==3) && ischar(x))
addOptional(p,'prestim', 0, @(x) length(varargin)==3 && x>=0 && isnumeric(x))
addOptional(p,'poststim', 0, @(x) length(varargin)==3 && x>=0 && isnumeric(x))

parse(p, dirname, arf_filename, groupname, varargin{:});
%converting p.Results into variables
for f = fieldnames(p.Results)' 
    evalc([f{1} '= p.Results.(f{1})']); 
end

if length(varargin) == 3
    separate_trials = true;
else
    separate_trials = false;
end


% rhd2(dirname, arf_filename, groupname) reads all of the intan .rhd files in the
% directory given by "dirname" and saves them in an arf file.
% The names of the channels should be in the following format:
% "DATASETNAME_DATATYPE_CHANNEL_NUMBER." "arf_filename" is the name
% of the arf file in which the data will be saved.  "groupname" is the
% name of the group in the arf file in which the data will be saved.

%% Getting file sizes and if needed stimulus times 
rhdfiles = dir([dirname '/*.rhd']);
% getting total length of recording

% size of tempurature array
temp_size = 0;
%file relative stimulus times
file_rel_stim = cell(length(rhdfiles),1);
dset_size = 0;
for  i = 1:length(rhdfiles)
    [amplifier_channels, amplifier_data, board_adc_channels,...
     board_adc_data,frequency_parameters,notes,spike_triggers,...,
     t_amplifier,t_board_adc,t_temp_sensor,temp_sensor_data]...
        = read_Intan_RHD2000_file([dirname '/' rhdfiles(i).name]);
          
    % ensures that channel names are consistent across rhd files
    if i == 1
        first_channels = [amplifier_channels, board_adc_channels];
        first_spike_trigger = spike_triggers;
        first_freq_params = frequency_parameters;
    else
        if ~isequal([amplifier_channels board_adc_channels],first_channels) || ...
                ~isequal(spike_triggers, first_spike_trigger) || ...
                ~isequal(first_freq_params, frequency_parameters)
            err = MException('rhd2arf:IncompatibleFiles',...
                'Files contain incompatible metadata');
        end
    end
    
    % getting pulse times and createing window to be saved
    if separate_trials
        pulse_ch_idx = find(any(strcmp(pulse_channel, {board_adc_channels.custom_channel_name})));
        if isempty(pulse_ch_idx)
            err = MException('rhd2arf:InvalidArgument',...
                'Pulse channel does not exist in all files');
            throw(err)
        end
        pulse_data = board_adc_data(pulse_ch_idx,:);           
        %Finds the sample of stimulus onset relative to the start of each
        %file.  Samples are found where signal crosses a threshold, and
        %for each crossing the maximum sample in a window following the 
        %threshold crossing is found.     
        thr_factor = 200; %factor number of standard deviations used as threshold to detect pulses
        pulse_winsize = 1000; %number of samples after threshold crossing to look for max
        thr = thr_factor * median(abs(pulse_data)/.6745);%0.6745 converts to estimator of standard deviation
        [~, crossings] = find(pulse_data(1:end-1) < thr ...
                          & pulse_data(2:end) > thr);
        for c = crossings
            [~, pulse_max] = max(pulse_data(c:c+pulse_winsize));
            file_rel_stim{i} = [file_rel_stim{i}, c+pulse_max];
        end
 
    else
        dset_size = dset_size + size(amplifier_data, 2);
    end
    
    % getting size of datasets to pass to arfcreate
    % temp_size = temp_size + size(temp_sensor_data,2);
    
end
if separate_trials
    sr = frequency_parameters.amplifier_sample_rate;
    dset_size = sr*prestim + sr*poststim;
end

%% Creating Datasets


% datasetnames

nstim = sum(cellfun(@length,file_rel_stim));

if separate_trials
    % 2-d cell array for trial-separated arf.  groups along first
    % dimension, channels along second.
    datasetnames = cell(nstim,length(amplifier_channels));
else
    datasetnames = cell(length(amplifier_channels),1);
end
sample_counter = 0;

all_channels = [amplifier_channels, board_adc_channels];
all_data = [amplifier_data;board_adc_data];

% creating datasets and saving attributes
for i = 1:length(all_channels)
    % checking names of channels, which will provide attributes
    channel_args = ...
        strsplit(all_channels(i).custom_channel_name,'_');
    if ~any(length(channel_args) == [2 3])
        warning(['Invalid channel name "' ...
            all_channels(i).custom_channel_name ...
            '". Datatype attribute will be undefined.'])
        ch_name = all_channels(i).custom_channel_name;
        datatype = '0';
        ch_num = '';
    else
        [ch_name, datatype] = channel_args{1:2};
        valid_datatypes = [0:6 1000:1002 2000:2002];
        if ~any(str2double(datatype) == valid_datatypes)
            warning(['Invalid datatype ' datatype ' for channel ' ...
                all_channels(i).custom_channel_name])
        end
        if length(channel_args) == 3
            ch_num = channel_args{3};
        else
            ch_num = '';
        end
    end
    %create datasets, write stimulus start datasets
    if separate_trials
        for group_idx = 1:nstim
            if isempty(ch_num)
                datasetnames{group_idx,i} = ['/' groupname '_' sprintf('%04d',group_idx-1)...
                    '/' ch_name];
            else
                datasetnames{group_idx,i} = ['/' groupname '_' sprintf('%04d',group_idx-1)...
                    '/' ch_name '_' ch_num];
            end
            arfcreate(arf_filename, datasetnames{group_idx,i}, dset_size, ...
                'arf_datatype', str2double(datatype), 'sampling_rate', ...
                frequency_parameters.amplifier_sample_rate);
            %saving channel attributes
            channel_fields = fields(all_channels)';
            for att_name = channel_fields
                arfwriteatt(arf_filename, datasetnames{group_idx,i}, ...
                    att_name{1}, all_channels(i).(att_name{1}))
            end
        end
        
    else
        datasetnames{i} = ['/' groupname '/' ch_name '_' ch_num];
        arfcreate(arf_filename, datasetnames{i}, dset_size, ...
            'arf_datatype', str2double(datatype), 'sampling_rate', ...
            frequency_parameters.amplifier_sample_rate);
        
        %saving channel attributes
        channel_fields = fields(all_channels)';
        for att_name = channel_fields
            arfwriteatt(arf_filename, datasetnames{i}, ...
                att_name{1}, all_channels(i).(att_name{1}))
        end
    end
    
    % saving spike_triggers as attributes if channel is from
    % amplifier_channels
    if i <= length(spike_triggers)
        spkt_fields = fields(spike_triggers)';
        for att_name = spkt_fields
            if separate_trials
                for group_idx = 1:nstim
                    arfwriteatt(arf_filename, datasetnames{group_idx,i}, ...
                        att_name{1}, spike_triggers(i).(att_name{1}));
                end
            else
                arfwriteatt(arf_filename, datasetnames{i}, ...
                    att_name{1}, spike_triggers(i).(att_name{1}));
            end
        end
    end
end

%% Writing Data
temp_data = zeros(2,temp_size);
t_temp = zeros(1,temp_size);
temp_counter = 0;

%trial counter for seperate_trials = true
trial_counter = 1;
% indicates whether to it is necessary to write trial from previous file
previous_stim = false;
for i = 1:length(rhdfiles)
    [amplifier_channels, amplifier_data, board_adc_channels,...
     board_adc_data,frequency_parameters,notes,spike_triggers,...
     t_amplifier,t_board_adc,t_temp_sensor,temp_sensor_data]...
        = read_Intan_RHD2000_file([dirname '/' rhdfiles(i).name]);
    all_data = [amplifier_data;board_adc_data];
      
    % reading tempurature data
    %    temp_idx = temp_counter+1:temp_counter+size(temp_sensor_data,2); 
    %t_temp(temp_idx) = t_temp_sensor;
    %temp_data(:,temp_idx) = temp_sensor_data;
    %temp_counter = temp_counter + size(temp_sensor_data,2);
        
    % saving channel data into arf
    
    if separate_trials
        if previous_stim
            stim_start = file_rel_stim{i-1}(end);
            sr = frequency_parameters.amplifier_sample_rate;
            trial_indices = stim_start + [(-1*round(sr*prestim)-1):(round(sr*poststim)-2)];
            for ch_idx = 1:length(all_channels)
                    h5write(arf_filename, datasetnames{trial_counter,ch_idx}, ...
                        combined_data(ch_idx,trial_indices))
            end
            previous_stim = false;
        end
        
        for stim_start = file_rel_stim{i}
            sr = frequency_parameters.amplifier_sample_rate;
            trial_indices = stim_start + [(-1*round(sr*prestim)-1):(round(sr*poststim)-2)];
            if any(trial_indices < 1)
                %covers the case where part of a trial is in the
                %previous file
                if i == 1
                    err = MException('rhd2arf:InvalidArgument', ...
                        'Not enough pre-trial data on first trial');
                    throw(err)
                end
                combined_data = [previous_data all_data];
                trial_indices = trial_indices + size(previous_data,2);        
                for ch_idx = 1:length(all_channels)
                    h5write(arf_filename, datasetnames{trial_counter,ch_idx}, ...
                        combined_data(ch_idx,trial_indices))
                end
            elseif any(trial_indices > size(all_data,2))
                continue
                previous_stim = true;
            else
                for ch_idx = 1:length(all_channels)
                    h5write(arf_filename, datasetnames{trial_counter,ch_idx}, ...
                        all_data(ch_idx,trial_indices))
                end
            end
            trial_counter = trial_counter + 1;
        end

        
    else
        for ch_idx = 1:length(all_channels)
            h5write(arf_filename, datasetnames{ch_idx}, ...
                all_data(ch_idx,:),sample_counter+1, ...
                size(all_data,2))
        end
    end

sample_counter = sample_counter + size(amplifier_data,2);
previous_data = all_data;
end
% saving group level attributes
%     arfwriteatt(arf_filename, ['/' groupname], 'tempurature', temp_data);
%     arfwriteatt(arf_filename, ['/' groupname], 'time_tempurature', t_temp);

notes_cell = struct2cell(notes)';
if separate_trials
    for group_idx = 1:nstim
        entryname = ['/' groupname '_' sprintf('%04d',group_idx-1)];
        arfwriteatt(arf_filename, entryname, ...
            'notes', strjoin(notes_cell,'\n'));
        freq_param_fields = fields(frequency_parameters)';
        for att_name = freq_param_fields
            arfwriteatt(arf_filename,entryname,att_name{1}, ...
                frequency_parameters.(att_name{1}));
        end
    end
else
    arfwriteatt(arf_filename, ['/' groupname], ...
        'notes', strjoin(notes_cell,'\n'));
    
    freq_param_fields = fields(frequency_parameters)';
    for att_name = freq_param_fields
        arfwriteatt(arf_filename,['/' groupname], att_name{1},...
            frequency_parameters.(att_name{1}));
    end
end