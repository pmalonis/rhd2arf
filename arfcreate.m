function arfcreate(filename, datasetname, size, varargin)
        p = inputParser;
        version = '2.2';
        spec_version = '2.1';
        datatype_codes = [0 1 2 3 4 5 6 23 1000 1001 1002 2000 2001 2002];
        
        % separate parameters passed to h5create and other parameters
        h5params = {};
        arfparams = {};
        h5_param_names = {'Datatype', 'ChunkSize', 'Deflate', 'FillValue', ...
            'Fletcher32', 'Shuffle'};

        for i = 1:2:length(varargin) - 1
            if any(strcmp(varargin{i},h5_param_names))
                h5params = [h5params varargin{i} varargin{i+1}];
            else
                arfparams = [arfparams varargin{i} varargin{i+1}];
            end
        end
         
        if verLessThan('matlab','8.2')
            addParamFunc = @addParamValue;
        else
            addParamFunc = @addParameter;
        end
        
        addParamFunc(p,'timestamp',0,@isnumeric)
        addParamFunc(p,'units','', @isstr)
        addParamFunc(p,'sampling_rate',0,@(x) (ceil(x)==x) & (x>0));
        addParamFunc(p,'arf_datatype',0,@(x) any(x==datatype_codes))
        parse(p,arfparams{:})
        %converting p.Results into variables
        for f = fieldnames(p.Results)' 
            evalc([f{1} '= p.Results.(f{1})']); 
        end

        %verify dataset attributes
        if strcmp(units,'')
            if sampling_rate == 0
                err = MException('arf:InvalidAttribute',...
                    ['Unitless data assumed time series and ' ...
                    'requires sampling_rate attribute.']);
                throw(err)
            end
        elseif strcmp(units,'samples') & sampling_rate == 0
            
            err = MException('arf:InvalidAttribute',...
                ['Data with units of "samples" requires ' ...
                'sampling_rate attribute.']);
            throw(err)
        end

        new_file = false; %indicates if h5create makes new file
        try
            info = h5info(filename);
        catch err
            if strcmp(err.identifier,'MATLAB:imagesci:h5info:libraryError')
                new_file = true;
            else
                rethrow(err) 
            end
        end

        % obtaining list of new groups h5create will create
        new_groups = {}; 
        groups = strsplit(datasetname,'/');
        groups = {groups{1:end-1}};
        if ~new_file
            file_id = H5F.open(filename);
            for i = 1:length(groups)
                group_path = strjoin({groups{1:i}},'/');
                try
                    H5G.close(H5G.open(file_id, group_path));
                catch err
                    if strcmp(err.identifier,'MATLAB:imagesci:hdf5lib:libraryError')
                        new_groups = [new_groups, group_path];
                    else
                        rethrow(err)
                    end
                end
            end
            H5F.close(file_id)
        else
            for i = 1:length(groups)
                group_path = strjoin({groups{1:i}},'/');
                new_groups = [new_groups, group_path];                
            end
        end
        
        % create dataset
        h5create(filename,datasetname,size, h5params{:})
        
        % add attributes to root group if new file
        if new_file
            arfwriteatt(filename,'/','arf_library', 'matlab');
            arfwriteatt(filename,'/','arf_library_version',version);
            arfwriteatt(filename,'/','arf_version',spec_version);
        end
        %add timestamp and uuid to group attributes
        for i = 1:length(new_groups)
            [~,uuid] = system('uuidgen');
            if timestamp == 0
                [~,tstamp] = system('date +%s');
                tstamp = uint64([str2double(tstamp) 0]);
            else
                tstamp = uint64(timestamp);
            end
            arfwriteatt(filename,new_groups{i},'uuid',uuid)
            arfwriteatt(filename,new_groups{i},'timestamp',tstamp)
        end
        
        %add attributes to dataset
        if ~strcmp(units,'')
            arfwriteatt(filename,datasetname,'units',units)
        end
        if sampling_rate > 0
            arfwriteatt(filename,datasetname,'sampling_rate',...
                sampling_rate)
        end 
        arfwriteatt(filename,datasetname,'datatype',...
            uint16(arf_datatype))
        
        
end