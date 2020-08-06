% Plot EEG from DBS data files
% Imran Quraishi
% Yale Comprehensive Epilepsy Center

% Load JSON files. Save last directory accessed.
pathfile = fullfile(pwd, '.jsondir');
if exist(pathfile, 'file')
    jsondir = strip(fileread(pathfile));
else
    jsondir = [];
end
jsondir = uigetdir(jsondir,'Select folder with DBS data');

fid = fopen(pathfile,'wt');
fprintf(fid,  '%s', jsondir);
fclose(fid);

filelist = dir(fullfile(jsondir, '**\*.json'));
filelist = fullfile(string({filelist.folder}), string({filelist.name}));

%filelist = filelist(5);

close all

for jsonfile = filelist
    %disp(jsonfile)
    dbsdata = jsondecode(fileread(jsonfile));
    patient.name = ([dbsdata.PatientInformation.Final.PatientFirstName, ' ', dbsdata.PatientInformation.Final.PatientLastName]);
    
    % "LFP Montage Time Domain" - ~30s snapshots of EEG obtained during
    % programming session
    
    lfp = [];
    if isfield(dbsdata, 'LfpMontageTimeDomain')
        lfp = [lfp;dbsdata.LfpMontageTimeDomain];
    end
    if isfield(dbsdata, 'SenseChannelTests')
        lfp = [lfp;dbsdata.SenseChannelTests];
    end
    
    if numel(lfp)>0
        
        figure
        set(gcf,'renderer','painter')
        
        passuniquename = cellstr(string({lfp.Pass;lfp.FirstPacketDateTime}').join);
        [lfp(:).passuniquename] = deal(passuniquename{:});
        passes = unique({lfp.passuniquename});
        numpasses = length(unique(passuniquename));
        totalchannels = length({lfp.Channel});
        channels_per_pass = ceil(length({lfp.Channel})/numpasses);
        index = reshape(1:totalchannels, numpasses, channels_per_pass).';
        clf
        ax = [];
        
        for passnum = 1:numpasses
            passname = passes(passnum);
            lfp_onepass = lfp(strcmp({lfp.passuniquename},passname),:);
            start_time = datetime(lfp_onepass(1).FirstPacketDateTime, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
            start_time.TimeZone = 'local';
            fs = lfp_onepass(1).SampleRateInHz;
            dt = seconds(1/fs);
            ns = length(lfp_onepass(1).TimeDomainData); % assume all are same length
            t = start_time + (0:dt:((ns-1)*dt));
            
            for k = 1:length({lfp_onepass.Channel})
                subplot(length({lfp_onepass.Channel}),numpasses,index(k,passnum));
                axis;
                a=gca;
                
                % Label plot
                channel_label = replace(lfp_onepass(k).Channel,'_','-');
                if contains(channel_label, 'LEFT')
                    side = 'L';
                else
                    side = 'R';
                end
                channel_label = replace(channel_label, {'ZERO','ONE','TWO','THREE'}, {'0','1','2','3'});
                channel_label = erase(channel_label, 'AND-');
                channel_label_parts = split(channel_label,'-');
                channel_label = [side, channel_label_parts{1}, '-', side, channel_label_parts{2}];
                ylabel(a,channel_label);
                ylim([-200 200])
                
                hold on
                plot(a,t,lfp_onepass(k).TimeDomainData);
                hold off
                
                ax(k,passnum) = gca;
                
                % Optional: make graph larger
                ax1 = gca;
                outerpos = ax1.OuterPosition;
                ti = ax1.TightInset;
                left = outerpos(1) + ti(1);
                bottom = outerpos(2) + ti(2);
                ax_width = outerpos(3) - ti(1) - ti(3);
                ax_height = outerpos(4) - ti(2) - ti(4);
                ax1.Position = [left bottom ax_width ax_height];
                
            end
            
            linkaxes(ax(:,passnum), 'xy'); % this command is slow
            clear lfp_onepass
            
        end
        
        sgtitle([patient.name, ' - Signal check of all electrode combinations - ' datestr(start_time, 'mm/dd/yy')])
        clear lfp
        
    end
    
    % "LFP Trend Log" - magnitude of signal at specified frequency, tracked
    % over time
    
    if isfield(dbsdata.DiagnosticData, 'LFPTrendLogs')
        
        lognames = unique([fields(dbsdata.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left), fields(dbsdata.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left)]);
        numlogs = length(lognames);
        if numlogs > 0
            trendtable = timetable;
            
            for k = 1:numlogs
                t = datetime({dbsdata.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left.(lognames{k}).DateTime}, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
                t.TimeZone = 'local';
                magleft = [dbsdata.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Left.(lognames{k}).LFP];
                magright = [dbsdata.DiagnosticData.LFPTrendLogs.HemisphereLocationDef_Right.(lognames{k}).LFP];
                trendtable = [trendtable; timetable(t',magleft',magright','VariableNames',{'Left','Right'})];
            end
            figure
            trendtable = sortrows(trendtable);
            s = stackedplot(trendtable,'GridVisible',true);
            yl=max(s.AxesProperties(1).YLimits(2),s.AxesProperties(2).YLimits(2));
            s.AxesProperties(1).YLimits = [0 yl];
            s.AxesProperties(2).YLimits = [0 yl];
            title([patient.name ' Trend'])
            
        end
    end
    
    % "BrainSense Time Domain" - EEG streamed in clinic
    
    if isfield(dbsdata, 'BrainSenseTimeDomain')
        subplot(3,1,1)
        start_time = datetime(dbsdata.BrainSenseTimeDomain(1).FirstPacketDateTime, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
        start_time.TimeZone = 'local';
        fs = dbsdata.BrainSenseTimeDomain(1).SampleRateInHz;
        dt = seconds(1/fs);
        ns = length(dbsdata.BrainSenseTimeDomain(1).TimeDomainData); % assume both sides are same length
        t = start_time + (0:dt:((ns-1)*dt));
        trendtable = timetable(t', dbsdata.BrainSenseTimeDomain(1).TimeDomainData, dbsdata.BrainSenseTimeDomain(2).TimeDomainData,'VariableNames',{'Left','Right'});
        s = stackedplot(trendtable);
        s.AxesProperties(1).YLimits = [-100 100];
        s.AxesProperties(2).YLimits = [-100 100];
        s.XLimits = [start_time, start_time + seconds(10)];

        subplot(3,1,2)
        s = stackedplot(trendtable);
        s.AxesProperties(1).YLimits = [-100 100];
        s.AxesProperties(2).YLimits = [-100 100];
        a2=gca;
        
        subplot(6,1,5)
        spectrogram(dbsdata.BrainSenseTimeDomain(1).TimeDomainData,kaiser(256,5),50,0:1:125,250,'yaxis')
        a5=gca;
        title('Left')
        
        subplot(6,1,6)
        spectrogram(dbsdata.BrainSenseTimeDomain(2).TimeDomainData,kaiser(256,5),50,0:1:125,250,'yaxis')
        title('Right')
        a6=gca;
        sgtitle([patient.name ' EEG streamed in clinic'])
        
        linkaxes([a5,a6],'x')
    end
    
    
end