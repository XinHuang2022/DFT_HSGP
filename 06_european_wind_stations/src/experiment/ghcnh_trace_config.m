function C = ghcnh_trace_config()
%GHCNH_TRACE_CONFIG Prespecified data and Algorithm-2 settings.

    here = fileparts(mfilename('fullpath'));
    C.version = "0.1.0";
    projectRoot=fileparts(fileparts(here));
    C.era5SourceRoot=fullfile(projectRoot,'src','model');
    C.adaptiveSourceRoot=fullfile(projectRoot,'src','adaptive');
    C.year = 2024;
    C.startTime = datetime(2024,6,1,0,0,0,'TimeZone','UTC');
    C.nDays = 90;
    C.synopticHours = [0 6 12 18];
    C.timeToleranceMinutes = 60;
    C.expectedTimePoints = C.nDays*numel(C.synopticHours);

    % Broad European screening window, fixed before inspecting responses.
    C.latitudeLimits = [35 72];
    C.longitudeLimits = [-12 35];
    C.targetStationRange = [500 1000];
    C.minimumSummerRecords = 2000;
    C.downloadCandidateCount = 1200;
    C.pilotDownloadCount = 100;
    C.downloadWorkers = 6;
    C.maximumNearCompleteGaps = 3;

    C.dataRawDir = fullfile(projectRoot,'input','raw_not_distributed');
    C.dataProcessedDir = fullfile(projectRoot,'input','data');
    C.resultsDir = getenv('EUROPEAN_WIND_REPRO_OUTPUT_DIR');
    if isempty(C.resultsDir)
        C.resultsDir = fullfile(projectRoot,'results','current');
    end
    C.stationListFile = fullfile(C.dataRawDir,'ghcnh-station-list.txt');
    C.inventoryFile = fullfile(C.dataRawDir,'ghcnh-inventory.txt');
    C.stationDataDir = fullfile(C.dataRawDir,sprintf('station_year_%d',C.year));
    C.progressLog = fullfile(C.resultsDir,'progress.log');

    C.stationUrl = ['https://www.ncei.noaa.gov/oa/global-historical-', ...
        'climatology-network/hourly/access/by-year/%d/psv/', ...
        'GHCNh_%s_%d.psv'];
    C.stationParquetUrl = ['https://www.ncei.noaa.gov/oa/global-historical-', ...
        'climatology-network/hourly/access/by-year/%d/parquet/', ...
        'GHCNh_%s_%d.parquet'];

    % Algorithm 2 settings. J denotes modes per coordinate and M=J^2.
    C.initialJ = 10;
    C.candidateJ = [8 10 12 15 18 20 24 28 32 36 40 45];
    C.referenceJ = 48;
    C.epsilonKlPerSite = 0.01;
    C.epsilonSensitivity = 0.05;
    C.maxCorrectionFits = 2;
    C.nFolds = 5;
    C.foldSeed = 20261113;
    C.progressEveryStations = 25;
    C.progressEveryFrequencies = 20;
end
