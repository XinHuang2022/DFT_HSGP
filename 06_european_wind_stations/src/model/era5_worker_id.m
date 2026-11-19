function id = era5_worker_id()
%ERA5_WORKER_ID Stable diagnostic worker identifier (zero when serial).
    id = 0;
    try
        item = getCurrentTask();
        if ~isempty(item), id = item.ID; end
    catch
    end
end
