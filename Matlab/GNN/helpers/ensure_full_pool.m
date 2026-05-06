function pool = ensure_full_pool()
    c = parcluster('local');
    targetWorkers = min(6, c.NumWorkers);  % if you dont want to use
    %everything the computer has
    %targetWorkers = c.NumWorkers;

    pool = gcp('nocreate');

    if ~isempty(pool) && pool.NumWorkers == targetWorkers
        return;
    end

    if ~isempty(pool)
        delete(pool);
    end

    try
        pool = parpool(c, targetWorkers);
    catch ME
        warning('Could not start pool with %d workers (%s). Falling back to default pool size.', ...
            targetWorkers, ME.message);
        pool = parpool(c);
    end
end