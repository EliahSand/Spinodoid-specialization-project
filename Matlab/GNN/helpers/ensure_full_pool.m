function pool = ensure_full_pool()
c = parcluster('local');
targetWorkers = c.NumWorkers;
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
