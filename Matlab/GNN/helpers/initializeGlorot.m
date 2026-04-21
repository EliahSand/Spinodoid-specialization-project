function w = initializeGlorot(sz)
    scale = 2 * sqrt(6 / (sz(1) + sz(2)));
    w = scale * (rand(sz) - 0.5);
end
