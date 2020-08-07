function SetupZeroPotential()
    global nx ny nz;
    global V;
    V = zeros(nx, ny, nz, 'gpuArray');
end