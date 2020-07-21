function gfx = SetupGraphicsVariables(sp)
  
    % Only initialise xScale, yscale, zScale if used
    if sp.gfxSteps > 0
        gfx.xScale = gpuArray(sp.dx*(0:sp.nx-1));
        gfx.yScale = gpuArray(sp.dy*(0:sp.ny-1));
        gfx.zScale = gpuArray(sp.dz*(0:sp.nz-1));
        
        gfx.kxScale = (2*pi/sp.lx)*[0:sp.nx/2-1, 0, -sp.nx/2+1:-1];
        gfx.kyScale = (2*pi/sp.ly)*[0:sp.ny/2-1, 0, -sp.ny/2+1:-1];    
        gfx.kzScale = (2*pi/sp.lz)*[0:sp.nz/2-1, 0, -sp.nz/2+1:-1];

    else
        % Not worth effort creating arrays if not used
        % TODO: examine the design logic of the code block below
        gfx.xScale = 0;
        gfx.yScale = 0;
        gfx.zScale = 0;
        
        gfx.kxScale = 0;
        gfx.kyScale = 0;
        gfx.kzScale = 0;
    end
end