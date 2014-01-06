function kernel = ADFKernel( numscales )

kernel = legion.Kernel;
kernel.addEnviro( @pitt.Depend.CCAadd );
kernel.add( @pitt.exp.StationaryScales.runKPSS, 'X', numscales );

end
