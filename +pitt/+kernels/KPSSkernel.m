function kernel = KPSSkernel( numscales )

kernel = legion.Kernel;
kernel.addEnviro( @pitt.Depend.CCAadd );
kernel.add( @pitt.exp.StationaryScales.runADF, 'X', numscales );

end
