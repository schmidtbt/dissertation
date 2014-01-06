function kernel = GCKernel( numScales, numlags )

kernel = legion.Kernel;
kernel.addEnviro( @pitt.Depend.CCAadd );
kernel.add( @pitt.exp.GCScales.runGC, 'X', numScales, numlags );
%kernel.add( @pitt.exp.GCScales.genGCPipe, 'X', numlags );

end