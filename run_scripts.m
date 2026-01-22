clear all
clc

 %% run scripts
    server      = 'leonhard-19';
    jobName     = 'job';    
    workdir     = '/Home/ode/kantner/multi-valley-envelope';
    script_arg  = 1;
    
    cmd = sprintf('ssh %s "screen -dmS %s bash -c ''cd %s && matlab -nodisplay -nosplash -batch \\"main_func(%d)\\"''"', ...
                  server, jobName, workdir, script_arg);    

    disp(cmd)      % always debug first
    
    [status, output] = system(cmd);
    
      

 %% run scripts
    server      = 'leonhard-20';
    jobName     = 'job';    
    workdir     = '/Home/ode/kantner/multi-valley-envelope';
    script_arg  = 2;
    
    cmd = sprintf('ssh %s "screen -dmS %s bash -c ''cd %s && matlab -nodisplay -nosplash -batch \\"main_func(%d)\\"''"', ...
                  server, jobName, workdir, script_arg);    

    disp(cmd)      % always debug first
    
    [status, output] = system(cmd);
    

 %% run scripts
    server      = 'leonhard-21';
    jobName     = 'job';    
    workdir     = '/Home/ode/kantner/multi-valley-envelope';
    script_arg  = 3;
    
    cmd = sprintf('ssh %s "screen -dmS %s bash -c ''cd %s && matlab -nodisplay -nosplash -batch \\"main_func(%d)\\"''"', ...
                  server, jobName, workdir, script_arg);    

    disp(cmd)      % always debug first
    
    [status, output] = system(cmd);
    

 %% run scripts
    server      = 'leonhard-22';
    jobName     = 'job';    
    workdir     = '/Home/ode/kantner/multi-valley-envelope';
    script_arg  = 4;
    
    cmd = sprintf('ssh %s "screen -dmS %s bash -c ''cd %s && matlab -nodisplay -nosplash -batch \\"main_func(%d)\\"''"', ...
                  server, jobName, workdir, script_arg);    

    disp(cmd)      % always debug first
    
    [status, output] = system(cmd);
    

 %% run scripts
    server      = 'leonhard-23';
    jobName     = 'job';    
    workdir     = '/Home/ode/kantner/multi-valley-envelope';
    script_arg  = 5;
    
    cmd = sprintf('ssh %s "screen -dmS %s bash -c ''cd %s && matlab -nodisplay -nosplash -batch \\"main_func(%d)\\"''"', ...
                  server, jobName, workdir, script_arg);    

    disp(cmd)      % always debug first
    
    [status, output] = system(cmd);
    

 %% run scripts
    server      = 'leonhard-24';
    jobName     = 'job';    
    workdir     = '/Home/ode/kantner/multi-valley-envelope';
    script_arg  = 6;
    
    cmd = sprintf('ssh %s "screen -dmS %s bash -c ''cd %s && matlab -nodisplay -nosplash -batch \\"main_func(%d)\\"''"', ...
                  server, jobName, workdir, script_arg);    

    disp(cmd)      % always debug first
    
    [status, output] = system(cmd);
        