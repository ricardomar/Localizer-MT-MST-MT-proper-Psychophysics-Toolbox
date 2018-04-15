
%--------------%
% MT localizer %
%--------------%

% Syncobx?
syncbox = 0;

% Close (eventually) open connections and PTB screens
IOPort('CloseAll');
Screen('CloseAll');

% Trick suggested by the PTB authors to avoid synchronization/calibration
% problems
figure(1)
plot(sin(0:0.1:3.14));
% Close figure with sin plot (PTB authors trick for synchronization)
close Figure 1

% Synchronization tests procedure - PTB
% Do you want to skipsync tests (1) or not (0) ?
skipsynctests = 0; 

% KbName will switch its internal naming
% scheme from the operating system specific scheme (which was used in
% the old Psychtoolboxes on MacOS-9 and on Windows) to the MacOS-X
% naming scheme, thereby allowing to use one common naming scheme for
% all operating systems
KbName('UnifyKeyNames');
% Code to identify "escape" key
escapekeycode = KbName('ESCAPE');

% open syncbox connection
if syncbox
    syncbox_handle = IOPort('OpenSerialPort', 'COM2', 'BaudRate=57600 DataBits=8 Parity=None StopBits=1 FlowControl=None');
    IOPort('Flush',syncbox_handle);
end


%----------------------------------------------------%

AssertOpenGL;

try

    % ------------------------
    % set dot field parameters
    % ------------------------

    nframes     = 21360; % number of animation frames in loop
    mon_width   = 39;   % horizontal dimension of viewable screen (cm)
    v_dist      = 60;   % viewing distance (cm)
    dot_speed   = 5; %7;    % dot speed (deg/sec)
    f_kill      = 0.00; % fraction of dots to kill each frame (limited lifetime)
    ndots       = 100; %2000; % number of dots
    max_d       = 5;%15;   % maximum radius of  annulus (degrees)
    min_d       = 0; %1;    % minumum (degrees)
    dot_w       = 0.1;  % width of dot (deg)
    fix_r       = 0.09; %0.15; % radius of fixation point (deg)
    waitframes = 1; % Show new dot-images at each waitframes'th monitor refresh.  'waitframes' Number of video refresh intervals to show each image before updating the dot field. Defaults to 1 if omitted.
    
  
    
    % ---------------
    % My Parameters
    % ---------------

    offset_left=10;  % (degrees)
    offset_right=10; % (degrees)
    offset_center=0;  % (degrees)
    
    efr=60; % estimated-target frame rate (Hz)
    TR=2;   % MRI TR (seconds)
    
    my_protocol=zeros(1,nframes);
    
    % rest
    my_protocol(1,1 : 5*TR*efr)=0;
    my_protocol(1, 173*TR*efr+1 : 178*TR*efr)=0;
    
    % movC
    my_protocol(1, 5*TR*efr+1 : 14*TR*efr)=1;
    my_protocol(1, 75*TR*efr+1 : 84*TR*efr)=1;
    my_protocol(1, 117*TR*efr+1 : 126*TR*efr)=1;
    my_protocol(1, 159*TR*efr+1 : 168*TR*efr)=1;
        
    % staC
    my_protocol(1, 14*TR*efr+1 : 19*TR*efr)=-1;
    my_protocol(1, 84*TR*efr+1 : 89*TR*efr)=-1;
    my_protocol(1, 126*TR*efr+1 : 131*TR*efr)=-1;
    my_protocol(1, 168*TR*efr+1 : 173*TR*efr)=-1;  
    
    % movL
    my_protocol(1, 33*TR*efr+1 : 42*TR*efr)=2;
    my_protocol(1, 61*TR*efr+1 : 70*TR*efr)=2;
    my_protocol(1, 89*TR*efr+1 : 98*TR*efr)=2;
    my_protocol(1, 145*TR*efr+1 : 154*TR*efr)=2;
        
    % staL
    my_protocol(1, 42*TR*efr+1 : 47*TR*efr)=-2;
    my_protocol(1, 70*TR*efr+1 : 75*TR*efr)=-2;
    my_protocol(1, 98*TR*efr+1 : 103*TR*efr)=-2;
    my_protocol(1, 154*TR*efr+1 : 159*TR*efr)=-2;
    
    % movR
    my_protocol(1, 19*TR*efr+1 : 28*TR*efr)=3;
    my_protocol(1, 47*TR*efr+1 : 56*TR*efr)=3;
    my_protocol(1, 103*TR*efr+1 : 112*TR*efr)=3;
    my_protocol(1, 131*TR*efr+1 : 140*TR*efr)=3;

    % staR
    my_protocol(1, 28*TR*efr+1 : 33*TR*efr)=-3;
    my_protocol(1, 56*TR*efr+1 : 61*TR*efr)=-3;
    my_protocol(1, 112*TR*efr+1 : 117*TR*efr)=-3;
    my_protocol(1, 141*TR*efr+1 : 145*TR*efr)=-3;     

    
    % ---------------
    % open the screen
    % ---------------

    screens=Screen('Screens');
	screenNumber=max(screens);
    [w, rect] = Screen('OpenWindow', screenNumber, 0);
    
    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [center(1), center(2)] = RectCenter(rect);
	fps=Screen('FrameRate',w);      % frames per second
    ifi=Screen('GetFlipInterval', w);
    if fps==0
       fps=1/ifi;
    end;
    
    white = WhiteIndex(w);
    HideCursor;	% Hide the mouse cursor
    Priority(MaxPriority(w));
    
    % Do initial flip...
    vbl=Screen('Flip', w);
    
    % ---------------------------------------
    % initialize dot positions and velocities
    % ---------------------------------------

    ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
    pfs = dot_speed * ppd / fps;                            % dot speed (pixels/frame)
    s = dot_w * ppd; % dot size (pixels)
    smin=s;
    smax=s+s;
    fix_cord = [center-fix_r*ppd center+fix_r*ppd];

    rmax = max_d * ppd;	% maximum radius of annulus (pixels from center)
    rmin = min_d * ppd; % minimum
    
    ms=(smax-smin)/(rmax-rmin);
    
    % IN #####################
    r_IN = rmax * sqrt(rand(ndots,1));	% r
    r_IN(r_IN<rmin) = rmin;
    t_IN = 2*pi*rand(ndots,1);                     % theta polar coordinate
    cs_IN = [cos(t_IN), sin(t_IN)];
    xy_IN = [r_IN r_IN] .* cs_IN;   % dot positions in Cartesian coordinates (pixels from center)

    %mdir = 2 * floor(rand(ndots,1)+0.5) - 1;    % motion direction (in or out) for each dot
    mdirIN=ones(ndots,1)-2;
    drIN = pfs * mdirIN;                            % change in radius per frame (pixels)
    dxdyIN = [drIN drIN] .* cs_IN;                       % change in x and y per frame (pixels)
    
    
    
    % OUT #####################
    r_EXT = rmax * sqrt(rand(ndots,1));	% r
    r_EXT(r_EXT<rmin) = rmin;
    t_EXT = 2*pi*rand(ndots,1);                     % theta polar coordinate
    cs_EXT = [cos(t_EXT), sin(t_EXT)];
    xy_EXT = [r_EXT r_EXT] .* cs_EXT;   % dot positions in Cartesian coordinates (pixels from center)

    %mdir = 2 * floor(rand(ndots,1)+0.5) - 1;    % motion direction (in or out) for each dot
    mdirEXT=ones(ndots,1);
    drEXT = pfs * mdirEXT;                            % change in radius per frame (pixels)
    dxdyEXT = [drEXT drEXT] .* cs_EXT;                       % change in x and y per frame (pixels)
    

    % Create a vector with different colors for each single dot, if
    % requested:
    colvect=white;

    %%
    Screen('FillOval', w, uint8(white), fix_cord);	% draw fixation dot (flip erases it)
    vbl=Screen('Flip', w);
     
    % Trigger
    % syncbox is 'on' or 'off'?
    if syncbox==1
        [gotTrigger, logdata.triggerTimeStamp] = waitForTrigger(syncbox_handle,1,300);
        
        if gotTrigger
            disp('Trigger OK! Starting stimulation...');
            
            % Clean syncbox buffer
            IOPort('Flush', syncbox_handle);
            IOPort('Purge', syncbox_handle);
            IOPort('Close', syncbox_handle);
            
        else
            disp('Absent trigger. Aborting...');
            throw(me)
        end
        
        [keyIsDown, secs, keyCode] = KbCheck;
        
        % The user asked to exit the program
        if keyIsDown==1 && keyCode(escapekeycode)
            
            % Close PTB screen and connections
            Screen('CloseAll');
            IOPort('CloseAll');
            ShowCursor;
            Priority(0);
            
            % Launch window with warning of early end of program
            warndlg('The run was terminated with ''Esc'' before the end!','Warning','modal')
            
            return % abort program
        end
        
    else
        WaitSecs(1)
        KbWait;
    end
    %%
    
    % --------------
    % animation loop
    % --------------    

    for i = 1:nframes
        if (i>1)
            Screen('FillOval', w, uint8(white), fix_cord);	% draw fixation dot (flip erases it)
               
            if my_protocol(i)==0
                % do nothing
                
            elseif my_protocol(i)==1
                my_xymatrix(1,:)=xymatrix(1,:)+offset_center*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN + dxdyIN;						% move dots
                r_IN = r_IN + drIN;							% update polar coordinates too
                xy_EXT = xy_EXT + dxdyEXT;				% move dots
                r_EXT = r_EXT + drEXT;					% update polar coordinates too
                
            
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif my_protocol(i)==-1
                my_xymatrix(1,:)=xymatrix(1,:)+offset_center*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN;  						% move dots
                r_IN = r_IN;							% update polar coordinates too
                xy_EXT = xy_EXT;      				% move dots
                r_EXT = r_EXT;   		    			% update polar coordinates too         
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif my_protocol(i)==2
                my_xymatrix(1,:)=xymatrix(1,:)-offset_left*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN + dxdyIN;						% move dots
                r_IN = r_IN + drIN;							% update polar coordinates too
                xy_EXT = xy_EXT + dxdyEXT;				% move dots
                r_EXT = r_EXT + drEXT;					% update polar coordinates too  
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif my_protocol(i)==-2
                my_xymatrix(1,:)=xymatrix(1,:)-offset_left*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN;						% move dots
                r_IN = r_IN;						% update polar coordinates too
                xy_EXT = xy_EXT;	     			% move dots
                r_EXT = r_EXT;					% update polar coordinates too     
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif my_protocol(i)==3
                my_xymatrix(1,:)=xymatrix(1,:)+offset_right*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN + dxdyIN;						% move dots
                r_IN = r_IN + drIN;							% update polar coordinates too
                xy_EXT = xy_EXT + dxdyEXT;				% move dots
                r_EXT = r_EXT + drEXT;					% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            elseif my_protocol(i)==-3
                my_xymatrix(1,:)=xymatrix(1,:)+offset_right*ppd;
                my_xymatrix(2,:)=xymatrix(2,:);
                xy_IN = xy_IN;						% move dots
                r_IN = r_IN;							% update polar coordinates too
                xy_EXT = xy_EXT;				% move dots
                r_EXT = r_EXT;					% update polar coordinates too
                Screen('DrawDots', w, my_xymatrix, my_s, colvect, center,1);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            end

        end
        
        
        [keyIsDown, secs, keyCode] = KbCheck;
        
        % The user asked to exit the program
        if keyIsDown==1 && keyCode(escapekeycode)
            
            % Close PTB screen and connections
            Screen('CloseAll');
            IOPort('CloseAll');
            ShowCursor;
            Priority(0);
            
            % Launch window with warning of early end of program
            warndlg('The run was terminated with ''Esc'' before the end!','Warning','modal')
            
            return % abort program
        end
                
        
        % check to see which dots have gone beyond the borders of the annuli
        r_out = find(r_EXT > rmax | r_EXT < rmin | rand(ndots,1) < f_kill);	% dots to reposition
        nout = length(r_out);
        if nout

            % choose new coordinates
            r_EXT(r_out) = rmin; %* sqrt(rand(nout,1));
            r_EXT(r_EXT<rmin) = rmin;
            t_EXT(r_out) = 2*pi*(rand(nout,1));

            % now convert the polar coordinates to Cartesian
            cs_EXT(r_out,:) = [cos(t_EXT(r_out)), sin(t_EXT(r_out))];
            xy_EXT(r_out,:) = [r_EXT(r_out) r_EXT(r_out)] .* cs_EXT(r_out,:);

            % compute the new cartesian velocities
            dxdyEXT(r_out,:) = [drEXT(r_out) drEXT(r_out)] .* cs_EXT(r_out,:);

        end;
        xymatrix_EXT = transpose(xy_EXT);
        
        
        
        
        
        % check to see which dots have gone beyond the borders of the annuli
        r_out = find(r_IN > rmax | r_IN < rmin | rand(ndots,1) < f_kill);	% dots to reposition
        nout = length(r_out);
        if nout

            % choose new coordinates
            r_IN(r_out) = rmax; %* sqrt(rand(nout,1));
            r_IN(r_IN<rmin) = rmax;
            t_IN(r_out) = 2*pi*(rand(nout,1));

            % now convert the polar coordinates to Cartesian
            cs_IN(r_out,:) = [cos(t_IN(r_out)), sin(t_IN(r_out))];
            xy_IN(r_out,:) = [r_IN(r_out) r_IN(r_out)] .* cs_IN(r_out,:);

            % compute the new cartesian velocities
            dxdyIN(r_out,:) = [drIN(r_out) drIN(r_out)] .* cs_IN(r_out,:);
        end;
        xymatrix_IN = transpose(xy_IN);
        
        
        
        my_s(1:ndots)=ms*r_IN'+smin;
        my_s(ndots+1:2*ndots)=ms*r_EXT'+smin;
        
        % RICARDO - record video
        %imageArrayN=Screen('GetImage', w, [], [], [], []);
        %imwrite(imageArrayN,[num2str(1000+i),'.png'],'png');
        % RICARDO - record video 
        
         % merge IN and OUT dots
        xymatrix=[xymatrix_IN xymatrix_EXT];
        vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
        


    end;
    
    Priority(0);
    ShowCursor
    Screen('CloseAll');
    
catch
    
    Priority(0);
    ShowCursor
    Screen('CloseAll');

end



