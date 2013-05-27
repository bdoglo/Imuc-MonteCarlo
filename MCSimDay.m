function [eventHistogram,simBySimCummulativeEventCountsActualTime,...
    multiSimEventCountsActualTime,multiSimEventCountsPatientTime,...
    rejectH0Count,pValues] = MCSimDay(arms, enrollment, survival, N, E, adj)
%
% calculate the event rates for the trial arms by performing N trial simulations
% 
%   arms:   vector of the rates at which enrolled patients are randomized to various
%       trial arms
%   enrollment: vector with the period by period estimated enrollment
%   survival: each column is the survival rates period by period for each arm respectively
%   pts: the number of patients in simulation

    daysPerMonth = 365.25/12;
    if nargin<5,
        E = 32;
    end;
    if nargin<4,
        N = 5000;
    end;
    if nargin<6,
        adj = 0;
    end;
    
    figure(1); plot(survival); grid; 
    title('OS Curves'); 
    legend('placebo','active');
    xlabel('days'); ylabel('probability');
    T = 124; % the number of pts we know were actually on trial

    simBySimCummulativeEventCountsActualTime = zeros(N,length(arms));
    eventHistogram = zeros(N,1);
    multiSimEventCountsActualTime = 0;
    multiSimEventCountsPatientTime = 0;
    rejectH0Count = 0;
    pValues=zeros(N,1);
    
    % repeat the MC simulation N times
    for i=1:N,
        % perform MC simulation of survival of all the arm assigned
        % enrolled patients
        [eventCountsActualTime,tp,eventCountsPatientTime,tPatient,...
            tActual,armPatient] = SimSurvival(arms,enrollment,survival);
        
        multiSimEventCountsPatientTime = multiSimEventCountsPatientTime + eventCountsPatientTime;
        multiSimEventCountsActualTime = multiSimEventCountsActualTime + eventCountsActualTime;
        % get the cummulative event count for each arm by period
        cummulativeEventCountsActualTime = cumsum(eventCountsActualTime);
        % get the total event count acroos all arms by period
        f1 = sum(cummulativeEventCountsActualTime,2);

        % adjust the number of target events 
        % proportional to the simulated number that ended up on trial
        targetEvents = E; 
        if adj
            % I'm thinking it is better to relax this normalization and
            % skip it and just let the final randomized count fall where it
            % may even though it doesn't match what we know to be true.
            % The mean will match and the normalization introduces
            % artifacts.  We should still use the same target threshold
            % though
            targetEvents = E*tp/T;
        end;
        % find the period when the number of total events exceeds the
        % target number of events (this is the interim analysis period)
        idx = find(f1>targetEvents,1);
        
        % save the interim analysis trigger time for this simulation
        eventHistogram(i) = idx;
        % save the event counts at this trigger time for each arm
        simBySimCummulativeEventCountsActualTime(i,:) = cummulativeEventCountsActualTime(idx,:);

        
        % p-value calculation
        % find patients surviving beyond the E-th event
        sele = tActual>idx;
        tp = tPatient;
        ta = tActual;
        % set their survival to the E-th event time
        tp(sele) = tPatient(sele)-(tActual(sele)-idx);
        ta(sele) = idx;

        % do we want to limit the results to those with over a certain time
        % of evaluation?
        tStart = tActual-tPatient;
        sel = tStart<idx-2*365.25; % must have started 2 or more years before

        A = length(arms);   % number of trial arms
        y = tp(armPatient~=A & sel);    % control patient survival times (up to E-th event)
        x = tp(armPatient==A & sel);    % active patient survival times (up to E-th event)
    
        %figure(3);
        %survivalBar(tPatient,tActual,armPatient,idx);
        
        ctrl_mean = mean(y);
        ctrl_sigma = std(y);

        pValue = rightSidedZtest(x,ctrl_mean,ctrl_sigma);
        rejectH0 = sum(pValue<=0.05);

        % the above two lines should produce the same as the statistics tool
        % box using the below function
        %[rejectH0,pValue] = ztest(x,ctrl_mean,ctrl_sigma,'Tail','right');

        rejectH0Count = rejectH0Count + rejectH0;
        pValues(i) = pValue;

        % show a counter once in a while so the user knows something is
        % happening
        if (mod(i,500)==0),
            fprintf('simulation %d\n',i);
        end;
    end;
    figure(2);
    hist(eventHistogram,min(eventHistogram):max(eventHistogram)); grid;
    title(sprintf('Histogram of %d event analysis period from %d simulations',targetEvents,N));
    ylabel('count');
    xlabel('days since trial start');
    mu = mean(eventHistogram);
    sigma = std(eventHistogram);
    med = median(eventHistogram);
    
    fprintf('mean = %g  (%d mo %d days)\n',mu,floor(mu/daysPerMonth),round(rem(mu,daysPerMonth)));
    fprintf('median = %g  (%d mo %d days)\n',med,floor(med/daysPerMonth),round(rem(med,daysPerMonth)));
    fprintf('std-dev = %g  (%d mo %d days)\n',sigma,floor(sigma/daysPerMonth),round(rem(sigma,daysPerMonth)));  
end

function [eventCountsActualTime,total,eventCountsPatientTime,tPatient,tActual,armPatient]...
    = SimSurvival(arms, enrollment, survival)
%
% calculate the event rates for the trial arms
% 
%   arms:   vector of the rates at which enrolled patients are randomized to various
%       trial arms
%   enrollment: the month by month estimated enrollment
%   survival: each column is the survival rates day by day for each arm respectively
%   pts: the number of patients in simulation
    T = size(survival,1); % number of days in survival curves
    A = length(arms);   % number of trial arms
    E = length(enrollment); % number of enrollment months
    daysPerMonth = 365.25/12;
    
    % some counters
    eventCountsActualTime = zeros(ceil(E*daysPerMonth)+T+1,A);  % counter of events for each day in the survival curve time domain
    eventCountsPatientTime = eventCountsActualTime;
    
    cummulativeArmAssignment = [0 cumsum(arms)]; % get the cummulative arm assignment
    mortality = cumsum(diff(1-survival)); % get the mortality probability from the survival curve
    total = 0;
    
    tPatient = zeros(1,0);
    armPatient = zeros(1,0);
    tActual = zeros(1,0);
    for i=0:E-1, % for each enrollment month
        % randomize enrollment
        r = rand(enrollment(i+1),1);
        for a=1:A,
            % take all the enrolled people for this month and see if they  
            % were randomized to this arm --> count how many were
            cnt = sum(r>=cummulativeArmAssignment(a) & r<cummulativeArmAssignment(a+1));
            % for each pt assigned to arm
            if cnt>0,
                % this part should help speed up the loop a little 
                tPatient(total+cnt) = 0; % resize x
                armPatient(total+cnt) = 0; % resize y
            end;
            for k=1:cnt,
                % get the start time in days that they were added
                rr = rand(1); % get randomly (uniform) the fraction of the month when they were added
                actualStartTime = (i+rr)*daysPerMonth; % start time is the current month + fraction --> converted to days
                actualStartTimeRoundedUp = ceil(actualStartTime); % then rounded up to the next day.
                
                % get the random survival time for the patient
                rt = rand(1); % random number for 0 to 1
                patientSurvivalTime = find (rt<mortality(:,a),1); % find the first occurance of mortality at that random time from 
                % there is a chance that the patient dies beyond the
                % mobidity curves, here we check that
                if isempty(patientSurvivalTime),
                    % it is unclear what to do at this point - my current
                    % solution is to pick the time at the end of the curve
                    % data.
                    % The issue is that the curve is about 5 years long and
                    % most of the time there won't be data that far out for
                    % some time after the trial is over.
                    patientSurvivalTime = length(mortality);  
                end;
                % count the deaths in the event counts in actual time and
                % in patient survivalTime (from patient start)
                if actualStartTimeRoundedUp+patientSurvivalTime <= length(eventCountsActualTime),
                    %sprintf('month %d, arm %d, ==> %d/%d   [%d,%d]\n',i,a,k,cnt,idx,idx+i)
                    eventCountsActualTime(actualStartTimeRoundedUp+patientSurvivalTime,a) = ...
                        eventCountsActualTime(actualStartTimeRoundedUp+patientSurvivalTime,a) + 1;
                    eventCountsPatientTime(patientSurvivalTime,a) = eventCountsPatientTime(patientSurvivalTime,a) + 1;
                end;
                armPatient(total+k) = a;
                tPatient(total+k) = patientSurvivalTime;
                tActual(total+k) = actualStartTimeRoundedUp+patientSurvivalTime;
            end;
            total = cnt + total; % count the total randomized patients
        end;
    end;
end

function z= rightSidedZtest(x, mu, sigma)
    % one-sided (right)
    z = 1-(erf((mean(x)-mu)/(sqrt(2/length(x))*sigma))/2 + 0.5);
end    
    

