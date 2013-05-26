function [eventHistogram,simBySimCummulativeEventCountsActualTime,...
    survival,multiSimEventCountsActualTime,multiSimEventCountsPatientTime,...
    hCount,pValues]=MCImucDay(armSurvivals, NumSimulations, targetEvents)
%function [ ]=MCImuc(armSurvivals, NumSimulations, targetEvents);
% inputs:
%   armSurvivals:   a vector of the median survivals for the treatment arms
%                   i.e. [placebo_med_OS_months, active_med_OS_months]
%   NumSimulations: the number of Monte Carlo simulations to run
%   TargetEvents: the target number of events to use to predict time of
%                   occurance

    daysPerMonth = 365.25/12;

    if nargin<3,
        targetEvents = 32;
    end;
    
    [survivalReferenceRate, survivalReferenceTimes] = getReferenceSurvival();

    % convert the reference curve to a daily OS curve
    survRef = StretchSurvival(survivalReferenceRate, daysPerMonth, survivalReferenceTimes);
    % get median OS of reference curve in months
    m0 = find(survRef<=0.5,1)/daysPerMonth; 
    
    mA = armSurvivals(1);
    mB = armSurvivals(2);
    
    % stretch the survival curves by (median survival)*daysPerMonth/m0
    survA = StretchSurvival(survivalReferenceRate,mA*daysPerMonth/m0,survivalReferenceTimes); % control
    survB = StretchSurvival(survivalReferenceRate,mB*daysPerMonth/m0,survivalReferenceTimes); % active
    
    survival = [survA,survB(1:length(survA))]; % make the survival curves the same length of time

    %pts = 278;
    % now get the enrollment numbers and the probabilities of assignment to
    % each arm
    enrollment = [2,2,4,1,4,9,14,12,18,17,12,24,25,20,23,25,24,25,17];
    adj=0; % adjust the target events to fit actual target enrollment
    arms = 0.4424*[1 2]/3; % placebo and active rates
    
    [eventHistogram,simBySimCummulativeEventCountsActualTime,...
        multiSimEventCountsActualTime,multiSimEventCountsPatientTime,...
        hCount,pValues] = ...
        MCSimDay(arms, enrollment, survival,NumSimulations, targetEvents, adj);
end

function [survivalReferenceRate, survivalReferenceTimes] = getReferenceSurvival()
    % this function provides the survival curve control points for
    % interpolation
    if 0, % this is my old chart eyeballed from graph
        survivalReferenceRate = [0.9800,0.9600,0.9400,0.9200,0.9000,0.8600,0.8200,...
        0.7800,0.7400,0.7000,0.6600,0.6200,0.5800,0.5400,0.5000,0.4640,...
        0.4280,0.3920,0.3560,0.3200,0.3040,0.2880,0.2720,0.2560,0.2400,...
        0.2240,0.2080,0.1920,0.1760,0.1600,0.1560,0.1520,0.1480,0.1440,...
        0.1400,0.1360,0.1320,0.1280,0.1240,0.1200,0.1160,0.1120,0.1080,...
        0.1040,0.1000,0.0960,0.0960,0.0960,0.0960,0.0960,0.0960,0.0960,...
        0.0960,0.0960,0.0960,0.0960,0.0960,0.0960,0.0960,0.0960]'; % 15 mo OS curve from 2005 Strupp
        survivalReferenceTimes = 1:length(survivalReferenceRate);
    else
        % this model fixes the observation that the original numbers
        % appeared to be piecewise linearly interpolated and that the
        % interpolation control points might be better as curve reference
        % points and then fit with a smooth i.e. curve like a cubic spline.
        survivalReferenceRate = [.9,.5,.32,.16,.12,.096]'; 
        survivalReferenceTimes = [5 15 20 30 40 46];
    end;
end

function s = StretchSurvival(survival, scale, x)
    % this function uses splines to interpolate the survival curves to a
    % specified scale
    L = max(x);
    s = interp1([0,x]'*scale,[1;survival],[0:L*scale]','spline');
    s = s(2:length(s));
end