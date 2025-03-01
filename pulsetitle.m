function [ tL ] = pulsetitle( varname,PULSE,time,t,titlerun,FREQ )
%pulsetitle returns a figure title based on the variable, value of PULSE,
%and timestep.
%   If PULSE is true (eruption flow is pulsating), the title includes the
%   frequency of pulsing. If PULSE is false, the title indicates steady
%   flow.

    if strcmp(PULSE,'T') == 1
        str = 'Unsteady flow';
        tL =  ({sprintf('%s, t=%d s',varname,time(t));sprintf('%s: %s, %g Hz',titlerun,str,FREQ)});
    elseif strcmp(PULSE,'F') == 1
        str = 'Steady flow';
        tL = ({sprintf('%s, t=%d s',varname,time(t));sprintf('%s: %s',titlerun,str)});
    end

end

