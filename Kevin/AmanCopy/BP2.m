%% Brake Profile 1

function [W_TB]=BP2(t,amp,bias)
% Creates custom brake profile that attempts to replicate braking scenario
% right before a turn:
% initial coast --> brake --> corner(coast) --> coast
%
% Inputs:
%   t; time array
%   amp; total brake torque
%   bias; brake bias 
%
% Outputs:
%   W_TB; brake torque per wheel, 4x1 array
%
%
% K. Barreto
% Version 10/10/21

x = length(t);
i_coast = x/5;
i_coast = round(i_coast);
i_decel = i_coast + round(x/5);
%period = i_decel-i_coast;
W_TB = zeros(4,x);

for i = 1:x
    if i <= i_coast
        %W_TB(:,i ) = [1;1;1;1];
        W_TB(:,i ) = [0;0;0;0];
    elseif (i>i_coast) && (i<=i_decel) 
        famp = -amp*bias*sin((pi)*(i/(i_coast+1)));
        ramp = -amp*(1-bias)*sin((pi)*(i/(i_coast+1)));
        W_TB(:,i ) = [famp;famp;ramp;ramp];
    else
        W_TB(:,i )=[0;0;0;0];
    end
end
end