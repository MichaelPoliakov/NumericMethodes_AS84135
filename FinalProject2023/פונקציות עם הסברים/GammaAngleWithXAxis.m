function gamma = GammaAngleWithXAxis(t)
% calculate Gamma angle from the skecth
% INPUT t- time
% OUTPUT Gamma angle

t90 =       2;          % [sec]
gamma_end = 10;         % [deg]
beta =      0.6;        % [1/sec]

gamma =  (t<=t90) .* ( 90                                                ) ... % [deg]
       + (t> t90) .* ( (90-gamma_end) .* exp(-beta.*(t-t90)) + gamma_end );

end