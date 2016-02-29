function x = pr_response(input,ada)
%
% predicted photoreceptor response 
%
% Normann & Perlman 1979

x = input./(input + ada);