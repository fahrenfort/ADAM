function output = cond_string(varargin)
% cond_string creates a comma separated string from the intersection of integers input arrays.
% Useful to convert trigger definitions from the experimental design to class specifications.
%
% Use as:
%   cond_string(varargin);
% 
% Example:
% Two experimental factors, each with three levels: (1) faces and (2) stimulus repetition in the
% trial.
%
% famous_faces = [5 6 7];                    % specifies triggers of all famous faces
% nonfamous_faces = [13 14 15];              % specifies triggers of all non-famous faces
% scrambled_faces = [17 18 19];              % specifies triggers of all scrambled faces
% first_presentation = [5 13 17];            % specifies triggers of all first presentations
% immediate_repeat = [6 14 18];              % specifies triggers of all immediate repeats
% delayed_repeat = [7 15 19];                % specifies triggers of all delayed repeats
%
% e.g.
% class_spec{1} = cond_string(nonfamous_faces,first_presentation);  % the first stimulus class
% class_spec{2} = cond_string(scrambled_faces,first_presentation);  % the second stimulus class
%
% or
%
% class_spec{1} = cond_string(nonfamous_faces,immediate_repeat);  % the first stimulus class
% class_spec{2} = cond_string(scrambled_faces,immediate_repeat);  % the second stimulus class
%
% User function of the ADAM toolbox, by J.J.Fahrenfort, VU, 2015
%
% See also: ADAM_MVPA_FIRSTLEVEL

output = varargin{1};
for c = 2:numel(varargin)
    output = intersect(output,varargin{c});
end
output = vec2str(output);