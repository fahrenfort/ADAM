function output = cond_string(varargin)
% cond_string creates a comma separated string from the intersection of integers input arrays
% containing event values.
% Useful to convert event definitions from the experimental design to class specifications.
% Optionally you can separate train and test events using a string that contains a semicolon ';',
% but pay close attention to the syntax of this function.
%
% Use as:
%   cond_string(events1,events2,....); OR cond_string([events1,...],';',[events2,...]);
%
% To prevent errors using the syntax of this function, it is recommended to closely follow the
% examples below, and modify these examples to match with your own experimental conditions.
% 
% Example event code definition:
% Two experimental factors, each with three levels: (1) stimulus type and (2) stimulus repetition in
% the trial.
%
% % FACTOR 1: stimulus type
% famous_faces = [5 6 7];                    % specifies triggers of all famous faces
% nonfamous_faces = [13 14 15];              % specifies triggers of all non-famous faces
% scrambled_faces = [17 18 19];              % specifies triggers of all scrambled faces
%
% % FACTOR 2: stimulus repetition
% first_presentation = [5 13 17];            % specifies triggers of all first presentations
% immediate_repeat = [6 14 18];              % specifies triggers of all immediate repeats
% delayed_repeat = [7 15 19];                % specifies triggers of all delayed repeats
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example 1:
%
% class_spec{1} = cond_string(nonfamous_faces);  % the first stimulus class
% class_spec{2} = cond_string(scrambled_faces);  % the second stimulus class
%
% output: 
% class_spec{1} 
%                   '13,14,15'
% class_spec{2}
%                   '17,18,19'
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example 2:
%
% class_spec{1} = cond_string(nonfamous_faces,immediate_repeat);  % the first stimulus class
% class_spec{2} = cond_string(scrambled_faces,immediate_repeat);  % the second stimulus class
%
% output: 
% class_spec{1} 
%                   '14'
% class_spec{2}
%                   '18'
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% It is also possible to create a definition in which separate event types are defined for training
% and testing, by inserting a semicolon in the definition, like this:
% 
% Example 3:
% class_spec{1} = cond_string(famous_faces,first_presentation,';',famous_faces,delayed_repeat);
% -> train on famous first presentations, test on famous delayed repeats
% class_spec{2} = cond_string(nonfamous_faces,first_presentation,';',nonfamous_faces,delayed_repeat);  % the second stimulus class
% -> train on nonfamous first presentations, test on nonfamous delayed repeats
%
% output: 
% class_spec{1} 
%                   '5;7'
% class_spec{2}
%                   '13;15'
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% User function of the ADAM toolbox, by J.J.Fahrenfort, VU, 2015/2018
%
% See also: ADAM_MVPA_FIRSTLEVEL

splitme = find(strcmpi(varargin,';'));
if isempty(splitme)
    output = get_intersection(varargin{:});
else
    output = [ get_intersection(varargin{1:splitme-1}) ';' get_intersection(varargin{splitme+1:end}) ];
end

function output = get_intersection(varargin)
output = varargin{1};
for c = 2:numel(varargin)
    output = intersect(output,varargin{c});
end
output = vec2str(output);