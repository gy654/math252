% Calculate your midterm grade in NA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the following lines, please modify the values of:
%
% Q1,Q2,...: to be your quiz **COMPLETION** score 
% which is just 1 if competed 0 if not
%
% H1, H2: to be your **RAW** (e.g not percentage) scores
% on HW1, and HW2 
%
% Mid_exam: to be your **RAW** (e.g not percentage) score
% on the midterm (out of 80)
%
% The values that I have for e.g Q1, H1, ... are just 
% examples that give 100 percent in the class. Note that 
% it is possible to have *more* than 100% in the class
% if you have done Extra Credit (EC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q1 = 1; %0 or 1 depending on if you completed 'Newton's Method Quiz'
Q2 = 1; %0 or 1 depending on if you completed 'Secant Method (and friends) Quiz'
Q3 = 1; %0 or 1 depending on if you completed 'LU factorization quiz'
Q4 = 1; %0 or 1 depending on if you completed 'Norms and Machine Numbers'
Q5 = 1; %0 or 1 depending on if you completed 'Power Method Quiz'

H1 = 25; %numerical score (NOT percentage) on HW1 (max [including EC is] is: 35)
H2 = 37; %numerical score (NOT percentage) on HW2 (max [including EC is] is: 48)

Mid_exam = 80; %numerical score (NOT percentage) on Midterm (max is: 80) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below this line computes your grade, modify 
% at your own risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quiz_avg = mean([Q1 Q2 Q3 Q4 Q5]); % Average quiz scores


H1_max = 25; % max score not includeing Extra Credit (EC)
H2_max = 37; % max score not includeing Extra Credit (EC)

HW_avg = mean([(H1/H1_max) (H2/H2_max)]);

Mid_max = 80; % max on midterm

Mid_avg = Mid_exam/Mid_max; %midterm percentage

Score = 0.1*quiz_avg + 0.4*HW_avg + 0.2*Mid_avg; % Midterm score in class
Max_Score = 0.1 + 0.4 + 0.2; % max score is not 1 because we haven't taken the final
Score_percent = 100*(Score/Max_Score);

disp('Your grade in the class is: ')
if Score_percent > 93
    disp('A')
elseif Score_percent > 90
    disp('A-')
elseif Score_percent > 87
    disp('B+')
elseif Score_percent > 84
    disp('B')
elseif Score_percent > 80
    disp('B-')
elseif Score_percent > 77
    disp('C+')
elseif Score_percent > 74
    disp('C')
elseif Score_percent > 70
    disp('C-')
elseif Score_percent > 67
    disp('D+')
elseif Score_percent > 65
    disp('D')
else 
    disp('F')
end
    
