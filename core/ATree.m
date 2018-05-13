
%The data structure used in the code
%Atree         : the combined formula which is a data struct has fields as
               %follows:
%               -tree: the structure of the formula
%               -rule: the sequence of rule used to generate the formula
%                  which uses index to denotes the rule as follows:
%                      1:A-->  \square A
%                      2:A-->   \diamond A
%                      3:A-->A v B
%                      4:A--> A ^ B
%                      5:A--> \mu
%               -time: the time bound sequence in the formula
%               -dir:  indicate < />= in predicate 
%                            -1---  <
%                            1---  >=
%               -mag:  the magnitude sequence of the predicate
%               -operator: the temporal operator sequence in the formula which use             
%                          index to denote  the operator
%                         1---  \square
%                         2---- \diamond

% By Gang Chen
% date: 1/8/2018






