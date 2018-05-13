% clear all
close all
clear all


results_inner_class  = Innerfault('class');
results_inner_robust = Innerfault('robust');

results_outer_class  = Outerfault('class');
results_outer_robust = Outerfault('robust');

results_ball_class  = Ballfault('class');
results_ball_robust = Ballfault('robust');

results = [{results_inner_class};{results_inner_robust};
          {results_outer_class};{results_outer_robust};
          {results_ball_class};{results_ball_robust}];

time  = clock;
name = [num2str(time(1)),'-',num2str(time(2)),'-',num2str(time(3)),'-',num2str(time(4)),'-', num2str(time(5)),'-', num2str(time(6)),'_results.mat'];
save(name,'results');
