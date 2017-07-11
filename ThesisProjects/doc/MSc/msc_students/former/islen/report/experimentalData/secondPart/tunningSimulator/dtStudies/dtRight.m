clear all;
close all;

%%%% TEST
%  %  matrix = [
%  %  %      2e-3  -2.8900198 4e-4;
%  %  %      4e-3 -2.8901421 3e-4;
%  %  %      6e-3 -2.8902146 2.4e-4];
%     8e-3 -2.8902043 2.e-4
%     1e-2 -2.8901204 1.8e-4
%     ];



%%%%%%%%%%%%%% DT STUDY FOR HE ON TITAN %%%%%%%%%%%%%%%%
%%%% Parameters: Monte Carlo cycles: 	10M
%%%%             Processors: 					4 
%%%%						 minblock:					 	1
%%%%						 maxblock:					 20k 
%%%%						 blocksize:					 500

% % matrix = [
% % 1e-4 		-2.8879435			1.6e-3;
% % 2e-4		-2.8883530			1.2e-3;
% % 3e-4		-2.8886752			1.0e-3;
% % 4e-4		-2.8891601			9.0e-4;
% % 5e-4		-2.8894674			8.0e-4;
% % 6e-4		-2.8894726			7.0e-4;
% % 7e-4		-2.8894645			7.0e-4;
% % 8e-4		-2.8894698			6.0e-4;
% % 9e-4		-2.8895527			6.0e-4;
% % 1e-3		-2.8896462			5.0e-4;
% % ];
% % 
% % %%%% parameter for plotting errorbars
% % hodet = 0.000008;
% % xstick = [1.0e-4:1.0e-4:1e-3];
% % ystick = [-2.8904:0.001:-2.8862];


% % matrix = [
% % 1e-3		-2.8896462			5.5e-4;   
% % 2e-3		-2.8900198			4.0e-4;   
% % 3e-3		-2.8901157			3.5e-4;   
% % 4e-3		-2.8901421			3.0e-4;   
% % 5e-3		-2.8901340			2.6e-4;   
% % 6e-3		-2.8902146			2.4e-4;   
% % 7e-3		-2.8902377			2.2e-4;   
% % 8e-3		-2.8902043			2.0e-4;   
% % 9e-3		-2.8901123			2.0e-4; 
% % 1e-2		-2.8901204			1.8e-4;  
% % ];
% % hodet = 0.00008;
% % xstick = [1.0e-3:2.0e-3:1e-2];
% % ystick = [-2.8904:0.0005:-2.8862];

  
  
% % matrix = [
% % 1e-2		-2.8901204			1.8e-4;  
% % 2e-2		-2.8901928			1.4e-4;  
% % 3e-2		-2.8903048			1.3e-4;  
% % 4e-2		-2.8902038			1.2e-4;  
% % 5e-2		-2.8902791			1.2e-4;  
% % 6e-2		-2.8900703			1.2e-4;  
% % 7e-2		-2.8901753			1.3e-4;  
% % 8e-2		-2.8903643			1.3e-4;  
% % 9e-2		-2.8901008			1.3e-4
% % 1e-1		-2.8900953			1.3e-4;  
% % ];
% % 
% % 
% % hodet = 0.0008;
% % xstick = [1.0e-2:1.0e-2:1e-1];
% % ystick = [-2.8904:0.00008:-2.8862];
% % 


% matrix = [
% 1e-1		-2.8900953			1.3e-4;  
% 2e-1		-2.8900165			2.0e-4;  
% 3e-1		-2.8902094			2.4e-4;  
% 4e-1		-2.8902177			3.0e-4;  
% 5e-1		-2.8900412			4.0e-4;  
% 6e-1		-2.8908504			5.0e-4;  
% 7e-1		-2.8888867			6.0e-4;  
% 8e-1		-2.8895206			7.0e-4;  
% %%%%9e-1		-2.8999409			9.0e-4
% ];
% 
% 
% hodet = 0.00008;
% xstick = [1.0e-3:1.0e-3:1e-2];
% ystick = [-2.8904:0.001:-2.8862];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix = [
%  1e-4 		-2.8879435			1.6e-3;
%  2e-4		-2.8883530			1.2e-3;
%  3e-4		-2.8886752			1.0e-3;
%  4e-4		-2.8891601			9.0e-4;
%  5e-4		-2.8894674			8.0e-4;
%  6e-4		-2.8894726			7.0e-4;
%  7e-4		-2.8894645			7.0e-4;
%  8e-4		-2.8894698			6.0e-4;
%  9e-4		-2.8895527			6.0e-4;
%  1e-3		-2.8896462			5.0e-4;
% % 
% % %%%% parameter for plotting errorbars
% % hodet = 0.000008;
% % xstick = [1.0e-4:1.0e-4:1e-3];
% % ystick = [-2.8904:0.001:-2.8862];

%  2e-3		-2.8900198			4.0e-4;   
%  3e-3		-2.8901157			3.5e-4;   
%  4e-3		-2.8901421			3.0e-4;   
%  5e-3		-2.8901340			2.6e-4;   
%  6e-3		-2.8902146			2.4e-4;   
%  7e-3		-2.8902377			2.2e-4;   
%  8e-3		-2.8902043			2.0e-4;   
%  9e-3		-2.8901123			2.0e-4; 
%  1e-2		-2.8901204			1.8e-4;  

% % hodet = 0.00008;
% % xstick = [1.0e-3:2.0e-3:1e-2];
% % ystick = [-2.8904:0.0005:-2.8862];

  
%  2e-2		-2.8901928			1.4e-4;  
%  3e-2		-2.8903048			1.3e-4;  
%  4e-2		-2.8902038			1.2e-4;  
%  5e-2		-2.8902791			1.2e-4;  
%  6e-2		-2.8900703			1.2e-4;  
%  7e-2		-2.8901753			1.3e-4;  
%  8e-2		-2.8903643			1.3e-4;  
%  9e-2		-2.8901008			1.3e-4
1e-1		-2.8900953			1.3e-4;  

% % hodet = 0.0008;
% % xstick = [1.0e-2:1.0e-2:1e-1];
% % ystick = [-2.8904:0.00008:-2.8862];
% % 

2e-1		-2.8900165			2.0e-4;  
3e-1		-2.8902094			2.4e-4;  
4e-1		-2.8902177			3.0e-4;  
5e-1		-2.8900412			4.0e-4;  
6e-1		-2.8908504			5.0e-4;  
7e-1		-2.8888867			6.0e-4;  
8e-1		-2.8895206			7.0e-4;  
% %%%%9e-1		-2.8999409			9.0e-4
];



dt = matrix(:,1);
E =  matrix(:,2);
err =  matrix(:,3);

%%%hE= ploterr(dt,E,[],err)
%%%%hE=plot(dt,E)

plot(dt,E, 'r', 'LineWidth', 2);
hold on;
hE = errorbar(dt, E, err);
  
 %%% Ajust line properties (functional)
set(hE(1)                            , ...
  'LineStyle'       , 'none'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.3 .3 .3]  );

% 
% %%%% Adjust line properties (esthetics)
set(hE(1)                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 .7 .7]  );
% 
% 
% % % % adjust error bar width
% hE_c  = ...
%     get(hE     , 'Children'    );
% errorbarXData          = ...
%     get(hE_c(2), 'XData'       );
% errorbarXData(4:9:end) = errorbarXData(1:9:end) - hodet;
% errorbarXData(7:9:end) = errorbarXData(1:9:end) - hodet;
% errorbarXData(5:9:end) = errorbarXData(1:9:end) + hodet;
% errorbarXData(8:9:end) = errorbarXData(1:9:end) + hodet;
% set(hE_c(2), 'XData', errorbarXData);



%%% Add legend and labels
hXLabel = xlabel('Time step, dt');
hYLabel = ylabel('Mean energy, E (au)');


%%% Add fonts and axis properties
%  set( gca                       , ...
%     'FontName'   , 'Helvetica' );
%  set([hXLabel, hYLabel], ...
%     'FontName'   , 'AvantGarde');
%  set([hXLabel, hYLabel]  , ...
%     'FontSize'   , 14 );
%  
% set(gca, ...
%  'Box'         , 'on'     , ...
%  'TickDir'     , 'in'     , ...
%  'TickLength'  , [.02 .02] , ...
%  'XMinorTick'  , 'on'      , ...
%  'YMinorTick'  , 'on'      , ...
%  'YGrid'       , 'on'      , ...
%  'XColor'      , [.3 .3 .3], ...
%  'YColor'      , [.3 .3 .3], ...
%  'XTick'       , xstick, ...
%  'YTick'       , ystick);
%   'LineWidth'   , 1         );
%  

%ylim([-2.8907 -2.8895])

