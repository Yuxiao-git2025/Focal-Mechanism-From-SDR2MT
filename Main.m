% =========================================================================
%                            -- Part I --
% =========================================================================
%% Initial Parameters
clc;
s0=120;
d0=45;
r0=120;
fprintf('# Initial Value: \n');
fprintf('# S/D/R= %.1f/%.1f/%.1f \n\n',s0,d0,r0);
%% SDR to TPB
[t,p,b]=MEC_SDR2TPB(s0,d0,r0);
[s1,d1,r1]=MEC_TPB2SDR(t,p,b);
fprintf('# From TPB Axis: \n');
fprintf('# S/D/R= %.1f/%.1f/%.1f \n\n',s1,d1,r1);

%% SDR to MT
[MT]=MEC_SDR2MT(s0,d0,r0);
[s2,d2,r2]=MEC_MT2SDR(MT);
fprintf('# From Harvard moment tensor: \n');
fprintf('# S/D/R= %.1f/%.1f/%.1f \n\n',s2,d2,r2);



% =========================================================================
%                            -- Part II --
% =========================================================================
%% MT Decomposition
% We present:
% The first input form is Harvard CMT
load('Data/MT_Harvard.mat');
ResultMT1=MEC_MTDecomposition(MT_Example,'harvard');

% The second input form is NED (e.g., Vavrycuk)
load('Data/MT_Bohemia.mat');
ResultMT2=MEC_MTDecomposition(MT_Bohemia,'ned');

%% Plot BeachBall
% The SDR value obtained through MT decomposition does not indicate which
% nodal plane it corresponds to, so the color zones here are meaningless
figure; tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
ax1=nexttile;
MEC_PlotBall2D(ResultMT2.strike2(1),ResultMT2.dip2(1),ResultMT2.rake2(1),'full',600,ax1);
title('Show only one');
ax2=nexttile;
MEC_PlotBall2DMore(ResultMT2.strike2,ResultMT2.dip2,ResultMT2.rake2,600,ax2);
title('Total');
set(gcf,'Position',[400,85,900,600]);

%% Plot Radiation (in NEZ axis)
figure; tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
for i=1:2
    axi=nexttile;
    MEC_PlotPRadiation(MT_Bohemia(i,:),30,axi);
end
set(gcf,'Position',[400,85,900,600]);

%% Plot Diamond
MEC_PlotDiamond(ResultMT2.clvd,ResultMT2.iso);

%% Plot Triangle
load("Data\SDR_Test.mat");
Result=MEC_PlotTriangleSDR(SDR_Test,'PlotOdd',1);


% =========================================================================
%                            -- Part III --
% =========================================================================
%% Identification (EXPL, CLAP, EQs)
% Columns are:
% [MNN, MEE, MUU, MEU, MNU, MNE] or [M11, M22, M33, M23, M13, M12]
% *The first six events are North Korean nuclear tests (2006–2017),
% *The seventh is a cavity collapse following the 2017 test,
% *The last two are South Korean earthquakes (2016 and 2017)
load('Data/MT_Test.mat');
IDResult=MEC_MTIdentify(MT_Test,'ned');

%% Plots
ResultMT3=MEC_MTDecomposition(MT_Test,'ned');
figure; tiledlayout(3,3,"TileSpacing","compact","Padding","compact");
for i=1:size(MT_Test,1)
    ax=nexttile;
    MEC_PlotBall2D(ResultMT3.strike2(i),ResultMT3.dip2(i),ResultMT3.rake2(i),'full',600,ax);
    title([IDResult.LikeEvent{i}]);
end
