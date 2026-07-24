% =========================================================================
% Plot source focal mechanisms into Frohlich triangle diagram
% Usage:
%   Result = MEC_PlotTriangleSDR(sdr)
%   Result = MEC_PlotTriangleSDR('example_sdr.txt')
%   Result = MEC_PlotTriangleSDR(sdr, 'PlotGrid', true, 'PlotSecGrid', true)
%   Result = MEC_PlotTriangleSDR(sdr, 'PlotBeachball', true, 'SavePath', 'out.png')
%
% Input:
%   sdrInput
%       N x 3 matrix:
%           [strike, dip, rake] in degrees
%       or text file containing at least three numeric columns.
%
% Optional name-value pairs:
%   'PlotGrid'
%       Plot major classification boundaries. Default: true.
%
%   'PlotSecGrid'
%       Plot secondary dip-angle grid. Default: true.
%
%   'PlotOdd'
%       Plot additional odd focal mechanism examples. Default: true.
%
%   'PlotBeachball'
%       Plot representative beachballs at triangle vertices. Default: true.
%
%   'Colors'
%       4 x 3 RGB matrix:
%           row 1: Odd
%           row 2: Strike-slip
%           row 3: Normal
%           row 4: Reverse
%
%   'MarkerSize'
%       Event marker size.
%
%   'BeachballSize'
%       Representative beachball diameter.
%
%   'SavePath'
%       Output image path. Empty means do not save. Default: ''.
%
% Output:
%   Result
%       Struct containing SDR, class, P/T/B-axis dips and triangle coordinates.
%
% Classification:
%   mClass = 0: Odd
%   mClass = 1: Strike-slip
%   mClass = 2: Normal
%   mClass = 3: Reverse
%
% Reference:
%   Frohlich, C. (1992): Triangle diagrams: ternary graphs to display
%   similarity and diversity of earthquake focal mechanisms.
% =========================================================================
function Result=MEC_PlotTriangleSDR(sdrInput, varargin)
p = inputParser;
addParameter(p, 'PlotGrid', true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p, 'PlotSecGrid', true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p, 'PlotOdd', true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p, 'PlotBeachball', true, @(x)islogical(x)||ismember(x,[0,1]));
addParameter(p, 'Colors', defaultColors(), @(x)isnumeric(x)&&isequal(size(x),[4,3]));
addParameter(p, 'MarkerSize', 70, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'BeachballSize', 0.18, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'SavePath', '', @(x)ischar(x)||isstring(x));
parse(p, varargin{:});
opt = p.Results;
coloF = opt.Colors;
sdr = localReadSDR(sdrInput);
if size(sdr,2) < 3
    error('MEC_PlotTriangleSDR:badInput', ...
        'Input SDR must have at least three columns: [strike, dip, rake].');
end
sdr = sdr(:,1:3);
[mClass, dP, dT, dB] = localMechClass(sdr(:,1), sdr(:,2), sdr(:,3));
[h, v] = localPTB2Triangle(dP, dT, dB);
figure('Color', 'w');
ax = axes;
hold(ax, 'on');
vertex = localTriangleVertices();
plot(ax, [vertex.h vertex.h(1)], [vertex.v vertex.v(1)], '-k', 'LineWidth', 1.2);
if opt.PlotSecGrid
    localPlotSecondaryGrid(ax);
end
if opt.PlotGrid
    localPlotMajorGrid(ax);
end

for k = 0:3
    idx = mClass == k;
    if any(idx)
        scatter(ax, h(idx), v(idx), opt.MarkerSize, ...
            'Marker', 'o', ...
            'MarkerFaceColor', coloF(k+1,:), ...
            'MarkerEdgeColor', coloF(k+1,:), ...
            'LineWidth', 0.6);
    end
end

if opt.PlotBeachball
    localPlotRepresentativeBeachballs(ax, vertex, coloF, opt.BeachballSize);
end

if opt.PlotOdd
    localPlotOddBeachballs(ax, coloF, opt.BeachballSize * 0.7);
end

axis(ax, 'equal');
axis(ax, 'off');
hold(ax, 'off');

Result = struct();
Result.SDR = sdr;
Result.Class = mClass;
Result.ClassName = localClassNames(mClass);
Result.dP = dP;
Result.dT = dT;
Result.dB = dB;
Result.TriangleX = h;
Result.TriangleY = v;
Result.Axes = ax;
if strlength(string(opt.SavePath)) > 0
    exportgraphics(gcf, char(opt.SavePath), 'Resolution', 300);
end
end







function colors = defaultColors()
colors = [
    0.20, 0.20, 0.20
    0.85, 0.20, 0.20
    0.20, 0.65, 0.30
    0.20, 0.45, 0.85
];
end


function sdr = localReadSDR(sdrInput)
if isnumeric(sdrInput)
    sdr = sdrInput;
    return;
end
fileName = char(sdrInput);
if exist(fileName, 'file') ~= 2
    error('MEC_PlotTriangleSDR:fileNotFound', ...
        'Cannot find SDR file: %s', fileName);
end
fid = fopen(fileName, 'r');
if fid < 0
    error('MEC_PlotTriangleSDR:fileOpenFailed', ...
        'Cannot open SDR file: %s', fileName);
end
cleaner = onCleanup(@() fclose(fid));
data = textscan(fid, '%f %f %f %[^\n]', 'CommentStyle', '#');
if isempty(data{1})
    error('MEC_PlotTriangleSDR:emptyFile', ...
        'No SDR data found in file: %s', fileName);
end
sdr = [data{1}, data{2}, data{3}];
end


function [mClass, dP, dT, dB] = localMechClass(strike, dip, rake)
% Classification of source focal mechanisms by P/T/B-axis dip angles
strike = strike(:);
dip = dip(:);
rake = rake(:);

n0(:,1) = -sind(dip) .* sind(strike);
n0(:,2) =  sind(dip) .* cosd(strike);
n0(:,3) = -cosd(dip);

u0(:,1) =  cosd(rake) .* cosd(strike) + cosd(dip) .* sind(rake) .* sind(strike);
u0(:,2) =  cosd(rake) .* sind(strike) - cosd(dip) .* sind(rake) .* cosd(strike);
u0(:,3) = -sind(rake) .* sind(dip);

P_osa = localNormalizeRows(n0 - u0);
T_osa = localNormalizeRows(n0 + u0);

P_osa(P_osa(:,3) > 0,:) = -P_osa(P_osa(:,3) > 0,:);
T_osa(T_osa(:,3) > 0,:) = -T_osa(T_osa(:,3) > 0,:);

P_theta = acosd(abs(P_osa(:,3)));
T_theta = acosd(abs(T_osa(:,3)));

dP = 90 - P_theta;
dT = 90 - T_theta;

tmp = 1 - sind(dT).^2 - sind(dP).^2;
tmp(tmp < 0 & tmp > -1e-10) = 0;
dB = asind(real(sqrt(max(tmp, 0))));

mClass = zeros(size(strike));

idxSS = sind(dB).^2 > 0.75;
idxN  = ~idxSS & sind(dP).^2 > 0.75;
idxR  = ~idxSS & ~idxN & sind(dT).^2 > 0.59;

mClass(idxSS) = 1;
mClass(idxN) = 2;
mClass(idxR) = 3;
end

function className = localClassNames(mClass)
className = strings(size(mClass));
className(mClass == 0) = "Odd";
className(mClass == 1) = "Strike-slip";
className(mClass == 2) = "Normal";
className(mClass == 3) = "Reverse";
end


function A = localNormalizeRows(A)
n = sqrt(sum(A.^2, 2));
n(n == 0) = NaN;
A = A ./ n;
end

function [h, v] = localPTB2Triangle(dP, dT, dB)
% Frohlich triangle projection.
dN = 35.26;

z = atan2d(sind(dT), sind(dP)) - 45;

den = sind(dN).*sind(dB) + cosd(dN).*cosd(dB).*cosd(z);

h = cosd(dB).*sind(z) ./ den;
v = (cosd(dN).*sind(dB) - sind(dN).*cosd(dB).*cosd(z)) ./ den;
end


function vertex = localTriangleVertices()
[dP, dT, dB] = deal(90, 0, 0);
[h1, v1] = localPTB2Triangle(dP, dT, dB);

[dP, dT, dB] = deal(0, 90, 0);
[h2, v2] = localPTB2Triangle(dP, dT, dB);

[dP, dT, dB] = deal(0, 0, 90);
[h3, v3] = localPTB2Triangle(dP, dT, dB);

vertex.h = [h1 h2 h3];
vertex.v = [v1 v2 v3];
vertex.normal = [h1 v1];
vertex.reverse = [h2 v2];
vertex.strikeSlip = [h3 v3];
end


function localPlotSecondaryGrid(ax)
Np = 120;
for dd = 20:10:80
    [h1, v1] = localGridCurve('P', dd, Np);
    plot(ax, h1, v1, '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.7);
    text(ax, h1(end), v1(end)-0.08, num2str(dd), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', ...
        'Color', 'k','FontSize',12);

    [h2, v2] = localGridCurve('B', dd, Np);
    [h2, ind] = sort(h2);
    v2 = v2(ind);
    plot(ax, h2, v2, '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.7);
    text(ax, h2(1)-0.1, v2(1), num2str(dd), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right', ...
        'Color', 'k','FontSize',12);

    [h3, v3] = localGridCurve('T', dd, Np);
    plot(ax, h3, v3, '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.7);
    text(ax, h3(1)+0.1, v3(1), num2str(dd), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'Color', 'k','FontSize',12);
end

text(ax, 0, -1.08, '$\delta{P}$', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'Color', 'k','FontSize',18,'Interpreter','latex');

text(ax, -0.9, 0.5, '$\delta{B}$', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'Color', 'k','FontSize',18,'Interpreter','latex');

text(ax, 0.9, 0.5, '$\delta{T}$', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'Color', 'k','FontSize',18,'Interpreter','latex');
Fun_Decorat;
end

function localPlotMajorGrid(ax)
Np = 140;
plot(ax, 0, 0, '+', 'Color', [0.9 0.9 0.9]);
[h1, v1] = localGridCurve('Pmajor', 60, Np);
plot(ax, h1, v1, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
[h2, v2] = localGridCurve('Bmajor', 60, Np);
[h2, ind] = sort(h2);
v2 = v2(ind);
plot(ax, h2, v2, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
[h3, v3] = localGridCurve('Tmajor', 50, Np);
plot(ax, h3, v3, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
end

function [h, v] = localGridCurve(mode, dd, Np)
dP = zeros(1, Np);
dT = zeros(1, Np);
dB = zeros(1, Np);
switch lower(mode)
    case 'p'
        for i = 1:Np
            dP(i) = dd;
            dT(i) = (i-1) * ((90 - dd) / (Np - 1));
            dB(i) = asind(sqrt(max(1 - sind(dP(i))^2 - sind(dT(i))^2, 0)));
        end
    case 'b'
        for i = 1:Np
            dB(i) = dd;
            dT(i) = (i-1) * ((90 - dd) / (Np - 1));
            dP(i) = asind(sqrt(max(1 - sind(dB(i))^2 - sind(dT(i))^2, 0)));
        end
    case 't'
        for i = 1:Np
            dT(i) = dd;
            dP(i) = (i-1) * ((90 - dd) / (Np - 1));
            dB(i) = asind(sqrt(max(1 - sind(dP(i))^2 - sind(dT(i))^2, 0)));
        end
    case 'pmajor'
        for i = 1:Np
            dP(i) = 60;
            dT(i) = (i-1) * (30 / (Np - 1));
            dB(i) = asind(sqrt(max(1 - sind(dP(i))^2 - sind(dT(i))^2, 0)));
        end
    case 'bmajor'
        for i = 1:Np
            dB(i) = 60;
            dT(i) = (i-1) * (30 / (Np - 1));
            dP(i) = asind(sqrt(max(1 - sind(dB(i))^2 - sind(dT(i))^2, 0)));
        end
    case 'tmajor'
        for i = 1:Np
            dT(i) = 50;
            dP(i) = (i-1) * (40 / (Np - 1));
            dB(i) = asind(sqrt(max(1 - sind(dP(i))^2 - sind(dT(i))^2, 0)));
        end
end
[h, v] = localPTB2Triangle(dP, dT, dB);
end


function localPlotRepresentativeBeachballs(ax, vertex, coloF, diam)
% Representative focal mechanisms:
%   Normal      : [0, 45, -90]
%   Reverse     : [0, 45,  90]
%   Strike-slip : [45, 90, 180]
MEC_PlotBall2DFilled([0,45,-90], vertex.normal(1), vertex.normal(2), ...
    diam, coloF(3,:), ax);

MEC_PlotBall2DFilled([0,45,90], vertex.reverse(1), vertex.reverse(2), ...
    diam, coloF(4,:), ax);

MEC_PlotBall2DFilled([45,90,180], vertex.strikeSlip(1), vertex.strikeSlip(2), ...
    diam, coloF(2,:), ax);
end


function localPlotOddBeachballs(ax, coloF, diam)
SDRodd = [
    0,   90, -90
    144, 60,  35
    216, 60, -35
];
[~, dP, dT, dB] = localMechClass(SDRodd(:,1), SDRodd(:,2), SDRodd(:,3));
[h, v] = localPTB2Triangle(dP, dT, dB);
for i = 1:size(SDRodd, 1)
    MEC_PlotBall2DFilled(SDRodd(i,:), h(i), v(i), diam, coloF(1,:), ax);
end
end


function h = MEC_PlotBall2DFilled(fm, centerX, centerY, radius, color, ax)
% =========================================================================
% MEC_PlotBall2DFilled
% Plot a filled lower-hemisphere focal-mechanism beachball.
%
% This function is designed to be consistent with MEC_PlotBall2D:
%   - SDR convention follows MEC_StrikeDip2Norm / MEC_SDR2Slip
%   - NEU coordinates: [North, East, Up]
%   - lower hemisphere: U <= 0
%   - Wulff (stereographic/equal-angle) projection
%   - tension quadrants: radiation amplitude amp > 0
%
% Input:
%   fm       = [strike, dip, rake], degree
%   centerX  = beachball center x-coordinate in parent axes
%   centerY  = beachball center y-coordinate in parent axes
%   radius   = radius of beachball in parent axes coordinates
%   color    = RGB color used for tension quadrants
%   ax       = axes handle
%
% Output:
%   h        = structure containing graphics handles
%
% Important:
%   This function intentionally uses the same SDR-to-vector routines and
%   Wulff projection convention as MEC_PlotBall2D.
% =========================================================================

if nargin < 6 || isempty(ax)
    ax = gca;
end

if numel(fm) ~= 3
    error('MEC_PlotBall2DFilled:badFM', ...
        'fm must be a 1x3 vector: [strike, dip, rake].');
end

if numel(color) ~= 3
    error('MEC_PlotBall2DFilled:badColor', ...
        'color must be an RGB vector [R G B].');
end

strike = fm(1);
dip    = fm(2);
rake   = fm(3);

color = reshape(color, 1, 1, 3);

% -------------------------------------------------------------------------
% Use exactly the same SDR-to-vector convention as MEC_PlotBall2D.
%
% Returned coordinate system:
%   [N, E, U]
% -------------------------------------------------------------------------
n1 = MEC_StrikeDip2Norm(strike, dip);
d1 = MEC_SDR2Slip(strike, dip, rake);

n1 = n1(:);
d1 = d1(:);

n1 = n1 ./ norm(n1);
d1 = d1 ./ norm(d1);

% Double-couple moment tensor:
% M = n*d' + d*n'
M = n1 * d1.' + d1 * n1.';

% -------------------------------------------------------------------------
% Construct a regular grid on the Wulff projection disk.
%
% Wulff lower-hemisphere inverse projection:
%
%   E = 2*x / (1 + r^2)
%   N = 2*y / (1 + r^2)
%   U = (r^2 - 1) / (1 + r^2)
%
% where:
%   x = East coordinate on the projection plane
%   y = North coordinate on the projection plane
%
% At disk center:
%   [N,E,U] = [0,0,-1]
%
% At disk boundary:
%   U = 0
% -------------------------------------------------------------------------
Ngrid = 401;

e = linspace(-1, 1, Ngrid);
n = linspace(-1, 1, Ngrid);

[Eplot, Nplot] = meshgrid(e, n);

R2 = Eplot.^2 + Nplot.^2;
inside = R2 <= 1;

den = 1 + R2;

Vn = 2 .* Nplot ./ den;
Ve = 2 .* Eplot ./ den;
Vu = (R2 - 1) ./ den;

% -------------------------------------------------------------------------
% Radiation amplitude:
%
%   amp = v' * M * v
%
% This definition is identical to MEC_PlotBall2D:
%
%   amp = sum((sphere*M).*sphere, 2)
%
% amp > 0: tension quadrant
% amp < 0: compression quadrant
% -------------------------------------------------------------------------
amp = ...
      M(1,1) .* Vn.^2 ...
    + M(2,2) .* Ve.^2 ...
    + M(3,3) .* Vu.^2 ...
    + 2 .* M(1,2) .* Vn .* Ve ...
    + 2 .* M(1,3) .* Vn .* Vu ...
    + 2 .* M(2,3) .* Ve .* Vu;

tensionMask = inside & (amp > 0);

% Beachball coordinates in the parent triangle axes
X = centerX + radius .* Eplot;
Y = centerY + radius .* Nplot;

hold(ax, 'on');

% -------------------------------------------------------------------------
% Background white disk
% -------------------------------------------------------------------------
theta = linspace(0, 2*pi, 600);

% Keep the same orientation as MEC_PlotBall2D:
% x = sin(theta), y = cos(theta)
circleX = centerX + radius .* sin(theta);
circleY = centerY + radius .* cos(theta);

h.background = fill(ax, circleX, circleY, 'w', ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

% -------------------------------------------------------------------------
% Filled tension quadrants
%
% Use RGB true-color surface + alpha mask.
% This avoids modifying the axes colormap, so red/green/blue/gray
% beachballs can coexist in the same triangle figure.
% -------------------------------------------------------------------------
rgb = repmat(color, Ngrid, Ngrid);

h.tension = surface(ax, X, Y, zeros(size(X)), rgb, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'AlphaData', double(tensionMask), ...
    'AlphaDataMapping', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'HandleVisibility', 'off');

% -------------------------------------------------------------------------
% Outer circle
% -------------------------------------------------------------------------
h.circle = plot(ax, circleX, circleY, 'k-', ...
    'LineWidth', 0.45, ...
    'HandleVisibility', 'off');

% -------------------------------------------------------------------------
% Nodal planes
%
% IMPORTANT:
% These use the exact same Wulff-projection routine as MEC_PlotBall2D.
% The auxiliary plane has normal vector equal to the slip vector d1.
% -------------------------------------------------------------------------
[planeX1, planeY1] = MEC_CircleWulffNEU(n1.', 800);
[planeX2, planeY2] = MEC_CircleWulffNEU(d1.', 800);

h.plane1 = plot(ax, ...
    centerX + radius .* planeX1, ...
    centerY + radius .* planeY1, ...
    'k-', ...
    'LineWidth', 0.65, ...
    'HandleVisibility', 'off');

h.plane2 = plot(ax, ...
    centerX + radius .* planeX2, ...
    centerY + radius .* planeY2, ...
    'k-', ...
    'LineWidth', 0.65, ...
    'HandleVisibility', 'off');

end


function n = localSDR2NormalNEU(strike, dip)
% Fault-plane normal vector in NEU coordinates

n = [
    -sind(dip) .* sind(strike)
     sind(dip) .* cosd(strike)
    -cosd(dip)
].';
end

function s = localSDR2SlipNEU(strike, dip, rake)
% Slip vector in NEU coordinates
s = [
     cosd(rake) .* cosd(strike) + cosd(dip) .* sind(rake) .* sind(strike)
     cosd(rake) .* sind(strike) - cosd(dip) .* sind(rake) .* cosd(strike)
    -sind(rake) .* sind(dip)
].';
end


function [xPlot, yPlot] = localProjectNodalPlane(normalVec, centerX, centerY, diam)
% Great circle of a plane projected to lower hemisphere
normalVec = normalVec(:) ./ norm(normalVec);
tmp = [0; 0; 1];
if abs(dot(tmp, normalVec)) > 0.95
    tmp = [1; 0; 0];
end
e1 = cross(normalVec, tmp);
e1 = e1 ./ norm(e1);
e2 = cross(normalVec, e1);
e2 = e2 ./ norm(e2);

t = linspace(0, 2*pi, 800);
v = e1*cos(t) + e2*sin(t);

% Keep lower hemisphere in NEU, where U <= 0.
v(:, v(3,:) > 0) = NaN;
xPlot = centerX + diam * v(2,:);
yPlot = centerY + diam * v(1,:);
end
