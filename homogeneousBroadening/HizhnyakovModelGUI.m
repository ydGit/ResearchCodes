function varargout = HizhnyakovModelGUI(varargin)
% HIZHNYAKOVMODELGUI MATLAB code for HizhnyakovModelGUI.fig
%      HIZHNYAKOVMODELGUI, by itself, creates a new HIZHNYAKOVMODELGUI or raises the existing
%      singleton*.
%
%      H = HIZHNYAKOVMODELGUI returns the handle to a new HIZHNYAKOVMODELGUI or the handle to
%      the existing singleton*.
%
%      HIZHNYAKOVMODELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HIZHNYAKOVMODELGUI.M with the given input arguments.
%
%      HIZHNYAKOVMODELGUI('Property','Value',...) creates a new HIZHNYAKOVMODELGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HizhnyakovModelGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HizhnyakovModelGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HizhnyakovModelGUI

% Last Modified by GUIDE v2.5 06-Nov-2013 16:28:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HizhnyakovModelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HizhnyakovModelGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function [w, D, G, T, P] = updateAllParams(handles)
%---Updating ALL parameters from the control values---------
w = get( handles.WControl, 'Value' );
D = get( handles.DControl, 'Value' );
G = get( handles.GControl, 'Value' );
T = get( handles.TControl, 'Value' );
P = get( handles.PositionControl, 'Value' );
paramLine = strcat('w = ', num2str(w, handles.formatSpec),' cm-1   |  ',' D = ', num2str(D, handles.formatSpec),' cm-1   |   ', ' G = ', num2str(G, handles.formatSpec), ' cm-1   |   ', ' T = ', num2str(T, handles.formatSpec), '  K ');
set(handles.textParameters,'String', paramLine);
%===========================================================

function [freqAllT, IAllT, freqLowT, ILowT] = calculateAll(handles)
%--readind current values of all parameters----
w = get( handles.WControl, 'Value' );
D = get( handles.DControl, 'Value' );
G = get( handles.GControl, 'Value' );
T = get( handles.TControl, 'Value' );
P = get( handles.PositionControl, 'Value' );
%--setting default values for the -------------
freqAllT = zeros(1, 10);
IAllT = zeros(1, 10);
freqLowT = zeros(1, 10);
ILowT = zeros(1, 10);
%---checking what we need to calculate
if (get(handles.checkAllT, 'Value') == 1)
    [freqAllT, IAllT] = HizhAll(w, D, G, T, P);    
end
if (get(handles.checkLowT, 'Value') == 1)
    [freqLowT, ILowT] = HizhLow(w, D, G, T, P);    
end
% Estimating critical temperature, where the Low-T model is supposedly
% fails
Tcrit = 2*1.44*w/( log(G/abs(D)+abs(D)/G) );
set( handles.textTcritical, 'String', strcat('Tc = ', num2str(Tcrit, handles.formatSpec), ' K') );


function h = plotAll(freqAll, IAll, freqLow, ILow, freqExp, IExp, handles)
%calculating current plotting limits
z = find(IAll > max(IAll)/2.0);
fw = freqAll(z(end))- freqAll(z(1));
xcenter = 0.5*(freqAll(z(end))+freqAll(z(1)));
span = 20;
xlower = xcenter - span*fw;
xupper = xcenter + span*fw;
xlimits = [xlower xupper];
ylimits = [0 1];

Gamma = fwhm(freqAll, IAll);
set(handles.textFWHMAllT, 'String', strcat('FWHM = ', num2str(Gamma)));

cla;
% plot All-T Model, Low-T, and Experiment
if (get(handles.checkAllT, 'Value') == 1) 
    plot(freqAll, IAll, 'red');    
    hold on;
end

if (get(handles.checkLowT, 'Value') == 1) 
    plot(freqLow, ILow, 'blue');
    hold on;
end

% Show experimental broadening: inhomogeneous and instrumental effects
%plot(handles.broadX, handles.broadY, 'Color', [.9, .9, .2]);
set(gca,'Color',[0.6,0.6,0.6])

if (get(handles.checkExperiment, 'Value') == 1)
    expT = get(handles.menuT, 'Value');
    switch expT
        case 1.0
            freqExp = handles.fExp25;
            IExp = handles.IExp25;
        case 2.0
            freqExp = handles.fExp30;
            IExp = handles.IExp30;
        case 3.0
            freqExp = handles.fExp35;
            IExp = handles.IExp35;
        case 4.0
            freqExp = handles.fExp40;
            IExp = handles.IExp40;
        case 5.0
            freqExp = handles.fExp45;
            IExp = handles.IExp45;
        case 6.0
            freqExp = handles.fExp60;
            IExp = handles.IExp60;
        case 7.0
            freqExp = handles.fExp70;
            IExp = handles.IExp70;
        case 8.0
            freqExp = handles.fExp90;
            IExp = handles.IExp90;
        case 9.0
            freqExp = handles.fExp120;
            IExp = handles.IExp120;
        case 10.0
            freqExp = handles.fExp150;
            IExp = handles.IExp150;
        case 11.0
            freqExp = handles.fExp180;
            IExp = handles.IExp180;
        case 12.0
            freqExp = handles.fExp240;
            IExp = handles.IExp240;
        case 13.0
            freqExp = handles.fExp300;
            IExp = handles.IExp300;
        otherwise
            freqExp = handles.fExp25;
            IExp = handles.IExp25;    
    end
    plot(freqExp, IExp, 'black');
    hold on;
end
hold off;
xlim(xlimits);
ylim(ylimits);
    

function [f, I] = HizhAll(w, D, G, T, P)
w = cast(w, 'double'); % QLM mode frequency, cm-1
D = cast(D, 'double'); % QLM frequency change, cm-1
G = 0.5*cast(G, 'double'); % QLM mode decay width, cm-1, !! IT IS HALF OF THE LOW T-Model Gamma !!
T = cast(T, 'double'); % Model temperature, K
P = cast(P, 'double'); % Emission line position, cm-1
%----------------------------------------------------------------

nbar = 1/(exp(1.44*w/T)-1); % QLM population; 1.44 is for cm-1 -> K conversion
b = sqrt( G - D^2/4 - 1i*G*D*(2*nbar+1) );
if ( real(b) > 0 )
    beta = b;
else
    beta = -b;
end

lambda = G - 1i*D/2 + beta;
alpha = -1i*nbar*D/lambda;

%==========Natural linewidth: gamma_0 FWHM ===============
g0 = 0.003; % FWHM, cm-1 or 0.90 GHz;

% generates Spectrum according to the All-Temperatures model
Fs = 12000; %GHz sampling frequency
tmax = 200; %ns integration time
numPoints = 1 + tmax*Fs; % number of data points in time domain
t = linspace(0, tmax, numPoints);
k = 30; % cm-1 to GHz=1/ns conversion
%--Spectrum using new, all-T model:
P0 = (Fs/2-k*D/2); % compensating for the initial shift at T = 0 K
y = exp( +1i*2*pi*t.*P0-pi*k*t.*g0 + 2*pi*1i*t.*P*k+2*pi*t.*1i*k*D*(nbar+1.0/2)-2*pi*t.*1i*k*D*(nbar+1)*alpha/(1+alpha)+(1i*D*(nbar+1)/(lambda*(1+alpha))*log(1+alpha-alpha*exp(-2*pi*t.*k*lambda))));

NFFT = 2^nextpow2(numPoints);
Y = fft(y,NFFT);
Y = Y/max(Y);

%f = 1/k*Fs/2*linspace(0,1,NFFT/2 + 1); % frequency domain spectrum
f_full = 1/k*Fs/2*linspace(-1,1,NFFT); % frequency domain spectrum
%I = Y(1:NFFT/2 + 1).*conj(Y(1:NFFT/2 + 1));
I_full = Y(1:NFFT).*conj(Y(1:NFFT));

% reducing the tremendous number of points in the Fourier Transform
% taking every 'nth = stride' point
stride = 100;
i = 1:stride:length(f_full);
f = zeros(1, length(i));
I = zeros(1, length(i));
k = 1;
for j = i
    f(k) = f_full(j);
    I(k) = I_full(j);
    k = k + 1;
end
    


function [f, I] = HizhLow(w, D, G, T, P)
w = cast(w, 'double'); % QLM mode frequency, cm-1
D = cast(D, 'double'); % QLM frequency change, cm-1
G = cast(G, 'double'); % QLM mode decay width, cm-1
T = cast(T, 'double'); % Model temperature, K
P = cast(P, 'double'); % Emission line position, cm-1
%----------------------------------------------------------------

nbar = 1/(exp(1.44*w/T)-1); % QLM population; 1.44 is for cm-1 -> K conversion
gamma = 2*D^2*G*nbar*(nbar + 1)/( D^2*(nbar+1)^2 + G^2 );
delta = D/2+2*D^2*G*nbar/( D^2*(nbar+1)^2 + G^2 );
zeta = 1/( 1 + 1i*G/(D*(nbar + 1)) );

%==========Natural linewidth: gamma_0 FWHM ===============
g0 = 0.003; % FWHM, cm-1 or 0.90 GHz;

% generates Spectrum according to the All-Temperatures model
Fs = 6000; %GHz sampling frequency
tmax = 200; %ns integration time
numPoints = 1 + tmax*Fs; % number of data points in time domain
t = linspace(0, tmax, numPoints);
k = 30; % cm-1 to GHz=1/ns conversion
%--Spectrum using new, all-T model:
%y = exp( -pi*k*t.*g0+2*pi*1i*t.*P*k);
P0 = (Fs/2-k*D/2); % compensating for the initial shift at T = 0 K
y = exp( +1i*2*pi*t.*P0-pi*k*t.*g0+2*pi*1i*t.*P*k+1i*delta*k*t*2*pi - 2*pi*t.*k*gamma/2 - zeta*log(1-zeta*exp(-1.44*w/T+2*pi*1i*t.*k*D - G*k*t.*2*pi)) );

NFFT = 2^nextpow2(numPoints);
Y = fft(y,NFFT);
Y = Y/max(Y);
%f = 1/k*Fs/2*linspace(0,1,NFFT/2 + 1); % frequency domain spectrum
f_full = 1/k*Fs/2*linspace(-1,1,NFFT); % frequency domain spectrum
%I = Y(1:NFFT/2 + 1).*conj(Y(1:NFFT/2 + 1));
I_full = Y(1:NFFT).*conj(Y(1:NFFT));

% reducing the tremendous number of points in the Fourier Transform
% taking every 'nth = stride' point
stride = 100;
i = 1:stride:length(f_full);
f = zeros(1, length(i));
I = zeros(1, length(i));
k = 1;
for j = i
    f(k) = f_full(j);
    I(k) = I_full(j);
    k = k + 1;
end


function [x, y] = readFile(fileName)
fileID = fopen(fileName);
A = fscanf(fileID, '%f\t%f\n', [2 inf]);
fclose(fileID);
x = A(1, 1:end);
y = A(2, 1:end);

function [x ,y] = lorentz(amplitude, center, width, xaxis)
x = xaxis;
y = amplitude./( 1+(2*(x-center)/width).^2 );

function [x ,y] = gauz(amplitude, center, width, xaxis)
x = xaxis;
y = amplitude*exp( -(2.40*(x-center)/width).^2 );

function g = fwhm( x, y )
z = find(y > max(y)/2.0);
g = x(z(end))- x(z(1));

function [TAll, gAll, TLow, gLow] = generateFWHM(handles)
%overSample = 1;
TAll = handles.e11T;%linspace(handles.e11T(1), handles.e11T(end), overSample*length(handles.e11T) );
gAll = zeros(1, length(TAll));
TLow = handles.e11T;
gLow = zeros(1, length(TLow));
%------------Updating ALL parameters------------------------
w = get( handles.WControl, 'Value' );
D = get( handles.DControl, 'Value' );
G = get( handles.GControl, 'Value' );
P = get( handles.PositionControl, 'Value' );
j = 1;
if (get(handles.checkAllT, 'Value') == 1)
    for i = TAll    
        [freq, I] = HizhAll(w, D, G, i, P);
        gAll(j) = fwhm(freq, I);
        j = j + 1;
    end
end
j = 1;    
if (get(handles.checkLowT, 'Value') == 1)
    for i = TLow    
        [freq, I] = HizhLow(w, D, G, i, P);
        gLow(j) = fwhm(freq, I);
        j = j + 1;
    end
end


% --- Executes just before HizhnyakovModelGUI is made visible.
function HizhnyakovModelGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HizhnyakovModelGUI (see VARARGIN)

% Reading experimental spectra from file
%[handles.x, handles.y] = readFile('C:\Documents and Settings\yury.dziashko\My Documents\MATLAB\Xe1e11Hizh\forAnalysis\Xe1e11_40Knorm.txt');

handles.e11T = [25.0 30.0 35.0 40.0 45.0 60.0 70.0 80.0 90.0 100.0 120.0 150.0 180.0 240.0 300.0];
handles.e11FWHM = [0.04 0.16 0.21 0.44 0.79 2.23 3.76 5.05 7.78 10.9 17.2 27.2 45.2 77.9 125.8];

handles.formatSpec = '%10.1f';

%------------Choose What to Display--------------------------------
handles.showAllT = 1.0; % All-T model should be calculated and displayed
handles.showLowT = 0.0; % Low-T model should NOT be calculated OR displayed
handles.showExperiment = 1.0; % Show experimental spectrum
set(handles.checkAllT, 'Value', handles.showAllT);
set(handles.checkLowT, 'Value', handles.showLowT);
set(handles.checkExperiment, 'Value', handles.showExperiment);
%==================================================================

%-----------QLM Frequency Related Controls and  Properties-----------
handles.Positionmin = -20.0;
handles.Positionmax = 20.0;
handles.Position = 0.0; % Quazi-Local Frequency
set(handles.PositionControl,'Min', handles.Positionmin);
set(handles.PositionControl,'Max', handles.Positionmax);
set(handles.PositionControl,'Value', handles.Position);
set(handles.editPositionmin,'String', num2str(handles.Positionmin, handles.formatSpec));
set(handles.editPositioncurrent,'String', num2str(handles.Position, handles.formatSpec));
set(handles.editPositionmax,'String', num2str(handles.Positionmax, handles.formatSpec));
%======================================================================


%-----------QLM Frequency Related Controls and  Properties-----------
handles.wmin = 1.0;
handles.wmax = 600.0;
handles.w = 115.0; % Quazi-Local Frequency
set(handles.WControl,'Min', handles.wmin);
set(handles.WControl,'Max', handles.wmax);
set(handles.WControl,'Value', handles.w);
set(handles.editWmin,'String', num2str(handles.wmin, handles.formatSpec));
set(handles.editWcurrent,'String', num2str(handles.w, handles.formatSpec));
set(handles.editWmax,'String', num2str(handles.wmax, handles.formatSpec));
%======================================================================

%---------QLM Frequency Change Related Controls and  Properties---------
handles.Dmin = -100.0;
handles.Dmax = 100.0;
handles.D = -40.0; % Change of the Frequency
set(handles.DControl,'Min', handles.Dmin);
set(handles.DControl,'Max', handles.Dmax);
set(handles.DControl,'Value', handles.D);
set(handles.editDmin,'String', num2str(handles.Dmin, handles.formatSpec));
set(handles.editDcurrent,'String', num2str(handles.D, handles.formatSpec));
set(handles.editDmax,'String', num2str(handles.Dmax, handles.formatSpec));
%======================================================================

%---------QLM Frequency Change Related Controls and  Properties---------
handles.Gmin = 1.0;
handles.Gmax = 100.0;
handles.G = 33.0; % Quazi-Local Mode Decay Rate
set(handles.GControl,'Min', handles.Gmin);
set(handles.GControl,'Max', handles.Gmax);
set(handles.GControl,'Value', handles.G);
set(handles.editGmin,'String', num2str(handles.Gmin, handles.formatSpec));
set(handles.editGcurrent,'String', num2str(handles.G, handles.formatSpec));
set(handles.editGmax,'String', num2str(handles.Gmax, handles.formatSpec));
%======================================================================

%-----------Temperature Related Controls and  Properties-----------
handles.Tmin = 0.1; %minumum model temperature, K
handles.Tmax = 300.0; %maximum model temperature, K
handles.T = 120.0; % Temperature of the Model
handles.Tcrit = 2*1.4388*handles.w/log(abs(handles.G/handles.D)+abs(handles.D/handles.G));
set(handles.TControl,'Min', handles.Tmin);
set(handles.TControl,'Max', handles.Tmax);
set(handles.TControl,'Value', handles.T);
set(handles.editTmin,'String', num2str(handles.Tmin, handles.formatSpec));
set(handles.editTcurrent,'String', num2str(handles.T, handles.formatSpec));
set(handles.editTmax,'String', num2str(handles.Tmax, handles.formatSpec));
%====================================================================

%---------------------Setting Parameters Summary Line--------------------
[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
%========================================================================

%-------------------- Loading Experimental Data----------------------
[handles.fExp25, handles.IExp25] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_25Knorm.txt');
[handles.fExp30, handles.IExp30] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_30Knorm.txt');
[handles.fExp35, handles.IExp35] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_35Knorm.txt');
[handles.fExp40, handles.IExp40] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_40Knorm.txt');
[handles.fExp45, handles.IExp45] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_45Knorm.txt');
[handles.fExp60, handles.IExp60] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_60Knorm.txt');
[handles.fExp70, handles.IExp70] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_70Knorm.txt');
[handles.fExp90, handles.IExp90] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_90Knorm.txt');
[handles.fExp120, handles.IExp120] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_120Knorm.txt');
[handles.fExp150, handles.IExp150] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_150Knorm.txt');
[handles.fExp180, handles.IExp180] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_180Knorm.txt');
[handles.fExp240, handles.IExp240] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_240Knorm.txt');
[handles.fExp300, handles.IExp300] = readFile('Xe1e11Hizh\forAnalysis\Xe1e11_300Knorm.txt');
% Flipping the experimental spectra
handles.IExp25 = handles.IExp25(end:-1:1);
handles.IExp30 = handles.IExp30(end:-1:1);
handles.IExp35 = handles.IExp35(end:-1:1);
handles.IExp40 = handles.IExp40(end:-1:1);
handles.IExp45 = handles.IExp45(end:-1:1);
handles.IExp60 = handles.IExp60(end:-1:1);
handles.IExp70 = handles.IExp70(end:-1:1);
handles.IExp90 = handles.IExp90(end:-1:1);
handles.IExp120 = handles.IExp120(end:-1:1);
handles.IExp150 = handles.IExp150(end:-1:1);
handles.IExp180 = handles.IExp180(end:-1:1);
handles.IExp240 = handles.IExp240(end:-1:1);
handles.IExp300 = handles.IExp300(end:-1:1);

set(handles.menuT, 'Value', 1.0);
handles.freqExp = handles.fExp120;
handles.IExp = handles.IExp120;

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);

% ------ Experimental effects --------------
set(handles.editLorentz, 'String', num2str(0.35, handles.formatSpec));
%[handles.inhomX, handles.inhomY] = lorentz(1.0, 0, 0.35, handles.freq);
[handles.inhomX, handles.inhomY] = lorentz(1.0, 0, 7.0, handles.freq);
[handles.instrumX, handles.instrumY] = gauz(1.0, 0, 0.40, handles.freq);
handles.broadX = handles.freq;
handles.broadY = conv(handles.inhomY, handles.instrumY, 'same');
handles.broadY = handles.broadY/max(handles.broadY);
%-------------------------------------------------------------------

plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);


%handles.x = conv(handles.I, handles.inhomY, 'same');
%handles.y = handles.freq;

%[handles.x, handles.y] = lorentz(1.0, 0, 0.2, handles.freq);
%[handles.x, handles.y] = readFile('C:\Documents and Settings\yury.dziashko\My Documents\MATLAB\Xe1e11Hizh\forAnalysis\Xe1e11_40Knorm.txt');

%plot(handles.freq, handles.I, handles.x, handles.y);
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));

% Choose default command line output for HizhnyakovModelGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HizhnyakovModelGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HizhnyakovModelGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function TControl_Callback(hObject, eventdata, handles)
% hObject    handle to TControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

set( handles.editTcurrent, 'String', num2str(handles.T, handles.formatSpec) );

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)

%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));


% --- Executes during object creation, after setting all properties.
function TControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns
% called

function editTmin_Callback(hObject, eventdata, handles)
% hObject    handle to editTmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% setting the minimum temperature limit
handles.Tmin = str2num(get(hObject, 'String'));
set(handles.TControl, 'Min', handles.Tmin);

% if the current temperature is lower than minimum allowed -- change it!
T = get(handles.TControl, 'Value');
if ( handles.Tmin > T )
    set(handles.TControl, 'Value', handles.Tmin);
    set(handles.editTcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editTmin as text
%        str2double(get(hObject,'String')) returns contents of editTmin as a double


% --- Executes during object creation, after setting all properties.
function editTmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function editTcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to editTcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.T = str2num(get(hObject, 'String'));
set(handles.TControl, 'Value', handles.T);

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));

% Hints: get(hObject,'String') returns contents of editTcurrent as text
%        str2double(get(hObject,'String')) returns contents of editTcurrent as a double


% --- Executes during object creation, after setting all properties.
function editTcurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editTmax_Callback(hObject, eventdata, handles)
% hObject    handle to editTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Tmax = str2num(get(hObject, 'String'));
set(handles.TControl, 'Max', handles.Tmax);

% if the current temperature is lower than minimum allowed -- change it!
T = get(handles.TControl, 'Value');
if ( handles.Tmax < T )
    set(handles.TControl, 'Value', handles.Tmax);
    set(handles.editTcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editTmax as text
%        str2double(get(hObject,'String')) returns contents of editTmax as a double


% --- Executes during object creation, after setting all properties.
function editTmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on slider movement.
function GControl_Callback(hObject, eventdata, handles)
% hObject    handle to GControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

set(handles.editGcurrent, 'String', num2str(handles.G, handles.formatSpec));

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function GControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editGmin_Callback(hObject, eventdata, handles)
% hObject    handle to editGmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Gmin = str2num(get(hObject, 'String'));
set(handles.GControl, 'Min', handles.Gmin);

% if the current temperature is lower than minimum allowed -- change it!
g = get(handles.GControl, 'Value');
if ( handles.Gmin > g )
    set(handles.GControl, 'Value', handles.Gmin);
    set(handles.editGcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editGmin as text
%        str2double(get(hObject,'String')) returns contents of editGmin as a double


% --- Executes during object creation, after setting all properties.
function editGmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function editGcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to editGcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.G = str2num(get(hObject, 'String'));
set(handles.GControl, 'Value', handles.G);

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'String') returns contents of editGcurrent as text
%        str2double(get(hObject,'String')) returns contents of editGcurrent as a double


% --- Executes during object creation, after setting all properties.
function editGcurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editGmax_Callback(hObject, eventdata, handles)
% hObject    handle to editGmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Gmax = str2num(get(hObject, 'String'));
set(handles.GControl, 'Max', handles.Gmax);

% if the current temperature is lower than minimum allowed -- change it!
g = get(handles.GControl, 'Value');
if ( handles.Gmax < g )
    set(handles.GControl, 'Value', handles.Gmax);
    set(handles.editGcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editGmax as text
%        str2double(get(hObject,'String')) returns contents of editGmax as a double


% --- Executes during object creation, after setting all properties.
function editGmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on slider movement.
function DControl_Callback(hObject, eventdata, handles)
% hObject    handle to DControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%------------Updating ALL parameters------------------------
[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
%===========================================================

set(handles.editDcurrent, 'String', num2str(handles.D, handles.formatSpec));
[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function DControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function editDmin_Callback(hObject, eventdata, handles)
% hObject    handle to editDmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Dmin = str2num(get(hObject, 'String'));
set(handles.DControl, 'Min', handles.Dmin);

% if the current Delta is lower than minimum allowed -- change it!
d = get(handles.DControl, 'Value');
if ( handles.Dmin > d )
    set(handles.DControl, 'Value', handles.Dmin);
    set(handles.editDcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editDmin as text
%        str2double(get(hObject,'String')) returns contents of editDmin as a double


% --- Executes during object creation, after setting all properties.
function editDmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editDcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to editDcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.D = str2num(get(hObject, 'String'));
set(handles.DControl, 'Value', handles.D);

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'String') returns contents of editDcurrent as text
%        str2double(get(hObject,'String')) returns contents of editDcurrent as a double


% --- Executes during object creation, after setting all properties.
function editDcurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editDmax_Callback(hObject, eventdata, handles)
% hObject    handle to editDmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Dmax = str2num(get(hObject, 'String'));
set(handles.DControl, 'Max', handles.Dmax);

% if the current temperature is lower than minimum allowed -- change it!
d = get(handles.DControl, 'Value');
if ( handles.Dmax < d )
    set(handles.DControl, 'Value', handles.Dmax);
    set(handles.editDcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editDmax as text
%        str2double(get(hObject,'String')) returns contents of editDmax as a double


% --- Executes during object creation, after setting all properties.
function editDmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function WControl_Callback(hObject, eventdata, handles)
% hObject    handle to WControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%------------Updating ALL parameters------------------------
[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
%===========================================================

set( handles.editWcurrent, 'String', num2str(handles.w, handles.formatSpec) );

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function WControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function editWmin_Callback(hObject, eventdata, handles)
% hObject    handle to editWmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Wmin = str2num(get(hObject, 'String'));
set(handles.WControl, 'Min', handles.Wmin);

% if the current temperature is lower than minimum allowed -- change it!
w = get(handles.WControl, 'Value');
if ( handles.Wmin > w )
    set(handles.WControl, 'Value', handles.Wmin);
    set(handles.editWcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editWmin as text
%        str2double(get(hObject,'String')) returns contents of editWmin as a double


% --- Executes during object creation, after setting all properties.
function editWmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function editWcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to editWcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.w = str2num(get(hObject, 'String'));
set(handles.WControl, 'Value', handles.w);

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

paramLine = strcat('Parameters:  w = ', num2str(handles.w, handles.formatSpec),'; D = ', num2str(handles.D, handles.formatSpec), '; G = ', num2str(handles.G, handles.formatSpec), '; T = ', num2str(handles.T, handles.formatSpec));
set(handles.textParameters,'String', paramLine);
[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'String') returns contents of editWcurrent as text
%        str2double(get(hObject,'String')) returns contents of editWcurrent as a double


% --- Executes during object creation, after setting all properties.
function editWcurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editWmax_Callback(hObject, eventdata, handles)
% hObject    handle to editWmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Wmax = str2num(get(hObject, 'String'));
set(handles.WControl, 'Max', handles.Wmax);

% if the current temperature is lower than minimum allowed -- change it!
w = get(handles.WControl, 'Value');
if ( handles.Wmax < w )
    set(handles.WControl, 'Value', handles.Wmax);
    set(handles.editWcurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editWmax as text
%        str2double(get(hObject,'String')) returns contents of editWmax as a double


% --- Executes during object creation, after setting all properties.
function editWmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function textParameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function PositionControl_Callback(hObject, eventdata, handles)
% hObject    handle to PositionControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%------------Updating ALL parameters------------------------
[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
%===========================================================

set( handles.editPositioncurrent, 'String', num2str(handles.Position, handles.formatSpec) );
[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
%set(handles.textFWHMAllT, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function PositionControl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PositionControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function editPositionmin_Callback(hObject, eventdata, handles)
% hObject    handle to editPositionmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Positionmin = str2num(get(hObject, 'String'));
set(handles.PositionControl, 'Min', handles.Positionmin);

% if the current temperature is lower than minimum allowed -- change it!
p = get(handles.PositionControl, 'Value');
if ( handles.Positionmin > p )
    set(handles.PositionControl, 'Value', handles.Positionmin);
    set(handles.editPositioncurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editPositionmin as text
%        str2double(get(hObject,'String')) returns contents of editPositionmin as a double


% --- Executes during object creation, after setting all properties.
function editPositionmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPositionmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function editPositioncurrent_Callback(hObject, eventdata, handles)
% hObject    handle to editPositioncurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Position = str2num(get(hObject, 'String'));
set(handles.PositionControl, 'Value', handles.Position);

[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);

[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)

% Hints: get(hObject,'String') returns contents of editPositioncurrent as text
%        str2double(get(hObject,'String')) returns contents of editPositioncurrent as a double


% --- Executes during object creation, after setting all properties.
function editPositioncurrent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPositioncurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editPositionmax_Callback(hObject, eventdata, handles)
% hObject    handle to editPositionmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Positionmax = str2num(get(hObject, 'String'));
set(handles.PositionControl, 'Max', handles.Positionmax);

% if the current temperature is lower than minimum allowed -- change it!
p = get(handles.PositionControl, 'Value');
if ( handles.Positionmax < p )
    set(handles.PositionControl, 'Value', handles.Positionmax);
    set(handles.editPositioncurrent, 'String', get(hObject, 'String'));
    [handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
    [handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);
    plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles);
end
% Hints: get(hObject,'String') returns contents of editPositionmax as text
%        str2double(get(hObject,'String')) returns contents of editPositionmax as a double


% --- Executes during object creation, after setting all properties.
function editPositionmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPositionmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btnGenerateFWHM.
function btnGenerateFWHM_Callback(hObject, eventdata, handles)
% hObject    handle to btnGenerateFWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[TAll, gAll, TLow, gLow] = generateFWHM(handles);
figure();
loglog(handles.e11T, handles.e11FWHM, 'o');
hold on;
if (get(handles.checkAllT, 'Value') == 1) 
    loglog(TAll, gAll, 'red');
    hold on;
end

if (get(handles.checkLowT, 'Value') == 1) 
    loglog(TLow, gLow, 'blue');
    hold on;
end
hold off;



% --- Executes during object creation, after setting all properties.
function btnGenerateFWHM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btnGenerateFWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkAllT.
function checkAllT_Callback(hObject, eventdata, handles)
% hObject    handle to checkAllT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showAllT = 1;
% Hint: get(hObject,'Value') returns toggle state of checkAllT


% --- Executes on button press in checkLowT.
function checkLowT_Callback(hObject, eventdata, handles)
% hObject    handle to checkLowT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkLowT


% --- Executes on button press in checkExperiment.
function checkExperiment_Callback(hObject, eventdata, handles)
% hObject    handle to checkExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkExperiment


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over checkAllT.
function checkAllT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to checkAllT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in menuT.
function menuT_Callback(hObject, eventdata, handles)
% hObject    handle to menuT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuT


% --- Executes during object creation, after setting all properties.
function menuT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes on button press in btnSave.

% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[w, D, G, T, P] = updateAllParams(handles);
[freqAllT, IAllT, freqLowT, ILowT] = calculateAll(handles);
% plot All-T Model, Low-T, and Experiment

%calculating current plotting limits
z = find(IAllT > max(IAllT)/2.0);
fw = freqAllT(z(end))- freqAllT(z(1));
xcenter = 0.5*(freqAllT(z(end))+freqAllT(z(1)));
span = 5;
xlower = xcenter(1) - span*fw;
xupper = xcenter(1) + span*fw;
xlimits = [xlower xupper];
ylimits = [0 1];

h = figure;
set(h, 'Visible', 'off');
if (get(handles.checkAllT, 'Value') == 1) 
    plot(freqAllT, IAllT, 'red');    
    hold on;
end

if (get(handles.checkLowT, 'Value') == 1) 
    plot(freqLowT, ILowT, 'blue');
    hold on;
end

if (get(handles.checkExperiment, 'Value') == 1)
    if (get(handles.menuT, 'Value') == 1.0)
        freqExp = handles.fExp25;
        IExp = handles.IExp25;
    elseif (get(handles.menuT, 'Value') == 2.0)
        freqExp = handles.fExp30;
        IExp = handles.IExp30;
    elseif (get(handles.menuT, 'Value') == 3.0)
        freqExp = handles.fExp35;
        IExp = handles.IExp35;
    elseif (get(handles.menuT, 'Value') == 4.0)
        freqExp = handles.fExp40;
        IExp = handles.IExp40;
    elseif (get(handles.menuT, 'Value') == 5.0)
        freqExp = handles.fExp45;
        IExp = handles.IExp45;
    elseif (get(handles.menuT, 'Value') == 6.0)
        freqExp = handles.fExp60;
        IExp = handles.IExp60;
    elseif (get(handles.menuT, 'Value') == 7.0)
        freqExp = handles.fExp70;
        IExp = handles.IExp70;
    elseif (get(handles.menuT, 'Value') == 8.0)
        freqExp = handles.fExp90;
        IExp = handles.IExp90;
    elseif (get(handles.menuT, 'Value') == 9.0)
        freqExp = handles.fExp120;
        IExp = handles.IExp120;
    end
    plot(freqExp, IExp, 'black');
    hold on;
end
hold off;
xlim(xlimits);
ylim(ylimits);

fileName = strcat('AllTRed_LowTBlue_w_', num2str(w, handles.formatSpec),'_D_', num2str(D,handles.formatSpec), '_G_', num2str(G, handles.formatSpec), '_T_', num2str(T,handles.formatSpec), '.png');
labelName = strcat('AllT is Red; LowT is Blue; w =', num2str(w, handles.formatSpec),' cm-1; D = ', num2str(D,handles.formatSpec), ' cm-1; G =', num2str(G, handles.formatSpec), ' cm-1; T=', num2str(T,handles.formatSpec));
text(xlower, 1.02, labelName);
print(h, '-dpng', fileName);


% --- Executes on button press in btnConvolve.
function btnConvolve_Callback(hObject, eventdata, handles)
% hObject    handle to btnConvolve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);

% ------ Experimental effects --------------
gL = str2num(get(handles.editLorentz, 'String')); % Lorentzian broadening
gG = 0.40; % Gaussian broadening
[handles.inhomX, handles.inhomY] = lorentz(1.0, 0, gL, handles.freq);
[handles.instrumX, handles.instrumY] = gauz(1.0, 0, gG, handles.freq);
handles.broadX = handles.freq;
handles.broadY = conv(handles.inhomY, handles.instrumY, 'same');
handles.broadY = handles.broadY/max(handles.broadY);
%-------------------------------------------------------------------

handles.I = conv(handles.broadY, handles.I, 'same');
handles.I = handles.I/max(handles.I);

handles.IL = conv(handles.broadY, handles.IL, 'same'); 
handles.IL = handles.IL/max(handles.IL);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)



function editFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to editFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrequency as text
%        str2double(get(hObject,'String')) returns contents of editFrequency as a double


% --- Executes during object creation, after setting all properties.
function editFrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editIntTime_Callback(hObject, eventdata, handles)
% hObject    handle to editIntTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIntTime as text
%        str2double(get(hObject,'String')) returns contents of editIntTime as a double


% --- Executes during object creation, after setting all properties.
function editIntTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIntTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLorentz_Callback(hObject, eventdata, handles)
% hObject    handle to editLorentz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLorentz as text
%        str2double(get(hObject,'String')) returns contents of editLorentz as a double


% --- Executes during object creation, after setting all properties.
function editLorentz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLorentz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnSConvolve.
function btnSConvolve_Callback(hObject, eventdata, handles)
% hObject    handle to btnSConvolve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.w, handles.D, handles.G, handles.T, handles.Position] = updateAllParams(handles);
[handles.freq, handles.I, handles.freqL, handles.IL] = calculateAll(handles);

% ------ Experimental effects --------------
gL = str2num(get(handles.editLorentz, 'String')); % Lorentzian broadening
gG = 0.40; % Gaussian broadening
[handles.inhomX, handles.inhomY] = lorentz(1.0, 0, gL, handles.freq);
[handles.instrumX, handles.instrumY] = gauz(1.0, 0, gG, handles.freq);
handles.broadX = handles.freq;
handles.broadY = conv(handles.inhomY, handles.instrumY, 'same');
handles.broadY = handles.broadY/max(handles.broadY);
%-------------------------------------------------------------------

handles.I = conv(handles.I, handles.I, 'same');
handles.I = conv(handles.I, handles.I, 'same');
handles.I = conv(handles.I, handles.broadY, 'same');
handles.I = handles.I/max(handles.I);

handles.IL = conv(handles.IL, handles.IL, 'same');
handles.IL = conv(handles.IL, handles.IL, 'same');
handles.IL = conv(handles.IL, handles.broadY, 'same'); 
handles.IL = handles.IL/max(handles.IL);
plotAll(handles.freq, handles.I, handles.freqL, handles.IL, handles.freqExp ,handles.IExp, handles)
