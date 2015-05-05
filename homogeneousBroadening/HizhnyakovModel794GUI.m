function varargout = HizhnyakovModel794GUI(varargin)
% HIZHNYAKOVMODEL794GUI MATLAB code for HizhnyakovModel794GUI.fig
%      HIZHNYAKOVMODEL794GUI, by itself, creates a new HIZHNYAKOVMODEL794GUI or raises the existing
%      singleton*.
%
%      H = HIZHNYAKOVMODEL794GUI returns the handle to a new HIZHNYAKOVMODEL794GUI or the handle to
%      the existing singleton*.
%
%      HIZHNYAKOVMODEL794GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HIZHNYAKOVMODEL794GUI.M with the given input arguments.
%
%      HIZHNYAKOVMODEL794GUI('Property','Value',...) creates a new HIZHNYAKOVMODEL794GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HizhnyakovModel794GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HizhnyakovModel794GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HizhnyakovModel794GUI

% Last Modified by GUIDE v2.5 10-Dec-2013 13:31:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HizhnyakovModel794GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HizhnyakovModel794GUI_OutputFcn, ...
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

function [w1, D1, G1, w2, D2, G2, w3, D3, G3, w4, D4, G4, T, P] = updateAllParams(handles)
%  Updating ALL parameters
w1str = get(handles.editW1, 'String');
D1str = get(handles.editD1, 'String');
G1str = get(handles.editG1, 'String');
w1 = str2num( w1str );
D1 = str2num( D1str );
G1 = str2num( G1str );
%
if ( get(handles.cbAll, 'Value') == 1.0 )
    set( handles.editW2, 'String', w1str);
    set( handles.editW3, 'String', w1str);
    set( handles.editW4, 'String', w1str);
    set( handles.editD2, 'String', D1str);
    set( handles.editD3, 'String', D1str);
    set( handles.editD4, 'String', D1str);
    set( handles.editG2, 'String', G1str);
    set( handles.editG3, 'String', G1str);
    set( handles.editG4, 'String', G1str);
end
    
w2 = str2num( get(handles.editW2, 'String') );
D2 = str2num( get(handles.editD2, 'String') );
G2 = str2num( get(handles.editG2, 'String') );
%
w3 = str2num( get(handles.editW3, 'String') );
D3 = str2num( get(handles.editD3, 'String') );
G3 = str2num( get(handles.editG3, 'String') );
%
w4 = str2num( get(handles.editW4, 'String') );
D4 = str2num( get(handles.editD4, 'String') );
G4 = str2num( get(handles.editG4, 'String') );

T = get( handles.ControlT, 'Value' );
P = get( handles.PositionControl, 'Value' );
%===========================================================

function [freq, I] = calculateAll(handles)

[w1, D1, G1, w2, D2, G2, w3, D3, G3, w4, D4, G4, T, P] = updateAllParams(handles);

Fs = str2num(get(handles.editFrequency, 'String')); %GHz sampling frequency
tmax = str2num(get(handles.editIntTime, 'String')); %ns integration time
numPoints = 1 + tmax*Fs; % number of data points in time domain
empty = true;
g = str2num( get(handles.editLorentz, 'String') );

% checking what we need to calculate
if (get(handles.cb1, 'Value') == 1)
    [freq, I] = HizhAll(w1, D1, G1, T, P, Fs, tmax);    
    empty = false;    
end
if (get(handles.cb2, 'Value') == 1)
    [freq2, I2] = HizhAll(w2, D2, G2, T, P, Fs, tmax);
    if empty
        freq = freq2;
        I = I2;
        empty = false;
    else
        I = conv(I, I2, 'same');
        I = I/max(I);
    end
end
if (get(handles.cb3, 'Value') == 1)
    [freq3, I3] = HizhAll(w3, D3, G3, T, P, Fs, tmax);
    if empty
        freq = freq3;
        I = I3;
        empty = false;
    else
        I = conv(I, I3, 'same');
        I = I/max(I);
    end
end
if (get(handles.cb4, 'Value') == 1)
    [freq4, I4] = HizhAll(w4, D4, G4, T, P, Fs, tmax);
    if empty
        freq = freq4;
        I = I4;
        empty = false;
    else
        I = conv(I, I4,'same');
        I = I/max(I);
    end
end


function [freqExp, IExp] = whichExp(handles)
expT = get(handles.menuT, 'Value');
switch expT
    case 1.0
        freqExp = handles.fExp80;
        IExp = handles.IExp80;
    case 2.0
        freqExp = handles.fExp100;
        IExp = handles.IExp100;
    case 3.0
        freqExp = handles.fExp120;
        IExp = handles.IExp120;
    case 4.0
        freqExp = handles.fExp150;
        IExp = handles.IExp150;    
    otherwise
        freqExp = handles.fExp80;
        IExp = handles.IExp80;    
end

function h = plotAll(handles)

[freq, I] = calculateAll(handles);

Gamma = fwhm(freq, I);
%set(handles.textFWHM, 'String', strcat('FWHM = ', num2str(Gamma), ' cm-1'));
set(handles.textFWHM, 'String', num2str(Gamma));

%if ( get(handles.cbLorentz, 'Value') == 1.0 )
%    gL = str2num(get(handles.editLorentz, 'String')); % Lorentzian broadening
%    gG = 0.40; % Gaussian broadening
%    [inhomX, inhomY] = lorentz(1.0, 0, gL, freq);
%    [instrumX, instrumY] = gauz(1.0, 0, gG, freq);
%    broadX = freq;
%    broadY = conv(inhomY, instrumY, 'same');
%    broadY = broadY/max(broadY);
%    %-------------------------------------------------------------------
%    I = conv(broadY, I, 'same');
%    I = I/max(I);
%end
% Clearing  figure
cla;
% Calculating asymmetrical lorentzian
if ( get(handles.cbLorentz, 'Value') == 1.0 )
    [af, aI] = alorentz(handles);
    plot(af, aI, 'green');    
    % convolving with the Hizhnyakov spectrum
    I = conv(aI, I, 'same');
    I = I/max(I);
end
hold on;

set(gca,'Color',[0.6,0.6,0.6])
if (get(handles.checkExperiment, 'Value') == 1)
    [freqExp, IExp] = whichExp(handles);    
    plot(freqExp, IExp, 'black');    
    hold on;
end

% Plotting Phonon Sideband
T = get(handles.ControlT, 'Value');
A = str2num( get(handles.editPhonon, 'String') );
% check where the maximum of the calculated spectrum is
z = find(I == max(I));
fmax = freq(z(1));
Ip = 0.00001*A*phonon(T, freq-fmax);
plot(freq, Ip);
hold on;
% see if we need to include the phonon sidebands 
% in the final spectrum
if (get(handles.cbPhonon, 'Value') == 1)
    I = I + Ip;
end

plot(freq, I, 'red');

hold off;

%calculating current plotting limits
z = find(I > max(I)/2.0);
fw = freq(z(end))- freq(z(1));
xcenter = 0.5*(freq(z(end))+freq(z(1)));
span = 5;
xlower = xcenter - span*fw;
xupper = xcenter + span*fw;
xlimits = [xlower xupper];
ylimits = [0 1];
xlim(xlimits);
ylim(ylimits);    

function [f, I] = HizhAll(w, D, G, T, P, Fs, tmax)

w = cast(w, 'double'); % QLM mode frequency, cm-1
D = cast(D, 'double'); % QLM frequency change, cm-1
G = 0.5*cast(G, 'double'); % QLM mode decay width, cm-1, !! IT IS HALF OF THE LOW T-Model Gamma !!
T = cast(T, 'double'); % Model temperature, K
P = cast(P, 'double'); % Emission line position, cm-1
Fs = cast(Fs, 'double');
tmax = cast(tmax, 'double');

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


function [x, y] = readFile(fileName)
fileID = fopen(fileName);
A = fscanf(fileID, '%f\t%f\n', [2 inf]);
fclose(fileID);
x = A(1, 1:end);
y = A(2, 1:end);


function [f, I] = alorentz(handles)
% generates asymmetric "Lorentzian"
T = get( handles.ControlT, 'Value' ); % current Temperature
Fs = str2num( get(handles.editFrequency, 'String') ); %GHz, sampling frequency
tmax = str2num( get(handles.editIntTime, 'String') ); %ns, integration time
numPoints = 1 + tmax*Fs; % number of data points in time domain
t = linspace(0, tmax, numPoints);
k = 30; % cm-1 to GHz=1/ns conversion

g1 = str2num( get(handles.editLorentz, 'String') ); % "Lorentzian" FWHM in cm-1
g2 = str2num( get(handles.editAsymmetry, 'String') ); %

y = exp( -g1*k*pi*abs(t) - 1i*1.44*g2*sign(t)/(4*T) );

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
I = fftshift(I);
I = I/max(I);

% Phonon sidebands: Stokes and anti-Stokes
function n = nbar(w, T)
b = 1.44/T;
n = ( exp( w*b )- 1 ).^(-1);

function Iphon = phonon(T, f)
Iphon = (f.^2).*( nbar(f, T).*heaviside(f)+(nbar(-f, T)+1).*heaviside(-f) );



function [x ,y] = lorentz(amplitude, center, width, xaxis)
x = xaxis;
y = amplitude./( 1+(2*(x-center)/width).^2 );

function [x ,y] = gauz(amplitude, center, width, xaxis)
x = xaxis;
y = amplitude*exp( -(2.40*(x-center)/width).^2 );

function g = fwhm( x, y )
z = find(y > max(y)/2.0);
g = x(z(end))- x(z(1));

function [T, g ] = generateFWHM(N, handles)
% Generates FWHM vs T for N temperature points
%overSample = 1;
T = handles.e12T(1:N);%linspace(handles.e11T(1), handles.e11T(end), overSample*length(handles.e11T) );
g = zeros(1, length(T));
j = 1;
for i = T
    set(handles.editTcurrent, 'String', num2str(i, handles.formatSpec));
    set(handles.ControlT, 'Value', i);
    [freq, I] = calculateAll(handles);
    g(j) = fwhm(freq, I);
    j = j + 1;
end


% --- Executes just before HizhnyakovModel794GUI is made visible.
function HizhnyakovModel794GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HizhnyakovModel794GUI (see VARARGIN)

% Reading experimental spectra from file
%[handles.x, handles.y] = readFile('C:\Documents and Settings\yury.dziashko\My Documents\MATLAB\Xe1e11Hizh\forAnalysis\Xe1e11_40Knorm.txt');
set(handles.editPhonon, 'String', '0.5');


handles.e12T = [80.0 100.0 120.0 150.0 180.0 240.0];
handles.e12FWHM = [15 17 20 26 32 48];

handles.formatSpec = '%10.1f';
set(handles.editFrequency, 'String', num2str(6000, handles.formatSpec));
set(handles.editIntTime, 'String', num2str(200, handles.formatSpec));

%------------Want to display experimental data?--------------------------------
handles.showExperiment = 1.0; % Show experimental spectrum
set(handles.checkExperiment, 'Value', handles.showExperiment);
%==================================================================

%----------Convinience Shift of the calculated spectra -----------
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

%=======================================================================
%-----------QLM 1 Parameters-----------
handles.W1 = 114.0;
handles.D1 = -15.0;
handles.G1 = 5.0;
set( handles.editW1,'String', num2str(handles.W1, handles.formatSpec) );
set( handles.editD1,'String', num2str(handles.D1, handles.formatSpec) );
set( handles.editG1,'String', num2str(handles.G1, handles.formatSpec) );
%-----------QLM 2 Parameters-----------
handles.W2 = 114.0;
handles.D2 = -15.0;
handles.G2 = 5.0;
set( handles.editW2,'String', num2str(handles.W2, handles.formatSpec) );
set( handles.editD2,'String', num2str(handles.D2, handles.formatSpec) );
set( handles.editG2,'String', num2str(handles.G2, handles.formatSpec) );
%-----------QLM 3 Parameters-----------
handles.W3 = 114.0;
handles.D3 = -15.0;
handles.G3 = 5.0;
set( handles.editW3,'String', num2str(handles.W3, handles.formatSpec) );
set( handles.editD3,'String', num2str(handles.D3, handles.formatSpec) );
set( handles.editG3,'String', num2str(handles.G3, handles.formatSpec) );
%-----------QLM 4 Parameters-----------
handles.W4 = 114.0;
handles.D4 = -15.0;
handles.G4 = 5.0;
set( handles.editW4,'String', num2str(handles.W4, handles.formatSpec) );
set( handles.editD4,'String', num2str(handles.D4, handles.formatSpec) );
set( handles.editG4,'String', num2str(handles.G4, handles.formatSpec) );
%======================================================================
% which modes are enabled by default?
set(handles.cb1, 'Value', 1.0); % 1st QLM ON
set(handles.cb2, 'Value', 0.0); % 2nd QLM OFF
set(handles.cb3, 'Value', 0.0); % 3rd QLM OFF
set(handles.cb4, 'Value', 0.0); % 4th QLM OFF

%-----------Temperature Related Controls and  Properties-----------
handles.Tmin = 20.0; %minumum model temperature, K
handles.Tmax = 300.0; %maximum model temperature, K
handles.T = 120.0; % Temperature of the Model
set( handles.ControlT,'Min', handles.Tmin );
set( handles.ControlT,'Max', handles.Tmax );
set( handles.ControlT,'Value', handles.T );
set( handles.editTmin, 'String', num2str(handles.Tmin, handles.formatSpec) );
set( handles.editTcurrent, 'String', num2str(handles.T, handles.formatSpec) );
set( handles.editTmax, 'String', num2str(handles.Tmax, handles.formatSpec) );
%====================================================================

%-------------------- Loading Experimental Data----------------------
% 1
[handles.fExp80, handles.IExp80] = readFile('Xe5e12_794analysis\Xe5e12_T80K_900_slit200_left.txt');
% 2
[handles.fExp100, handles.IExp100] = readFile('Xe5e12_794analysis\Xe5e12_T100K_900_slit200_left.txt');
% 3
[handles.fExp120, handles.IExp120] = readFile('Xe5e12_794analysis\Xe5e12_T120K_900_slit200_left.txt');
% 4
[handles.fExp150, handles.IExp150] = readFile('Xe5e12_794analysis\Xe5e12_T150K_900_slit200_left.txt');

% Flipping the experimental spectra
%handles.IExp80 = handles.IExp80(end:-1:1);
%handles.IExp100 = handles.IExp100(end:-1:1);
%handles.IExp120 = handles.IExp120(end:-1:1);
%andles.IExp150 = handles.IExp150(end:-1:1);

set(handles.menuT, 'Value', 1.0);
handles.freqExp = handles.fExp80;
handles.IExp = handles.IExp80;

% ------ Experimental effects --------------
lor = 0.35;
gau = 0.40;
set(handles.editLorentz, 'String', num2str(lor, handles.formatSpec));
%-------------------------------------------------------------------
set(handles.editAsymmetry, 'String', num2str(lor, handles.formatSpec));


[handles.freq, handles.I] = calculateAll(handles);

[handles.inhomX, handles.inhomY] = lorentz(1.0, 0, lor, handles.freq);
[handles.instrumX, handles.instrumY] = gauz(1.0, 0, gau, handles.freq);
handles.broadX = handles.freq;
handles.broadY = conv(handles.inhomY, handles.instrumY, 'same');
handles.broadY = handles.broadY/max(handles.broadY);

plotAll(handles);
% Choose default command line output for HizhnyakovModel794GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HizhnyakovModel794GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HizhnyakovModel794GUI_OutputFcn(hObject, eventdata, handles) 
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
function ControlT_Callback(hObject, eventdata, handles)
% hObject    handle to ControlT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set( handles.editTcurrent, 'String', num2str(get(hObject, 'Value'), handles.formatSpec) );
plotAll(handles);

function editW1_Callback(hObject, eventdata, handles)
% hObject    handle to editW1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editW1 as text
%        str2double(get(hObject,'String')) returns contents of editW1 as a double


% --- Executes during object creation, after setting all properties.
function editW1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editW1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editD1_Callback(hObject, eventdata, handles)
% hObject    handle to editD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
%set(handles.textFWHM, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));

% Hints: get(hObject,'String') returns contents of editD1 as text
%        str2double(get(hObject,'String')) returns contents of editD1 as a double


% --- Executes during object creation, after setting all properties.
function editD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function editG1_Callback(hObject, eventdata, handles)
% hObject    handle to editG1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editG1 as text
%        str2double(get(hObject,'String')) returns contents of editG1 as a double


% --- Executes during object creation, after setting all properties.
function editG1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editG1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function PositionControl_Callback(hObject, eventdata, handles)
% hObject    handle to PositionControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set( handles.editPositioncurrent, 'String', num2str(get(hObject, 'Value'), handles.formatSpec) );
plotAll(handles);
%set(handles.textFWHM, 'String', strcat('FWHM =  ', num2str(fwhm(handles.freq, handles.I), handles.formatSpec)));
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
plotAll(handles);

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
N = 6; % Number of experimental data points to calculate
[T, g] = generateFWHM(N, handles);
figure();
loglog(handles.e12T, handles.e12FWHM, 'o');
hold on;
loglog(T, g, 'red');
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
    [freqExp, IExp] = whichExp(handles);    
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


function editFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to editFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
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
plotAll(handles);
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
plotAll(handles);
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

% --- Executes during object creation, after setting all properties.
function ControlT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ControlT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function editW4_Callback(hObject, eventdata, handles)
% hObject    handle to editW4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editW4 as text
%        str2double(get(hObject,'String')) returns contents of editW4 as a double


% --- Executes during object creation, after setting all properties.
function editW4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editW4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editD4_Callback(hObject, eventdata, handles)
% hObject    handle to editD4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editD4 as text
%        str2double(get(hObject,'String')) returns contents of editD4 as a double


% --- Executes during object creation, after setting all properties.
function editD4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editG4_Callback(hObject, eventdata, handles)
% hObject    handle to editG4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editG4 as text
%        str2double(get(hObject,'String')) returns contents of editG4 as a double


% --- Executes during object creation, after setting all properties.
function editG4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editG4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb4.
function cb4_Callback(hObject, eventdata, handles)
% hObject    handle to cb4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hint: get(hObject,'Value') returns toggle state of cb4



function editW3_Callback(hObject, eventdata, handles)
% hObject    handle to editW3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editW3 as text
%        str2double(get(hObject,'String')) returns contents of editW3 as a double


% --- Executes during object creation, after setting all properties.
function editW3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editW3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editD3_Callback(hObject, eventdata, handles)
% hObject    handle to editD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editD3 as text
%        str2double(get(hObject,'String')) returns contents of editD3 as a double


% --- Executes during object creation, after setting all properties.
function editD3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editG3_Callback(hObject, eventdata, handles)
% hObject    handle to editG3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editG3 as text
%        str2double(get(hObject,'String')) returns contents of editG3 as a double


% --- Executes during object creation, after setting all properties.
function editG3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editG3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb3.
function cb3_Callback(hObject, eventdata, handles)
% hObject    handle to cb3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hint: get(hObject,'Value') returns toggle state of cb3



function editW2_Callback(hObject, eventdata, handles)
% hObject    handle to editW2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editW2 as text
%        str2double(get(hObject,'String')) returns contents of editW2 as a double


% --- Executes during object creation, after setting all properties.
function editW2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editW2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editD2_Callback(hObject, eventdata, handles)
% hObject    handle to editD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editD2 as text
%        str2double(get(hObject,'String')) returns contents of editD2 as a double


% --- Executes during object creation, after setting all properties.
function editD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editG2_Callback(hObject, eventdata, handles)
% hObject    handle to editG2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editG2 as text
%        str2double(get(hObject,'String')) returns contents of editG2 as a double


% --- Executes during object creation, after setting all properties.
function editG2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editG2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb2.
function cb2_Callback(hObject, eventdata, handles)
% hObject    handle to cb2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hint: get(hObject,'Value') returns toggle state of cb2


% --- Executes on button press in cb1.
function cb1_Callback(hObject, eventdata, handles)
% hObject    handle to cb1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hint: get(hObject,'Value') returns toggle state of cb1



function editTcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to editTcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Tcurrent = str2num(get(hObject, 'String'));
set(handles.ControlT, 'Value', handles.Tcurrent);
plotAll(handles);


% Hints: get(hObject,'String') returns contents of editTcurrent as text
%        str2double(get(hObject,'String')) returns contents of editTcurrent as a double


% --- Executes on button press in cbLorentz.
function cbLorentz_Callback(hObject, eventdata, handles)
% hObject    handle to cbLorentz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hint: get(hObject,'Value') returns toggle state of cbLorentz


% --- Executes on button press in cbAll.
function cbAll_Callback(hObject, eventdata, handles)
% hObject    handle to cbAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbAll



function editPhonon_Callback(hObject, eventdata, handles)
% hObject    handle to editPhonon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editPhonon as text
%        str2double(get(hObject,'String')) returns contents of editPhonon as a double


% --- Executes during object creation, after setting all properties.
function editPhonon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPhonon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbPhonon.
function cbPhonon_Callback(hObject, eventdata, handles)
% hObject    handle to cbPhonon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbPhonon
plotAll(handles);



function editAsymmetry_Callback(hObject, eventdata, handles)
% hObject    handle to editAsymmetry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAll(handles);
% Hints: get(hObject,'String') returns contents of editAsymmetry as text
%        str2double(get(hObject,'String')) returns contents of editAsymmetry as a double


% --- Executes during object creation, after setting all properties.
function editAsymmetry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAsymmetry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
