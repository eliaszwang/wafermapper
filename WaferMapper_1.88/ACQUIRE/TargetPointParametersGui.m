function varargout = TargetPointParametersGui(varargin)
% TARGETPOINTPARAMETERSGUI MATLAB code for TargetPointParametersGui.fig
%      TARGETPOINTPARAMETERSGUI, by itself, creates a new TARGETPOINTPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = TARGETPOINTPARAMETERSGUI returns the handle to a new TARGETPOINTPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      TARGETPOINTPARAMETERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TARGETPOINTPARAMETERSGUI.M with the given input arguments.
%
%      TARGETPOINTPARAMETERSGUI('Property','Value',...) creates a new TARGETPOINTPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TargetPointParametersGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TargetPointParametersGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TargetPointParametersGui

% Last Modified by GUIDE v2.5 05-Nov-2013 10:44:56

% Begin initialization code - DO NOT EDIT



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TargetPointParametersGui_OpeningFcn, ...
                   'gui_OutputFcn',  @TargetPointParametersGui_OutputFcn, ...
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


% --- Executes just before TargetPointParametersGui is made visible.
function TargetPointParametersGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TargetPointParametersGui (see VARARGIN)

% Choose default command line output for TargetPointParametersGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TargetPointParametersGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TargetPointParametersGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function LowResForAlignWidthInMicrons_editBox_Callback(hObject, eventdata, handles)

global GuiGlobalsStruct



MyAnswer = get(handles.LowResForAlignWidthInMicrons_editBox,'String');
MyAnswer2Num = str2num(MyAnswer);

if ~( isnan(MyAnswer2Num) || (MyAnswer2Num < 1) || (MyAnswer2Num > 4000) )
    GuiGlobalsStruct.MontageTarget.LowResForAlignWidthInMicrons = MyAnswer2Num;
    GuiGlobalsStruct.MontageTarget.LowResForAlignHeightInMicrons = MyAnswer2Num;
else
    disp('Illegal value. Not updating.');
    set(handles.LowResForAlignWidthInMicrons_editBox,'String',GuiGlobalsStruct.MontageTarget.LowResForAlignWidthInMicrons) ;
end
disp(GuiGlobalsStruct.MontageTarget.LowResForAlignHeightInMicrons)

% --- Executes during object creation, after setting all properties.
function LowResForAlignWidthInMicrons_editBox_CreateFcn(hObject, eventdata, handles)
global GuiGlobalsStruct

% hObject    handle to LowResForAlignWidthInMicrons_editBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',num2str(GuiGlobalsStruct.MontageTarget.LowResForAlignWidthInMicrons)) ;



function useCrossCorrelation_Callback(hObject, eventdata, handles)
global GuiGlobalsStruct

useCC = get(hObject,'Value');
GuiGlobalsStruct.TargetPointParameters.useCrossCorrelation = useCC;


function useCrossCorrelation_CreateFcn(hObject, eventdata, handles)
global GuiGlobalsStruct

%% Check GuiGlobalsStruct for necessary default values

if ~isfield(GuiGlobalsStruct,'TargetPointParameters')
    GuiGlobalsStruct.TargetPointParameters.useCrossCorrelation = 1;
else
    if ~isfield(GuiGlobalsStruct.TargetPointParameters,'useCrossCorrelation')
        GuiGlobalsStruct.TargetPointParameters.useCrossCorrelation =1;
    end
end
    

 set(hObject,'Value',GuiGlobalsStruct.TargetPointParameters.useCrossCorrelation);
% hObject    handle to useCrossCorrelation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
