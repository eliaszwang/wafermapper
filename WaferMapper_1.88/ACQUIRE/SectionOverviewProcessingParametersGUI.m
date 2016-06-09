function varargout = SectionOverviewProcessingParametersGUI(varargin)
% SECTIONOVERVIEWPROCESSINGPARAMETERSGUI MATLAB code for SectionOverviewProcessingParametersGUI.fig
%      SECTIONOVERVIEWPROCESSINGPARAMETERSGUI, by itself, creates a new SECTIONOVERVIEWPROCESSINGPARAMETERSGUI or raises the existing
%      singleton*.
%
%      H = SECTIONOVERVIEWPROCESSINGPARAMETERSGUI returns the handle to a new SECTIONOVERVIEWPROCESSINGPARAMETERSGUI or the handle to
%      the existing singleton*.
%
%      SECTIONOVERVIEWPROCESSINGPARAMETERSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SECTIONOVERVIEWPROCESSINGPARAMETERSGUI.M with the given input arguments.
%
%      SECTIONOVERVIEWPROCESSINGPARAMETERSGUI('Property','Value',...) creates a new SECTIONOVERVIEWPROCESSINGPARAMETERSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SectionOverviewProcessingParametersGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SectionOverviewProcessingParametersGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SectionOverviewProcessingParametersGUI

% Last Modified by GUIDE v2.5 27-Nov-2013 16:55:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SectionOverviewProcessingParametersGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SectionOverviewProcessingParametersGUI_OutputFcn, ...
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


% --- Executes just before SectionOverviewProcessingParametersGUI is made visible.
function SectionOverviewProcessingParametersGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SectionOverviewProcessingParametersGUI (see VARARGIN)

% Choose default command line output for SectionOverviewProcessingParametersGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SectionOverviewProcessingParametersGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
UpdateAllFields(handles);



function UpdateAllFields(handles)
global GuiGlobalsStruct;

set(handles.CenterAngle_edit,'String',num2str(GuiGlobalsStruct.SectionOverviewProcessingParameters.CenterAngle));
set(handles.AngleIncrement_edit,'String',num2str(GuiGlobalsStruct.SectionOverviewProcessingParameters.AngleIncrement));
set(handles.NumMultiResSteps_edit,'String',num2str(GuiGlobalsStruct.SectionOverviewProcessingParameters.NumMultiResSteps));

CenterAngle = GuiGlobalsStruct.SectionOverviewProcessingParameters.CenterAngle;
AngleIncrement = GuiGlobalsStruct.SectionOverviewProcessingParameters.AngleIncrement;
NumMultiResSteps = GuiGlobalsStruct.SectionOverviewProcessingParameters.NumMultiResSteps;
AnglesInDegreesToTryArray = [CenterAngle-2*AngleIncrement, CenterAngle-AngleIncrement,  CenterAngle,   CenterAngle+AngleIncrement ,   CenterAngle+2*AngleIncrement];






ScaleFactor = 1/(2^(NumMultiResSteps-1));
MyAccumString = '';
%Note: this mimics how the angles to try are computer in alignment function
for MultiResStepIndex = 1:NumMultiResSteps
    ScaleFactor = 1/(2^(NumMultiResSteps+-MultiResStepIndex));
    AnglesInDegreesToTryArray = [CenterAngle-2*AngleIncrement, CenterAngle-AngleIncrement,  CenterAngle,   CenterAngle+AngleIncrement ,   CenterAngle+2*AngleIncrement];
    
    MyAccumString = sprintf('%s\n%s',MyAccumString,...
        num2str(AnglesInDegreesToTryArray));
    
    %Prepare for next cycle centered on the found angle
    ScaleFactor = 2*ScaleFactor;
    AngleIncrement = AngleIncrement/2;
    CenterAngle = CenterAngle;
end

set(handles.StartingAnglesToTry_text,'String',MyAccumString);

% --- Outputs from this function are returned to the command line.
function varargout = SectionOverviewProcessingParametersGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function CenterAngle_edit_Callback(hObject, eventdata, handles)
% hObject    handle to CenterAngle_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CenterAngle_edit as text
%        str2double(get(hObject,'String')) returns contents of CenterAngle_edit as a double
global GuiGlobalsStruct;

ValueString = get(handles.CenterAngle_edit,'String');
Value = str2num(ValueString);

if ~isempty(Value) && (Value >= -180) && (Value <= 180)
    GuiGlobalsStruct.SectionOverviewProcessingParameters.CenterAngle = Value;
else
    uiwait(msgbox('Illegal value. Not updating.'));
end

UpdateAllFields(handles);

% --- Executes during object creation, after setting all properties.
function CenterAngle_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CenterAngle_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AngleIncrement_edit_Callback(hObject, eventdata, handles)
% hObject    handle to AngleIncrement_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AngleIncrement_edit as text
%        str2double(get(hObject,'String')) returns contents of AngleIncrement_edit as a double
global GuiGlobalsStruct;

ValueString = get(handles.AngleIncrement_edit,'String');
Value = str2num(ValueString);

if ~isempty(Value) && (Value >= 0) && (Value <= 90)
    GuiGlobalsStruct.SectionOverviewProcessingParameters.AngleIncrement = Value;
else
    uiwait(msgbox('Illegal value. Not updating.'));
end

UpdateAllFields(handles);

% --- Executes during object creation, after setting all properties.
function AngleIncrement_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AngleIncrement_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NumMultiResSteps_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NumMultiResSteps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumMultiResSteps_edit as text
%        str2double(get(hObject,'String')) returns contents of NumMultiResSteps_edit as a double
global GuiGlobalsStruct;

ValueString = get(handles.NumMultiResSteps_edit,'String');
Value = str2num(ValueString);

if ~isempty(Value) && (Value >= 1) && (Value <= 7)
    GuiGlobalsStruct.SectionOverviewProcessingParameters.NumMultiResSteps = Value;
else
    uiwait(msgbox('Illegal value. Not updating.'));
end

UpdateAllFields(handles);

% --- Executes during object creation, after setting all properties.
function NumMultiResSteps_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumMultiResSteps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ResetToDefault_button.
function ResetToDefault_button_Callback(hObject, eventdata, handles)
% hObject    handle to ResetToDefault_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SetSectionOverviewProcessingParametersDefaults();
UpdateAllFields(handles);



function StartingAnglesToTry_edit_Callback(hObject, eventdata, handles)
% hObject    handle to StartingAnglesToTry_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartingAnglesToTry_edit as text
%        str2double(get(hObject,'String')) returns contents of StartingAnglesToTry_edit as a double


% --- Executes during object creation, after setting all properties.
function StartingAnglesToTry_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartingAnglesToTry_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SmallestAngleResolution_edit_Callback(hObject, eventdata, handles)
% hObject    handle to SmallestAngleResolution_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmallestAngleResolution_edit as text
%        str2double(get(hObject,'String')) returns contents of SmallestAngleResolution_edit as a double


% --- Executes during object creation, after setting all properties.
function SmallestAngleResolution_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmallestAngleResolution_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function SO_AlignmentType_uigroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SO_AlignmentType_uigroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%function SO_AlignmentType_uigroup_SelectionChangeFnc(hObject, eventdata, handles)
function SO_AlignmentType_uigroup_SelectionChangeFcn(hObject, eventdata, handles)

disp('bark')
% 
% %function selcbk(source,eventdata)
% disp(source);
% disp([eventdata.EventName,'  ',... 
%      get(eventdata.OldValue,'String'),'  ', ...
%      get(eventdata.NewValue,'String')]);
% disp(get(get(source,'SelectedObject'),'String'))



% 
% 
% function AutofocusMethod_uipanel_SelectionChangeFcn(hObject, eventdata, handles)
% % hObject    handle to the selected object in AutofocusMethod_uipanel 
% % eventdata  structure with the following fields (see UIBUTTONGROUP)
% %	EventName: string 'SelectionChanged' (read only)
% %	OldValue: handle of the previously selected object or empty if none was selected
% %	NewValue: handle of the currently selected object
% % handles    structure with handles and user data (see GUIDATA)
% global GuiGlobalsStruct;
% 
% %set(handles.IsSingle_AF_ForWholeMontage_radiobutton, 'Value', GuiGlobalsStruct.MontageParameters.IsSingle_AF_ForWholeMontage);
% Value = get(handles.IsSingle_AF_ForWholeMontage_radiobutton,'Value');
% GuiGlobalsStruct.MontageParameters.IsSingle_AF_ForWholeMontage = Value;
% 
% Value = get(handles.IsSingle_AFASAF_ForWholeMontage_radiobutton, 'Value');
% GuiGlobalsStruct.MontageParameters.IsSingle_AFASAF_ForWholeMontage = Value;
% 
% Value = get(handles.IsAFOnEveryTile_radiobutton, 'Value');
% GuiGlobalsStruct.MontageParameters.IsAFOnEveryTile = Value;
% 
% Value = get(handles.IsAFASAFOnEveryTile_radiobutton, 'Value');
% GuiGlobalsStruct.MontageParameters.IsAFASAFOnEveryTile = Value;
% 
% Value = get(handles.IsPlaneFit_radiobutton, 'Value');
% GuiGlobalsStruct.MontageParameters.IsPlaneFit = Value;
% 
% Value = get(handles.Is4square_radiobutton, 'Value');
% GuiGlobalsStruct.MontageParameters.Is4square = Value;
% 
% Value = get(handles.IsXFit_radiobutton, 'Value');
% GuiGlobalsStruct.MontageParameters.IsXFit = Value;
% 
% UpdateAllFields(handles);
% 
% 
% 
















% 




function MaxRes_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MaxRes_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxRes_edit as text
%        str2double(get(hObject,'String')) returns contents of MaxRes_edit as a double
global GuiGlobalsStruct;

ValueString = get(handles.MaxRes_edit,'String');
Value = str2num(ValueString);

if ~isempty(Value) && (Value >= 64)
    GuiGlobalsStruct.SectionOverviewProcessingParameters.MaxRes = Value;
else
    uiwait(msgbox(sprintf('%s is not valid.  Try 64 or above.',ValueString)));
end

% --- Executes during object creation, after setting all properties.
function MaxRes_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxRes_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global GuiGlobalsStruct;

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


try
    MaxRes = GuiGlobalsStruct.SectionOverviewProcessingParameters.MaxRes;
catch err

    MaxRes = 512;
      GuiGlobalsStruct.SectionOverviewProcessingParameters.MaxRes =MaxRes;

end

set(hObject,'String',num2str(MaxRes))


% --- Executes on button press in checkbox_SURF.
function checkbox_SURF_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SURF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SURF


global GuiGlobalsStruct;


Value = get(hObject,'Value');
GuiGlobalsStruct.SectionOverviewProcessingParameters.SURF = Value;

if ~Value
    set(handles.SURFall_checkbox,'Value',Value);
    GuiGlobalsStruct.SectionOverviewProcessingParameters.SURFall = Value;

end


% --- Executes during object creation, after setting all properties.
function checkbox_SURF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_SURF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global GuiGlobalsStruct;

try
Value = GuiGlobalsStruct.SectionOverviewProcessingParameters.SURF;
catch err
    Value = 0;
    GuiGlobalsStruct.SectionOverviewProcessingParameters.SURF = Value;
end

set(hObject,'Value',Value)


function SURFall_checkbox_Callback(hObject, eventdata, handles)

global GuiGlobalsStruct;

Value = get(hObject,'Value');
GuiGlobalsStruct.SectionOverviewProcessingParameters.SURFall = Value;

if Value %turn Surf on in Surf all is on
set(handles.checkbox_SURF,'Value',Value);
GuiGlobalsStruct.SectionOverviewProcessingParameters.SURF = Value;
end


function SURFall_checkbox_CreateFcn(hObject, eventdata, handles)

global GuiGlobalsStruct;

try
Value = GuiGlobalsStruct.SectionOverviewProcessingParameters.SURFall;
catch err
    Value = 0;
    GuiGlobalsStruct.SectionOverviewProcessingParameters.SURFall = Value;
end

set(hObject,'Value',Value)
