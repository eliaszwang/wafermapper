function varargout = ProcessOverviewsRobustRigidParameters2(varargin)
% PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2 MATLAB code for ProcessOverviewsRobustRigidParameters2.fig
%      PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2, by itself, creates a new PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2 or raises the existing
%      singleton*.
%
%      H = PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2 returns the handle to a new PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2 or the handle to
%      the existing singleton*.
%
%      PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2.M with the given input arguments.
%
%      PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2('Property','Value',...) creates a new PROCESSOVERVIEWSROBUSTRIGIDPARAMETERS2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcessOverviewsRobustRigidParameters2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcessOverviewsRobustRigidParameters2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcessOverviewsRobustRigidParameters2

% Last Modified by GUIDE v2.5 01-Aug-2013 15:07:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessOverviewsRobustRigidParameters2_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessOverviewsRobustRigidParameters2_OutputFcn, ...
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


% --- Executes just before ProcessOverviewsRobustRigidParameters2 is made visible.
function ProcessOverviewsRobustRigidParameters2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcessOverviewsRobustRigidParameters2 (see VARARGIN)

% Choose default command line output for ProcessOverviewsRobustRigidParameters2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProcessOverviewsRobustRigidParameters2 wait for user response (see UIRESUME)
UpdateAllFields(handles);
if ~isempty(varargin)
    SetSURFOptionsFromOptions(varargin{1},handles)
    SetRanSacOptionsFromOptions(varargin{2},handles)
    SetMetaAlignmentOptionsFromOptions(varargin{3},handles)
end
uiwait(handles.figure1);



function UpdateAllFields(handles)
global GuiGlobalsStruct;
filestruct = dir([GuiGlobalsStruct.SectionOverviewsDirectory filesep '*.tif']);
numsections=length(filestruct);
set(handles.slider1,'Max',numsections);
set(handles.slider1,'SliderStep',[1/numsections 5/numsections]);
set(handles.slider2,'Max',numsections);
set(handles.slider2,'SliderStep',[1/numsections 5/numsections]);

% --- Outputs from this function are returned to the command line.
function varargout = ProcessOverviewsRobustRigidParameters2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = GetSURFOptionsFromHandles(handles);
varargout{2} = GetRanSacOptionsFromHandles(handles);
varargout{3} = GetMetaAlignmentOptionsFromHandles(handles);
delete(handles.figure1);


% --- Executes on button press in MexicanHat_CheckBox.
function MexicanHat_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to MexicanHat_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MexicanHat_CheckBox



function InitialDownsample_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to InitialDownsample_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InitialDownsample_EditBox as text
%        str2double(get(hObject,'String')) returns contents of InitialDownsample_EditBox as a double


% --- Executes during object creation, after setting all properties.
function InitialDownsample_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InitialDownsample_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function Octaves_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOctaves_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Octaves_EditBox_Callback(hObject,eventdata,handles)

function Threshold_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Threshold_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Threshold_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CenterFraction_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to CenterFraction_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CenterFraction_EditBox as text
%        str2double(get(hObject,'String')) returns contents of CenterFraction_EditBox as a double


% --- Executes during object creation, after setting all properties.
function CenterFraction_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CenterFraction_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ReNormalize_checkbox.
function ReNormalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ReNormalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ReNormalize_checkbox



function Inner_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Inner_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inner_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Inner_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Inner_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inner_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Outer_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Outer_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Outer_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Outer_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Outer_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Outer_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sigma_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Sigma_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sigma_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Sigma_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Sigma_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sigma_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DistThresh_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to DistThresh_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DistThresh_EditBox as text
%        str2double(get(hObject,'String')) returns contents of DistThresh_EditBox as a double


% --- Executes during object creation, after setting all properties.
function DistThresh_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DistThresh_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NBest_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to NBest_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NBest_EditBox as text
%        str2double(get(hObject,'String')) returns contents of NBest_EditBox as a double


% --- Executes during object creation, after setting all properties.
function NBest_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NBest_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinInliers_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to MinInliers_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinInliers_EditBox as text
%        str2double(get(hObject,'String')) returns contents of MinInliers_EditBox as a double


% --- Executes during object creation, after setting all properties.
function MinInliers_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinInliers_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Delta_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Delta_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Delta_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Delta_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Delta_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Delta_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lambda_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Lambda_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lambda_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Lambda_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Lambda_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lambda_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TimeSteps_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSteps_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeSteps_EditBox as text
%        str2double(get(hObject,'String')) returns contents of TimeSteps_EditBox as a double


% --- Executes during object creation, after setting all properties.
function TimeSteps_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeSteps_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxDetChange_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDetChange_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxDetChange_EditBox as text
%        str2double(get(hObject,'String')) returns contents of MaxDetChange_EditBox as a double


% --- Executes during object creation, after setting all properties.
function MaxDetChange_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDetChange_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value=get(hObject,'Value');
set(hObject,'Value',round(value));
set(handles.Section1_EditBox,'String',num2str(int8(round(value))));


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value=get(hObject,'Value');
set(hObject,'Value',round(value));
set(handles.Section2_EditBox,'String',num2str(int8(round(value))));


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Section1_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Section1_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Section1_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Section1_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Section1_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Section1_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Section2_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Section2_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Section2_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Section2_EditBox as a double


% --- Executes during object creation, after setting all properties.
function Section2_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Section2_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [I1,I2,sect1,sect2]=readUserSelectedImages(handles)
global GuiGlobalsStruct;
%get the list of files and directories
SURFOptions=GetSURFOptionsFromHandles(handles);
[Files,MatFiles,labels]=GetSortedImagesAndMatfiles(GuiGlobalsStruct.SectionOverviewsDirectory);

%pull out which sections the user has chosen
sect1=round(str2double(get(handles.Section1_EditBox,'String')));
sect2=round(str2double(get(handles.Section2_EditBox,'String')));

%calculate the Pixel region corresponding to the given center_fraction
PixelRegion=DefinePixelRegionFromCenterFrac(Files{1},SURFOptions.center_frac);

I1=imread(Files{sect1},'PixelRegion',PixelRegion);
I2=imread(Files{sect2},'PixelRegion',PixelRegion);

% --- Executes on button press in TestSURF_Button.
function TestSURF_Button_Callback(hObject, eventdata, handles)
% hObject    handle to TestSURF_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GuiGlobalsStruct;

%pull out the SURF options from the GUI
SURFOptions=GetSURFOptionsFromHandles(handles);

%read in images user selected (pulling out center_frac inside function)
[I1,I2,sect1,sect2]=readUserSelectedImages(handles);
I1 = PreFilterImage(I1,SURFOptions);
I2 = PreFilterImage(I2,SURFOptions);

%extract the SURF points
p1=OpenSurf(I1,SURFOptions);
p2=OpenSurf(I2,SURFOptions);

%update GUI with the number of points
set(handles.NumPtsIn1_Text,'String',sprintf('# pts in 1: %d',length(p1)));
set(handles.NumPtsIn2_Text,'String',sprintf('# pts in 2: %d',length(p2)));

%plot out the results
%shift the points for image 2 right by the width of image 1
for j=1:length(p2)
    p2(j).x=p2(j).x+size(I1,2);
end

%make a concatenated image
I = zeros([size(I1,1) size(I1,2)*2 size(I1,3)],'double');
I(:,1:size(I1,2),:)=double(I1);
I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=double(I2);
I=I/max(I(:));
%select main axes and make sure whatever comes next overwrites what is
%there
axes(handles.MainAxes);
hold off;
%user a modified version of OpenSurf's function to plot points
PaintSURF_modified(I,[p1 p2]);
hold on;
%draw a vertical line to seperate two images
line([size(I1,2) size(I1,2)],[1 size(I1,1)]);

function SetSURFOptionsFromOptions(Options,handles)
% handles    structure with handles and user data (see GUIDATA)
if isfield(Options,'init_sample'),set(handles.InitialDownsample_EditBox,'String',num2str(Options.init_sample));end
if isfield(Options,'octaves'),set(handles.Octaves_EditBox,'String',num2str(Options.octaves));end %how many octaves of upsampling to try
if isfield(Options,'tresh'),set(handles.Threshold_EditBox,'String',num2str(Options.tresh));end %what threshold to use (lower for more features, raise for fewer)
if isfield(Options,'center_frac'),set(handles.CenterFraction_EditBox,'String',num2str(Options.center_frac));end
if isfield(Options,'MexicanHat'),set(handles.MexicanHat_CheckBox,'Value',Options.MexicanHat);end
if isfield(Options,'Outer'),set(handles.Outer_EditBox,'String',num2str(Options.Outer));end
if isfield(Options,'Inner'),set(handles.Inner_EditBox,'String',num2str(Options.Inner));end
if isfield(Options,'ReNormalize'),set(handles.ReNormalize_checkbox,'Value',Options.ReNormalize);end
if isfield(Options,'Sigma'),set(handles.Sigma_EditBox,'String',num2str(Options.Sigma));end

function SetRanSacOptionsFromOptions(Options,handles)
% handles    structure with handles and user data (see GUIDATA)
set(handles.DistThresh_EditBox,'String',num2str(Options.dist_thresh));
set(handles.NBest_EditBox,'String',num2str(Options.NBest)); %how many octaves of upsampling to try
eval(['set(handles.Model_Panel,''SelectedObject'',handles.' Options.Model ');']); %what threshold to use (lower for more features, raise for fewer)
set(handles.MaxDetChange_EditBox,'String',num2str(Options.MaxDetChange));

function SetMetaAlignmentOptionsFromOptions(Options,handles)
% handles    structure with handles and user data (see GUIDATA)


set(handles.MinInliers_EditBox,'String',num2str(Options.MinInliers));
switch Options.Method
    case{'GradDescent'}
        set(handles.OptMethod_Panel,'SelectedObject',handles.GradDescent_radiobutton)
        set(handles.LocalDepth_EditBox,'String',num2str(Options.LocalDepth));
        set(handles.LongRangeDelta_EditBox,'String',num2str(Options.LongRangeDelta));
        set(handles.LongRangeDepth_EditBox,'String',num2str(Options.LongRangeDepth));  
        set(handles.TimeSteps_EditBox,'String',num2str(Options.TimeSteps));
        set(handles.Delta_EditBox,'String',num2str(Options.Delta));
        
        ToggleMatrixInversionOptions(handles,'off');
        ToggleShortBridgeOptions(handles,'off');
        ToggleGradDescentOptions(handles,'on');
        
    case{'MatrixInversion'}
        set(handles.OptMethod_Panel,'SelectedObject',handles.MatrixInversion_radiobutton)
        set(handles.LocalDepth_EditBox,'String',num2str(Options.LocalDepth));
        set(handles.LongRangeDelta_EditBox,'String',num2str(Options.LongRangeDelta));
        set(handles.LongRangeDepth_EditBox,'String',num2str(Options.LongRangeDepth));
        set(handles.Lambda_EditBox,'String',num2str(Options.Lambda));
        
        ToggleShortBridgeOptions(handles,'off');
        ToggleGradDescentOptions(handles,'off');
        ToggleMatrixInversionOptions(handles,'on');
    case{'ShortBridge'}
        set(handles.OptMethod_Panel,'SelectedObject',handles.ShortBridge_radiobutton)
        set(handles.MaxDelta_EditBox,'String',num2str(Options.MaxDelta));
        
        ToggleMatrixInversionOptions(handles,'off');
        ToggleGradDescentOptions(handles,'off');
        ToggleShortBridgeOptions(handles,'on');
    otherwise
        disp('Error!');
end



function Options=GetSURFOptionsFromHandles(handles)
% handles    structure with handles and user data (see GUIDATA)
Options.verbose=0; %whether to produce a verbose output for debugging
Options.init_sample=round(str2double(get(handles.InitialDownsample_EditBox,'String')));
Options.octaves=round(str2double(get(handles.Octaves_EditBox,'String'))); %how many octaves of upsampling to try
Options.tresh=str2double(get(handles.Threshold_EditBox,'String')); %what threshold to use (lower for more features, raise for fewer)
Options.center_frac=str2double(get(handles.CenterFraction_EditBox,'String'));
Options.MexicanHat=get(handles.MexicanHat_CheckBox,'Value');
Options.Outer=round(str2double(get(handles.Outer_EditBox,'String')));
Options.Inner=round(str2double(get(handles.Inner_EditBox,'String')));
Options.ReNormalize=get(handles.ReNormalize_checkbox,'Value');
Options.Sigma=str2double(get(handles.Sigma_EditBox,'String'));

function Options=GetRanSacOptionsFromHandles(handles)
% handles    structure with handles and user data (see GUIDATA)
Options.verbose=0; %whether to produce a verbose output for debugging
Options.dist_thresh=str2double(get(handles.DistThresh_EditBox,'String'));
Options.NBest=round(str2double(get(handles.NBest_EditBox,'String'))); %how many octaves of upsampling to try
Options.Model=get(get(handles.Model_Panel,'SelectedObject'),'Tag'); %what threshold to use (lower for more features, raise for fewer)
Options.MaxDetChange=str2double(get(handles.MaxDetChange_EditBox,'String'));

function Options=GetMetaAlignmentOptionsFromHandles(handles)
% handles    structure with handles and user data (see GUIDATA)

selectedTag=get(get(handles.OptMethod_Panel,'SelectedObject'),'Tag');
Options.MinInliers=round(str2double(get(handles.MinInliers_EditBox,'String')));
switch selectedTag
    case{'GradDescent_radiobutton'}
        Options.Method='GradDescent';
        
        Options.LocalDepth=round(str2double(get(handles.LocalDepth_EditBox,'String')));
        Options.LongRangeDelta=round(str2double(get(handles.LongRangeDelta_EditBox,'String')));
        Options.LongRangeDepth=round(str2double(get(handles.LongRangeDepth_EditBox,'String')));
        
        Options.TimeSteps=round(str2double(get(handles.TimeSteps_EditBox,'String')));
        Options.Delta=str2double(get(handles.Delta_EditBox,'String'));

    case{'MatrixInversion_radiobutton'}
        Options.Method='MatrixInversion';
    
        Options.LocalDepth=round(str2double(get(handles.LocalDepth_EditBox,'String')));
        Options.LongRangeDelta=round(str2double(get(handles.LongRangeDelta_EditBox,'String')));
        Options.LongRangeDepth=round(str2double(get(handles.LongRangeDepth_EditBox,'String')));
        
        Options.Lambda=str2double(get(handles.Lambda_EditBox,'String'));
      
    case{'ShortBridge_radiobutton'}
        Options.Method='ShortBridge';
        Options.MaxDelta=round(str2double(get(handles.MaxDelta_EditBox,'String')));
        
    otherwise
        disp('Error!');
end


% --- Executes on button press in TestRanSac_Button.
function TestRanSac_Button_Callback(hObject, eventdata, handles)
% hObject    handle to TestRanSac_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GuiGlobalsStruct;


%get the SURF and RANSAC options from the user
SURFOptions=GetSURFOptionsFromHandles(handles);
RanSacOptions=GetRanSacOptionsFromHandles(handles);
tic
%read in images user selected (pulling out center_frac inside function)
[I1,I2,sect1,sect2]=readUserSelectedImages(handles);

I1 = PreFilterImage(I1,SURFOptions);
I2 = PreFilterImage(I2,SURFOptions);
%extract the SURF features
p1=OpenSurf(I1,SURFOptions);
p2=OpenSurf(I2,SURFOptions);
t1=toc; %time after reading both images and calculating surf points
%find the best correspondances between the two images
p1=shiftSURFpoints(p1,-size(I1,2)/2,-size(I1,1)/2);
p2=shiftSURFpoints(p2,-size(I2,2)/2,-size(I2,1)/2);


[Files,MatFiles,labels]=GetSortedImagesAndMatfiles(GuiGlobalsStruct.SectionOverviewsDirectory);
PixelRegion=DefinePixelRegionFromCenterFrac(Files{sect1},SURFOptions.center_frac);
[p1,p2]=RemovePointsFromOverlapZone(p1,p2,PixelRegion,MatFiles{sect1},MatFiles{sect2});
 
          
[Pos1,Pos2]=find_best_SURF_match(p1,p2,RanSacOptions.NBest);

%setup the ransac options
ransacCoef.minPtNum=2;   %for a rigid or similarity transform 2 is the number needed
ransacCoef.iterNum= round(.25*nchoosek(size(Pos1,1),2)); %run through a number of iterations equal to the number of points
ransacCoef.thDist= RanSacOptions.dist_thresh*size(I1,1); %should be within 20 microns to be right
ransacCoef.thInlrRatio=.02; %at least 5 percent should be right
ransacCoef.thDet=RanSacOptions.MaxDetChange;

%run Ransac
[f inlierIdx] = ransac1( flipdim(Pos1',1),flipdim(Pos2',1),ransacCoef,@fit_rigid,@EuclideanDistance,@getDeterminantOfAffineTransform);
M=f.tdata.T(1:2,1:2);
t2=toc;%time after running ransac

%update GUI with information on SURF points and inliers after ransac
set(handles.NumInliers_Text,'String',sprintf('NumInliers: %d',length(inlierIdx)));
set(handles.NumPtsIn1_Text,'String',sprintf('# pts in 1: %d',length(p1)));
set(handles.NumPtsIn2_Text,'String',sprintf('# pts in 2: %d',length(p2)));
set(handles.Determinant_Text','String',sprintf('Determinant^.5: %3.2f',sqrt(det(M))));
readSURFtime=t1/2;
ransactime=t2-t1;
numsections=get(handles.slider1,'Max');
numadjacent=round(str2double(get(handles.LocalDepth_EditBox,'String')));
numlong=round(str2double(get(handles.LongRangeDepth_EditBox,'String')));

totaltime_sec=readSURFtime*numsections+numadjacent*numsections+numlong*numsections;
totaltime_min=totaltime_sec/60;
if totaltime_min<120
    set(handles.EstTime_Text,'String',sprintf('Estimated Total Time:%3.2f min',totaltime_min));
else
    set(handles.EstTime_Text,'String',sprintf('Estimated Total Time:%3.2f hrs',totaltime_min/60));
end

%plot the correspondances in the MainAxes
axes(handles.MainAxes);
hold off;
PlotCorrespondances(I1,I2,Pos1,Pos2,inlierIdx);


function LocalDepth_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to LocalDepth_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LocalDepth_EditBox as text
%        str2double(get(hObject,'String')) returns contents of LocalDepth_EditBox as a double


% --- Executes during object creation, after setting all properties.
function LocalDepth_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LocalDepth_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LongRangeDelta_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to LongRangeDelta_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LongRangeDelta_EditBox as text
%        str2double(get(hObject,'String')) returns contents of LongRangeDelta_EditBox as a double


% --- Executes during object creation, after setting all properties.
function LongRangeDelta_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LongRangeDelta_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LongRangeDepth_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to LongRangeDepth_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LongRangeDepth_EditBox as text
%        str2double(get(hObject,'String')) returns contents of LongRangeDepth_EditBox as a double


% --- Executes during object creation, after setting all properties.
function LongRangeDepth_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LongRangeDepth_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end



function MaxDelta_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDelta_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxDelta_EditBox as text
%        str2double(get(hObject,'String')) returns contents of MaxDelta_EditBox as a double


% --- Executes during object creation, after setting all properties.
function MaxDelta_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDelta_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ToggleShortBridgeOptions(handles,onoff)
   set(handles.MaxDelta_EditBox,'Visible',onoff);
   set(handles.MaxDelta_Text,'Visible',onoff);
 
   
function ToggleGradDescentOptions(handles,onoff)
   set(handles.LocalDepth_EditBox,'Visible',onoff);
   set(handles.LocalDepth_Text,'Visible',onoff);
   set(handles.LongRangeDelta_EditBox,'Visible',onoff);
   set(handles.LongRangeDepth_EditBox,'Visible',onoff);
   set(handles.LongRange_Text,'Visible',onoff);
   set(handles.TimeSteps_EditBox,'Visible',onoff);
   set(handles.TimeSteps_Text,'Visible',onoff);
   set(handles.Delta_EditBox,'Visible',onoff);
   set(handles.Delta_Text,'Visible',onoff);

function ToggleMatrixInversionOptions(handles,onoff)
   set(handles.LocalDepth_EditBox,'Visible',onoff);
   set(handles.LocalDepth_Text,'Visible',onoff);
   set(handles.LongRangeDelta_EditBox,'Visible',onoff);
   set(handles.LongRangeDepth_EditBox,'Visible',onoff);
   set(handles.LongRange_Text,'Visible',onoff);
   set(handles.Lambda_EditBox,'Visible',onoff);
   set(handles.Lambda_Text,'Visible',onoff);


% --- Executes when selected object is changed in OptMethod_Panel.
function OptMethod_Panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in OptMethod_Panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

selectedTag=get(hObject,'Tag');

switch selectedTag
    case{'GradDescent_radiobutton'}     
        ToggleMatrixInversionOptions(handles,'off');
        ToggleShortBridgeOptions(handles,'off');
        ToggleGradDescentOptions(handles,'on'); 
    case{'MatrixInversion_radiobutton'} 
        ToggleShortBridgeOptions(handles,'off');
        ToggleGradDescentOptions(handles,'off');
        ToggleMatrixInversionOptions(handles,'on');
    case{'ShortBridge_radiobutton'}
        ToggleGradDescentOptions(handles,'off');
        ToggleMatrixInversionOptions(handles,'off');
        ToggleShortBridgeOptions(handles,'on');
               
    otherwise
        disp('Error!')
end
