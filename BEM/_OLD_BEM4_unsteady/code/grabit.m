function grabit(fname)

%GRABIT Extracts data points from an image file.
%
% GRABIT starts a GUI program for extracting data from an image file.
% It is capable of reading in BMP, JPG, TIF, GIF, and PNG files (anything
% that is readable by IMREAD).  Multiple data sets can be extracted from a
% single image file, and the data is saved as an n-by-2 matrix variable in
% the workspace.  It can also be renamed and saved as a MAT file.
%
% Following steps should be taken:
%   1. Load the image file.
%   2. Calibrate axes dimensions.  You will be prompted to select 4 points
%      on the image. Zoom and pan enabled.
%   3. Grab points by clicking on points.  Right-click to delete a point.
%      Image can be zoomed  and panned.
%   4. Multiple data sets will remain in memory so long as the GUI is open.
%      Variables can be renamed, saved to file, or edited in Array Editor.
% 
% Panning is achieved by clicking and dragging on the image. Double-click
% to center view. Right click and drag to zoom in and out. In addition,
% there are keyboard shortcuts for zooming:
%   <a>     - zoom in
%   <b>     - zoom out
%   <space> - reset view
%
% This code will also work for extracting data points from a tilted or a
% skewed image (even upside-down or mirrored).  The calibration stage
% ensures that the imperfect orientation or quality of the image is
% accounted for.
%
% The types of files that will most likely work are BMP, JPG, TIF, GIF (up
% to 8-bit), and PNG files.  Basically,  any format supported by the IMREAD
% is accepted.
%
% GRABIT(FILENAME) will start the GUI program and open the image file
% FILENAME.
%
% Type GRABIT('-sample') to load a sample image.
%
% 
% VERSIONS:
%   v1.0 - first version
%   v1.1 - use imshow instead of image (takes care of colormap)
%   v1.5 - convert to a GUI version
%   v1.6 - added functionality to open a file from the command window and
%          embedded a sample image file to the function
%   v1.6.1 - changed cross cursor to crosshair
%   v1.6.2 - brought back 'image' in case the user doesn't have Image Toolbox
%   v1.6.5 - fixed zoom problem in R14
%   v2.0 - major code change. added zoom feature during calibration. added
%          panning feature. (March 3, 2006)
%   v2.1 - store sample image as HEX to reduce file size. (March 6, 2006)
%   v2.1.1 - animate view change and zooming for a better visual perception
%            (March 11, 2006)
%   v2.1.5 - added features: double-click to center view. right-click and
%            drag to zoom. other minor code changes. (March 16, 2006)
%   v2.2 - fixed loadImageFcn bug. (May 4, 2006)
%   v2.3 - fixed bug to work with grayscale image (Jan, 2007)
%

% Created in Matlab R13. Tested up to R2006b
%
% Copyright 2003
% Jiro Doke
%
% To Do: Capability to deal with logarithmic axes
%

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
sh = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');

% close existing windows
im  =  findobj('type',  'figure',  'tag',  'GrabitGUI');
if ishandle(im)
  close(im);
end

set(0, 'ShowHiddenHandles', sh);

% background colors
bgcolor1 = [.8, .8, .8];
bgcolor2 = [ 1,  1,  1];
bgcolor3 = [.7, .7, .7];
bgcolor4 = [ 1, .5, .5];

%--------------------------------------------------------------------------
% Custom cursor pointers
%--------------------------------------------------------------------------
zoomPointer = [
     2     2     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     1     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2   NaN   NaN   NaN   NaN     2     2     2     2   NaN   NaN   NaN   NaN   NaN
     2     1     2   NaN   NaN     2     2     1     1     1     1     2     2   NaN   NaN   NaN
     2     1     2   NaN     2     1     1     2     2     2     2     1     1     2   NaN   NaN
   NaN     2   NaN   NaN     2     1     2     2     1     1     2     2     1     2   NaN   NaN
   NaN   NaN   NaN     2     1     2     2     2     1     1     2     2     2     1     2   NaN
   NaN   NaN   NaN     2     1     2     1     1     1     1     1     1     2     1     2   NaN
   NaN   NaN   NaN     2     1     2     1     1     1     1     1     1     2     1     2   NaN
   NaN   NaN   NaN     2     1     2     2     2     1     1     2     2     2     1     2   NaN
   NaN   NaN   NaN   NaN     2     1     2     2     1     1     2     2     1     2   NaN   NaN
   NaN   NaN   NaN   NaN     2     1     1     2     2     2     2     1     1     1     2   NaN
   NaN   NaN   NaN   NaN   NaN     2     2     1     1     1     1     2     1     1     1     2
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     2     2     2   NaN     2     1     1     1
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     2
 ];

zoomInOutPointer = [
   NaN   NaN   NaN     2     2   NaN   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN   NaN
   NaN   NaN     2     1     1     2   NaN   NaN   NaN   NaN     1     2     2     2     1   NaN
   NaN     2     1     1     1     1     2   NaN   NaN     1     2     2     1     2     2     1
     2     1     1     1     1     1     1     2   NaN     1     2     1     1     1     2     1
     2     1     2     1     1     2     1     2   NaN     1     2     2     1     2     2     1
   NaN     2     2     1     1     2     2   NaN   NaN     1     1     2     2     2     1   NaN
   NaN   NaN     2     1     1     2   NaN   NaN     1     1     1     1     1     1   NaN   NaN
   NaN   NaN   NaN     2     2   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     2     2   NaN   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN   NaN
   NaN   NaN     2     1     1     2   NaN   NaN   NaN   NaN     1     2     2     2     1   NaN
   NaN     2     2     1     1     2     2   NaN   NaN     1     2     2     2     2     2     1
     2     1     2     1     1     2     1     2   NaN     1     2     1     1     1     2     1
     2     1     1     1     1     1     1     2   NaN     1     2     2     2     2     2     1
   NaN     2     1     1     1     1     2   NaN   NaN     1     1     2     2     2     1   NaN
   NaN   NaN     2     1     1     2   NaN   NaN     1     1     1     1     1     1   NaN   NaN
   NaN   NaN   NaN     2     2   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
 ];
 
% closed hand pointer (from Jérôme Briot)
closedHandPointer = [
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN     2     2   NaN     2     2   NaN     2     2   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     2     1     1     2     1     1     2     1     1     2     2   NaN   NaN
   NaN   NaN     2     1     2     2     1     2     2     1     2     2     1     1     2   NaN
   NaN   NaN     2     1     2     2     2     2     2     2     2     2     1     2     1     2
   NaN   NaN   NaN     2     1     2     2     2     2     2     2     2     2     2     1     2
   NaN   NaN     2     1     1     2     2     2     2     2     2     2     2     2     1     2
   NaN     2     1     2     2     2     2     2     2     2     2     2     2     2     1     2
   NaN     2     1     2     2     2     2     2     2     2     2     2     2     2     1     2
   NaN     2     1     2     2     2     2     2     2     2     2     2     2     1     2   NaN
   NaN   NaN     2     1     2     2     2     2     2     2     2     2     2     1     2   NaN
   NaN   NaN   NaN     2     1     2     2     2     2     2     2     2     1     2   NaN   NaN
   NaN   NaN   NaN   NaN     2     1     2     2     2     2     2     2     1     2   NaN   NaN
   NaN   NaN   NaN   NaN     2     1     2     2     2     2     2     2     1     2   NaN   NaN
 ];

% X-axis Origin pointer
xoPointer = [
     2     2     2     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     1     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     2     2   NaN   NaN     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1     1   NaN     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN     1     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN   NaN   NaN   NaN     1     1     1   NaN
   NaN   NaN   NaN   NaN     1     1     1     1     1   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN     1     1     1   NaN     1     1     1   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN     1     1   NaN   NaN   NaN     1     1   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN
 ];

% X-axis Max pointer
xmPointer = [
     2     2     2     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     1     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     2     2   NaN   NaN     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1     1   NaN     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN     1     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN   NaN     1     1     1     1     1   NaN   NaN     1     1   NaN     1     1
   NaN   NaN   NaN     1     1     1   NaN     1     1     1   NaN     1   NaN     1   NaN     1
   NaN   NaN   NaN     1     1   NaN   NaN   NaN     1     1   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN     1
 ];
 
% Y-axis Origin pointer
yoPointer = [
     2     2     2     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     1     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     2     2   NaN   NaN     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1     1   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN     1     1     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN     1     1     1     1   NaN   NaN   NaN     1     1     1   NaN
   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN     1     1     1     1     1   NaN   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN     1     1     1   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN
 ];
 
% Y-axis Max pointer
ymPointer = [
     2     2     2     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     1     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     1     2     2     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
     2     2     2   NaN   NaN     2     2   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN     1     1     1   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN     1     1     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN     1     1     1     1   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN   NaN   NaN     1     1     1   NaN   NaN   NaN     1     1   NaN     1     1
   NaN   NaN   NaN     1     1     1     1     1   NaN   NaN   NaN     1   NaN     1   NaN     1
   NaN   NaN   NaN     1     1     1   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN     1
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN     1
 ];

% zoom button icon
zoomIcon = [
     1     1     1     1     1     0     0     0     0     1     1     1     1     1     1     1
     1     1     1     0     0     0     0     0     0     0     0     1     1     1     1     1
     1     1     0     0     0     1     1     1     1     0     0     0     1     1     1     1
     1     0     0     1     1     1     0     0     1     1     1     0     0     1     1     1
     1     0     0     1     1     1     0     0     1     1     1     0     0     1     1     1
     0     0     1     1     1     1     0     0     1     1     1     1     0     0     1     1
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     1     1
     0     0     1     0     0     0     0     0     0     0     0     1     0     0     1     1
     0     0     1     1     1     1     0     0     1     1     1     1     0     0     1     1
     1     0     0     1     1     1     0     0     1     1     1     0     0     1     1     1
     1     0     0     1     1     1     0     0     1     1     1     0     0     1     1     1
     1     1     0     0     0     1     1     1     1     0     0     0     0     0     1     1
     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     1
     1     1     1     1     1     0     0     0     0     1     1     0     0     0     0     0
     1     1     1     1     1     1     1     1     1     1     1     1     0     0     0     0
     1     1     1     1     1     1     1     1     1     1     1     1     1     0     0     0
   ];  

% create zoom icon (RGB matrix)
fgID = zoomIcon==0;
bgID = ~fgID;
zI1  = zeros([size(zoomIcon), 3]);
zI2  = zeros([size(zoomIcon), 3]);
for id = 1:3
  tmp         = zoomIcon;
  tmp(fgID)   = 0;
  tmp(bgID)   = bgcolor3(id);
  zI1(:,:,id) = tmp;
  tmp(bgID)   = bgcolor4(id);
  zI2(:,:,id) = tmp;
end

% get screen size in pixels
un          = get(0, 'units');
set(0, 'units', 'pixels');
screenSize  = get(0, 'ScreenSize');
sW          = screenSize(3);
sH          = screenSize(4);
set(0, 'units', un);

% figure width and height (in pixels)
fW          = sW-200;
fH          = sH-100;

im = figure(...
  'units'                           , 'pixels', ...
  'position'                        , [100, 50, fW, fH], ...
  'backingstore'                    , 'off', ...
  'doublebuffer'                    , 'on', ...
  'name'                            , 'Grabit', ...
  'numbertitle'                     , 'off', ...
  'menubar'                         , 'none', ...
  'color'                           , bgcolor1, ...
  'pointer'                         , 'arrow', ...
  'visible'                         , 'off', ...
  'interruptible'                   , 'off', ...
  'busyaction'                      , 'cancel', ...
  'resizefcn'                       , @figResizeFcn, ...
  'windowbuttonupfcn'               , @winBtnUpFcn, ...
  'keypressfcn'                     , @keyPressFcn, ...
  'deletefcn'                       , 'delete(timerfind(''name'', ''BtnUpTimer''));', ...
  'tag'                             , 'GrabitGUI', ...
  'defaultUicontrolUnits'           , 'pixels', ...
  'defaultUicontrolBackgroundColor' , bgcolor3, ...
  'defaultUicontrolFontname'        , 'Verdana', ...
  'defaultUicontrolFontUnits'       , 'pixels', ...
  'defaultUicontrolFontsize'        , 10, ...
  'defaultUicontrolInterruptible'   , 'off', ...
  'defaultUicontrolBusyAction'      , 'cancel', ...
  'defaultAxesFontName'             , 'Verdana', ...
  'defaultAxesUnits'                , 'pixels', ...
  'defaultAxesFontSize'             , 8);

panelW = 0.8*fW;
uicontrol(...
  'style'                           , 'frame', ...
  'position'                        , [10, fH-110, panelW, 100], ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'pushbutton', ...
  'string'                          , 'Load Image...', ...
  'callback'                        , {@loadImageFcn, []}, ...
  'position'                        , [15, fH-40, 200, 25], ...
  'tag'                             , 'LoadImageBtn');
uicontrol(...
  'style'                           , 'edit', ...
  'backgroundcolor'                 , bgcolor2, ...
  'string'                          , '', ...
  'position'                        , [220, fH-40, panelW-215, 25], ...
  'horizontalalignment'             , 'left', ...
  'enable'                          , 'inactive', ...
  'tag'                             , 'ImageFileLoc');

uicontrol(...
  'style'                           , 'frame', ...
  'position'                        , [15, fH-75, 3*(panelW-20)/4-5, 30], ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'togglebutton', ...
  'string'                          , 'Calibrate', ...
  'buttondownfcn'                   , @calibrateImageFcn, ...
  'position'                        , [20, fH-72, panelW/4-25, 25], ...
  'enable'                          , 'off', ...
  'tag'                             , 'CalibrateImageBtn');
calibFrameX = panelW/4-5;
calibFrameW = panelW/2-10;
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5, fH-70, calibFrameW/8, 15], ...
  'string'                          , 'Xo:', ...
  'fontweight'                      , 'bold', ...
  'horizontalalignment'             , 'right', ...
  'tooltipstring'                   , sprintf('Calibration:\nX-Axis Origin'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+calibFrameW/8, fH-70, calibFrameW/8, 15], ...
  'string'                          , ' NaN', ...
  'horizontalalignment'             , 'left', ...
  'tag'                             , 'hXoValue', ...
  'tooltipstring'                   , sprintf('Calibration:\nX-Axis Origin'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+calibFrameW/4, fH-70, calibFrameW/8, 15], ...
  'string'                          , 'Xm:', ...
  'fontweight'                      , 'bold', ...
  'horizontalalignment'             , 'right', ...
  'tooltipstring'                   , sprintf('Calibration:\nX-Axis Max'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+3*calibFrameW/8, fH-70, calibFrameW/8, 15], ...
  'string'                          , ' NaN', ...
  'horizontalalignment'             , 'left', ...
  'tag'                             , 'hXmValue', ...
  'tooltipstring'                   , sprintf('Calibration:\nX-Axis Max'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+calibFrameW/2, fH-70, calibFrameW/8, 15], ...
  'string'                          , 'Yo:', ...
  'fontweight'                      , 'bold', ...
  'horizontalalignment'             , 'right', ...
  'tooltipstring'                   , sprintf('Calibration:\nY-Axis Origin'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+5*calibFrameW/8, fH-70, calibFrameW/8, 15], ...
  'string'                          , ' NaN', ...
  'horizontalalignment'             , 'left', ...
  'tag'                             , 'hYoValue', ...
  'tooltipstring'                   , sprintf('Calibration:\nY-Axis Origin'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+3*calibFrameW/4, fH-70, calibFrameW/8, 15], ...
  'string'                          , 'Ym:', ...
  'fontweight'                      , 'bold', ...
  'horizontalalignment'             , 'right', ...
  'tooltipstring'                   , sprintf('Calibration:\nY-Axis Max'), ...
  'backgroundcolor'                 , bgcolor2);
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [calibFrameX+5+7*calibFrameW/8, fH-70, calibFrameW/8, 15], ...
  'string'                          , ' NaN', ...
  'horizontalalignment'             , 'left', ...
  'tag'                             , 'hYmValue', ...
  'tooltipstring'                   , sprintf('Calibration:\nY-Axis Max'), ...
  'backgroundcolor'                 , bgcolor2);

uicontrol(...
  'style'                           , 'togglebutton', ...
  'string'                          , 'Grab Points', ...
  'buttondownfcn'                   , @grabPointsFcn, ...
  'position'                        , [3*(panelW)/4, fH-72, panelW/4+5, 25], ...
  'enable'                          , 'off', ...
  'tag'                             , 'GrabPointsBtn');
uicontrol(...
  'style'                           , 'togglebutton', ...
  'cdata'                           , zI1, ...
  'buttondownfcn'                   , @zoomBtnFcn, ...
  'position'                        , [15, fH-105, 25, 25], ...
  'enable'                          , 'inactive', ...
  'tag'                             , 'ZoomStateBtn');
uicontrol(...
  'style'                           , 'pushbutton', ...
  'string'                          , 'Reset View ', ...
  'callback'                        , @resetViewFcn, ...
  'position'                        , [45, fH-105, 100, 25], ...
  'tag'                             , 'ResetViewBtn');
uicontrol(...
  'style'                           , 'text', ...
  'position'                        , [150, fH-108, panelW-145, 31], ...
  'horizontalalignment'             , 'center', ...
  'foregroundcolor'                 , [0 0 .5], ...
  'backgroundcolor'                 , bgcolor2, ...
  'fontsize'                        , 10, ...
  'string'                          , {'Click and drag to pan. Double-click to center. Right-click and drag to zoom.', ...
    'Keyboard Shortcuts: <a> - zoom in, <z> - zoom out, <space> - reset view'});

rPanelX = 0.82 * fW;
rPanelW = fW - rPanelX - 10;
uicontrol(...
  'style'                           , 'listbox', ...
  'callback'                        , @selectVariableFcn, ...
  'position'                        , [rPanelX+10, 100, rPanelW-20, 0.6*fH], ...
  'backgroundcolor'                 , bgcolor2, ...
  'tooltipstring'                   , sprintf('Double-click to edit\nvariable in Array Editor'), ...
  'tag'                             , 'VariableList');
uicontrol(...
  'style'                           , 'text', ...
  'string'                          , 'Data in Memory', ...
  'position'                        , [rPanelX+10, 0.6*fH+100, rPanelW-20, 20], ...
  'backgroundcolor'                 , bgcolor1, ...
  'fontweight'                      , 'bold');
uicontrol(...
  'style'                           , 'pushbutton', ...
  'string'                          , 'Save to file as...', ...
  'position'                        , [rPanelX+10, 70, rPanelW-20, 25], ...
  'callback'                        , @variableManipulationFcn, ...
  'tooltipstring'                   , sprintf('Save variable as a MAT file or\nDouble-precision, tab-delimited TXT file'), ...
  'tag'                             , 'SaveAs');
uicontrol(...
  'style'                           , 'pushbutton', ...
  'string'                          , 'Rename...', ...
  'position'                        , [rPanelX+10, 45, rPanelW-20, 25], ...
  'callback'                        , @variableManipulationFcn, ...
  'tooltipstring'                   , 'Rename variable in workspace', ...
  'tag'                             , 'Rename');
uicontrol(...
  'style'                           , 'pushbutton', ...
  'string'                          , 'Delete...', ...
  'position'                        , [rPanelX+10, 20, rPanelW-20, 25], ...
  'callback'                        , @variableManipulationFcn, ...
  'tooltipstring'                   , 'Delete variable from workspace', ...
  'tag'                             , 'Delete');

axes(...
  'position'                        , [10, 10, 0.8*fW, fH-200], ...
  'visible'                         , 'on', ...
  'handlevisibility'                , 'on', ...
  'box'                             , 'on', ...
  'drawmode'                        , 'fast', ...
  'xtick'                           , [], ...
  'ytick'                           , [], ...
  'interruptible'                   , 'off', ...
  'buttondownfcn'                   , @winBtnDownFcn, ...
  'tag'                             , 'ImageAxis');
imAxRatio = (fH - 200) / (0.8 * fW);
title('', 'fontunits', 'pixels', 'fontsize', 24, 'color', 'red');

axes(...
  'position'                        , [rPanelX+20, 0.6*fH+150, rPanelW-30, 0.4*fH-190], ...
  'box'                             , 'on', ...
  'tag'                             , 'PreviewAxis');

uicontrol(...
  'style'                           , 'text', ...
  'string'                          , 'Preview Box', ...
  'fontweight'                      , 'bold', ...
  'position'                        , [rPanelX+20, fH-40, rPanelW-30, 25], ...
  'backgroundcolor'                 , bgcolor1);

uicontrol(...
  'style'                           , 'edit', ...
  'backgroundcolor'                 , bgcolor4, ...
  'string'                          , 'Enter Value', ...
  'position'                        , [0 0 100 25], ...
  'horizontalalignment'             , 'left', ...
  'visible'                         , 'off', ...
  'callback'                        , @coordinateEditFcn, ...
  'tag'                             , 'CoordinateEdit');

set(findobj(im, 'type', 'uicontrol'), 'units', 'normalized');
set(findobj(im, 'type', 'axes'), 'units', 'normalized');

% create handles structure
handles                           = guihandles(im);
handles.zoomPointer               = zoomPointer;
handles.zoomPointerHotSpot        = [2 2];
handles.zoomInOutPointer          = zoomInOutPointer;
handles.zoomInOutPointerHotSpot   = [9 9];
handles.closedHandPointer         = closedHandPointer;
handles.closedHandPointerHotSpot  = [9 9];
handles.xoPointer                 = xoPointer;
handles.xoPointerHotSpot          = [2 2];
handles.xmPointer                 = xmPointer;
handles.xmPointerHotSpot          = [2 2];
handles.yoPointer                 = yoPointer;
handles.yoPointerHotSpot          = [2 2];
handles.ymPointer                 = ymPointer;
handles.ymPointerHotSpot          = [2 2];
handles.zoomIconUp                = zI1;
handles.zoomIconDown              = zI2;
handles.curPointer                = 'arrow';
handles.curPointerData.CData      = zoomPointer;
handles.curPointerData.HotSpot    = [1 1];
handles.state                     = 'normal';
handles.bgcolor1                  = bgcolor1;
handles.bgcolor2                  = bgcolor2;
handles.bgcolor3                  = bgcolor3;
handles.bgcolor4                  = bgcolor4;
handles.imAxRatio                 = imAxRatio;
handles.I                         = [];
handles.map                       = [];
handles.savedVars                 = struct;
handles.isPanning                 = false;
handles.curTitle                  = '';
handles.CurrentPointAxes          = [];
handles.CurrentPointFig           = [];
handles.timer                     = timer(...
                                          'Name', 'BtnUpTimer', ...
                                          'StartDelay', 0.2, ...
                                          'TimerFcn', {@btnUpTimerFcn, im});

guidata(im, handles);

set(im, 'handlevisibility', 'callback');

if nargin == 1 && ischar(fname)
  switch lower(fname)
    case '-sample'
      fname = createSampleImageFcn;
      if ~isempty(fname)
        loadImageFcn(handles.LoadImageBtn, [], fname);
        try
          pause(0.5);
          delete(fname);
        catch
          errordlg(lasterr);
        end
      end
      
    otherwise
      filename = which(fname);
      if ~isempty(filename)
        loadImageFcn(handles.LoadImageBtn, [], filename);
      else
        errordlg(sprintf('%s\nnot found.', fname));
        return;
      end
  end
end

set(im, 'visible', 'on');


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% variableManipulationFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function variableManipulationFcn(varargin)
% This callback is called when one of the three buttons under the variable
% listbox is pressed.

obj = varargin{1};

handles = guidata(obj);

vars = fieldnames(handles.savedVars);
if isempty(vars)
  return;
end
listboxVal = get(handles.VariableList, 'value');
varName = vars{listboxVal};
switch get(obj, 'Tag')
  
  case 'SaveAs' % save the variable as a .mat file
    [fname, pname] = uiputfile(...
      {'*.mat', 'MAT files (*.mat)'; ...
      '*.txt', 'TXT files (*.txt)'}, ...
      sprintf('Save ''%s'' as:', varName), ...
      varName);
    if ~isequal(fname, 0) && ~isequal(pname, 0)
      saveDatFcn(pname, fname, handles.savedVars.(varName))
    end      
    
  case 'Rename'  % rename the variable in the base workspace
    answer = inputdlg({sprintf('Rename ''%s'' as:', varName)}, ...
      'Rename...', 1, {varName});
    if ~isempty(answer) && ~strcmp(varName, answer{1})
      if ~(evalin('base', sprintf('exist(''%s'', ''var'')', answer{1}))) && ...
          isempty(strmatch(answer{1}, vars, 'exact'))
        newVarNames = vars;
        newVarNames{listboxVal} = answer{1};
        for id = 1:length(vars)
          tmp.(newVarNames{id}) = handles.savedVars.(vars{id});
        end
        handles.savedVars = tmp;
        showAllVarsFcn(handles);
        set(handles.VariableList, 'value', listboxVal);
        assignin('base', answer{1}, handles.savedVars.(answer{1}));
        evalin('base', sprintf('clear %s;', varName));
      else
        errordlg('That name is already in use.');
        return;
      end
    end
    
  case 'Delete' % delete the variable from the base workspace
    btn = questdlg(sprintf('Delete ''%s'' from the workspace?', varName), ...
      'Delete from workspace...', 'Yes', 'No', 'No');
    switch btn
      case 'Yes'
        handles.savedVars = rmfield(handles.savedVars, varName);
        showAllVarsFcn(handles);
        evalin('base', sprintf('clear %s;', varName));
      case 'No'
    end
    
end

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% saveDatFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function saveDatFcn(pname, fname, var)
% this function saves the variable to file

[p, fname, ext] = fileparts(fname);
switch lower(ext)
  case '.mat'
    eval(sprintf('%s = var;', fname));
    save(fullfile(pname, [fname, ext]), fname, '-v6');
  case '.txt'
    eval(sprintf('%s = var;', fname));
    save(fullfile(pname, [fname, ext]), fname, '-ascii', '-double', '-tabs');
  otherwise
    errordlg('Unknown extension.');
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% selectVariableFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function selectVariableFcn(varargin)
% this callback is called when a variable is selected from the variable
% listbox.  When the variable is clicked, it will plot the data in the
% Preview Box above.  If the variable is double-clicked, it will be opened
% in the Array Editor.  It also ensures that the variable in the base
% workspace is the same copy as the variable stored in the GRABIT
% workspace.  This means that if you change the contents of a variable in
% the base workspace (via other functions), then it will allow you to
% update the variable in the GRABIT workspace.

obj = varargin{1};

handles = guidata(obj);

sType = get(handles.GrabitGUI, 'SelectionType');
switch sType
  case {'normal', 'open'} % single or double click
    vars = fieldnames(handles.savedVars);
    if isempty(vars)
      return;
    end
    listVal = get(obj, 'value');
    varName = vars{listVal};
    
    % check to see if the stored var is the same as the var in the
    % base workspace
    if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
      
      % the copy in the base workspace must be a n-by-2 DOUBLE array
      if strcmp(class(evalin('base', varName)), 'double') && ...
          ndims(evalin('base', varName)) == 2 && ...
          size(evalin('base', varName), 2) == 2
        
        if ~isequal(evalin('base', varName), handles.savedVars.(varName))
          % may have been modified in the base workspace
          btn = questdlg(sprintf('''%s'' may have been modified in the base workspace (e.g. Array Editor).\nUpdate the variable?', varName), ...
            'Modified Variable', 'Update Base Workspace', 'Update Grabit', 'Neither', 'Update Grabit');
        
          switch btn
            case 'Update Base Workspace'
              assignin('base', varName, handles.savedVars.(varName));
              
            case 'Update Grabit'
              handles.savedVars.(varName) = evalin('base', varName);
              showAllVarsFcn(handles);
              set(obj, 'value', listVal);
              
          end
          
        end
        
      else
        btn = questdlg(sprintf('''%s'' in the base workspace is different from the one in Grabit.\nUpdate the base workspace variable?', varName), ...
          'Wrong Variable', 'Update Base Workspace', 'Leave Untouched', 'Update Base Workspace');
        
        switch btn
          case 'Update Base Workspace'
            assignin('base', varName, handles.savedVars.(varName));
            
        end
        
      end
      
    else % the variable does not exist in base workspace, so make a copy
      assignin('base', varName, handles.savedVars.(varName));
    end
    
    switch sType
      case 'normal' % single click
        axes(handles.PreviewAxis);
        set(handles.PreviewLine, ...
          'xdata', handles.savedVars.(varName)(:, 1), ...
          'ydata', handles.savedVars.(varName)(:, 2));
        axis auto;
        
      case 'open' % double click
        try
          openvar(varName);
        catch
          errordlg(lasterr);
        end
        
    end
end

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% loadImageFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function loadImageFcn(varargin)
% this function loads an image file

[obj, filename] = splitvar(varargin([1, 3]));

handles = guidata(obj);

if isempty(filename)
  [fname, pathname] = uigetfile(...
    {'*.bmp;*.jpg;*.jpeg;*.tif;*.tiff;*.gif;*.png', ...
    'Image Files (*.bmp,  *.jpg,  *.jpeg,  *.tif,  *.tiff,  *.gif,  *.png)';
    '*.bpm', 'BPM files (*.bpm)';
    '*.jpg;*jpeg', 'JPG files (*.jpg,  *.jpeg)';
    '*.tif;*tiff', 'TIFF files (*.tif,  *.tiff)';
    '*.gif', 'GIF files (*.gif)';
    '*.png', 'PNG files (*.png)';
    '*.*', 'All files (*.*)'}, 'Select an image file');
  if ischar(fname)
    filename = fullfile(pathname, fname);
  else
    return;
  end
end

set(handles.ImageFileLoc, 'string', filename);

try
  %warning off;
  [A, map] = imread(filename);
  %warning on;
catch
  errordlg(lasterr);
  return;
end

if ndims(A) == 3   %some TIFF files have wrong size matrices.
  if size(A, 3)>3
    A = A(:, :, 1:3);
  elseif size(A, 3)<3
    errordlg('This is a weird image file...possibly a bad TIFF file...');
    return;
  end
end

% Need this so that image shows up when not called by a CALLBACK
set(0, 'ShowHiddenHandles', 'on');

cla(handles.PreviewAxis);
handles.PreviewLine = line(NaN, NaN, ...
  'color'     , 'blue', ...
  'linestyle' , '-', ...
  'marker'    , '.', ...
  'tag'       , 'PreviewLine', ...
  'parent'    , handles.PreviewAxis);
axes(handles.ImageAxis);
cla(handles.ImageAxis);
NP = get(handles.GrabitGUI, 'nextplot');  % for compatibility with R14

if isempty(map)
  imageInfo = imfinfo(filename);
  if strcmpi(imageInfo.ColorType, 'grayscale')
    colormap(gray(2^imageInfo.BitDepth));
  end
else
  colormap(map);
end
iH = image(A);
set(iH, 'HitTest', 'off', 'EraseMode', 'normal');
set(handles.ImageAxis, 'xtick', [], 'ytick', []);
axis equal;

set(handles.GrabitGUI, 'nextplot', NP);
set(handles.ImageAxis, ...
  'drawmode'        , 'fast', ...
  'tag'             , 'ImageAxis', ...
  'handlevisibility', 'callback', ...
  'buttondownfcn'   , @winBtnDownFcn);
set(get(handles.ImageAxis, 'title'), ...
  'string'          , '', ...
  'fontunits'       , 'pixels', ...
  'fontsize'        , 24, ...
  'color'           , 'red');
handles.ImageLine = line(NaN, NaN, ...
  'color'           , 'red', ...
  'linestyle'       , 'none', ...
  'marker'          , '.', ...
  'tag'             , 'ImageLine', ...
  'hittest'         , 'off');
handles.CalibPtsH(:,1) = line(repmat(NaN, 2, 4), repmat(NaN, 2, 4));
handles.CalibPtsH(:,2) = line(repmat(NaN, 2, 4), repmat(NaN, 2, 4));
set(handles.CalibPtsH(:,1), ...
  {'marker','color'}, {'o','r';'o','r';'o','b';'o','b'}, ...
  'markersize'      , 10, ...
  'hittest'         , 'off');
set(handles.CalibPtsH(:,2), ...
  {'marker','color'}, {'+','r';'x','r';'+','b';'x','b'}, ...
  'markersize'      , 20, ...
  'hittest'         , 'off');

% initialize image data
handles.I           = A;
handles.map         = map;
handles.CalibVals   = [];
handles.CalibPts    = [NaN, NaN, NaN, NaN];
handles.CalibPtsIm  = repmat(NaN, 4, 2);
handles.ImLimits    = [get(handles.ImageAxis, 'xlim'); ...
                       get(handles.ImageAxis, 'ylim')];

set(handles.CalibrateImageBtn, ...
  'string', 'Calibrate', ...
  'enable', 'inactive', ...
  'value' , 0);
set(handles.GrabPointsBtn, ...
  'enable', 'off', ...
  'value' , 0);
set(handles.ZoomStateBtn, ...
  'enable', 'inactive', ...
  'value' , 0);
zoom off;

set(0, 'ShowHiddenHandles', 'off');

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% calibrateImageFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function calibrateImageFcn(varargin)
% this function performs calibration of the image by prompting the user to
% select 4 points on the image as reference points.

obj = varargin{1};

handles = guidata(obj);

if isempty(handles.I)
  set(obj, 'enable', 'off');
else
  switch get(obj, 'value')
    case 0  % start calibration
      set(obj, 'value', 1, 'backgroundcolor', handles.bgcolor4);
            
      handles.curPointer              = 'custom';
      handles.curPointerData.CData    = handles.xoPointer;
      handles.curPointerData.HotSpot  = handles.xoPointerHotSpot;
      set(handles.GrabitGUI, ...
        'PointerShapeCData'     , handles.curPointerData.CData, ...
        'PointerShapeHotSpot'   , handles.curPointerData.HotSpot, ...
        'WindowButtonMotionFcn' , {@pointerFcn, handles, handles.curPointer});
      
      set([handles.LoadImageBtn, ...
          handles.SaveAs, ...
          handles.Rename, ...
          handles.Delete], ...
        'enable', 'off');
      
      set(handles.CalibPtsH, 'xdata', NaN, 'ydata', NaN);
      handles.CalibVals = [];
      handles.CalibPts = [NaN, NaN, NaN, NaN];
      handles.CalibPtsIm = repmat(NaN, 4, 2);
      set([handles.hXoValue, handles.hXmValue, ...
          handles.hYoValue, handles.hYmValue], ...
        'string', ' NaN');
      
      % this is the first calibration point: X-Axis Origin
      handles.curTitle = 'Click on the ORIGIN of x-axis';
      set(get(handles.ImageAxis, 'title'), ...
        'string', handles.curTitle);
      
      % change state to CALIBRATION
      handles.state = 'calibration';

    case 1  % stop (prematurely) calibration
      set(obj, ...
        'value'           , 0, ...
        'backgroundcolor' , handles.bgcolor3, ...
        'string'          , 'Calibrate');
      handles.curTitle = '';
      set(get(handles.ImageAxis, 'title'), 'string', '');
            
      handles.curPointer          = 'arrow';
      set(handles.GrabitGUI, ...
        'WindowButtonMotionFcn' , {@pointerFcn, handles, handles.curPointer});
           
      set([handles.LoadImageBtn, ...
          handles.SaveAs, ...
          handles.Rename, ...
          handles.Delete], ...
        'enable', 'on');
      
      % calibration was prematurely stopped, so reset all values
      set(handles.CalibPtsH, 'xdata', NaN, 'ydata', NaN);
      handles.CalibVals = [];
      handles.CalibPts = [NaN, NaN, NaN, NaN];
      handles.CalibPtsIm = repmat(NaN, 4, 2);
      set([handles.hXoValue, handles.hXmValue, ...
          handles.hYoValue, handles.hYmValue], ...
        'string', ' NaN');
      
      set(handles.GrabPointsBtn, 'enable', 'off');
      set(handles.CoordinateEdit, 'visible', 'off');
      
      % change state to NORMAL
      handles.state = 'normal';
  end
end

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% grabPointsFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function grabPointsFcn(varargin)
% this function is used to extract data points by prompting the user to
% select points on the image.

obj = varargin{1};

handles = guidata(obj);

switch get(obj, 'value')
  case 0  % initiate point grabbing
    calib = handles.CalibVals;
    
    axes(handles.ImageAxis);
    set(handles.ImageAxis, ...
      'xlim', handles.ImLimits(1,:), ...
      'ylim', handles.ImLimits(2,:));
    set(handles.ImageLine, 'xdata', NaN, 'ydata', NaN);
    handles.curTitle = {'Grab points by clicking on data points.', ...
    '<BACKSPACE> or <DEL> to delete previous point. <ENTER> to finish.'};
    set(get(handles.ImageAxis, 'title'), ...
      'string', handles.curTitle);
    
    set(handles.PreviewLine, 'xdata', NaN, 'ydata', NaN);
    set(handles.PreviewAxis, ...
      'xlim', [min([calib.Xo calib.Xm]) max([calib.Xo calib.Xm])], ...
      'ylim', [min([calib.Yo calib.Ym]) max([calib.Yo calib.Ym])]);
    
    handles.ImDat   = [];
    handles.TrueDat = [];
    
    handles.curPointer          = 'crosshair';
    set(handles.GrabitGUI, ...
      'WindowButtonMotionFcn', {@pointerFcn, handles, handles.curPointer});
        
    set(obj, ...
      'value'           , 1, ...
      'string'          , 'Grabbing Points (0)', ...
      'backgroundcolor' , handles.bgcolor4);
    set([handles.CalibrateImageBtn, ...
        handles.LoadImageBtn, ...
        handles.SaveAs, ...
        handles.Rename, ...
        handles.Delete], ...
      'enable', 'off');

    % change state to GRAB
    handles.state = 'grab';
    
  case 1  % finish point grabbing
    set(obj, ...
      'value'           , 0, ...
      'string'          , 'Grab Points', ...
      'backgroundcolor' , handles.bgcolor1);

    handles.curTitle = '';
    set(get(handles.ImageAxis, 'title'), 'string', '');
    handles.curPointer            = 'arrow';
    set(handles.GrabitGUI, 'WindowButtonMotionFcn', {@pointerFcn, handles, handles.curPointer});
    set(handles.CalibrateImageBtn, 'enable', 'inactive');
    set([handles.LoadImageBtn, ...
        handles.SaveAs, ...
        handles.Rename, ...
        handles.Delete], ...
      'enable', 'on');
    if ~isempty(handles.TrueDat) % some points were grabbed
      varNames                          = fieldnames(handles.savedVars);
      varNames{end + 1}                 = findNextVarNameFcn(varNames);
      handles.savedVars.(varNames{end}) = handles.TrueDat;
      assignin('base', varNames{end}, handles.savedVars.(varNames{end}));
      showAllVarsFcn(handles);
    end
    
    % change to NORMAL state
    handles.state = 'normal';
end  

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% showAllVarsFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function showAllVarsFcn(handles)
% this function shows all data set variables that exist in the GRABIT
% workspace in the variable listbox.

varNames = fieldnames(handles.savedVars);

if isempty(varNames)
  listboxStr = {''};
else
  for id = 1:length(varNames)
    [m, n] = size(handles.savedVars.(varNames{id}));
    listboxStr{id} = sprintf('%s [%dx%d]', varNames{id}, m, n);
  end
end

set(handles.VariableList, 'string', listboxStr, 'value', length(listboxStr));


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% findNextVarNameFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function newVarName = findNextVarNameFcn(varNames)
% this helper function determines the next available variable name by
% checking the existing variable names in the base workspace and GRABIT
% workspace.

wsVarNames = evalin('base', 'who');
vars = unique([wsVarNames(:); varNames(:)]);
idx = 1;
while 1
  if isempty(strmatch(sprintf('Data%03d', idx), vars, 'exact'))
    newVarName = sprintf('Data%03d', idx);
    return;
  else
    idx = idx + 1;
  end
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% zoomBtnFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function zoomBtnFcn(varargin)
% this function toggles the zoom state

obj = varargin{1};

handles = guidata(obj);

switch get(obj, 'value')
  case 0
    set(obj, ...
      'value'          , 1, ...
      'backgroundcolor', handles.bgcolor4, ...
      'cdata'          , handles.zoomIconDown);
    udata.titlestring = get(get(handles.ImageAxis, 'Title'), 'string');
    udata.btnstate = get([handles.LoadImageBtn, ...
        handles.CalibrateImageBtn, ...
        handles.GrabPointsBtn, ...
        handles.VariableList, ...
        handles.SaveAs, ...
        handles.Rename, ...
        handles.Delete],'enable');
    udata.imgstate = get(handles.ImageAxis, 'ButtonDownFcn');
    udata.handlevisibility = get(handles.ImageAxis, 'handlevisibility');
    set(obj, 'userdata', udata);
    
    set(get(handles.ImageAxis, 'Title'), 'string', 'Zoom ON');
    set([handles.LoadImageBtn, ...
        handles.CalibrateImageBtn, ...
        handles.GrabPointsBtn, ...
        handles.VariableList, ...
        handles.SaveAs, ...
        handles.Rename, ...
        handles.Delete], 'enable', 'off');
    set(handles.ImageAxis, 'ButtonDownFcn', '');
    zoom('on');
    handles.curPointerData.CData   = get(handles.GrabitGUI, 'PointerShapeCData');
    handles.curPointerData.HotSpot = get(handles.GrabitGUI, 'PointerShapeHotSpot');
    set(handles.GrabitGUI, ...
      'PointerShapeCData'     , handles.zoomPointer, ...
      'PointerShapeHotSpot'   , handles.zoomPointerHotSpot, ...
      'WindowButtonMotionFcn' , {@pointerFcn, handles, 'custom'}, ...
      'keypressfcn'           , ';'); % this prevents switching to command window

    % this seems necessary in some versions of Matlab
    set(handles.ImageAxis, 'handlevisibility', 'on');
  
  case 1
    set(obj, ...
      'value', 0, ...
      'backgroundcolor', handles.bgcolor3, ...
      'cdata', handles.zoomIconUp);
    zoom('off');
    
    %----------------------------------------------------------------------
    % If zoom created a compact view window, expand it to fill the whole
    % axes.
    %----------------------------------------------------------------------
    xl = get(handles.ImageAxis, 'xlim'); xrng = diff(xl);
    yl = get(handles.ImageAxis, 'ylim'); yrng = diff(yl);
    if abs(yrng/xrng - handles.imAxRatio) > .01 % wrong axes ratio
      if yrng/xrng > handles.imAxRatio
        xrng = yrng / handles.imAxRatio;
        xl   = mean(xl) + [-0.5, 0.5] * xrng;
      else
        yrng = xrng * handles.imAxRatio;
        yl   = mean(yl) + [-0.5, 0.5] * yrng;
      end
      set(handles.ImageAxis, 'xlim', xl, 'ylim', yl);
    end
    %----------------------------------------------------------------------
    
    udata = get(obj, 'userdata');
    set(get(handles.ImageAxis, 'Title'), 'string', udata.titlestring);
    set(handles.ImageAxis, 'handlevisibility', udata.handlevisibility);
    set([handles.LoadImageBtn, ...
        handles.CalibrateImageBtn, ...
        handles.GrabPointsBtn, ...
        handles.VariableList, ...
        handles.SaveAs, ...
        handles.Rename, ...
        handles.Delete], {'enable'}, udata.btnstate);
    set(handles.ImageAxis, 'ButtonDownFcn', udata.imgstate);
    set(handles.GrabitGUI, ...
      'PointerShapeCData'     , handles.curPointerData.CData, ...
      'PointerShapeHotSpot'   , handles.curPointerData.HotSpot, ...
      'WindowButtonMotionFcn' , {@pointerFcn, handles, handles.curPointer}, ...
      'keypressfcn'           , @keyPressFcn);
end

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% pointerFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function pointerFcn(varargin)
% this changes the pointer based on whether the cursor is in the image
% axes.

[handles, ptr] = splitvar(varargin(3:4));

pt = get(handles.ImageAxis, 'CurrentPoint');
xl = get(handles.ImageAxis, 'xlim');
yl = get(handles.ImageAxis, 'ylim');
if pt(1,1) > xl(1) && pt(1,1) < xl(2) && pt(1,2) > yl(1) && pt(1,2) < yl(2)
  set(handles.GrabitGUI, 'pointer', ptr);
else
  set(handles.GrabitGUI, 'pointer', 'arrow');
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% figResize
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function figResizeFcn(varargin)
% this function makes sure the axis fills the whole axes extent

obj = varargin{1};

handles = guidata(obj);

axis(handles.ImageAxis, 'equal');
handles.imAxRatio = diff(get(handles.ImageAxis, 'ylim')) / ...
  diff(get(handles.ImageAxis, 'xlim'));
handles.ImLimits = [get(handles.ImageAxis, 'xlim'); ...
    get(handles.ImageAxis, 'ylim')];

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% keyPressFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function keyPressFcn(varargin)
% this is for the keyboard shortcuts. During 'grab' mode, <backspace>
% deletes the last point clicked and <return> ends the 'grab' mode.

obj = varargin{1};

handles = guidata(obj);

if ~isempty(handles.I)
  k = lower(get(obj, 'CurrentKey'));
  
  switch k
    case 'a'         % zoom in
      xl = get(handles.ImageAxis, 'xlim'); xrng = diff(xl);
      yl = get(handles.ImageAxis, 'ylim'); yrng = diff(yl);
      
      % prevent zooming in too much.
      % set the limit to 64x zoom.
      if xrng >= size(handles.I, 2)/64*2
        % animate zoom
        for id = 0:0.2:1
          set(handles.ImageAxis, ...
            'xlim', xl + id * xrng / 4 * [1, -1], ...
            'ylim', yl + id * yrng / 4 * [1, -1]);
          drawnow;
        end
      end
      
    case 'z'        % zoom out
      xl = get(handles.ImageAxis, 'xlim'); xrng = diff(xl);
      yl = get(handles.ImageAxis, 'ylim'); yrng = diff(yl);
      % animate zoom
      for id = 0:0.2:1
        set(handles.ImageAxis, ...
          'xlim', xl + id * xrng / 2 * [-1, 1], ...
          'ylim', yl + id * yrng / 2 * [-1, 1]);
        drawnow;
      end
      
    case 'space'    % reset view
      resetViewFcn(handles.ResetViewBtn);
      
    case {'backspace', 'delete'}
      
      switch handles.state
        case 'grab'
          if ~handles.isPanning
            
            if isempty(handles.ImDat)
              return;
            else
              handles.ImDat(end, :)   = [];
              handles.TrueDat(end, :) = [];
            end
            
            set(handles.PreviewLine, ...
              'xdata', handles.TrueDat(:, 1), ...
              'ydata', handles.TrueDat(:, 2));
            set(handles.ImageLine, ...
              'xdata', handles.ImDat(:, 1), ...
              'ydata', handles.ImDat(:, 2));
            
            set(handles.GrabPointsBtn, 'string', sprintf('Grabbing Points (%d)', size(handles.ImDat, 1)));
            
            guidata(obj, handles);
          end
          
      end
      
    case {'return', 'enter'}
      
      switch handles.state
        case 'grab'
          grabPointsFcn(handles.GrabPointsBtn);
          
      end

  end
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% resetViewFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function resetViewFcn(varargin)
% this function resets the view

obj = varargin{1};

handles = guidata(obj);

if ~isempty(handles.I)
  xl = get(handles.ImageAxis, 'xlim');
  yl = get(handles.ImageAxis, 'ylim');
  xd = (handles.ImLimits(1, :) - xl) / 10;
  yd = (handles.ImLimits(2, :) - yl) / 10;
  % animate zoom
  for id = 0:10
    set(handles.ImageAxis, ...
      'xlim', xl + id * xd, ...
      'ylim', yl + id * yd);
    drawnow;
  end
end

% take focus away
loseFocusFcn(handles)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% loseFocusFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function loseFocusFcn(handles)
% attempt to take focus away by setting the ENABLE property to off and then
% back to the original setting

settings = get([handles.LoadImageBtn, ...
    handles.ResetViewBtn, ...
    handles.VariableList, ...
    handles.SaveAs, ...
    handles.Rename, ...
    handles.Delete], ...
  'enable');
set([handles.LoadImageBtn, ...
    handles.ResetViewBtn, ...
    handles.VariableList, ...
    handles.SaveAs, ...
    handles.Rename, ...
    handles.Delete], ...
  'enable', 'off');
drawnow;
set([handles.LoadImageBtn, ...
    handles.ResetViewBtn, ...
    handles.VariableList, ...
    handles.SaveAs, ...
    handles.Rename, ...
    handles.Delete], ...
 {'enable'}, settings);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% winBtnDownFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function winBtnDownFcn(varargin)
% this function is called when the mouse click initiates

obj = varargin{1};

handles = guidata(obj);

if strcmpi(get(handles.timer, 'Running'), 'on') || isempty(handles.I)
  return;
end

handles.CurrentPointAxes = get(handles.ImageAxis, 'CurrentPoint');
handles.CurrentPointFig  = get(handles.GrabitGUI, 'CurrentPoint');

switch get(handles.GrabitGUI, 'SelectionType')
  case 'normal'
    
    set(handles.CoordinateEdit, 'visible', 'off');
    id = find(isnan(handles.CalibPts));
    if ~isempty(id)
      set(handles.CalibPtsH(id(1),:), 'xdata', NaN, 'ydata', NaN);
    end
    
    % first call winBtnMotionPauseFcn to prevent immediate click-n-drag
    set(handles.GrabitGUI, ...
      'WindowButtonMotionFcn', ...
      {@winBtnMotionPauseFcn, handles, handles.CurrentPointAxes(1,1:2), clock});
    
  case 'alt'
    
    set(handles.CoordinateEdit, 'visible', 'off');
    id = find(isnan(handles.CalibPts));
    if ~isempty(id)
      set(handles.CalibPtsH(id(1),:), 'xdata', NaN, 'ydata', NaN);
    end
    
    xl = get(handles.ImageAxis, 'XLim');midX = (xl(1)+xl(2))/2;
    yl = get(handles.ImageAxis, 'YLim');midY = (yl(1)+yl(2))/2;
    figPos = get(handles.GrabitGUI, 'Position');
    handles.curTitle = get(get(handles.ImageAxis, 'title'), 'string');
    set(handles.GrabitGUI, ...
      'Pointer', 'custom', ...
      'PointerShapeCData'     , handles.zoomInOutPointer, ...
      'PointerShapeHotSpot'   , handles.zoomInOutPointerHotSpot, ...
      'WindowButtonMotionFcn' , {@zoomMotionFcn, handles, ...
                                 get(handles.GrabitGUI, 'CurrentPoint'), ...
                                 figPos(4), ...
                                 size(handles.I, 2), ...
                                 midX, ...
                                 midY, ...
                                 diff(xl)/2, ...
                                 diff(yl)/2});      
    
    set(get(handles.ImageAxis, 'title'), 'string', 'Zooming...');
                             
end

guidata(obj, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% zoomMotionFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function zoomMotionFcn(varargin)
% this performs the click-n-drag zooming function. The pointer location
% relative to the initial point determines the amount of zoom (in or out).

[obj, handles, initPt, figHt, horizPx, midX, ...
  midY, rngXhalf, rngYhalf] = splitvar(varargin([1, 3:end]));

pt = get(obj, 'CurrentPoint');

% get relative pointer location (y-coord only).
% power law allows for the inverse to work:
%   C^(x) * C^(-x) = 1
% Choose C to get "appropriate" zoom factor
r = 30 ^ ((initPt(2) - pt(2)) / figHt);

% make sure it doesn't zoom in too much.
% the limit is based on the size of the original image.
% set limit to 64x zoom.
if r < horizPx/64/rngXhalf/2  % stop zoom
  set(get(handles.ImageAxis, 'title'), 'string', 'Max Zoom Reached');
  set(obj, ...
    'Pointer'               , 'arrow', ...
    'WindowButtonMotionFcn' , '');
else
  set(handles.ImageAxis, ...
    'XLim', [midX - r * rngXhalf, midX + r * rngXhalf], ...
    'YLim', [midY - r * rngYhalf, midY + r * rngYhalf]);
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% winBtnMotionPauseFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function winBtnMotionPauseFcn(varargin)
% this prevents click-n-drag from happening for X seconds. This is useful
% because users may move the mouse as they are clicking.

[obj, handles, xy, c] = splitvar(varargin([1, 3:end]));

if etime(clock, c) > .15  % waits .15 seconds before dragging occurs
  set(obj, ...
    'Pointer'             , 'custom', ...
    'PointerShapeCData'   , handles.closedHandPointer, ...
    'PointerShapeHotSpot' , handles.closedHandPointerHotSpot);
  set(obj, 'WindowButtonMotionFcn', {@winBtnMotionFcn, handles.ImageAxis, xy});
  handles.curTitle = get(get(handles.ImageAxis, 'title'), 'string');
  set(get(handles.ImageAxis, 'title'), 'string', 'Panning...');
  handles.isPanning = true;
  guidata(obj, handles);
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% winBtnMotionFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function winBtnMotionFcn(varargin)
% this function is called when click-n-drag (panning) is happening

[axH, xy] = splitvar(varargin(3:4));

pt = get(axH, 'CurrentPoint');
set(axH, ...
  'xlim', get(axH, 'xlim')+(xy(1)-pt(1,1)), ...
  'ylim', get(axH, 'ylim')+(xy(2)-pt(1,2)));


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% winBtnUpFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function winBtnUpFcn(varargin)
% this is called when the mouse button is released. It initiates the button
% down timer object (see btnUpTimerFcn). This timer waits for 0.2 seconds
% before any action is taken. This is in order to detect a double click. If
% a double click is detected, the single click action is NOT executed.

obj = varargin{1};

handles = guidata(obj);
if ~isempty(handles.I)  % there is an image displayed
  switch get(obj, 'SelectionType')
    case 'normal'
      if strcmpi(get(handles.timer, 'Running'), 'off')
        
        % start the timer which waits some time to see if double-clicking
        % occurs
        start(handles.timer);
        
        set(obj, ...
          'Pointer'               , handles.curPointer, ...
          'PointerShapeCData'     , handles.curPointerData.CData, ...
          'PointerShapeHotSpot'   , handles.curPointerData.HotSpot, ...
          'WindowButtonMotionFcn' , {@pointerFcn, handles, handles.curPointer});
        
        set(get(handles.ImageAxis, 'title'), 'string', handles.curTitle);
      end
      
    case 'alt'
      set(obj, ...
        'Pointer'               , handles.curPointer, ...
        'PointerShapeCData'     , handles.curPointerData.CData, ...
        'PointerShapeHotSpot'   , handles.curPointerData.HotSpot, ...
        'WindowButtonMotionFcn' , {@pointerFcn, handles, handles.curPointer});
      
      set(get(handles.ImageAxis, 'title'), 'string', handles.curTitle);
  end
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% btnUpTimerFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function btnUpTimerFcn(varargin)

figH = varargin{3};

handles = guidata(figH);

switch get(handles.GrabitGUI, 'SelectionType')
  case 'open' % double click
    % this will center the view
    
    % get current units
    un = get([0, handles.GrabitGUI, handles.ImageAxis], 'units');
    
    set([0, handles.GrabitGUI, handles.ImageAxis], 'units', 'pixels');
    pt     = get(0, 'PointerLocation');
    figPos = get(handles.GrabitGUI, 'Position');
    pt2    = pt - figPos(1:2);
    axPos  = get(handles.ImageAxis, 'position');
    
    % check to see if pointer is inside the image axes
    if pt2(1) > axPos(1) && pt2(1) < axPos(1)+axPos(3) && ...
        pt2(2) > axPos(2) && pt2(2) < axPos(2)+axPos(4)
      xl = get(handles.ImageAxis, 'xlim'); xrng = diff(xl);
      yl = get(handles.ImageAxis, 'ylim'); yrng = diff(yl);
      x  = (pt2(1)-axPos(1))/axPos(3)*xrng+xl(1);
      y  = (axPos(2)+axPos(4)-pt2(2))/axPos(4)*yrng+yl(1);
      
      % animate (slide) to the new location. this may give a better
      % visual perception of the view change
      
      % determine how fast to animate, depending on the location of the
      % pointer
      interval = ceil(norm((figPos(1:2)+axPos(1:2)+axPos(3:4)/2) - pt)/30);
      
      if interval
        ld = (([x, y] - [xrng, yrng]/2) - [xl(1), yl(1)])/interval;
        pd = ((figPos(1:2)+axPos(1:2)+axPos(3:4)/2) - pt)/interval;
      else
        % interval == 0, meaning it's the same point
        ld = [0, 0];
        pd = [0, 0];
      end
      
      for id = 0:interval
        set(handles.ImageAxis, ...
          'xlim', xl + id * ld(1), ...
          'ylim', yl + id * ld(2));
        
        % center the pointer location
        set(0, ...
          'PointerLocation', pt + id * pd);
        drawnow;
      end
    end
    
    % reset UNITS
    set([0, handles.GrabitGUI, handles.ImageAxis], {'units'}, un);
    
  case 'normal'  % single click
    
    switch handles.state
      case 'grab'
        if ~handles.isPanning && ~isempty(handles.CurrentPointAxes)
          calib = handles.CalibVals;
          
          X = handles.CurrentPointAxes(1, 1);
          Y = handles.CurrentPointAxes(1, 2);
          
          % the point on the image
          handles.ImDat(end + 1, :) = [X, Y];
          
          % the true data point coordinates
          handles.TrueDat(end + 1, :) = ([calib.Xo; calib.Yo] + ...
            [calib.e1 calib.e2] \ [X - calib.Xxo; Y - calib.Yyo])';
          set(handles.PreviewLine, ...
            'xdata', handles.TrueDat(:, 1), ...
            'ydata', handles.TrueDat(:, 2));
          set(handles.ImageLine, ...
            'xdata', handles.ImDat(:, 1), ...
            'ydata', handles.ImDat(:, 2));
          
          set(handles.GrabPointsBtn, 'string', sprintf('Grabbing Points (%d)', size(handles.ImDat, 1)));
          
        end
        
      case 'calibration'
        if ~handles.isPanning && ~isempty(handles.CurrentPointAxes)
          id = find(isnan(handles.CalibPts));
          pt = handles.CurrentPointAxes;
          set(handles.CalibPtsH(id(1),:), 'xdata', pt(1,1), 'ydata', pt(1,2));
          handles.CalibPtsIm(id(1),:) = pt(1,1:2);
          set(handles.CoordinateEdit, ...
            'units'   , 'pixels', ...
            'position', [handles.CurrentPointFig+[5, 5], 100, 25], ...
            'visible' , 'on', ...
            'string'  , 'Enter Value');
        end
        
    end
end

handles.CurrentPointAxes = [];
handles.CurrentPointFig  = [];
handles.isPanning        = false;

guidata(handles.GrabitGUI, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% coordinateEditFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function coordinateEditFcn(varargin)
% this function is called when a value is entered in the edit box for
% calibration values

obj = varargin{1};

handles = guidata(obj);

strs = {'Click on the MAXIMUM value of x-axis', ...
    'Click on the ORIGIN of y-axis', ...
    'Click on the MAXIMUM value of y-axis', ...
    ''};
ptrs = {'xmPointer', 'yoPointer', 'ymPointer'};
hs = {'hXoValue', 'hXmValue', 'hYoValue', 'hYmValue'};

val = str2double(get(obj, 'string'));
if ~isnan(val)
  id = find(isnan(handles.CalibPts));
  if ~isempty(id)
    handles.CalibPts(id(1)) = val;
    set(handles.(hs{id(1)}), 'string', sprintf(' %g', val));
    set(obj, 'visible', 'off');
    handles.curTitle = strs{id(1)};
    set(get(handles.ImageAxis, 'title'), 'string', strs{id(1)});
    if id(1) == 4 % last calibration point
      handles = calculateCalibrationFcn(handles);
    else
      handles.curPointerData.CData    = handles.(ptrs{id(1)});
      handles.curPointerData.HotSpot  = handles.([ptrs{id(1)} 'HotSpot']);
      set(handles.GrabitGUI, ...
        'PointerShapeCData'   , handles.curPointerData.CData, ...
        'PointerShapeHotSpot' , handles.curPointerData.HotSpot);
    end
  else
    set(obj, 'string', '');
  end
else
  set(obj, 'string', '');
end

guidata(handles.GrabitGUI, handles);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% calculateCalibrationFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function handles = calculateCalibrationFcn(handles)
% this function calculates the calibration values

Xxo = handles.CalibPtsIm(1,1);
Yxo = handles.CalibPtsIm(1,2);
Xxm = handles.CalibPtsIm(2,1);
Yxm = handles.CalibPtsIm(2,2);
Xyo = handles.CalibPtsIm(3,1);
Yyo = handles.CalibPtsIm(3,2);
Xym = handles.CalibPtsIm(4,1);
Yym = handles.CalibPtsIm(4,2);

v1  = [Xxm - Xxo; Yxm - Yxo];
v2  = [Xym - Xyo; Yym - Yyo];

Xo  = handles.CalibPts(1);
Xm  = handles.CalibPts(2);
Yo  = handles.CalibPts(3);
Ym  = handles.CalibPts(4);

Xlength = Xm - Xo;
Ylength = Ym - Yo;

% the basis vectors in the MATLAB domain
e1 = v1 / Xlength;
e2 = v2 / Ylength;

% rearrage axes
C     = [e1 e2] \ [Xyo - Xxo; Yyo - Yxo];
blahX = [Xxo; Yxo] + C(2) * e2; Xxo = blahX(1); Yxo = blahX(2);
blahY = [Xyo; Yyo] - C(1) * e1; Xyo = blahY(1); Yyo = blahY(2);

calib.Xo  = Xo;
calib.Xm  = Xm;
calib.Yo  = Yo;
calib.Ym  = Ym;
calib.e1  = e1;
calib.e2  = e2;
calib.Xxo = Xxo;
calib.Yyo = Yyo;

handles.CalibVals = calib;
set(handles.GrabPointsBtn, ...
  'enable', 'inactive');
set(handles.CalibrateImageBtn, ...
  'value'           , 0, ...
  'backgroundcolor' , handles.bgcolor3, ...
  'enable'          , 'inactive', ...
  'string'          , 'Re-Calibrate');

handles.state = 'normal';

handles.curPointer = 'arrow';
set(handles.GrabitGUI     , ...
  'WindowButtonMotionFcn' , {@pointerFcn, handles, handles.curPointer});

set([handles.LoadImageBtn, ...
    handles.SaveAs, ...
    handles.Rename, ...
    handles.Delete], ...
  'enable', 'on');


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% splitvar
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function varargout = splitvar(varargout)
% this function splits input arguments into individual variables.


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% createSampleImageFcn
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function fname = createSampleImageFcn
% this function creates a temporary image file in the temp directory by
% loading in a sample binary image data.

str = [
  '4749463839611002CA01F70000000000800000008000808000000080800080008080C0'
  'C0C0C0DCC0A6CAF0402000602000802000A02000C02000E02000004000204000404000'
  '604000804000A04000C04000E04000006000206000406000606000806000A06000C060'
  '00E06000008000208000408000608000808000A08000C08000E0800000A00020A00040'
  'A00060A00080A000A0A000C0A000E0A00000C00020C00040C00060C00080C000A0C000'
  'C0C000E0C00000E00020E00040E00060E00080E000A0E000C0E000E0E0000000402000'
  '40400040600040800040A00040C00040E00040002040202040402040602040802040A0'
  '2040C02040E02040004040204040404040604040804040A04040C04040E04040006040'
  '206040406040606040806040A06040C06040E060400080402080404080406080408080'
  '40A08040C08040E0804000A04020A04040A04060A04080A040A0A040C0A040E0A04000'
  'C04020C04040C04060C04080C040A0C040C0C040E0C04000E04020E04040E04060E040'
  '80E040A0E040C0E040E0E040000080200080400080600080800080A00080C00080E000'
  '80002080202080402080602080802080A02080C02080E0208000408020408040408060'
  '4080804080A04080C04080E04080006080206080406080606080806080A06080C06080'
  'E06080008080208080408080608080808080A08080C08080E0808000A08020A08040A0'
  '8060A08080A080A0A080C0A080E0A08000C08020C08040C08060C08080C080A0C080C0'
  'C080E0C08000E08020E08040E08060E08080E080A0E080C0E080E0E0800000C02000C0'
  '4000C06000C08000C0A000C0C000C0E000C00020C02020C04020C06020C08020C0A020'
  'C0C020C0E020C00040C02040C04040C06040C08040C0A040C0C040C0E040C00060C020'
  '60C04060C06060C08060C0A060C0C060C0E060C00080C02080C04080C06080C08080C0'
  'A080C0C080C0E080C000A0C020A0C040A0C060A0C080A0C0A0A0C0C0A0C0E0A0C000C0'
  'C020C0C040C0C060C0C080C0C0A0C0C0FFFBF0A0A0A4808080FF000000FF00FFFF0000'
  '00FFFF00FF00FFFFFFFFFF2C000000001002CA010008FF00FF091C48B0A0C18308132A'
  '5CC8B0A1C38710234A9C48B1A2C58B18336ADCC8B1A3C78F20438A1C49B2A4C9932853'
  'AA5CC9B2A5CB973063CA9C49B3A6CD9B3873EADCC9B3A7CF9F40830A1D4AB4A8D1A348'
  '932A5DCAB4A9D3A750A34A9D4AB5AAD5AB58B36ADDCAB5ABD7AF60C38A1D4BB6ACD986'
  '00D29E5DCBB66DD9B470010894EBB6AEDDBB4FE3D235B817AFDFBF806FEAED9B9070E0'
  'C388137F1C3CD1B0E2C79023178E8BD1B1E4CB98FF32E66839B3E7CF613787EC0CBAB4'
  '69A97A53923ECDBAB550D12B57BB9E4D5B266C98B26BEBDE3D3275CEDCBC830B6FECBB'
  '27F0E1C84916F07B1BE8F1E4D01B2748686FF9C0B4E8C81647FA3CBA77EA01E0DAFF0B'
  '400FA13D02D37FC9FD85BE6BF3A4FCF87D9F5FB13A41F6D701A01B4F77FC747B854C67'
  'D57B4EC5271F7D083E649F4000A6B71E7A0BE6D75E54DB5165E08509D2176142D63168'
  '4202FEFDD36084E389182053045A78A18119CEB7E1411BDAF3A184205A57DD742F4C48'
  '548A59AD88618BD1D9C79F8E0C7668E27F011C79E33FF80140CF8CCE51D69D8A3EC607'
  '2472FF2DD72093446EF80B382216802384050000A20953B6D45C9A525569E595BA5567'
  '4B5A62DAF7A5870216495000E911F08F931FCAE8138F7399E5269C81A1891D4366DE92'
  '1692D9D9F3023D765648A28032B2191B6511691A55958802664F76229A406A417202A0'
  '5E79772EF9CFA47682FF39509E228607D74F841EF663A87FBD9067A5D881485E985986'
  '3910AC5A9E0820AA8356F8D9AEBCDEB52CAA7C725963B1C6BA5AA78DB6D2AA53AE9E41'
  '1BAD5B32A2B3A84123DAA86EAD9162DBECADC2893B2EB949F2472A7F0090A967B605CC'
  'E9294AE0BA76E1BFF37675E392E9EE2BA49182390BDD8A0417BC958C09E0F7AABEE631'
  '1C93C3F46118B1C453592CD02D60660A28432F6ECA29AF1E834CAE09E0C025EBBBF082'
  'DCB2CB6EBD80EBCA380BF4E3C73DD71670C13F070DE4D03D176DF4771C2F9DD0AE403B'
  '7D18D252130475D5AE518DF541576FED596A51570D6DD85E73B566D918B1781DDA80E5'
  '4AF6D263B35D57D385C93DD18A04BD6D77515AEFFF7D91B87AFBCD53DF826B0478E155'
  'D18DF849872F8EA2E28EA3246FE091934478E5234D8E39CD946F9EB6DA0575EE39713C'
  '8F9E13D4A59BAE1AE4AAD7A474A1AD9BC47AEC37C9FB27ED205D8EBB4CA087BE7B65B3'
  'FFCE93EDB70B0F51F0C6FB84F741A2A3AD7BF23B115F3CF4CF43FF93F4CDE38CBCF54A'
  '2DCF7CEBD5736F94F4D33B7E5BF6E2D7DE3B5F882B8E7EFAAEAFEF3BDB48BF0F3FEFF2'
  'CF2FB5EEF6DF0F13F9E573D9F6FC371600F62F33E123E0570CC8AB042AB0806F52C801'
  'DB96BA07460680B0F38E032DC8160C0670371BE4A05D3C3841B00C508492F19E045B13'
  '4214EA2A7FDFFBDA095D681A0F7E702A699919672A48C387D9B084283159A6FFFC7491'
  '16F6F059308C61C89C2410F530F178333C62BC92A8BF90E5C989A79A4CCDA438AF1F6E'
  'C551FAE14B14B9F81D1BDE702AF82A4FDEC828352AB22F2B6024D219D948343756D12A'
  'B652E31AE998343BDE912AF6C8A1129F220001FCC3907CA4D2811E024495F047877B74'
  '4A212759C844B6C98F7F8C4AA620B9364252729296848A0A1DD2C893982C2EB42A2548'
  '3E49C95016C88C9D9C0A18D58210557E8495A07465F77C24115B6ECC93B844A42E8FE2'
  'A64592722CBEE448307339CCA214F378DA01E6329B499467421399925CE627A919146B'
  '32329A4B69E5210D394D6E5E0F54D7144B322FB24D8304D39CCA1BE531B199947622A4'
  '9CF0C4092C33E9FF1EA4D85321F8CC67FC22D84B7A1685950F09A840F147D04E197428'
  '0895082E0752C985BA649F9134A151FE19116D32D3A292C3642D1F0A148E4EC4A31505'
  '29E3443AC8AFAC3321AD346945B4A9D29460742D2F75A7423342D39A96449EDF24A94D'
  '3CEA919EFA34732C9D0B0FCDB613A32AF39D4725C93EFB9253353515AAB79C68544502'
  'D4A0AA53273B7D6A44B70A92A4B6D4A560C5EA2AB54AD68EDC949F5BA92A41C66A1299'
  'B6D5226FCD285AAFFAD19388F3AE8633EB1B43C393BFA6C4AE808D485EF5EA15B95234'
  'A52BE96B6229B2D8586A34AD903DAC61277B37C10EF6B23991AC6633CB598854368384'
  '0DAD685182D8D22AE4B4738C2B6663B259D7FF3604B68E55895C6B1B59DEDA1621B08D'
  'AD56764BDA96F8F6B706096E6E558393E3F6B6B8C8B59A67CFDAD8E6AE9625AD456E70'
  '858B95AA3A17BBD7D5EE76977B929C7E17BCE1B5ED76B97B15F34297B6E975ED7AC92B'
  'BBA1C617BDEFFDED7CBB9245F68E16B3F92DED7AFD9B94B3D1E4BC2F417062F75B1527'
  'EA91BA2EB96F82157C570653C5C12B9C0985232CE1B60E98BE0A01A3B70A82E1BA6938'
  'C007EEF0563F7C14542AC48970418FB7C6C853154F18C5159EEE488952A2173CEAC5FD'
  '5114247DB9619814B9A63A36B1517C6C2605858788A8BDB1307D92DD152779C74BFEB1'
  '43F80357951CD9C8388EEA80A33C14268FD83C4FEE326B6D0C5F365B74FFCC04CE8999'
  '1D824535FB35CCF6C5B34ABB5A44A4CCB9206A7CC54086F859E3EAB92655F6299CE37C'
  '933F8F2C2D6ACC5102E0D25F327BD9CD273EF442E10CE28438FA1F6054A31317524A4C'
  '67DAD4D4E4B49FB5FCE8EE34F2CB8866AB951BBA914ED7322E954E277AEB296B315F19'
  'CB42E172AD39ACE99B8475CFBF8630530D5DEC3CD355D1C956B66CF1CBEB5E43DB981E'
  'B17554601DEB675F1BDB1DD1F668881D4E6F2339DAD2C6CAA217C2ED6E9BFBCDE84E77'
  '55D0B9E6A69033D1DC5CB7A5D54DEF92B41BC0533E37ADC3BDC043D5B5D954FE3719F9'
  '9C1171BBD5E0245138C08FAA6F46BFB2DF2179B752346ECE8A3B7C23DE1409BE372AF1'
  '1E32BCE160FF81F85A477ED092D3B0E2167F8ACA3BA2D66DBB1C859CA631534ECECE9A'
  '4B53E0E0CE3604E33DD76397FBE62224BABCDBCB738A189D292CA723CCF71DD7A63BDD'
  'DAD20C78C795BEF401F132AB512739AA5FCEF54257D7EA274578CB91EEBFA99B852E8B'
  '667B4FC27E44B77F7C313E83B3DC0B4B77B20F1CEFF4FCFA46F63E77C25BCFED31A7D0'
  '40D0EE10C3F355ED1C44FCDD859EF7783BFEF1F9963C382B1FF4995E5EB56357A0E02D'
  'B779C633E4F3A00FFDFD302E92C97346BADB457DEAB52EF5D193FEABB0FFBB44656F5D'
  'D58B2FE4BD21A9E961EAFBA3F47DF5C06FBD50877F10DECFDE92C9573EEE73DFF98E42'
  'FEE8C54FDECC833FFDC533FFB1D787BAFFF331C77AE977DFFB3A1E3FE6A15F762567E5'
  'CC8EF9BE40D43FF14422DEB255113184E54FFFE7DBBFFDC02615F8C26A8CC5799E177E'
  'F6D67F88737FF847159FD680FF607B09957D1B777CBF078001E88004984184F17D0AE8'
  '7F6CC4805407154CC63C5B847ED5075014887DB46772BAC77D56F18023387C1F587F52'
  '247F7D761532185B12C86E2BC8820B8781EEA7814D667605C1781638151CF740220881'
  '51B183AB818408281549283C4D98784501850C617535C8773F383A5788854461666154'
  '804768435DE88553D83A61E87A0F612B2738829C9782E0D78266B184C8F78230884668'
  'B1851E94868F078870238443B857AFF5875F988058D776849881FDFF745B0024880027'
  '897DA4877B68880B717294686C3E4740612887DDE5103D588776F174E9D38885F88890'
  '98449BD87B55683A9F088AEDF5107CD68A9C8887D083831AE186BB685A12688BC6364E'
  'C0182DB128864EF11C2AF78A5B318C2C838A8EA88AABE83DCC6883796889B7777E9948'
  '6FD3385BA2E78CCFB86CA2A88D895815CA1839C5688C4D011CDEB48D00D68DD6788DA9'
  'E587E2B8865C518E8B738EBCA81210C78EED588DB6215487082AFC488DDC538C71088D'
  'B4883703D98FA7F88EB9434BF958199485210BC990B9E88D19168F1271211549905688'
  '91A9088E1B69201D6991C6738EB2388B13191FE3D81525191C289992038257F638162F'
  'C91B31FF898E2852112CD9922E7993BA91931129913C0994A9A77D42B97957478F6461'
  '94B391933AB9149EF2497498184ED91A50399416F12F940492547195A7A18B94A79169'
  '57485EA98460591A50199505D67366799654E89338B3965A591153394962C916698944'
  '55191B0039816F099771C99483E890E5A594BB574979B9167B79196BC9966DB99488F4'
  '9815489841F3987569970738877D79188D994282D975335996C2B49867519395D899BA'
  '85988DF751A6799A72498CA1299A89939899F59A65819A2E43999929419466673AD557'
  'B8999BB189289409993A016321497C1C358A90A19B756498CCD56090B69CCD5773C3D9'
  '94C57925C7E986690165065162DFFF789DD6969D62019DB2399BB4691399026AF97266'
  '18F69E6F444BCC696EE6799E9F0918EA399E3BC11F94121EB9C62471A12833469F2AB8'
  '88C7596DF97917099A65ED41680B610F8F049C75C85BF719168BC83609FA7169562AF2'
  '19A11DEA84F7F46ECEA918A688350DEA7013EAA1720423212A93F3975E17EA15272A36'
  'FBB99E33B1A243746624762E224A9E768882D2598A9DA8A1376A843DA1A368D2A2FF20'
  '6803B1A4483AA23836A3CB58A465D3A0C8694AE1212B43D26A91F6A06548A1313AA558'
  '0A5156BA3558AAA21D6A2F5EDA443CF41C1446A55AE154463AA4AB391492229F8E3221'
  'A1E6550D11A7654A1467DA46471AA53CE19F7FC2490E059841FF7A86812AA8E8D93185'
  '6AA887AA2828C7A8BE38A94058A7767AA76479A07A56A28F818BBB29A744F9A9F53992'
  '9DEA17A41A9D9A6A8622E9838766AA5FD1AACD28AAFFF897C1A969E57719B66A9CDB67'
  '13BDD918AD596CC12A19BF9A8032177D3531AC05F5A7CD76ACC81AA9EBB714CCDAACAC'
  '09A42B89AB264AADB79896D24A13CEBAA8B2DAA84F13AE91E1ADEEA6AE2B45AB3988AA'
  'BBBAADEE4AA3EC5A54E63AA6350A148F9AA59D029E2E911BB076AD9951AF62455A2895'
  'ACD1F3AA949A6D3AF77AA7676AE83AAD0B79B0144BB05CB5AFD9236C38211BFFD6AB9E'
  '81B0CC56B1220B1F0ABBB0BDE1AF5655AE9FB3AFFE64B11621B2CC54B1CE54B24A7590'
  '22FF3180003029B90AAA8145B39BAA86745A10E2746F207B5165BA17FD4332A0561E29'
  'C3121C5B7C2C6B7C2E6B7DF90AAADBD9B3AA3913499B1D7722322FB11A4817B562074A'
  'F77A69414B9A86B589F37AA928D12ACB5184B8A1B21C21B6633BB5407AB692895236B1'
  'B6A78A122F901D8E22AE3CEB1174BB7645AB1144659907BBB73E2BA6B529A584C9B7F4'
  '7AB813B8AB16ABB73421B619DB0E3B411A8427B9F598A1883BB2CF35A836A5B926612F'
  '612AB8A95A56851B143455B6688BB77726BA2B41B7EF331E9912A0290BB917BBAA7F41'
  'B1F64A546D66BB21D5B88EDB118E820E28EBB483EBBAC01BBCC23B78C42B13552B55C8'
  '9BBC1EF10B5ADBBA2101BAC4FF645A302BBB966BBC210B88855B42CBDBBDBE8B54D1CB'
  '15C4E3265DC90F8B5B74D055BFB377B5C085BB2BB1BBD8DABEBFFBBEFCB63CC5E42381'
  '49BFF43BBE158AB9404BBED06B1407B4BAEC1BAF3F95BD465BC0189CC1E3ABC03F41B9'
  '78F5BAF653221B0BC0D82BC0F396C106B72B06BCC1FF746F664A7FE0DBB723C1A6234C'
  '9EC79BB55D81C2544450F251251CAC88FA1B8116ACBDF5C1B90C7222AC2BB457FBBA26'
  '81C119B1C20C9C750E8CB52DC6125B52ADED6AC2A224B08A05C51EDCC1161BC3326C4A'
  '48BCB3F63BC54F3CC44DCCC57723C43D4995385CB768CC936A0CABC3A5ADFA58C701CC'
  'AD7F93C0F31BC7723CC7AA0AC8FFFBA9532BC6998B4E84FF3CC88139A7A887C8EF3AC3'
  '85B031F9725630CCC473CBC70F17B1406C998EBAC8494C123E06268A3AC368E29E7264'
  '18FD07C9170CBA059C155F1C905ADCBB27E1288AF29D2481A834DC800A88C91FECB115'
  'CCC69D2CC8D1189929610FA51C123EF6A04B9A4A14CC12AC8C54D407CA1F21CC2C48CC'
  'F2A814217CCA41B4A6CD9C37945183BE5CCCC0DCCAB35C787217CD5B694A6892CC21A1'
  'A41F1A655D38CEE7FACA3AA1C995B982EAAC99A93BC9FDBBA51E12CF0134CFFB9CC6F6'
  '3C3C052D72E61B8EE7FCAF5501CFFEDA17415CCFF4CC992C8BCF0A5A85094DAC8004D0'
  'B552C9912488E5DC4D071D4F1BBD72947BD2E42A806B0AA01985BE23ADAF9CAC3E2AFF'
  '3DBC0A57D1E32A10CB9C007B3A63673C13D6DC1341CDB8186DB81658D3BAA649E1511E'
  '00A0A87B01D3332DD4319DB0155D9F875BD539EDA19DF1D413CD35518DD048CDD07A5C'
  'AC293DD6263B9AC288CD9B5CD444FD26D4DCD60DED6CED16D6F3847B76BBBF433D5071'
  '0DD43FF2D6F9AB7655ADD63B99D6704DD7D37C49790DC60BFDC9406CD75D9D899C2914'
  '862DAF6C4D65693DC515BDD8A881AF60DDB8952DD9892D145147CFA65B9B779DA9345B'
  'D5623DD9A489DA7E6DBDA58D43C2784E09FAD9E313DAB0AB62AC0DB9B1CC37846DD239'
  '694CBB6DD0B6ADD87846DAB1BDD9FC58DC996C25C38DB5AA7DB7AB2D7EBDCD3715C9DC'
  '7DFCD55B3CD5527B5FCFFF7DB79EEC1C1D89DD941DDD9DC5DD635BCC8DAD4EE34DDE83'
  'FCDDD5ECDE58ECD566ADC4E1FD1A2589DED00DDFD05BDF0F7B5DFC6DDF0290D5644D4C'
  'F22D8FFA3DC0FE5DAE5A97E07CC5AF868B140E7EDEDADD23E6FDCC16FDDAB05D4904FE'
  'DFF79DC55CA7C36E71E0AE3865158E59360BCBD56DCE0386C20C1AE01E3ECD7B5DBA30'
  '6A73A7DDDCFA4D6BC0A7E1277CE1A4E8C6244E6E352EC51F0EE28CD7D725AD9F306EB5'
  '40BEE47FD5E1BCFDD884EBC442BC78429AE48A11E4A736995ABE6B3FFA95C92DD53A3C'
  'E65DEE153E8EAF8654E6663B65502EE037BED6643EE67CB9E0A9BAE49CFDE560BEE2FF'
  '13E72E5E43760EB9747EE7788E966F3EE52AC42272FFCE1A6A0E6678C9E31B1E706DEE'
  '4EE3F429559E908B3E74812EDDEB3DE82A89E9994ED59F2EE0D9A475916E9DF01BEA9D'
  '3DE331E1C7455EBB38BAD967E1E87F71E9B7FBC6AD1E71A255EAFC194AB47EBCD77B8B'
  'BB1E8A9BD6EBC17CD952EE962DA8EB67CD7E13AED7C67EEC9B19EC683DECC4DEDF68DE'
  '91C7A5EC441C82279EC807F2EB29F65EDA6EC7CD84DB0CF5ED611EEE413AEE9CCEECD5'
  'BED6675CE8B69991A02550584EDBC604EE8CBEEE5ED15FEC9E20F7DED992AED9C54BEF'
  '546160F066EED01C89FC8862592D9EAF2E754D9EBD18C48EE7F5F0D569EAF9F4EEBF0C'
  'ADD0EEF17D6814B7EC273C0AF1CB9E79818A86CC8860771748E5B1A7010A63DFFF9900'
  '059AE21B5FA649B5F261F671E3F1A02F7A1F412648EDBE751C6FE9687B94215F6620BD'
  'A328131EDB5EEE458FE047EF8A751D1410CDA383F6F3FF3E2E512FCB842CEFA25EF540'
  '71F5200AD2436FEFFBBC589B48616ADAD4018DF5035167E4FE6DAA7EDEC82ED834A767'
  'B6E6F6C206A10311688306A54F4FF428795A60FFDB49AF1303C8201EDDA5EEF94492F6'
  '9B279FF0620C5B92A8626DEFF365D8A702BA1A5B6F33AEDCC3265CE8B6FA71CBFCF21F'
  '0A462B2D607CDBEC964B6C6CC6F370E8F6912C5F22E8FA788C5F9AD6F226E0B0FAA5F0'
  'AB3D602EFBEB9F1F5D133FA3DDCEE4A3F5AB5B5FFC7309FCD988FBBCBD66040FE143E1'
  'FCA9999D015F63C5FF09EECD6FFC19DED0540E767AEEE647F6FDE03FFEC42DFD0C4E78'
  'FA8EFEE99FC1BF0CFD64FD8508D56C9F8FFDCF2FFFEF0D92E95EFF00F14FE04082050D'
  '0E047050E142860D1D3E8418F15F4289152D5EC49851E3468E1D3D7E8CC84F24BF8123'
  '4DFE2309D1E448901505BC8409B365CC971D29B6C49953E14D9D3D7DFE041A54284891'
  '25579E2458D4E851A54307D2A4F911EA479E4EAD46AC7A55EB56AE5DBD3664CA12E5C9'
  'B05EA14E15A0312AC8AC5FBBB6751B57EE5CBA1CC3DE15EBF6ECD98C6BA9D67D0B58F0'
  '60C28597E25D2917ED5E973273C235FC38F264CA95AF224E09B8E653BE69194EF50840'
  'F468CB3F2197469D5AB5CAB2A5F77A16B8F91FFF688CA349AF368D5BF76EDE8713ABE6'
  '3B3BA670C7126D9FEEFD37F972E6935BA77EFD1AE271E4CD955BC79EDD2DD3DDD1FD1A'
  'A4CE96807686D5C99F474FB429EFB4DE0B1EC769CF84E85B000824484FD07C7EFEFDFD'
  '1F7C2D3CAAC021C89E0000A0E71F030140C7BFFDFE8330C2F3A40349BED1E61B4FA017'
  'ECC3CFC20CF97B50421147DC6D34C640FA0541055354F19F0333F4103F1049A4B1C612'
  'E18B8DB68F7E91F11F161B5C90400531EC31BD106D4432C9AB0454A83DD87CAA8F4105'
  '0F1432C6FE8E54324B2DA9B26DB003550C52202BF3C3724B33CF748849C2EC11EDBE29'
  '01A8F2C52BD1A4B34EAC70B42CCA04DF84F1C00667B433500737C253B530FFC5DCD0CD'
  'FADC04545047D13B524DDD0E1D68C104851CF4514DB5434ED2E42CCC0AD4327B1B7553'
  '530DCBAAD0EC400567B422472CF5545901A32ED6D4A2BC2D495B67E5752B4F7B4D6E57'
  '60877D4C5562757B5520618F65F622639BC52DCA642782B65A9D7EB5D6D02F454B76D9'
  '6C9BC5362329BF956BC3B6BC2577D667350A37DDABCC9DF62674DD45B35D8B8A64714F'
  '7AE382F7A02EF765D65E8CA425285F80F9E5D6DF83875DD7A3055D2D58347DA75D38A8'
  '7EC1AB78538171BA582016472300008A330E0D530D13C698E47A1B1EAA631F2B9D0F4E'
  '953F7AB82D97119A594996B7BAB940397306AA676A819670E3AE8416F367A27B82775C'
  '9C973672E7BAFF907EF34FA85BDA563485AF16CC56A90BA33AC50FB976F8A179C90EAD'
  '22A32773F915372D44DBAAB3E3E6A8D3AF5573D95CD1ACA63BA8B9FB1617BC7FD1D3FB'
  'CF1401D7EA6FC42DAAB5BFAC155F7CBAC813BF1B3B7B26970B72CC951D7C739535077C'
  '6DCFE905FD6AD1475FB8F4992B473D67D5013EBD21510960B475706D9F8875875E3D30'
  'C128F9C69D61CF758788E04A39E43378E5358D5DA2AC8B4C544C229707F6758D89CFE8'
  'F27F3A3E744CEAD55D1AFB96B8A7527AE4BF073FE3E67D225FE621CF47FF54EB6D5CBF'
  '6594DFAC12C3F8D3B7567CAEDAEF93D3F6A7B12CEDAA7E6E0160870250BB0132AF808C'
  'F31F60F2D6BB7FD4A781B29A5F65CC1341FFC3B84C5426BBA00375E62F0E56A670950A'
  '61F51E7840D59430855BCA20AA3AE71F03BD305B31A40B0B6DB8437639C8853C04620F'
  'CFF3C320169150D6D161A0EC5100236A0987384922428067A3292A88894D1C616A88C8'
  'B95C29E905D35A2216B31899281AC46036920FDF6A6890308A515783296343CE28A117'
  'F0ED17B5BBA38F46838E04B4D18DF49BCB162532C79135E716E0E810390AB2C6ED25C8'
  '40F8B9631FAFF8C71A3D317733DCCAC7DA24B2F4C887405F2C481DC5548857D983947E'
  'A42489AC17C71D518455F9F9453D14694637FD02530FBB0F2A532922CD0932280B1A1B'
  '795E5000309AA0412FD8D3123B74CA49EEB26850F4A55580991F36E9FF6B227FBA4501'
  '1030C93C6E2F97CD742679F6F89EBA45B32B0BAA22766E81804208C482302385D54045'
  '0F1384F33C9A5C4899CC2917B151B39DB7A8A23D5254487B6A678E5B6B082BEBF28A7F'
  '74683EFCB985234959D02B15F2A0297BCF3E05F38290ED8D3F7E3C2445D3633C85180C'
  '7E4FD3A86158F41FF95853A4E4A9192717A2C998BDAA8B2FC5A96EA8B6484AE5D4A7CB'
  'D9E92295F653A2EE26A83EB34F5195AAD3FB39E4474B856A6A76CAD081E82FAA57A595'
  '6840E8323D69E86DC08B5952098A55B256E8400BE1AAC4B8A8B502A96851E92C6B5CED'
  '2753CEEDE9170B31D0DB862A57BECE755CB78848F4DEC7C0BE16362759C348F7A66758'
  'C6B6447B4B1B51EC491B3BD9AF4496B094C5AC562C9B59CE72659A55EB6C68A5F92205'
  'D255B4A7FD89B97C9754D4B6D6273503A16B655BA17ACED6B6B7C56D6E75BB5BDEF6D6'
  'B7BF056E70854BA78000003B                                              '];

% convert hex to dec
s = hex2dec(reshape(deblank(reshape(str', 1, [])), 2, [])');

fname = [tempname '.gif'];
fid = fopen(fname, 'w');
if fid > 0
  fwrite(fid, s);
  fclose(fid);
else
  fname = '';
  errordlg('Error trying to create sample image.');
end