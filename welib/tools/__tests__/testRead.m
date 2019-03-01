function test_suite=testRead
initTestSuite;
% 
setDefaultPath;
global VERSION
if isempty(VERSION)
    warning('Using default version...');
    VERSION='v01';
    global VERSIONNUM
    VERSIONNUM=str2num(VERSION(2:3));
end
require('WTlib',VERSION,1);



function testInit
global PATH
global VERSIONNUM

if(VERSIONNUM>1)
    WT=fInitWT('NREL5MW','hawc',PATH.DATA_WT);
    WT=fInitWT('SB1','xblade',PATH.DATA_WT);
    WT=fInitWT('Tjaere','flex',PATH.DATA_WT);
%     fprintf('S');
else
    fprintf('S');
    WT=fInitWT('SB1','xblade',PATH.DATA_IN);
%     fprintf('S');
%     fprintf('S');
end


function testFun
