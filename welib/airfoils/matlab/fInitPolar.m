function Polar = fInitPolar(varargin)
%  Assumes sets of key , value arguments

%% Polar Defintion
Polar.id=[];
Polar.Re=[];
Polar.alpha=[];
Polar.Cl=[];
Polar.Cd=[];
Polar.Cm=[];
% 
Polar.f_st=[];
Polar.Cl_inv=[];
Polar.Cl_fs=[];


%% Argument handling
format=[]; % because it already exists in matlab
keys_string  = {'filename=','format=delim-alpha-Cl-Cd-Cm','id='} ;  % keys that accept strings as value
keys_numb    = {'nlines_header=0','Re=[]'} ;  

key_value_handler('keys_string',keys_string,'keys_number',keys_numb,{varargin{:}});



%% 
if length(filename)==0
    error('TODO Computed Polars (naca, Karman Treftzs')
else


    % inputs for fully separated polar (in degrees) - TODO
    alpha_merge=40;
    delta_alpha=40;
    dclda_fs=pi*(pi/180); % deg
    
    % 
    filename
    Re
    format
    if isequal(format,'delim-alpha-Cl-Cd-Cm')
        Data=load(filename);
        Polar.alpha=Data(:,1);
        Polar.Cl=Data(:,2);
        Polar.Cd=Data(:,3);
        [Polar.Cl_inv Cl_inv_sin alpha0]      = fPolarInviscid(Polar.alpha,Polar.Cl)                                        ; 
        [Polar.Cl_fs,alpha1,alpha2,alr1,alr2] = fPolarFullySeparated(Polar.alpha,Polar.Cl,dclda_fs,alpha_merge,delta_alpha) ; 
        Polar.f_st=(Polar.Cl-Polar.Cl_fs)./(Polar.Cl_inv-Polar.Cl_fs);
    else
        error('TODO a good handling of various format by splitting the format string')
    end
end
