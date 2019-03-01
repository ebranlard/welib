function [Cl_fs,alpha1,alpha2,alr1,alr2]=fPolarFullySeparated(alpha,Cl,dclda_fs,alpha_merge,delta_alpha)
% returns a Fully separated polar based on input polar (extended to 360 if possible) in a way similar to hawc2/Oye.
% Note: Input should be consistent : if alpha is in degree then all the rest should be in degree (even dcl_dalpha!!!)

% dclda_fs: slope of fully separated airfoil in linear region (usually half the airfoil:  pi [1/rad] or pi*pi/180 [1/deg])
% alpha_merge: angle where the fully separated polar merges with the original polar
% delta_alpha: characterize the extent of the interpolation region between the linear region and the merging point


% Examples of input (in degrees!!!!)
% alpha_merge=40;
% delta_alpha=35;
% dclda_fs=pi*(pi/180); % deg
% Examples of input (in radians!!!!)
% alpha_merge=40*pi/180;
% delta_alpha=35*pi/180;
% dclda_fs=pi; % rad


% Find alpha0
[alpha0]=fAlpha0(alpha,Cl);


%% Everything below is indepent of the dimension used for alpha (deg or rad)

% Positive stall
Cl1=interp1(alpha,Cl,alpha_merge);
alpha1=alpha0+Cl1/dclda_fs; % deg
if(delta_alpha>alpha1-alpha0)
    alr1=alpha1-alpha0;
else
    alr1=delta_alpha
end
% Negative stall
Cl2=interp1(alpha,Cl,-alpha_merge);
alpha2=alpha0+Cl2/dclda_fs; % deg
if(delta_alpha>alpha0-alpha2)
    alr2=alpha0-alpha2;
else
    alr2=delta_alpha
end

% Create full separated profile Cl coef
Cl_fs=zeros(size(Cl));
for i=1:length(alpha)
    alfa=alpha(i); 
    if (alfa >= alpha0)  
        % Positive part values of Cl
        if (alfa < alpha1-alr1)  
            % linear region
            Cl_fs(i)=dclda_fs*(alfa-alpha0) ;
        elseif (alfa < alpha1+alr1) 
            % interpolation
            Cl_fs(i)=dclda_fs*(alfa-alpha0)-dclda_fs/4.0d0/alr1*(alfa-alpha1+alr1)^2;
        elseif (alfa < alpha_merge)  
            % possible constant zone
            Cl_fs(i)=Cl1 ;
        else
            % merged polar
            Cl_fs(i)=Cl(i);
        end
    else
        % Negative part values of Cl
        if (alfa > alpha2+alr2)  
            Cl_fs(i)=dclda_fs*(alfa-alpha0) ;
        elseif (alfa > alpha2-alr2) 
            Cl_fs(i)=dclda_fs*(alfa-alpha0)+dclda_fs/4.0d0/alr2*(alfa-alpha2-alr2)^2 ;
        elseif (alfa > alpha_merge)  
            Cl_fs(i)=Cl2 ;
        else
            Cl_fs(i)=Cl(i);
        end
    end
end
