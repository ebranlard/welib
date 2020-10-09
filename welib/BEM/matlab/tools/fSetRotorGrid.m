function [ WT ] = fSetRotorGrid(r_new,WT,varargin)

%% interpolation of Rotor radial positions if needed
if(r_new==0)
    % OK, we use the original one
    [ WT ] = fInterpRotor( WT.Sources.Rotor.r ,WT);    
elseif(length(r_new)>1)
    % We use the r vector the user give
    [ WT ] = fInterpRotor( r_new, WT);    
else
    % We respect the Algo.Ngrid parameter
    Ngrid=r_new; %that's the trick
    rfull = linspace(WT.Rotor.rhub,WT.Rotor.R,Ngrid+1)'; % could be change for more fancy cos
    r_mid = (rfull(1:end-1)+rfull(2:end))/2; 
    [ WT ] = fInterpRotor( r_mid,WT ); 
end
WT.Rotor.rfull=sort(unique([WT.Rotor.r(:)' WT.Rotor.rhub WT.Rotor.R]));

if(WT.Rotor.r(1)==WT.Rotor.rhub)
    disp('!Warning the first node is at the Rotor hub, don''t use AeroDyn integration method!')
end
if(WT.Rotor.r(end)==WT.Rotor.R)
    disp('!Warning the last node is the Rotor radius, don''t use AeroDyn integration method!')
end

WT.Rotor.e_ref_for_khi=find(WT.Rotor.r>0.7*WT.Rotor.R,1,'first');
WT.Rotor.ne=length(WT.Rotor.r);


%% Rotor dr if needed
WT.Rotor.dr=WT.Rotor.r*0;
WT.Rotor.dr(1)=2*(WT.Rotor.r(1)-WT.Rotor.rhub);
for i=2:length(WT.Rotor.r)
    WT.Rotor.dr(i)=2*(WT.Rotor.r(i)-WT.Rotor.r(i-1)-WT.Rotor.dr(i-1)*0.5);
end

if (WT.Rotor.dr(end)*0.5+WT.Rotor.r(end)-WT.Rotor.R)/WT.Rotor.R >0.005
    disp('!Warning Radial positions with dr invalid')
end

end



function [ WT  ] = fInterpRotor( r, WT )

if(isequal(WT.Sources.Format,'wtperf'))
    disp('!Warning: Format is WT_perf, no interpolation possible, if needed consideration interpolation of ProfileSet!')
else
    WT.Rotor.chord         = interp1(WT.Sources.Rotor.r,WT.Sources.Rotor.chord,r) ;
    WT.Rotor.thickness_rel = interp1(WT.Sources.Rotor.r,WT.Sources.Rotor.thickness_rel,r);
    WT.Rotor.twist         = interp1(WT.Sources.Rotor.r,WT.Sources.Rotor.twist,r);
    if sum(isnan(WT.Rotor.chord))>0
        %ok, we need to do some extrap...
        warning('I have to extrapolate chord twist and thickness...Is that normal?')
        I=1:length(r);
        IEnd=I(r>max(WT.Sources.Rotor.r));
        IStart=I(r<min(WT.Sources.Rotor.r));
        if ~isempty(IEnd) 
            WT.Rotor.chord(IEnd)=WT.Rotor.chord(IEnd(1)-1);
            WT.Rotor.thickness_rel(IEnd)=WT.Rotor.thickness_rel(IEnd(1)-1);
            WT.Rotor.twist(IEnd)=WT.Rotor.twist(IEnd(1)-1);
        end 
        if ~isempty(IStart) 
            WT.Rotor.chord(IStart)=WT.Rotor.chord(IStart(end)+1);
            WT.Rotor.thickness_rel(IStart)=WT.Rotor.thickness_rel(IStart(end)+1);
            WT.Rotor.twist(IStart)=WT.Rotor.twist(IStart(end)+1);
        end
    end

    WT.Rotor.r             = r(:)';
    WT.Rotor.chord=WT.Rotor.chord(:)';
    WT.Rotor.twist=WT.Rotor.twist(:)';
    WT.Rotor.thickness_rel=WT.Rotor.thickness_rel(:)';
    
    if(isequal(WT.Sources.Format,'flex'))
        Rotor.ProfileSet=zeros(3,length(r));
        for i=1:length(r);
           profilesup=find(WT.Profiles.r >= WT.Rotor.r(i));
           profileinf=max(find(WT.Profiles.r <= WT.Rotor.r(i)));
           if isempty(profilesup)
               profilesup=profileinf;
           elseif(isempty(profileinf))
               profileinf=1;
           end
           profilesup=profilesup(1);
           WT.Rotor.ProfileSet(:,i) = [(profileinf~=profilesup)+1,profileinf,profilesup];
        end
    else
        % Profile Sets
        WT.Rotor.ProfileSet=zeros(3,length(r));
        for i=1:length(r)
            profilesup=find(WT.Profiles.thickness_rel >= WT.Rotor.thickness_rel(i));
            profileinf=max(find(WT.Profiles.thickness_rel <= WT.Rotor.thickness_rel(i)));
            if isempty(profilesup)
                profilesup=profileinf;
            elseif(isempty(profileinf))
                profileinf=1;
            end
            if isempty(profilesup)
                warning('profilesup empty')
                keyboard
            end
            profilesup=profilesup(1);
            WT.Rotor.ProfileSet(:,i) = [(profileinf~=profilesup)+1,profileinf,profilesup];
        end
    end
    
    
end
end


    
