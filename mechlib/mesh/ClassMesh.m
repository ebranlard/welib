%%
classdef ClassMesh < handle

properties (SetAccess=private);
    nDim;
    n;
    nGridPoints;
    nCells;
    v1,v2,v3;
    dv1,dv2,dv3;
    dCell;
    xMesh_min;
    xMesh_max;
    bCartesian;
    bRegular;

    bGridPointsComputed=false;
    bCellCentersComputed=false;
    bCellVolumesComputed=false;
    bFakeCellVolumesComputed=false;

    GP; % GridPoints
    CC; % CellCenters
    CellVolumes;
    FakeCellVolumes;
end
properties
end
methods 
    % constructor 
function o=ClassMesh(v1,v2,v3);
    if nargin==1
        o.nDim = 1   ; 
        o.v1=v1;
        o.v2=[];
        o.v3=[];
        o.dCell(1)=v1(2)-v1(1);
        o.xMesh_min=[v1(1)];
        o.xMesh_max=[v1(end)];
        o.n=[length(v1)];
        o.bRegular=mean(diff(v1))==v1(2)-v1(1);
        if ~o.bRegular
            o.dCell=[]; % to trigger errors
        end
    elseif nargin==2
        o.nDim = 2   ; 
        o.v1=v1;
        o.v2=v2;
        o.v3=[];
        o.dCell(1)=v1(2)-v1(1);
        o.dCell(2)=v2(2)-v2(1);
        o.xMesh_min=[v1(1) v2(1)];
        o.xMesh_max=[v1(end) v2(end)];
        o.n=[length(v1) length(v2)];
        o.bRegular=sum(abs(diff(v1,2)))<1e-12 & sum(abs(diff(v2,2)))<1e-12 ;
        if ~o.bRegular
            o.dCell=[]; % to trigger errors
        end
    elseif nargin==3
        o.nDim = 3   ; 
        o.v1=v1;
        o.v2=v2;
        o.v3=v3;
        o.dCell(1)=v1(2)-v1(1);
        o.dCell(2)=v2(2)-v2(1);
        o.dCell(3)=v3(2)-v3(1);
        o.xMesh_min=[v1(1) v2(1) v3(1)];
        o.xMesh_max=[v1(end) v2(end) v3(end)];
        o.n=[length(v1) length(v2) length(v3)];
        o.bRegular=sum(abs(diff(v1,2)))<1e-12 & sum(abs(diff(v2,2)))<1e-12 & sum(abs(diff(v1,2)))<1e-12 ;
        if ~o.bRegular
            o.dCell=[]; % to trigger errors
        end
    end
    o.dv1=diff(o.v1);
    o.dv2=diff(o.v2);
    o.dv3=diff(o.v3);
    o.nGridPoints=prod(o.n);
    o.nCells=prod(o.n-1);

    o.bGridPointsComputed=false;
    o.bCellCentersComputed=false;
    o.bCellVolumesComputed=false;
    o.bFakeCellVolumesComputed=false;
end

    function Volume=getCellVolume(o)
        if o.bRegular
            Volume=prod(o.dCell);
        else
            error('Cell volume not defined')
        end
    end
function Volumes=getCellVolumes(o)
    if (~o.bCellVolumesComputed)
        o.CellVolumes=zeros(o.nCells,1);
        if o.bRegular
            o.CellVolumes(:)=prod(o.dCell);
        else
            switch(o.nDim)
                case(1)
                    o.CellVolumes=o.dv1;
                case(2)
                    p=0;
                    for i=1:o.n(1)-1
                        for j=1:o.n(2)-1
                            p=p+1;
                            o.CellVolumes(p)=o.dv1(i)*o.dv2(j);
                        end
                    end
                case(3)
                    p=0;
                    for i=1:o.n(1)-1
                        for j=1:o.n(2)-1
                            for k=1:o.n(3)-1
                                p=p+1;
                                o.CellVolumes(p)=o.dv1(i)*o.dv2(j)*o.dv3(k);
                            end
                        end
                    end
            end % switch ndim
        end % if not regular
        o.bCellVolumesComputed=true;
    end
    Volumes=o.CellVolumes;
end % getVolumes

    function Volumes=getFakeCellVolumes(o)
        if (~o.bFakeCellVolumesComputed)
            o.FakeCellVolumes=zeros(o.nGridPoints,1);
            if o.bRegular
                o.FakeCellVolumes(:)=prod(o.dCell);
            else
                switch(o.nDim)
                    case(1)
                        o.FakeCellVolumes = ([0 o.dv1]+[o.dv1 0])/2;
                    case(2)
                        X1=([0 o.dv1]+[o.dv1 0])/2;
                        X2=([0 o.dv2]+[o.dv2 0])/2;
                        p=0;
                        for i=1:o.n(1)
                            for j=1:o.n(2)
                                p=p+1;
                                o.FakeCellVolumes(p)=X1(i)*X2(j);
                            end
                        end
                        Zc=[];
                    case(3)
                        X1=([0 o.dv1]+[o.dv1 0])/2;
                        X2=([0 o.dv2]+[o.dv2 0])/2;
                        X3=([0 o.dv3]+[o.dv3 0])/3;
                        p=0;
                        for i=1:o.n(1)
                            for j=1:o.n(2)
                                for k=1:o.n(3)
                                    p=p+1;
                                    o.FakeCellVolumes(p)=X1(i)*X2(j)*X3(k);
                                end
                            end
                        end
                end % switch ndim
            end % if not regular
            o.bFakeCellVolumesComputed=true;
        end
        Volumes=o.FakeCellVolumes;
    end % getVolumes


function [CC] = getFlatCellCenters(o)
    if ~o.bCellCentersComputed
        Xc=zeros(o.nCells,1);
        Yc=zeros(o.nCells,1);
        Zc=zeros(o.nCells,1);
        switch(o.nDim)
            case(1)
                X1=o.v1(1:end-1)+o.dv1/2;
                Xc=X1; Yc=[]; Zc=[];
            case(2)
                X1=o.v1(1:end-1)+o.dv1/2;
                X2=o.v2(1:end-1)+o.dv2/2;
                p=0;
                for i=1:o.n(1)-1
                    for j=1:o.n(2)-1
                        p=p+1;
                        Xc(p)=X1(i);
                        Yc(p)=X2(j);
                    end
                end
                Zc=[];
            case(3)
                X1=o.v1(1:end-1)+o.dv1/2;
                X2=o.v2(1:end-1)+o.dv2/2;
                X3=o.v3(1:end-1)+o.dv3/2;
                p=0;
                for i=1:o.n(1)-1
                    for j=1:o.n(2)-1
                        for k=1:o.n(3)-1
                            p=p+1;
                            Xc(p)=X1(i);
                            Yc(p)=X2(j);
                            Zc(p)=X3(k);
                        end
                    end
                end
        end
        o.CC=[Xc Yc Zc];
        o.bCellCentersComputed=true;
    end % already computed
    CC=o.CC;
end

    function [GP] = getFlatGridPoints(o)
        if ~o.bGridPointsComputed
            Xc=zeros(o.nGridPoints,1);
            Yc=zeros(o.nGridPoints,1);
            Zc=zeros(o.nGridPoints,1);
            switch(o.nDim)
                case(1)
                    Xc=o.v1; Yc=[]; Zc=[];
                case(2)
                    p=0;
                    for i=1:o.n(1)
                        for j=1:o.n(2)
                            p=p+1;
                            Xc(p)=o.v1(i);
                            Yc(p)=o.v2(j);
                        end
                    end
                    Zc=[];
                case(3)
                    p=0;
                    for i=1:o.n(1)
                        for j=1:o.n(2)
                            for k=1:o.n(3)
                                p=p+1;
                                Xc(p)=o.v1(i);
                                Yc(p)=o.v2(j);
                                Zc(p)=o.v3(k);
                            end
                        end
                    end
            end
            o.GP=[Xc Yc Zc];
            o.bGridPointsComputed=true;
        end % already computed
        GP=o.GP;
    end


function v_flat=flattenValues(o,m)
    % Takes Mesh values m(1:nd,1:n1,1:n2,1:n3)
    nGridPoints=prod(o.n);
    s=size(m);
    if(length(s)==o.nDim)
        m=reshape(m,[1 s(:)']);
        nd=1;
    else
        nd=size(m,1);
    end
    nInput=prod(size(m))/nd;
    if(nInput~=nGridPoints) 
        warn('Wrong mesh input')
        kbd
    end
    v_flat=zeros(nd,nGridPoints);
    switch(o.nDim)
        case(1)
            v_flat=m(:,:);
        case(2)
            p=0;
            for i=1:o.n(1)
                for j=1:o.n(2)
                    p=p+1;
                    v_flat(:,p)=m(:,i,j);
                end
            end
        case(3)
            p=0;
            for i=1:o.n(1)
                for j=1:o.n(2)
                    for k=1:o.n(3)
                        p=p+1;
                        v_flat(:,p)=m(:,i,j);
                    end
                end
            end
    end
end % function flatten values

end
end
