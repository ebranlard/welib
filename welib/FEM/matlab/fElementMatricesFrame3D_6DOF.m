function [ke,me]=fElementMatricesFrame3D_6DOF(E,G,Kv,EA,EIx,EIy,EIz,L,A,Mass,R)
% Stiffness and mass matrices for Hermitian beam element with 6DOF per node
% Euler-Bernoulli beam model. The torsion is de-coupled to the rest, not like Timoshenko. 
%
% The beam coordinate system is such that the cross section is assumed to be in the y-z plane
% 
% (ux uy uz thetax thetay thetaz)
% (ux1 uy1 uz1 tx1 ty1 tz1 ux2 uy2 yz2 tx2 ty2 tz2)

% The torsional equation is fully decoupled and as follows:
%   Ipx/3A Mass txddot  + G Kv /L tx = Torsional Moment
%
% INPUTS
%     E : Young's (elastic) modulus
%     Gs: Shear modulus. For an isotropic material G = E/2(nu+1) with nu the Poission's ratio
%     Kv: Saint-Venant's torsion constant, Polar moment of i
%     L :    Element length
%     A :    Cross section area
%     Mass :    Element mass = rho * A * L [kg]
% NOTE: The following values, may not be equal to E*A and E*I to allow the user to tweak these value.
%     EA  : Young Modulus times Cross section. 
%     EIx : Young Modulus times Polar  second moment of area,local x-axis. Ix=\iint(y^2+z^2) dy dz [m4]
%     EIy : Young Modulus times Planar second moment of area,local y-axis. Iy=\iint z^2 dy dz [m4]
%     EIz : Young Modulus times Planar second moment of area,local z-axis. Iz=\iint y^2 dy dz [m4]
% OPTIONAL INPUTS
%     R   : Transformation matrix from global coord to element coord: x_e = R.x_g 
%          if provided, element matrix is provided in global coord
%
% OUTPUTS
%     ke: Element stiffness matrix (12x12)   
%     me: Element mass matrix (12x12)
%
% AUTHOR: E. Branlard
%
   
 
% --- Stiffness matrix
a =      EA  / L   ; % a1
b = 12 * EIz / L^3 ; % b1
c = 6  * EIz / L^2 ; % b2
d = 12 * EIy / L^3 ; % c1
e = 6  * EIy / L^2 ; % c2
f =     G*Kv / L   ; % a2 = G*J/L
g = 2  * EIy / L   ; % c3
h = 2  * EIz / L   ; % b3

% NOTE: OK with
%          - Serano beam3e function
%          - Matlab FEM Book
%          - frame3d_6j
%          - Panzer-Hubele
% NOTE: compatible with Timoshenko with shear offsets
%    ux1 uy1 uz1 tx1 ty1 tz1   ux2  uy2 yz2 tx2 ty2 tz2
ke = [ a   0  0   0   0   0    -a   0   0   0   0   0  ;
      0   b   0   0   0   c     0  -b   0   0   0   c  ; 
      0   0   d   0  -e   0     0   0  -d   0  -e   0  ; 
      0   0   0   f   0   0     0   0   0  -f   0   0  ; 
      0   0  -e   0  2*g  0     0   0   e   0   g   0  ; 
      0   c   0   0   0  2*h    0  -c   0   0   0   h  ; 

     -a   0   0   0   0   0     a   0   0   0   0   0  ; 
      0  -b   0   0   0  -c     0   b   0   0   0  -c  ; 
      0   0  -d   0   e   0     0   0   d   0   e   0  ; 
      0   0   0  -f   0   0     0   0   0   f   0   0  ; 
      0   0  -e   0   g   0     0   0   e   0  2*g  0  ; 
      0   c   0   0   0   h     0  -c   0   0   0  2*h]; 


% ---
% SOURCE: frame3d_6j.m -  Structural analysis of 3d frame by Ace_ventura  - https://se.mathworks.com/matlabcentral/fileexchange/49559-structural-analysis-of-3d-frames
% NOTE: compatible with above
% GIx=G*Kv;
% ke2=[EA/L      0             0             0      0            0            -EA/L 0             0             0      0            0
%      0         12*EIz/(L^3)  0             0      0            6*EIz/(L^2)  0     -12*EIz/(L^3) 0             0      0            6*EIz/(L^2)
%      0         0             12*EIy/(L^3)  0      -6*EIy/(L^2) 0            0     0             -12*EIy/(L^3) 0      -6*EIy/(L^2) 0
%      0         0             0             GIx/L  0            0            0     0             0             -GIx/L 0            0
%      0         0             -6*EIy/(L^2)  0      4*EIy/L      0            0     0             6*EIy/(L^2)   0      2*EIy/L      0
%      0         6*EIz/(L^2)   0             0      0            4*EIz/L      0     -6*EIz/(L^2)  0             0      0            2*EIz/L
%      -EA/L     0             0             0      0            0            EA/L  0             0             0      0            0
%      0         -12*EIz/(L^3) 0             0      0            -6*EIz/(L^2) 0     12*EIz/(L^3)  0             0      0            -6*EIz/(L^2)
%      0         0             -12*EIy/(L^3) 0      6*EIy/(L^2)  0            0     0             12*EIy/(L^3)  0      6*EIy/(L^2)  0
%      0         0             0             -GIx/L 0            0            0     0             0             GIx/L  0            0
%      0         0             -6*EIy/(L^2)  0      2*EIy/L      0            0     0             6*EIy/(L^2)   0      4*EIy/L      0
%      0         6*EIz/(L^2)   0             0      0            2*EIz/L      0     -6*EIz/(L^2)  0             0      0            4*EIz/L      ]; %formation of element stiffness matrix IN MEMBER AXIS

% ---
% SOURCE: Panzer-Hubele - Generating a Parametric Finite Element Model of a 3D Cantilever Timoshenko Beam Using Matlab
% NOTE: compatible with above
% Py=0; Pz=0; It=Kv; l=L; Iz=EIz/E; Iy=EIy/E;
% K11 = zeros(6,6);
% K11(1,1) = E * A/l ;
% K11(2,2) = 12 * E * Iz/(l^3 * (1+Py)) ;
% K11(2,6) = 6 * E * Iz/(l^2 * (1+Py)) ;
% K11(3,3) = 12 * E * Iy/(l^3 * (1+Pz)) ;
% K11(3,5) = -6 * E * Iy/(l^2 * (1+Pz)) ;
% K11(4,4) = G * It/l ;
% K11(5,5) = (4+Pz) * E * Iy/(l * (1+Pz)) ;
% K11(5,3) = K11(3,5) ;
% K11(6,6) = (4+Py) * E * Iz/(l * (1+Py)) ;
% K11(6,2) = K11(2,6) ;

% K22 = -K11 + 2 * diag(diag(K11));
% K21 = K11 - 2 * diag(diag(K11));
% K21(5,5) = (2-Pz) * E * Iy/(l * (1+Pz)) ;
% K21(6,6) = (2-Py) * E * Iz/(l * (1+Py)) ;
% K21(2,6) = -K21(6,2);
% K21(3,5) = -K21(5,3);
% ke3 = [K11, K21'; K21, K22];



% ---  Mass matrix
% SOURCE: What-When-How-FEM-For-Frames. NOTE: the sign was reveresed in front of 35*r2!!!, to be consistent with Panzer-Hubele with Iy and Iz=0
a  = L/2 ; a2 = a^2 ; r2 = EIx/E/A;
me = Mass/2/105 * [
%ux1  uy1     uz1   tx1     ty1   tz1     ux2    uy2    yz2    tx2    ty2    tz2
 70     0     0       0     0      0      35     0      0       0      0      0    ;  
  0    78     0       0     0   22*a       0    27      0       0      0  -13*a    ;  
  0     0    78       0 -22*a      0       0     0     27       0   13*a      0    ;  
  0     0     0   70*r2     0      0       0     0      0   35*r2      0      0    ;  
  0     0 -22*a       0  8*a2      0       0     0  -13*a       0  -6*a2     0     ;  
  0  22*a     0       0     0   8*a2       0  13*a      0       0      0  -6*a2    ;  
 35     0     0       0     0      0      70     0      0       0      0      0    ;  
  0    27     0       0     0   13*a       0    78      0       0      0  -22*a    ;  
  0     0    27       0 -13*a      0       0     0     78       0   22*a      0    ;  
  0     0     0   35*r2     0      0       0     0      0   70*r2      0      0    ;  
  0     0  13*a       0 -6*a2      0       0     0   22*a       0   8*a2      0    ;  
  0 -13*a     0       0     0  -6*a2       0 -22*a      0       0      0   8*a2   ];  



%% Element in global coord
if exist('R','var')
    RR=blkdiag(R,R,R,R);
    me = RR'*me*RR;
    ke = RR'*ke*RR;
end


% ---
% SOURCE: Panzer-Hubele - Generating a Parametric Finite Element Model of a 3D Cantilever Timoshenko Beam Using Matlab
% Iz=EIz/E; Iy=EIy/E; Ip=EIx/E; l=L;
% M11 = zeros(6,6);
% M11(1,1) = 1/3;
% M11(2,2) = 13/35 + 6 * Iz/(5 * A * l^2);
% M11(3,3) = 13/35 + 6 * Iy/(5 * A * l^2);
% M11(4,4) = Ip/(3 * A);
% M11(5,5) = l^2/105 + 2 * Iy/(15 * A);
% M11(6,6) = l^2/105 + 2 * Iz/(15 * A);
% M11(6,2) = 11 * l/210 + Iz/(10 * A * l);
% M11(2,6) = M11(6,2) ;
% M11(5,3) = -11 * l/210 - Iy/(10 * A * l);
% M11(3,5) = M11(5,3) ;
% M22 = -M11 + 2 * diag(diag(M11));
% M21 = zeros(6,6);
% M21(1,1) = 1/6;
% M21(2,2) = 9/70 - 6 * Iz/(5 * A * l^2);
% M21(3,3) = 9/70 - 6 * Iy/(5 * A * l^2);
% M21(4,4) = Ip/(6 * A);
% M21(5,5) = -l^2/140 - Iy/(30 * A);
% M21(6,6) = -l^2/140 - Iz/(30 * A);
% M21(6,2) = -13 * l/420 + Iz/(10 * A * l);
% M21(2,6) = -M21(6,2);
% M21(5,3) = 13 * l/420 - Iy/(10 * A * l);
% M21(3,5) = -M21(5,3);
% me= Mass * [M11, M21'; M21, M22];
% keyboard
 %
%   b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
%   L=sqrt(b'*b);  n1=b/L;
%   lc=sqrt(eo*eo'); n3=eo/lc;
%  % Elemenr load vector
%      if nargin==5;   eq=[0 0 0 0];  end 
%  eq = [qx qy qz qw];    distributed loads
%      qx=eq(1); qy=eq(2); qz=eq(3); qw=eq(4);
%    fle=L/2*[qx qy qz qw -1/6*qz*L 1/6*qy*L qx qy qz qw 1/6*qz*L -1/6*qy*L]';
% 
%  %
%     n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
%     n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
%     n2(3)=n3(1)*n1(2)-n1(1)*n3(2);
% %
%     An=[n1';
%         n2;
%         n3];
% %
%     Grot=[  An     zeros(3) zeros(3) zeros(3);
%        zeros(3)   An     zeros(3) zeros(3);
%        zeros(3) zeros(3)   An     zeros(3);
%        zeros(3) zeros(3) zeros(3)   An    ];
%  %
%     Ke1=Grot'*Kle*Grot;  fe1=Grot'*fle;

