!> 
subroutine fUi_CylinderTangSemiInfinite(Xcp, Ycp, Zcp, gamma_t, R, z0, ui)
implicit none 
! Use somemodule 
! Arguments declarations 
real*8, intent(in) :: gamma_t  !<  
real*8, intent(in) :: R  !<  
real*8, dimension(:,:), intent(out) :: ui  !< ! m2f:check dim(ncp, 
real*8, intent(in) :: Xcp  !<  
real*8, intent(in) :: Ycp  !<  
real*8, intent(in) :: z0  !<  
real*8, intent(in) :: Zcp  !<  
! Variable declarations 
real*8 :: EE  !<  
real*8 :: EPSILON  !<  
integer, dimension(:) :: Iz  !< !m2f: check dim(:) 
real*8 :: k  !<  
real*8 :: k0_2  !<  
real*8 :: k_2  !<  
real*8 :: KK  !<  
integer :: ncp  !<  
real*8, dimension(:) :: PI  !< ! m2f:check dim(1,si 
real*8, dimension(:) :: PI2  !< ! m2f:check dim(1,si 
real*8 :: r  !<  
real*8 :: z  !<  
! 
! Induced velocity from a semi infinite cylinder extending along the z axis, starting at z=z0 
! 
! INPUTS ARE FLAT 
! gamma_t: tangential vorticity sheet strength of the cylinder 
! R: cylinder radius 
! z0: origin of the cylinder 
! Came from the need of an optimized version, easy to use, takes flat vector and returns flat vector. 
! 
! If mesh needed, use the following before and after 
! --------------------------------------------------- 
! [Xcp Ycp Zcp ncp Xcp Ycp Zcp] = fgetControlPoints( v1,v2,v3,bPolarIn,1 ); 
! [ui]=fUi_Ring1_Opt(Xcp,Ycp,Zcp,Gamma,R,z) 
! ui=fReshape(ui,length(v1),length(v2),length(v3),3); 

! 
! it returns "polar" velocity so use the following to go to cartesian before reshaping 
! ------------------------------------------------- 
! phi=atan2(Ycp,Xcp); ! TODO should make use of input!!!! if polar in... 
! urad=ui(:,1); 
! upsi=ui(:,2); 
! ui(:,1)=upsi.*sin(phi)+urad.*cos(phi); 
! ui(:,2)=upsi.*cos(phi)-urad.*sin(phi); 

EPSILON=1e-7 ! threshold for using axis formula

ncp=size(Xcp)
!m2f: ui=zeros(ncp,3)
if (allocated(ui)) deallocate(ui)
allocate(ui(ncp,3))
ui = 0.0D0 


!! Vectorialization 
r=sqrt(Xcp.**2+Ycp.**2)
z=Zcp-z0 ! !!!! removing z0


! Eliptic integrals 
k_2 = 4*r*R./((R+r).**2+z.**2)
k = sqrt(k_2)
k0_2 = 4*r*R./((R+r).**2)
call ellipke( k_2 , KK,EE) 
! [PI] = ellipticPiManu(k0_2,k_2) ; 
!m2f: PI=zeros(1,size(k0_2))
if (allocated(PI)) deallocate(PI)
allocate(PI(1,size(k0_2)))
PI = 0.0D0 

!m2f: PI2=zeros(1,size(k0_2))
if (allocated(PI2)) deallocate(PI2)
allocate(PI2(1,size(k0_2)))
PI2 = 0.0D0 

do i=1,size(k0_2) 
disp([k0_2(i) k_2(i)])
PI(i) = ellipticPiManu(k0_2(i),k_2(i))
PI2(i) = ellipticPi(k0_2(i),k_2(i))
disp([PI(i) PI2(i) PI(i)-PI2(i)])
end do 


! ur 
ui(:,1)=-gamma_t/(2*pi)*sqrt(R./r).*( (2-k_2)./k.*KK - 2./k.*EE )
! ur 
ui(:,3)=gamma_t/2*( (R-r+abs(R-r))./(2*abs(R-r))+z.*k./(2*pi*sqrt(r*R)).*(KK + (R-r)./(R+r).*PI ) )

! T1=gamma_t/2*( (R-r+abs(R-r))./(2*abs(R-r))); 
! T2=gamma_t/2*( z.*k./(2*pi*sqrt(r*R)).*(KK )); 
! T3=gamma_t/2*( z.*k./(2*pi*sqrt(r*R)).*((R-r)./(R+r).*PI ) ); 


!! Using Axis formula : v_z=-Gamma/(2) *( 1 + z / sqrt(R^2+z^2)) 
Iz=r<EPSILON
ui(Iz,1)=0
ui(Iz,3)=gamma_t/(2)*(1 + z(Iz)./sqrt(z(Iz).**2+R**2)) ! no minus sign... it's a matter of convention




end subroutine fUi_CylinderTangSemiInfinite 
