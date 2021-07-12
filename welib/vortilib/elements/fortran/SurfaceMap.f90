!> 
module SurfaceMap
    use PrecisionMod, only: MK
    implicit none
            
contains
            
    ! function ue=fSurfaceIntegral(vX0,vY0,vZ0,Ux,Uy,Uz,n,XCP,YCP,ZCP)
    ! Compute the integral of the source and doublet from a given "Velocity on a Face" using point approximation
    subroutine fUi_PointSourceVortex_surface_mesh(P0,Sigmas,Gammas,Mesh,h)
        !use MathConstants, only: fourpi_inv
        use MeshTypes, only: T_Mesh
        use UIBlobs, only: fUi_Particle_11
        use UIPointSource, only: fUi_PointSource_11
        ! Arguments
        type(T_Mesh), intent(inout)         :: Mesh
        real(MK), dimension(:,:),intent(in) :: P0     !< Position of sources and doublet
        real(MK), dimension(:),intent(in)   :: Sigmas
        real(MK), dimension(:,:),intent(in) :: Gammas
        real(MK),intent(in)                 :: h
        ! Variables
        integer :: im1,im2,im3
        integer :: nm1,nm2,nm3
        integer :: p
        real(MK), dimension(3) :: CP
        real(MK), dimension(3) :: r
        real(MK), dimension(3) :: us
        real(MK), dimension(3) :: uv
        real(MK), dimension(3) :: Uind_cum
        real(MK), dimension(9) :: Grad_cum
        real(MK), dimension(3) :: Uind_tmp
        real(MK), dimension(9) :: Grad_tmp
        real(MK) :: r_norm
        real(MK) :: r3_inv
        real(MK) :: r2_inv
        real(MK) :: factor
        integer :: isd,nsd
        integer :: SmoothModel
        real(MK) :: SmoothParam

        ! Init
        nm1=Mesh%GP%n1
        nm2=Mesh%GP%n2
        nm3=Mesh%GP%n3
        nsd=size(P0,2)

        SmoothParam=1.5*h;

        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP DO PRIVATE(p,im1,im2,im3,CP,Uind_cum,Grad_cum,Uind_tmp,Grad_tmp,&
        !$OMP& isd,r,r_norm,r3_inv,r2_inv,factor,SmoothModel,us,uv&
        !$OMP& ) schedule(runtime)
        ! Loop on CPs is unrolled
        do p=1,Mesh%nCPs
            ! Mesh indexes
            im1=mod((p-1),nm1)+1
            im2=mod((p-1)/nm1,nm2)+1
            im3=(p-1)/(nm1*nm2)+1
            !
            CP(1:3)=(/Mesh%GP%v1(im1),Mesh%GP%v2(im2),Mesh%GP%v3(im3)/)
            Uind_cum(1:3)=0.0_MK
            Grad_cum(1:9)=0.0_MK

            ! The 'if' below was a safety to avoid the points next to the singularities
            !if(im1>1 .and. im2>1 .and. im3>1 .and. im1<nm1 .and. im2<nm2 .and. im3<nm3) then
                ! Loop on Sources, unrolled
                do isd=1,nsd
                    ! ue=-(un(i01,i02)*r+cross(squeeze(ut(i01,i02,1:3)),r))/(4*pi*r_norm^3)*dS;
                    r(1:3)      = CP(1:3)-P0(1:3,isd)
                    r_norm = sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))

                    !
                    SmoothModel=1;
                    Uind_tmp(1:3)=0.0_MK
                    Grad_tmp(1:9)=0.0_MK
                    call fUi_PointSource_11(r,Sigmas(isd), SmoothModel, SmoothParam,.true.,Uind_tmp, Grad_tmp)
                    Uind_cum(1:3)= Uind_cum(1:3)+Uind_tmp(1:3)
                    Grad_cum(1:9)= Grad_cum(1:9)+Grad_tmp(1:9)
                    ! 
                    SmoothModel=2;
                    Uind_tmp(1:3)=0.0_MK
                    Grad_tmp(1:9)=0.0_MK
                    call fUi_Particle_11(r, Gammas(1:3,isd), SmoothModel, SmoothParam, .true., Uind_tmp, Grad_tmp)
                    Uind_cum(1:3)= Uind_cum(1:3)+Uind_tmp(1:3)
                    Grad_cum(1:9)= Grad_cum(1:9)+Grad_tmp(1:9)
!                     if (r_norm>h) then
!                         us(1) = Sigmas(isd)*r(1)  ! 1/(4*pi*r_norm^3)
!                         us(2) = Sigmas(isd)*r(2)  ! 1/(4*pi*r_norm^3)
!                         us(3) = Sigmas(isd)*r(3)  ! 1/(4*pi*r_norm^3)
!                         uv(1) = Gammas(2,isd) * r(3) - Gammas(3,isd) * r(2)
!                         uv(2) = Gammas(3,isd) * r(1) - Gammas(1,isd) * r(3)
!                         uv(3) = Gammas(1,isd) * r(2) - Gammas(2,isd) * r(1)
!                         r3_inv=1._MK/(r_norm**3);
!                         factor=fourpi_inv*r3_inv !1/(4 pi r^3)
!                         Uind_cum(1:3)= Uind_cum(1:3)+(us(1:3)+uv(1:3))*factor
! 
!                         r2_inv=1._MK/(r_norm**2);
!                         Grad_cum(1) = Grad_cum(1)+factor*Sigmas(isd)*(1 -3*r(1)*r(1)*r2_inv) ! 11
!                         Grad_cum(2) = Grad_cum(2)+factor*Sigmas(isd)*(  -3*r(2)*r(1)*r2_inv) ! 12
!                         Grad_cum(3) = Grad_cum(3)+factor*Sigmas(isd)*(  -3*r(3)*r(1)*r2_inv) ! 13
!                                                         
!                         Grad_cum(4) = Grad_cum(4)+factor*Sigmas(isd)*(  -3*r(1)*r(2)*r2_inv) ! 21
!                         Grad_cum(5) = Grad_cum(5)+factor*Sigmas(isd)*(1 -3*r(2)*r(2)*r2_inv) ! 22
!                         Grad_cum(6) = Grad_cum(6)+factor*Sigmas(isd)*(  -3*r(3)*r(2)*r2_inv) ! 23
!                                                         
!                         Grad_cum(7) = Grad_cum(7)+factor*Sigmas(isd)*(  -3*r(1)*r(3)*r2_inv) ! 31
!                         Grad_cum(8) = Grad_cum(8)+factor*Sigmas(isd)*(  -3*r(2)*r(3)*r2_inv) ! 32
!                         Grad_cum(9) = Grad_cum(9)+factor*Sigmas(isd)*(1 -3*r(3)*r(3)*r2_inv) ! 33
! 
!                     endif
                enddo
                Mesh%Values(1:3,im1,im2,im3)  = Mesh%Values(1:3,im1,im2,im3)+Uind_cum(1:3)
                Mesh%Values(4:12,im1,im2,im3) = Mesh%Values(4:12,im1,im2,im3)+Grad_cum(1:9)
            !endif
        enddo
        !$OMP END DO 
        !$OMP END PARALLEL
    end subroutine
end module SurfaceMap
