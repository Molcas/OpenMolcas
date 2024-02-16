!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine susceptibility_pa(exch,nLoc,nCenter,nneq,neqv,neq,nss,nexch,nTempMagn,nT,Tmin,Tmax,XTexp,eso,dipso,s_so,w,dipexch, &
                             s_exch,T,R_LG,zJ,tinput,XLM,ZLM,XRM,ZRM,iopt,chiT_theta,doplot,mem)

use Constants, only: Zero, Three
Use Definitions, only: wp, u6

! chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]
implicit none
integer, intent(in) :: nLoc, nCenter, nTempMagn, nT, mem
integer, intent(in) :: exch, nneq, neqv, iopt
integer, intent(in) :: neq(nneq), nss(nneq), nexch(nneq)
real(kind=8), intent(in) :: T(nT+nTempMagn)
real(kind=8), intent(in) :: W(exch)
real(kind=8), intent(in) :: eso(nneq,nLoc)
real(kind=8), intent(in) :: zJ
real(kind=8), intent(in) :: Tmin, Tmax
real(kind=8), intent(in) :: XTexp(nT+nTempMagn)
real(kind=8), intent(in) :: R_LG(nneq,neqv,3,3)
real(kind=8), intent(out) :: chit_theta(nT+nTempMagn)
! contributions from local excited states, computed in the XT section:
real(kind=8), intent(out) :: XLM(nCenter,nTempMagn,3,3)
real(kind=8), intent(out) :: ZLM(nCenter,nTempMagn)
real(kind=8), intent(out) :: XRM(nCenter,nTempMagn,3,3)
real(kind=8), intent(out) :: ZRM(nCenter,nTempMagn)
logical, intent(in) :: tinput, doplot
! BIG matrices:
complex(kind=8), intent(in) :: dipexch(3,exch,exch)
complex(kind=8), intent(in) :: s_exch(3,exch,exch)
complex(kind=8), intent(in) :: dipso(nneq,3,nLoc,nLoc)
complex(kind=8), intent(in) :: s_so(nneq,3,nLoc,nLoc)
#include "stdalloc.fh"
! local variables
real(kind=8), allocatable :: chit_tens_l(:,:,:)       !chit_tens_l( nneq,3,3)
real(kind=8), allocatable :: smu_chit_tens_l(:,:,:)   !smu_chit_tens_l( nneq,3,3)
real(kind=8), allocatable :: ss_chit_tens_l(:,:,:)    !ss_chit_tens_l( nneq,3,3)
real(kind=8), allocatable :: chit_tens_lr(:,:,:)      !chit_tens_lr(nneq,3,3)
real(kind=8), allocatable :: smu_chit_tens_lr(:,:,:)  !smu_chit_tens_lr(nneq,3,3)
real(kind=8), allocatable :: ss_chit_tens_lr(:,:,:)   !ss_chit_tens_lr(nneq,3,3)
real(kind=8), allocatable :: chit_tens_ex(:,:)        !chit_tens_ex(3,3)
real(kind=8), allocatable :: smu_chit_tens_ex(:,:)    !smu_chit_tens_ex(3,3)
real(kind=8), allocatable :: ss_chit_tens_ex(:,:)     !ss_chit_tens_ex(3,3)
real(kind=8), allocatable :: chit_tens_tot(:,:,:)     !chit_tens_tot(nT+nTempMagn,3,3)
real(kind=8), allocatable :: smu_chit_tens_tot(:,:)   !smu_chit_tens_tot(3,3)
real(kind=8), allocatable :: ss_chit_tens_tot(:,:)    !ss_chit_tens_tot(3,3)
real(kind=8), allocatable :: chit_theta_tens(:,:,:)   !chit_theta_tens(nT+nTempMagn,3,3)
real(kind=8), allocatable :: zstat_l(:)               !zstat_l( nneq)
real(kind=8), allocatable :: zstat_lr(:)              !zstat_lr(nneq)
real(kind=8), allocatable :: zstat_tot(:)             !zstat_tot(nT+nTempMagn)
real(kind=8), allocatable :: chit(:)                  !chit(nT+nTempMagn)
real(kind=8), allocatable :: chi_theta_1(:)           !chi_theta_1(nT+nTempMagn)
real(kind=8), allocatable :: XL(:,:,:)                !XL(nCenter,3,3)
real(kind=8), allocatable :: ZL(:)                    !ZL(nCenter)
real(kind=8), allocatable :: XR(:,:,:)                !XR(nCenter,3,3)
real(kind=8), allocatable :: ZR(:)                    !ZR(nCenter)
real(kind=8), allocatable :: SMUR(:,:,:)              !SMUR(nCenter,3,3)
real(kind=8), allocatable :: SMUL(:,:,:)              !SMUL(nCenter,3,3)
real(kind=8), allocatable :: SSR(:,:,:)               !SSR( nCenter,3,3)
real(kind=8), allocatable :: SSL(:,:,:)               !SSL( nCenter,3,3)
real(kind=8), allocatable :: wt(:), zt(:,:)           !wt(3),zt(3,3)
real(kind=8), allocatable :: A_dir(:,:)               !A_dir(3,3)
real(kind=8), allocatable :: A_inv(:,:)               !A_inv(3,3)
real(kind=8), allocatable :: unity(:,:)               !unity(3,3)
real(kind=8) :: xxm
real(kind=8) :: zstat_ex
real(kind=8) :: boltz_k, coeff_chi
real(kind=8) :: det
real(kind=8) :: dev, Fa, Fb, Fc, Fd, Fe, Ff
external dev
integer :: i, iT, jT, ic, jc
integer :: j, n1, n2, im, jm
integer :: isite, info, mem_local, RtoB
logical :: dbg
character(len=50) :: label
real(wp), external :: dnrm2_

mem_local = 0
dbg = .false.
RtoB = 8
coeff_chi = 0.1250486120_wp*Three ! IFG = n_a*mu_bohr^2/(k_boltz) in cm^3*k/mol
boltz_k = 0.69503560_wp           ! IFG   in cm^-1*k-1
!-----------------------------------------------------------------------
if (dbg) then
  write(u6,*) 'Verification of input data on entrance to PA-SUSC:'
  write(u6,*) 'exch:         ',exch
  write(u6,*) 'nLoc:         ',nLoc
  write(u6,*) 'nCenter:      ',nCenter
  write(u6,*) 'nneq:         ',nneq
  write(u6,*) 'neqv:         ',neqv
  write(u6,*) 'neq():        ',(neq(i),i=1,nneq)
  write(u6,*) 'nss():        ',(nss(i),i=1,nneq)
  write(u6,*) 'nexch():      ',(nexch(i),i=1,nneq)
  write(u6,*) 'nT:           ',nT
  write(u6,*) 'iopt:         ',iopt
  write(u6,*) 'mem:          ',mem
  write(u6,*) 'nTempMagn:    ',nTempMagn
  write(u6,*) 'Tmin:         ',Tmin
  write(u6,*) 'Tmax:         ',Tmax
  write(u6,*) 'XTexp():      ',(XTexp(i),i=1,nT+nTempMagn)
  write(u6,*) 'T():          ',(T(i),i=1,nT+nTempMagn)
  write(u6,*) 'chit_theta(): ',(chit_theta(i),i=1,nT+nTempMagn)
  write(u6,*) 'W()           ',(W(i),i=1,exch)
  write(u6,*) 'zJ:           ',zJ
  write(u6,*) 'tinput:       ',tinput
  write(u6,*) 'doplot:       ',doplot
  !do i=1,nneq
  !  write(u6,*) 'eso()         ',(eso(i,j),j=1,nss(i))
  !  call prMom('SUSC: input  s_so(i,:,:,:):',s_so(i,:,:,:),nexch(i))
  !  call prMom('SUSC: input dipso(i,:,:,:):',dipso(i,:,:,:),nexch(i))
  !  call atens(s_so(i,:,:,:),nexch(i),gtens,maxes,2)
  !  call atens(dipso(i,:,:,:),nexch(i),gtens,maxes,2)
  !End Do
  !call prMom('SUSC: input  S_EXCH(l,i,j):',s_exch,exch)
  !call prMom('SUSC: input DIPEXCH(l,i,j):',dipexch,exch)
  !call atens(s_exch,exch,gtens,maxes,2)
  !call atens(dipexch,exch,gtens,maxes,2)
  ! change to pseudospin:
  !call zcopy_(exch*exch,[cZero],0,Z,1)
  !call zcopy_(3*exch*exch,[cZero],0,dipexch2,1)
  !call zcopy_(3*exch*exch,[cZero],0,s_exch2,1)
  !call zcopy_(3*exch*exch,dipexch,1,dipexch2,1)
  !call zcopy_(3*exch*exch,s_exch,1,s_exch2,1)
  !call pseudospin(dipexch2,exch,Z,3,1)
  !
  !call zcopy_(3*exch*exch,[cZero],0,dipexch,1)
  !call zcopy_(3*exch*exch,[cZero],0,s_exch,1)
  !call UTMU(exch,exch,Z,dipexch2,dipexch)
  !call UTMU(exch,exch,Z,s_exch2,s_exch)
  !
  !call prMom('SUSC: input  S_EXCH2(l,i,j):',s_exch,exch)
  !call prMom('SUSC: input DIPEXCH2(l,i,j):',dipexch,exch)
end if ! dbg
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(100A)') (('%'),J=1,95)
write(u6,'(35X,A)') 'CALCULATION OF THE MAGNETIC SUSCEPTIBILITY'
write(u6,'(100A)') (('%'),J=1,95)
write(u6,*)

if (tinput) then
  write(u6,'(2x,a)') 'Temperature dependence of the magnetic susceptibility and'
  write(u6,'(2x,a)') 'high-field magnetization will be calculated according to '
  write(u6,'(2x,a)') 'experimental values provided in the input.'
else
  write(u6,'(2x,a,i3,a)') 'Temperature dependence of the magnetic susceptibility will be calculated in',nT,' points, '
  write(u6,'(2x,a,f4.1,a,f6.1,a)') 'equally distributed in temperature range ',tmin,' ---',tmax,' K.'
end if
call XFlush(u6)
!-----------------------------------------------------------------------
call mma_allocate(chit_tens_tot,(nT+nTempMagn),3,3,'XT_tens_tot')
call mma_allocate(chit_theta_tens,(nT+nTempMagn),3,3,'XT_theta_t')
call mma_allocate(chit_tens_l,nneq,3,3,'XT_tens_l')
call mma_allocate(chit_tens_lr,nneq,3,3,'XT_tens_lr')
call mma_allocate(chit_tens_ex,3,3,'XT_tens_ex')
call mma_allocate(zstat_l,nneq,'Zstat_l')
call mma_allocate(zstat_lr,nneq,'Zstat_lr')
mem_local = mem_local+(2+2*3*3)*nneq*RtoB+9*RtoB

call mma_allocate(ZR,nCenter,'ZR')
call mma_allocate(ZL,nCenter,'ZL')
call mma_allocate(XL,nCenter,3,3,'XL')
call mma_allocate(XR,nCenter,3,3,'XR')
mem_local = mem_local+(5+2*3*3)*nCenter*RtoB

call mma_allocate(zstat_tot,(nT+nTempMagn),'Zstat_tot')
call mma_allocate(chiT,(nT+nTempMagn),'chiT')
call mma_allocate(chi_theta_1,(nT+nTempMagn),'chi_theta_1')
mem_local = mem_local+(3+2*3*3)*(nT+nTempMagn)*RtoB

call dcopy_(3*3*(nT+nTempMagn),[Zero],0,chit_tens_tot,1)
call dcopy_(3*3*(nT+nTempMagn),[Zero],0,chit_theta_tens,1)
call dcopy_((nT+nTempMagn),[Zero],0,zstat_tot,1)

if (zJ == Zero) then
  if (dbg) write(u6,*) 'SUSC:  memory allocated (local):'
  if (dbg) write(u6,*) 'mem_local=',mem_local
  if (dbg) write(u6,*) 'SUSC:  memory allocated (total):'
  if (dbg) write(u6,*) 'mem_total=',mem+mem_local

  do iT=1,nT+nTempMagn
    ! initialize temporary variables:
    call dcopy_(3*3,[Zero],0,chit_tens_ex,1)
    call dcopy_(nneq*3*3,[Zero],0,chit_tens_l,1)
    call dcopy_(nneq*3*3,[Zero],0,chit_tens_lr,1)
    call dcopy_(nneq,[Zero],0,zstat_l,1)
    call dcopy_(nneq,[Zero],0,zstat_lr,1)
    call dcopy_(nCenter,[Zero],0,ZR,1)
    call dcopy_(nCenter,[Zero],0,ZL,1)
    call dcopy_(3*3*nCenter,[Zero],0,XL,1)
    call dcopy_(3*3*nCenter,[Zero],0,XR,1)
    zstat_ex = Zero
    !-------------------------------------------------------------------
    ! local susceptibility= total susceptibility coming from individual magnetic centers
    do i=1,nneq
      if (dbg) write(u6,'(A,2I5)') 'nss(i)=',nss(i)
      if (dbg) write(u6,'(A,2I5)') 'nexch(i)=',nexch(i)
      if (dbg) write(u6,'(A,9F10.6)') 'eso(i,:)=',eso(i,1:nss(i))
      if (dbg) write(u6,'(A,9F10.6)') 'W(:)    =',W(1:exch)
      if (dbg) write(u6,'(A,9F10.6)') 'T(iT)=',T(iT)
      call chi(dipso(i,1:3,1:nss(i),1:nss(i)),dipso(i,1:3,1:nss(i),1:nss(i)),eso(i,1:nss(i)),nss(i),T(it),zstat_l(i), &
               chit_tens_l(i,1:3,1:3))

      call chi(dipso(i,1:3,1:nexch(i),1:nexch(i)),dipso(i,1:3,1:nexch(i),1:nexch(i)),eso(i,1:nexch(i)),nexch(i),T(iT),zstat_lr(i), &
               chit_tens_lr(i,1:3,1:3))

    end do ! i (nneq)
    call chi(dipexch,dipexch,W,exch,T(iT),zstat_ex,chit_tens_ex)

    Fa = dnrm2_(9*nneq,chit_tens_l,1)
    Fb = dnrm2_(9*nneq,chit_tens_lr,1)
    Fc = dnrm2_(9,chit_tens_ex,1)
    Fd = dnrm2_(nneq,zstat_l,1)
    Fe = dnrm2_(nneq,zstat_lr,1)
    call Add_Info('XT:  chit_tens_l',[Fa],1,6)
    call Add_Info('XT:  chit_tens_lr',[Fb],1,6)
    call Add_Info('XT:  chit_tens_exch',[Fc],1,6)
    call Add_Info('XT:  zstat_exch',[zstat_ex],1,6)
    call Add_Info('XT:  zstat_l',[Fd],1,6)
    call Add_Info('XT:  zstat_lr',[Fe],1,6)
    ! expand the basis and rotate local tensors to the general
    ! coordinate system:
    isite = 0
    do i=1,nneq
      do j=1,neq(i)
        isite = isite+1
        ZL(isite) = zstat_l(i)
        ZR(isite) = zstat_lr(i)
        do ic=1,3
          do jc=1,3
            do n1=1,3
              do n2=1,3
                XL(isite,ic,jc) = XL(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*chit_tens_l(i,n1,n2)
                XR(isite,ic,jc) = XR(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*chit_tens_lr(i,n1,n2)
              end do
            end do
          end do
        end do
      end do
    end do

    Fa = dnrm2_(9*nCenter,XL,1)
    Fb = dnrm2_(9*nCenter,XR,1)
    call Add_Info('XT:  XL',[Fa],1,6)
    call Add_Info('XT:  XR',[Fb],1,6)
    ! save some data:
    if (it <= nTempMagn) then
      call dscal_(3*3*nCenter,coeff_chi,XR,1)
      call dscal_(3*3*nCenter,coeff_chi,XL,1)
      call dcopy_(3*3*nCenter,XR,1,XRM(:,iT,:,:),1)
      call dcopy_(3*3*nCenter,XL,1,XLM(:,iT,:,:),1)
      call dcopy_(nCenter,ZR,1,ZRM(:,iT),1)
      call dcopy_(nCenter,ZL,1,ZLM(:,iT),1)
    end if
    call chi_sum(nCenter,chit_tens_ex,zstat_ex,XL,ZL,XR,ZR,iopt,chit_tens_tot(it,1:3,1:3),zstat_tot(it))

    Fa = dnrm2_(nCenter,ZL,1)
    Fb = dnrm2_(nCenter,ZR,1)
    call Add_Info('XT:  ZL',[Fa],1,6)
    call Add_Info('XT:  ZR',[Fb],1,6)
    call Add_Info('XT: ZEx',[zstat_ex],1,6)

    chit(it) = coeff_chi*(chit_tens_tot(iT,1,1)+chit_tens_tot(iT,2,2)+chit_tens_tot(iT,3,3))/Three
    chit_theta(iT) = chit(iT)

    if (abs(chit(iT)) < 1.0e-20_wp) then
      chit(iT) = 1.0e-20_wp
      chit_theta(iT) = 1.0e-20_wp
    end if
    chi_theta_1(iT) = T(iT)/chit(iT)
    ! add some verification data:
    Fa = dnrm2_(9,chit_theta_tens(iT,1:3,1:3),1)
    call Add_Info('XT: chit_theta_tens',[Fa],1,6)
  end do ! iT
  Fb = dnrm2_(nT+nTempMagn,T,1)
  call Add_Info('XT: T',[Fb],1,6)

else ! i.e. when (zJ /= 0)

  ! allocate memory for temporary arrays:
  call mma_allocate(smu_chit_tens_l,nneq,3,3,'smu_X_tens_l')
  call mma_allocate(smu_chit_tens_lr,nneq,3,3,'smu_X_tens_lr')
  call mma_allocate(smu_chit_tens_ex,3,3,'smu_chit_X_ex')
  call mma_allocate(smu_chit_tens_tot,3,3,'SM_XT')
  call mma_allocate(ss_chit_tens_l,nneq,3,3,'ss_chit_tens_l')
  call mma_allocate(ss_chit_tens_lr,nneq,3,3,'ss_chit_tens_lr')
  call mma_allocate(ss_chit_tens_ex,3,3,'ss_chit_tens_ex')
  call mma_allocate(ss_chit_tens_tot,3,3,'SS_XT')
  call mma_allocate(SMUR,nCenter,3,3,'SMUR')
  call mma_allocate(SMUL,nCenter,3,3,'SMUL')
  call mma_allocate(SSR,nCenter,3,3,'SSR')
  call mma_allocate(SSL,nCenter,3,3,'SSL')
  call mma_allocate(A_dir,3,3,'A_dir')
  call mma_allocate(A_inv,3,3,'A_inv')
  call mma_allocate(Unity,3,3,'Unity')
  mem_local = mem_local+4*3*3*(nCenter+nneq)*RtoB
  mem_local = mem_local+7*3*3*RtoB

  if (dbg) write(u6,*) 'SUSC:  memory allocated (local):'
  if (dbg) write(u6,*) 'mem_local=',mem_local
  if (dbg) write(u6,*) 'SUSC:  memory allocated (total):'
  if (dbg) write(u6,*) 'mem_total=',mem+mem_local

  do iT=1,nT+nTempMagn
    ! initialization:
    call dcopy_(3*3,[Zero],0,chit_tens_ex,1)
    call dcopy_(3*3,[Zero],0,smu_chit_tens_ex,1)
    call dcopy_(3*3,[Zero],0,ss_chit_tens_ex,1)
    call dcopy_(nneq*3*3,[Zero],0,chit_tens_l,1)
    call dcopy_(nneq*3*3,[Zero],0,smu_chit_tens_l,1)
    call dcopy_(nneq*3*3,[Zero],0,ss_chit_tens_l,1)
    call dcopy_(nneq*3*3,[Zero],0,chit_tens_lr,1)
    call dcopy_(nneq*3*3,[Zero],0,smu_chit_tens_lr,1)
    call dcopy_(nneq*3*3,[Zero],0,ss_chit_tens_lr,1)
    call dcopy_(3*3,[Zero],0,smu_chit_tens_tot,1)
    call dcopy_(3*3,[Zero],0,ss_chit_tens_tot,1)
    call dcopy_(nneq,[Zero],0,zstat_l,1)
    call dcopy_(nneq,[Zero],0,zstat_lr,1)
    call dcopy_(nCenter,[Zero],0,ZR,1)
    call dcopy_(nCenter,[Zero],0,ZL,1)
    call dcopy_(nCenter*3*3,[Zero],0,XL,1)
    call dcopy_(nCenter*3*3,[Zero],0,XR,1)
    call dcopy_(nCenter*3*3,[Zero],0,SMUR,1)
    call dcopy_(nCenter*3*3,[Zero],0,SMUL,1)
    call dcopy_(nCenter*3*3,[Zero],0,SSR,1)
    call dcopy_(nCenter*3*3,[Zero],0,SSL,1)
    call dcopy_(3*3,[Zero],0,A_dir,1)
    call dcopy_(3*3,[Zero],0,A_inv,1)
    zstat_ex = Zero
    det = Zero
    call unitmat(unity,3)

    ! compute local tensors  L, and LR:
    do i=1,nneq
      call chi(dipso(i,1:3,1:nss(i),1:nss(i)),dipso(i,1:3,1:nss(i),1:nss(i)),eso(i,1:nss(i)),nss(i),T(iT),zstat_l(i), &
               chit_tens_l(i,1:3,1:3))

      call chi(s_so(i,1:3,1:nss(i),1:nss(i)),dipso(i,1:3,1:nss(i),1:nss(i)),eso(i,1:nss(i)),nss(i),T(iT),zstat_l(i), &
               smu_chit_tens_l(i,1:3,1:3))

      call chi(s_so(i,1:3,1:nss(i),1:nss(i)),s_so(i,1:3,1:nss(i),1:nss(i)),eso(i,1:nss(i)),nss(i),T(iT),zstat_l(i), &
               ss_chit_tens_l(i,1:3,1:3))

      call chi(dipso(i,1:3,1:nexch(i),1:nexch(i)),dipso(i,1:3,1:nexch(i),1:nexch(i)),eso(i,1:nexch(i)),nexch(i),T(iT),zstat_lr(i), &
               chit_tens_lr(i,1:3,1:3))

      call chi(s_so(i,1:3,1:nexch(i),1:nexch(i)),dipso(i,1:3,1:nexch(i),1:nexch(i)),eso(i,1:nexch(i)),nexch(i),T(iT),zstat_lr(i), &
               smu_chit_tens_lr(i,1:3,1:3))

      call chi(s_so(i,1:3,1:nexch(i),1:nexch(i)),s_so(i,1:3,1:nexch(i),1:nexch(i)),eso(i,1:nexch(i)),nexch(i),T(iT),zstat_lr(i), &
               ss_chit_tens_lr(i,1:3,1:3))
    end do ! i (nneq)

    ! compute exchange tensors:
    call chi(dipexch,dipexch,W,exch,T(it),zstat_ex,chit_tens_ex)
    call chi(s_exch,dipexch,W,exch,T(iT),zstat_ex,smu_chit_tens_ex)
    call chi(s_exch,s_exch,W,exch,T(iT),zstat_ex,ss_chit_tens_ex)
    ! expand the basis and rotate local tensors to the general
    ! coordinate system:
    isite = 0
    do i=1,nneq
      do j=1,neq(i)
        isite = isite+1
        ZL(isite) = zstat_l(i)
        ZR(isite) = zstat_lr(i)
        ! use R_lg matrices, which have arbitrary determinant:  +1  or -1;
        ! reason:  X_ab is a bi-dimensional tensor.
        !          We need to rotate twice ==> the sign of R does not
        !          matter (+1 * +1) = (-1 * -1)
        !  >> R_rot matrices have determinant strict +1, and are used to
        !          rotate vectors
        do ic=1,3
          do jc=1,3
            do n1=1,3
              do n2=1,3
                XR(isite,ic,jc) = XR(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*chit_tens_lr(i,n1,n2)
                XL(isite,ic,jc) = XL(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*chit_tens_l(i,n1,n2)

                SMUL(isite,ic,jc) = SMUL(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*smu_chit_tens_l(i,n1,n2)
                SMUR(isite,ic,jc) = SMUR(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*smu_chit_tens_lr(i,n1,n2)

                SSL(isite,ic,jc) = SSL(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*ss_chit_tens_l(i,n1,n2)
                SSR(isite,ic,jc) = SSR(isite,ic,jc)+r_lg(i,j,ic,n1)*r_lg(i,j,jc,n2)*ss_chit_tens_lr(i,n1,n2)
              end do
            end do
          end do
        end do
      end do
    end do

    ! add verification data:
    Fa = dnrm2_(9*nCenter,XR,1)
    Fb = dnrm2_(9*nCenter,XL,1)
    Fc = dnrm2_(9*nCenter,SMUL,1)
    Fd = dnrm2_(9*nCenter,SMUR,1)
    Fe = dnrm2_(9*nCenter,SSL,1)
    Ff = dnrm2_(9*nCenter,SSR,1)
    call Add_Info('XT:    XR',[Fa],1,6)
    call Add_Info('XT:    XL',[Fb],1,6)
    call Add_Info('XT:  SMUL',[Fc],1,6)
    call Add_Info('XT:  SMUR',[Fd],1,6)
    call Add_Info('XT:   SSL',[Fe],1,6)
    call Add_Info('XT:   SSR',[Ff],1,6)

    ! save some data:
    if (iT <= nTempMagn) then
      call dscal_(3*3*nCenter,coeff_chi,XR,1)
      call dscal_(3*3*nCenter,coeff_chi,XL,1)
      call dcopy_(3*3*nCenter,XR,1,XRM(:,iT,:,:),1)
      call dcopy_(3*3*nCenter,XL,1,XLM(:,iT,:,:),1)
      call dcopy_(nCenter,ZR,1,ZRM(:,iT),1)
      call dcopy_(nCenter,ZL,1,ZLM(:,iT),1)
    end if

    ! compute total tensors:
    call chi_sum(nCenter,chit_tens_ex,zstat_ex,XL,ZL,XR,ZR,iopt,chit_tens_tot(iT,1:3,1:3),zstat_tot(iT))

    call chi_sum(nCenter,smu_chit_tens_ex,zstat_ex,SMUL,ZL,SMUR,ZR,iopt,smu_chit_tens_tot,zstat_tot(iT))

    call chi_sum(nCenter,ss_chit_tens_ex,zstat_ex,SSL,ZL,SSR,ZR,iopt,ss_chit_tens_tot,zstat_tot(iT))

    ! form the A_dir matrix:
    ! A_dir(:,:) = 1(:,:)* kB * T(iT) - zJ*XSS(:,:)
    call daxpy_(3*3,Boltz_k*T(iT),Unity,1,A_dir,1)
    call daxpy_(3*3,-zJ,ss_chit_tens_tot,1,A_dir,1)
    ! invert it:
    call REVERSE(A_dir,A_inv,DET)
    do ic=1,3
      do jc=1,3
        xxm = Zero
        do im=1,3
          do jm=1,3
            xxm = xxm+smu_chit_tens_tot(im,ic)*a_inv(im,jm)*smu_chit_tens_tot(jm,jc)
          end do
        end do
        chit_theta_tens(iT,ic,jc) = chit_tens_tot(iT,ic,jc)+zj*xxm
      end do ! jc
    end do ! ic

    chit(iT) = coeff_chi*(chit_tens_tot(iT,1,1)+chit_tens_tot(iT,2,2)+chit_tens_tot(iT,3,3))/Three

    chit_theta(iT) = coeff_chi*(chit_theta_tens(iT,1,1)+chit_theta_tens(iT,2,2)+chit_theta_tens(iT,3,3))/Three
    if (abs(chit(iT)) < 1.0e-20_wp) then
      chit(iT) = 1.0e-20_wp
      chit_theta(iT) = 1.0e-20_wp
    end if
    if (abs(chit_theta(iT)) < 1.0e-20_wp) then
      chit_theta(iT) = 1.0e-20_wp
    end if

    chi_theta_1(iT) = t(iT)/chit_theta(iT)

    ! add some verification data:
    Fa = dnrm2_(9,chit_tens_tot(it,1:3,1:3),1)
    Fb = dnrm2_(9,chit_theta_tens(it,1:3,1:3),1)
    Fc = dnrm2_(9,smu_chit_tens_tot(1:3,1:3),1)
    Fd = dnrm2_(9,ss_chit_tens_tot(1:3,1:3),1)

    call Add_Info('XT: chit_tens_tot',[Fa],1,6)
    call Add_Info('XT: chit_theta_tens',[Fb],1,6)
    call Add_Info('XT: smu_chit_tens_tot',[Fc],1,6)
    call Add_Info('XT:  ss_chit_tens_tot',[Fd],1,6)
  end do ! it
  Fb = dnrm2_(nT+nTempMagn,T,1)
  call Add_Info('XT: T',[Fb],1,6)

  call mma_deallocate(smu_chit_tens_l)
  call mma_deallocate(smu_chit_tens_lr)
  call mma_deallocate(smu_chit_tens_ex)
  call mma_deallocate(smu_chit_tens_tot)
  call mma_deallocate(ss_chit_tens_l)
  call mma_deallocate(ss_chit_tens_lr)
  call mma_deallocate(ss_chit_tens_ex)
  call mma_deallocate(ss_chit_tens_tot)
  call mma_deallocate(SMUR)
  call mma_deallocate(SMUL)
  call mma_deallocate(SSR)
  call mma_deallocate(SSL)
  call mma_deallocate(A_dir)
  call mma_deallocate(A_inv)
  call mma_deallocate(Unity)

end if  !zJ
!-----------------------------------------------------------------------
! printing the results
do iT=1,nT+nTempMagn
  do ic=1,3
    do jc=1,3
      chit_tens_tot(iT,ic,jc) = coeff_chi*chit_tens_tot(iT,ic,jc)
      if (zJ /= Zero) chit_theta_tens(iT,ic,jc) = coeff_chi*chit_theta_tens(iT,ic,jc)
    end do
  end do
end do
write(u6,*)
write(u6,'(A)') '----------------------------------------------------------------------------------------|'
write(u6,'(A)') '     |     T      | Statistical |    CHI*T    |    CHI*T    |     CHI     |    1/CHI    |'
write(u6,'(A)') '     |            |  Sum (Z)    |    (zJ=0)   |             |             |             |'
write(u6,'(A)') '-----|----------------------------------------------------------------------------------|'
write(u6,'(A)') 'Units|   kelvin   |    ---      |  cm3*K/mol  |  cm3*K/mol  |   cm3/mol   |   mol/cm3   |'
write(u6,'(A)') '-----|----------------------------------------------------------------------------------|'

do iT=1,nT
  jT = iT+nTempMagn
  write(u6,'(A,F11.6,A,ES12.5,A,F12.8,A,F12.8,A,ES12.5,A,ES12.5,A)') '     |',T(jT),' |',zstat_tot(jT),' |',chiT(jT),' |', &
                                                                     chiT_theta(jT),' |',chit_theta(jT)/T(jT),' |', &
                                                                     chi_theta_1(jT),' |'
end do
write(u6,'(A)') '-----|----------------------------------------------------------------------------------|'
Fb = dnrm2_(nT+nTempMagn,chiT,1)
call Add_Info('XT: T',[Fb],1,6)
Fa = dnrm2_(nT+nTempMagn,chiT_theta,1)
call Add_Info('XT: CHIT_THETA',[Fa],1,6)
Fa = dnrm2_(nT+nTempMagn,zstat_tot,1)
call Add_Info('XT: CHIT_THETA',[Fa],1,6)
!  calcualtion of the standard deviation:
if (tinput) &
  write(u6,'(a,5x, f20.14)') 'ST.DEV.CHIT:',dev(nT,chit_theta((1+nTempMagn):(nT+nTempMagn)),XTexp((1+nTempMagn):(nT+nTempMagn)))

!-------------------------  PLOTs -------------------------------------!
write(label,'(A)') 'no_field'
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T((1+nTempMagn):(nT+nTempMagn)),chit_theta((1+nTempMagn):(nT+nTempMagn)), &
                          XTexp((1+nTempMagn):(nT+nTempMagn)))
  else
    call plot_XT_no_Exp(label,nT,T((1+nTempMagn):(nT+nTempMagn)),chit_theta((1+nTempMagn):(nT+nTempMagn)))
  end if
end if
!---------------------- END PLOTs -------------------------------------!

! print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
call mma_allocate(wt,3,'wt')
call mma_allocate(zt,3,3,'zt')
if (zJ == Zero) then
  write(u6,'(/)')
  write(u6,'(111A)') ('-',i=1,110),'|'
  write(u6,'(31X,A,22x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ = 0,  in cm3*K/mol','|'
  write(u6,'(111A)') ('-',i=1,110),'|'
  write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
  do iT=1,nT
    jT = iT+nTempMagn
    info = 0
    call dcopy_(3,[Zero],0,wt,1)
    call dcopy_(3*3,[Zero],0,zt,1)
    call DIAG_R2(chit_tens_tot(jT,:,:),3,info,wt,zt)
    write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- '// &
                    'z ----|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(chit_tens_tot(jT,1,j),j=1,3),' |  X:',wt(1),'|', &
                                                    (zt(j,1),j=1,3),'|'
    write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(chiT_tens_tot(jT,2,j),j=1,3),' |  Y:',wt(2),'|', &
                                                             (zt(j,2),j=1,3),'|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(chiT_tens_tot(jT,3,j),j=1,3),' |  Z:',wt(3),'|', &
                                                    (zt(j,3),j=1,3),'|'
  end do
  write(u6,'(111A)') ('-',i=1,110),'|'

else ! zJ /= Zero

  write(u6,'(/)')
  write(u6,'(111A)') ('-',i=1,110),'|'

  write(u6,'(31X,A,F9.6,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ =',zJ,',  in cm3*K/mol','|'
  write(u6,'(111A)') ('-',i=1,110),'|'
  write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
  do iT=1,nT
    jT = iT+nTempMagn
    info = 0
    call dcopy_(3,[Zero],0,wt,1)
    call dcopy_(3*3,[Zero],0,zt,1)
    call DIAG_R2(chit_theta_tens(jT,:,:),3,info,wt,zt)
    write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- '// &
                    'z ----|'

    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(chit_theta_tens(jT,1,j),j=1,3),' |  X:',wt(1),'|', &
                                                    (zt(j,1),j=1,3),'|'
    write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(chiT_theta_tens(jT,2,j),j=1,3),' |  Y:',wt(2),'|', &
                                                             (zt(j,2),j=1,3),'|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(chiT_theta_tens(jT,3,j),j=1,3),' |  Z:',wt(3),'|', &
                                                    (zt(j,3),j=1,3),'|'
  end do
  write(u6,'(111A)') ('-',i=1,110),'|'
end if
call mma_deallocate(wt)
call mma_deallocate(zt)

call mma_deallocate(chit_tens_tot)
call mma_deallocate(chit_theta_tens)
call mma_deallocate(chit_tens_l)
call mma_deallocate(chit_tens_lr)
call mma_deallocate(chit_tens_ex)
call mma_deallocate(zstat_l)
call mma_deallocate(zstat_lr)
call mma_deallocate(zstat_tot)
call mma_deallocate(ZR)
call mma_deallocate(ZL)
call mma_deallocate(XL)
call mma_deallocate(XR)
call mma_deallocate(chiT)
call mma_deallocate(chi_theta_1)

Go To 190
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(5x,a)') 'on user request, the magnetic susceptibility and '
write(u6,'(5x,a)') 'the magnetic susceptibility tensor were not calculated.'
write(u6,*)
190 continue

return

end subroutine susceptibility_pa
