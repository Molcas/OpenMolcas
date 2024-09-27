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
! chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Three, Ten, cm_s, hPlanck, kBoltzmann, mBohr, rNAVO
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: exch, nLoc, nCenter, nneq, neqv, neq(nneq), nss(nneq), nexch(nneq), nTempMagn, nT, iopt, mem
real(kind=wp), intent(in) :: Tmin, Tmax, XTexp(nT+nTempMagn), eso(nneq,nLoc), W(exch), T(nT+nTempMagn), R_LG(nneq,neqv,3,3), zJ
complex(kind=wp), intent(in) :: dipso(nneq,3,nLoc,nLoc), s_so(nneq,3,nLoc,nLoc), dipexch(3,exch,exch), s_exch(3,exch,exch)
logical(kind=iwp), intent(in) :: tinput, doplot
real(kind=wp), intent(out) :: XLM(nCenter,nTempMagn,3,3), ZLM(nCenter,nTempMagn), XRM(nCenter,nTempMagn,3,3), &
                              ZRM(nCenter,nTempMagn), chit_theta(nT+nTempMagn)
integer(kind=iwp) :: i, ic, im, info, isite, iT, j, jc, jm, jT, mem_local, n1, n2
real(kind=wp) :: A_dir(3,3), A_inv(3,3), chit_tens_ex(3,3), det, Fa, Fb, Fc, Fd, Fe, Ff, smu_chit_tens_ex(3,3), &
                 smu_chit_tens_tot(3,3), ss_chit_tens_ex(3,3), ss_chit_tens_tot(3,3), unity(3,3), xxm, zstat_ex
character(len=50) :: label
real(kind=wp), allocatable :: chi_theta_1(:), chit(:), chit_tens_l(:,:,:), chit_tens_lr(:,:,:), chit_tens_tot(:,:,:), &
                              chit_theta_tens(:,:,:), eso_tmp(:), smu_chit_tens_l(:,:,:), smu_chit_tens_lr(:,:,:), SMUL(:,:,:), &
                              SMUR(:,:,:), ss_chit_tens_l(:,:,:), ss_chit_tens_lr(:,:,:), SSL(:,:,:), SSR(:,:,:), wt(:), zt(:,:), &
                              XL(:,:,:), XR(:,:,:), ZL(:), ZR(:), zstat_l(:), zstat_lr(:), zstat_tot(:)
complex(kind=wp), allocatable :: dipso_tmp(:,:,:), s_so_tmp(:,:,:)
real(kind=wp), parameter :: boltz_k = kBoltzmann/(cm_s*hPlanck), & ! in cm^-1*k-1
                            coeff_chi = rNAVO*mBohr**2/kBoltzmann/Ten ! = n_a*mu_bohr^2/(k_boltz) in cm^3*k/mol
real(kind=wp), external :: dev, dnrm2_

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(mem)
#endif

mem_local = 0
!-----------------------------------------------------------------------
#ifdef _DEBUGPRINT_
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
!call mma_allocate(s_so_tmp,3,nLoc,nLoc,label='s_so_tmp')
!call mma_allocate(dipso_tmp,3,nLoc,nLoc,label='dipso_tmp')
!do i=1,nneq
!  write(u6,*) 'eso()         ',(eso(i,j),j=1,nss(i))
!  s_so_tmp(:,:,:) = s_so(i,:,:,:)
!  dipso_tmp(:,:,:) = dipso(i,:,:,:)
!  call prMom('SUSC: input  s_so(i,:,:,:):',s_so_tmp,nexch(i))
!  call prMom('SUSC: input dipso(i,:,:,:):',dipso_tmp,nexch(i))
!  call atens(s_so_tmp,nexch(i),gtens,maxes,2)
!  call atens(dipso_tmp,nexch(i),gtens,maxes,2)
!end do
!call mma_deallocate(s_so_tmp)
!call mma_deallocate(dipso_tmp)
!call prMom('SUSC: input  S_EXCH(l,i,j):',s_exch,exch)
!call prMom('SUSC: input DIPEXCH(l,i,j):',dipexch,exch)
!call atens(s_exch,exch,gtens,maxes,2)
!call atens(dipexch,exch,gtens,maxes,2)
! change to pseudospin:
!dipexch2(:,:,:) = dipexch(:,:,:)
!exch2(:,:,:) = s_exch(:,:,:)
!call pseudospin(dipexch2,exch,Z,3,1)
!
!call UTMU(exch,exch,Z,dipexch2,dipexch)
!call UTMU(exch,exch,Z,s_exch2,s_exch)
!
!call prMom('SUSC: input  S_EXCH2(l,i,j):',s_exch,exch)
!call prMom('SUSC: input DIPEXCH2(l,i,j):',dipexch,exch)
#endif
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(A)') repeat('%',95)
write(u6,'(35X,A)') 'CALCULATION OF THE MAGNETIC SUSCEPTIBILITY'
write(u6,'(A)') repeat('%',95)
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
call mma_allocate(chit_tens_tot,3,3,(nT+nTempMagn),'XT_tens_tot')
call mma_allocate(chit_theta_tens,3,3,(nT+nTempMagn),'XT_theta_t')
call mma_allocate(chit_tens_l,3,3,nneq,'XT_tens_l')
call mma_allocate(chit_tens_lr,3,3,nneq,'XT_tens_lr')
call mma_allocate(zstat_l,nneq,'Zstat_l')
call mma_allocate(zstat_lr,nneq,'Zstat_lr')
mem_local = mem_local+size(chit_tens_tot)*RtoB
mem_local = mem_local+size(chit_theta_tens)*RtoB
mem_local = mem_local+size(chit_tens_l)*RtoB
mem_local = mem_local+size(chit_tens_lr)*RtoB
mem_local = mem_local+size(zstat_l)*RtoB
mem_local = mem_local+size(zstat_lr)*RtoB

call mma_allocate(ZR,nCenter,'ZR')
call mma_allocate(ZL,nCenter,'ZL')
call mma_allocate(XL,nCenter,3,3,'XL')
call mma_allocate(XR,nCenter,3,3,'XR')
mem_local = mem_local+size(ZR)*RtoB
mem_local = mem_local+size(ZL)*RtoB
mem_local = mem_local+size(XR)*RtoB
mem_local = mem_local+size(XL)*RtoB

call mma_allocate(zstat_tot,(nT+nTempMagn),'Zstat_tot')
call mma_allocate(chiT,(nT+nTempMagn),'chiT')
call mma_allocate(chi_theta_1,(nT+nTempMagn),'chi_theta_1')
mem_local = mem_local+size(zstat_tot)*RtoB
mem_local = mem_local+size(chiT)*RtoB
mem_local = mem_local+size(chi_theta_1)*RtoB

if (zJ == Zero) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'SUSC:  memory allocated (local):'
  write(u6,*) 'mem_local=',mem_local
  write(u6,*) 'SUSC:  memory allocated (total):'
  write(u6,*) 'mem_total=',mem+mem_local
#endif

  do iT=1,nT+nTempMagn
    ! initialize temporary variables:
    XL(:,:,:) = Zero
    XR(:,:,:) = Zero
    !-------------------------------------------------------------------
    ! local susceptibility= total susceptibility coming from individual magnetic centers
    do i=1,nneq
#     ifdef _DEBUGPRINT_
      write(u6,'(A,2I5)') 'nss(i)=',nss(i)
      write(u6,'(A,2I5)') 'nexch(i)=',nexch(i)
      write(u6,'(A,9F10.6)') 'eso(i,:)=',eso(i,1:nss(i))
      write(u6,'(A,9F10.6)') 'W(:)    =',W(:)
      write(u6,'(A,9F10.6)') 'T(iT)=',T(iT)
#     endif
      call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
      call mma_allocate(eso_tmp,nss(i),label='eso_tmp')
      dipso_tmp(:,:,:) = dipso(i,:,1:nss(i),1:nss(i))
      eso_tmp(:) = eso(i,1:nss(i))
      call chi(dipso_tmp,dipso_tmp,eso_tmp,nss(i),T(it),zstat_l(i),chit_tens_l(:,:,i))
      call mma_deallocate(dipso_tmp)
      call mma_deallocate(eso_tmp)

      call mma_allocate(dipso_tmp,3,nexch(i),nexch(i),label='dipso_tmp')
      call mma_allocate(eso_tmp,nexch(i),label='eso_tmp')
      dipso_tmp(:,:,:) = dipso(i,:,1:nexch(i),1:nexch(i))
      eso_tmp(:) = eso(i,1:nexch(i))
      call chi(dipso_tmp,dipso_tmp,eso_tmp,nexch(i),T(iT),zstat_lr(i),chit_tens_lr(:,:,i))
      call mma_deallocate(dipso_tmp)
      call mma_deallocate(eso_tmp)

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
        do jc=1,3
          do n1=1,3
            do n2=1,3
              XL(isite,:,jc) = XL(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*chit_tens_l(n1,n2,i)
              XR(isite,:,jc) = XR(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*chit_tens_lr(n1,n2,i)
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
      XR(:,:,:) = coeff_chi*XR(:,:,:)
      XRM(:,iT,:,:) = XR(:,:,:)
      XL(:,:,:) = coeff_chi*XL(:,:,:)
      XLM(:,iT,:,:) = XL(:,:,:)
      ZRM(:,iT) = ZR(:)
      ZLM(:,iT) = ZL(:)
    end if
    call chi_sum(nCenter,chit_tens_ex,zstat_ex,XL,ZL,XR,ZR,iopt,chit_tens_tot(:,:,it),zstat_tot(it))

    Fa = dnrm2_(nCenter,ZL,1)
    Fb = dnrm2_(nCenter,ZR,1)
    call Add_Info('XT:  ZL',[Fa],1,6)
    call Add_Info('XT:  ZR',[Fb],1,6)
    call Add_Info('XT: ZEx',[zstat_ex],1,6)

    chit(it) = coeff_chi*(chit_tens_tot(1,1,iT)+chit_tens_tot(2,2,iT)+chit_tens_tot(3,3,iT))/Three
    chit_theta(iT) = chit(iT)

    if (abs(chit(iT)) < 1.0e-20_wp) then
      chit(iT) = 1.0e-20_wp
      chit_theta(iT) = 1.0e-20_wp
    end if
    chi_theta_1(iT) = T(iT)/chit(iT)
    ! add some verification data:
    chit_theta_tens(:,:,iT) = Zero
    Fa = dnrm2_(9,chit_theta_tens(:,:,iT),1)
    call Add_Info('XT: chit_theta_tens',[Fa],1,6)
  end do ! iT
  Fb = dnrm2_(nT+nTempMagn,T,1)
  call Add_Info('XT: T',[Fb],1,6)

else ! i.e. when (zJ /= 0)

  ! allocate memory for temporary arrays:
  call mma_allocate(smu_chit_tens_l,3,3,nneq,'smu_X_tens_l')
  call mma_allocate(smu_chit_tens_lr,3,3,nneq,'smu_X_tens_lr')
  call mma_allocate(ss_chit_tens_l,3,3,nneq,'ss_chit_tens_l')
  call mma_allocate(ss_chit_tens_lr,3,3,nneq,'ss_chit_tens_lr')
  call mma_allocate(SMUR,nCenter,3,3,'SMUR')
  call mma_allocate(SMUL,nCenter,3,3,'SMUL')
  call mma_allocate(SSR,nCenter,3,3,'SSR')
  call mma_allocate(SSL,nCenter,3,3,'SSL')
  mem_local = mem_local+size(smu_chit_tens_l)*RtoB
  mem_local = mem_local+size(smu_chit_tens_lr)*RtoB
  mem_local = mem_local+size(ss_chit_tens_l)*RtoB
  mem_local = mem_local+size(ss_chit_tens_lr)*RtoB
  mem_local = mem_local+size(SMUR)*RtoB
  mem_local = mem_local+size(SMUL)*RtoB
  mem_local = mem_local+size(SSR)*RtoB
  mem_local = mem_local+size(SSL)*RtoB

# ifdef _DEBUGPRINT_
  write(u6,*) 'SUSC:  memory allocated (local):'
  write(u6,*) 'mem_local=',mem_local
  write(u6,*) 'SUSC:  memory allocated (total):'
  write(u6,*) 'mem_total=',mem+mem_local
# endif

  do iT=1,nT+nTempMagn
    ! initialization:
    XR(:,:,:) = Zero
    XL(:,:,:) = Zero
    SMUR(:,:,:) = Zero
    SMUL(:,:,:) = Zero
    SSR(:,:,:) = Zero
    SSL(:,:,:) = Zero
    call unitmat(unity,3)

    ! compute local tensors  L, and LR:
    do i=1,nneq
      call mma_allocate(eso_tmp,nss(i),label='eso_tmp')
      call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
      call mma_allocate(s_so_tmp,3,nss(i),nss(i),label='s_so_tmp')
      eso_tmp(:) = eso(i,1:nss(i))
      dipso_tmp(:,:,:) = dipso(i,:,1:nss(i),1:nss(i))
      s_so_tmp(:,:,:) = s_so(i,:,1:nss(i),1:nss(i))
      call chi(dipso_tmp,dipso_tmp,eso_tmp,nss(i),T(iT),zstat_l(i),chit_tens_l(:,:,i))
      call chi(s_so_tmp,dipso_tmp,eso_tmp,nss(i),T(iT),zstat_l(i),smu_chit_tens_l(:,:,i))
      call chi(s_so_tmp,s_so_tmp,eso_tmp,nss(i),T(iT),zstat_l(i),ss_chit_tens_l(:,:,i))
      call mma_deallocate(eso_tmp)
      call mma_deallocate(dipso_tmp)
      call mma_deallocate(s_so_tmp)

      call mma_allocate(eso_tmp,nexch(i),label='eso_tmp')
      call mma_allocate(dipso_tmp,3,nexch(i),nexch(i),label='dipso_tmp')
      call mma_allocate(s_so_tmp,3,nexch(i),nexch(i),label='s_so_tmp')
      eso_tmp(:) = eso(i,1:nexch(i))
      dipso_tmp(:,:,:) = dipso(i,:,1:nexch(i),1:nexch(i))
      s_so_tmp(:,:,:) = s_so(i,:,1:nexch(i),1:nexch(i))
      call chi(dipso_tmp,dipso_tmp,eso_tmp,nexch(i),T(iT),zstat_lr(i),chit_tens_lr(:,:,i))
      call chi(s_so_tmp,dipso_tmp,eso_tmp,nexch(i),T(iT),zstat_lr(i),smu_chit_tens_lr(:,:,i))
      call chi(s_so_tmp,s_so_tmp,eso_tmp,nexch(i),T(iT),zstat_lr(i),ss_chit_tens_lr(:,:,i))
      call mma_deallocate(eso_tmp)
      call mma_deallocate(dipso_tmp)
      call mma_deallocate(s_so_tmp)
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
        do jc=1,3
          do n1=1,3
            do n2=1,3
              XR(isite,:,jc) = XR(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*chit_tens_lr(n1,n2,i)
              XL(isite,:,jc) = XL(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*chit_tens_l(n1,n2,i)

              SMUL(isite,:,jc) = SMUL(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*smu_chit_tens_l(n1,n2,i)
              SMUR(isite,:,jc) = SMUR(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*smu_chit_tens_lr(n1,n2,i)

              SSL(isite,:,jc) = SSL(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*ss_chit_tens_l(n1,n2,i)
              SSR(isite,:,jc) = SSR(isite,:,jc)+r_lg(i,j,:,n1)*r_lg(i,j,jc,n2)*ss_chit_tens_lr(n1,n2,i)
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
      XR(:,:,:) = coeff_chi*XR(:,:,:)
      XRM(:,iT,:,:) = XR(:,:,:)
      XL(:,:,:) = coeff_chi*XL(:,:,:)
      XLM(:,iT,:,:) = XL(:,:,:)
      ZRM(:,iT) = ZR(:)
      ZLM(:,iT) = ZL(:)
    end if

    ! compute total tensors:
    call chi_sum(nCenter,chit_tens_ex,zstat_ex,XL,ZL,XR,ZR,iopt,chit_tens_tot(:,:,iT),zstat_tot(iT))

    call chi_sum(nCenter,smu_chit_tens_ex,zstat_ex,SMUL,ZL,SMUR,ZR,iopt,smu_chit_tens_tot,zstat_tot(iT))

    call chi_sum(nCenter,ss_chit_tens_ex,zstat_ex,SSL,ZL,SSR,ZR,iopt,ss_chit_tens_tot,zstat_tot(iT))

    ! form the A_dir matrix:
    ! A_dir(:,:) = 1(:,:)* kB * T(iT) - zJ*XSS(:,:)
    A_dir(:,:) = Boltz_k*T(iT)*Unity(:,:)-zJ*ss_chit_tens_tot(:,:)
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
        chit_theta_tens(ic,jc,iT) = chit_tens_tot(ic,jc,iT)+zj*xxm
      end do ! jc
    end do ! ic

    chit(iT) = coeff_chi*(chit_tens_tot(1,1,iT)+chit_tens_tot(2,2,iT)+chit_tens_tot(3,3,iT))/Three

    chit_theta(iT) = coeff_chi*(chit_theta_tens(1,1,iT)+chit_theta_tens(2,2,iT)+chit_theta_tens(3,3,iT))/Three
    if (abs(chit(iT)) < 1.0e-20_wp) then
      chit(iT) = 1.0e-20_wp
      chit_theta(iT) = 1.0e-20_wp
    end if
    if (abs(chit_theta(iT)) < 1.0e-20_wp) chit_theta(iT) = 1.0e-20_wp

    chi_theta_1(iT) = t(iT)/chit_theta(iT)

    ! add some verification data:
    Fa = dnrm2_(9,chit_tens_tot(:,:,it),1)
    Fb = dnrm2_(9,chit_theta_tens(:,:,it),1)
    Fc = dnrm2_(9,smu_chit_tens_tot,1)
    Fd = dnrm2_(9,ss_chit_tens_tot,1)

    call Add_Info('XT: chit_tens_tot',[Fa],1,6)
    call Add_Info('XT: chit_theta_tens',[Fb],1,6)
    call Add_Info('XT: smu_chit_tens_tot',[Fc],1,6)
    call Add_Info('XT:  ss_chit_tens_tot',[Fd],1,6)
  end do ! it
  Fb = dnrm2_(nT+nTempMagn,T,1)
  call Add_Info('XT: T',[Fb],1,6)

  call mma_deallocate(smu_chit_tens_l)
  call mma_deallocate(smu_chit_tens_lr)
  call mma_deallocate(ss_chit_tens_l)
  call mma_deallocate(ss_chit_tens_lr)
  call mma_deallocate(SMUR)
  call mma_deallocate(SMUL)
  call mma_deallocate(SSR)
  call mma_deallocate(SSL)

end if  !zJ
!-----------------------------------------------------------------------
! printing the results
do iT=1,nT+nTempMagn
  do ic=1,3
    do jc=1,3
      chit_tens_tot(ic,jc,iT) = coeff_chi*chit_tens_tot(ic,jc,iT)
      if (zJ /= Zero) chit_theta_tens(ic,jc,iT) = coeff_chi*chit_theta_tens(ic,jc,iT)
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
  write(u6,'(a,5x, f20.14)') 'ST.DEV.CHIT:',dev(nT,chit_theta(1+nTempMagn:),XTexp(1+nTempMagn:))

!-------------------------  PLOTs -------------------------------------!
write(label,'(A)') 'no_field'
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T(1+nTempMagn:),chit_theta(1+nTempMagn:),XTexp(1+nTempMagn:))
  else
    call plot_XT_no_Exp(label,nT,T(1+nTempMagn:),chit_theta(1+nTempMagn:))
  end if
end if
!---------------------- END PLOTs -------------------------------------!

! print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
call mma_allocate(wt,3,'wt')
call mma_allocate(zt,3,3,'zt')
if (zJ == Zero) then
  write(u6,'(/)')
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(31X,A,22x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ = 0,  in cm3*K/mol','|'
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
  do iT=1,nT
    jT = iT+nTempMagn
    info = 0
    call DIAG_R2(chit_tens_tot(:,:,jT),3,info,wt,zt)
    write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- '// &
                    'z ----|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(chit_tens_tot(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                    (zt(j,1),j=1,3),'|'
    write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(chiT_tens_tot(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                             (zt(j,2),j=1,3),'|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(chiT_tens_tot(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                    (zt(j,3),j=1,3),'|'
  end do
  write(u6,'(2A)') repeat('-',110),'|'

else ! zJ /= Zero

  write(u6,'(/)')
  write(u6,'(2A)') repeat('-',110),'|'

  write(u6,'(31X,A,F9.6,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR zJ =',zJ,',  in cm3*K/mol','|'
  write(u6,'(2A)') repeat('-',110),'|'
  write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
  do iT=1,nT
    jT = iT+nTempMagn
    info = 0
    call DIAG_R2(chit_theta_tens(:,:,jT),3,info,wt,zt)
    write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- '// &
                    'z ----|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(chit_theta_tens(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                    (zt(j,1),j=1,3),'|'
    write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(chiT_theta_tens(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                             (zt(j,2),j=1,3),'|'
    write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(chiT_theta_tens(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                    (zt(j,3),j=1,3),'|'
  end do
  write(u6,'(2A)') repeat('-',110),'|'
end if
call mma_deallocate(wt)
call mma_deallocate(zt)

call mma_deallocate(chit_tens_tot)
call mma_deallocate(chit_theta_tens)
call mma_deallocate(chit_tens_l)
call mma_deallocate(chit_tens_lr)
call mma_deallocate(zstat_l)
call mma_deallocate(zstat_lr)
call mma_deallocate(zstat_tot)
call mma_deallocate(ZR)
call mma_deallocate(ZL)
call mma_deallocate(XL)
call mma_deallocate(XR)
call mma_deallocate(chiT)
call mma_deallocate(chi_theta_1)

!-----------------------------------------------------------------------
!write(u6,*)
!write(u6,'(5x,a)') 'on user request, the magnetic susceptibility and '
!write(u6,'(5x,a)') 'the magnetic susceptibility tensor were not calculated.'
!write(u6,*)

return

end subroutine susceptibility_pa
