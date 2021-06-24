!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!               2018, Sijia S. Dong                                    *
!***********************************************************************

subroutine Prpt_(nIrrep,nBas,nDim,Occ,n2Tot,Vec,MaxScr,Scr,var,Short,iUHF,ifallorb)
!***********************************************************************
!                                                                      *
! Purpose: calculation of expectation values of different              *
!          operators as available on the 'ONEINT' file                 *
!                                                                      *
! Caution: before calling this subroutine one needs to                 *
!          open the ONEINT file                                        *
!                                                                      *
! Calling parameters:                                                  *
!   nIrRep            number of irreducible representations            *
!   nBas(0:nIrRep)    number of basis functions in each repre-         *
!                     sentation                                        *
!   ndim              total number of basis functions                  *
!   occ(1:ndim)       occupation number for all eigenvectors           *
!   maxscr            the maximum available size of the                *
!                     scratch area                                     *
!   scr(1:maxscr)     a scratch area whose size, nscr, can be          *
!                     calculated in the following way:                 *
!                                                                      *
!                     nscr =sum(i,i=0,nirrep-1)(nbas(i)*(nbas(i)+1)/2) *
!                           +3                                         *
!                           +3                                         *
!                           +nComp                                     *
!                           +nComp                                     *
!                           +2*nComp                                   *
!                           +ndim*(ndim+1)/2+4                         *
!                           +2*ntComp*(ntComp+1)                       *
!                                                                      *
!                     and should not exceed maxscr.                    *
!                     nComp is the number of cartresian components     *
!                     for the given operator                           *
!                     ntComp is the number of components of            *
!                     the l-th moment opartors which are transformed   *
!                     into l-pole moment. Currently,                   *
!                     ntComp=15 (hexadecapole moments)                 *
!     ifallorb        logical option for whether the property of       *
!                     all orbitals are printed (and not weighted by    *
!                     occupation number)in property calculation        *
!                     (S.S.Dong, 2018)                                 *
!                                                                      *
! 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden.          *
! Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
! - Enable properties to be printed for all orbitals                   *
! (including virtuals) and not weighted by occupation numbers          *
!***********************************************************************

use iso_c_binding, only: c_f_pointer, c_loc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas(0:nIrrep-1), nDim, n2Tot, MaxScr, iUHF
real(kind=wp), intent(in) :: Occ(nDim), Vec(n2Tot)
real(kind=wp), intent(out) :: Scr(MaxScr)
logical(kind=iwp), intent(in) :: var, Short, ifallorb
integer(kind=iwp) :: iadC1, iadC2, iadDen, iadDen_ab, iadEl, iadLab, iadNuc, iadopr, iadtmp, iadtmt, idum(1), iOcc_ab, iopt, &
                     ip_Scr, iPL, ipScr, ipVec, irc, iscr, iSmLbl, iTol, jRC, lpole, maxCen, maxGG, mBas(0:7), mDim, mInt, nbast, &
                     nblock, nCen, nComp, nfblock, nscr, tNUC
logical(kind=iwp) :: NxtOpr
character(len=8) :: label
real(kind=wp), allocatable :: D1ao(:), El_Work(:,:), ElSum(:), NucSum(:)
real(kind=wp), parameter :: Thrs = 1.0e-6_wp
integer(kind=iwp), external :: ip_of_Work, iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
#include "hfc_logical.fh"

call Prpt_Internal(Scr)

return

! This is to allow type punning without an explicit interface
contains

subroutine Prpt_Internal(Scr)

  real(kind=wp), target :: Scr(*)
  character, pointer :: cScr(:)
  integer(kind=iwp) :: i, iComp, iEF, iIrrep, iOcc, j, k

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  iPL = iPrintLevel(-1)
  if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

  if (iPL >= 2) then
    write(u6,*)
    call CollapseOutput(1,'   Molecular properties:')
    write(u6,'(3X,A)') '   ---------------------'
    write(u6,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call ICopy(nIrrep,nBas,1,mBas,1)
  if (Short) then
    mDim = 1
  else
    mDim = nDim
    if (iUHF == 1) then
      mDim = 2*mDim
      do iIrrep=0,nIrrep-1
        mBas(iIrrep) = 2*mBas(iIrrep)
      end do
    end if
  end if
  iadopr = -1   ! dummy initialize
  iadtmt = -1   ! dummy initialize
  iadtmp = -1   ! dummy initialize
  iadDen_ab = 1 ! dummy initialize
  iOcc_ab = 1   ! dummy initialize
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nblock = 0
  do i=0,nirrep-1
    nblock = nblock+nbas(i)*(nbas(i)+1)/2
  end do
  nfblock = ndim*(ndim+1)/2

  ! estimate the minimum scratch area needed for calculatig
  ! the average values of a 3 component operator

  nscr = nblock+nfblock
  iscr = nscr+76

  if (iscr > maxscr) then
    call Error()
    return
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Short) then
    iadDen = 1

    ! calculate the density matrix with all off-diagonal elements
    ! multipled by 2

    !if (iUHF ==  0) then
    call dcopy_(nblock,[Zero],0,Scr(iadDen),1)
    call mma_allocate(D1ao,nBlock,label='D1ao')
    if (var) then
      call Get_D1ao_Var(D1ao,nBlock)
    else
      call Get_D1ao(D1ao,nBlock)
    end if
    Scr(1:nBlock) = D1ao(:)
    call mma_deallocate(D1ao)
    !end if
    iadC1 = iadDen+nblock

  else

    ! Make iadDen to point at the MO vectors
    ipVec = ip_of_Work(Vec(1))
    ipScr = ip_of_Work(Scr(1))
    iadDen = ipVec-ipScr+1
    iadC1 = 1
    if (iUHF == 1) then
      iadDen_ab = iadDen+nDim**2
      iOcc_ab = 1+nDim
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scan the ONEINT file for multipole moment operators

  iadC2 = iadC1+3
  iadNuc = iadC2+3

  !write(u6,*) ' Starting scan of ONEINT for multipole moments'
  do i=0,99
    NxtOpr = .false.
    nComp = (i+1)*(i+2)/2

    iadEl = iadNuc+nComp
    iadLab = iadEl+nComp
    if (.not. Short) then
      call mma_allocate(El_Work,nDim,nComp,label='iadEl1')
      ip_Scr = ip_of_Work(Scr(1))
      iadEl = ip_of_Work(El_Work(1,1))-(ip_Scr-1)
    end if
    iscr = nscr+10+4*nComp
    if (iscr > maxscr) then
      if (.not. Short) call mma_deallocate(El_Work)
      call Error()
      return
    end if

    call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
    call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl),1)
    call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
    write(label,'(a,i2)') 'MLTPL ',i
    do iComp=1,nComp
      irc = -1
      iopt = 1
      call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
      if (irc == 0) mInt = idum(1)
      if (irc /= 0) cycle
      NxtOpr = .true.
      irc = -1
      iopt = 0
      iadOpr = iadLab+2*nComp
      call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
      if (irc /= 0) cycle
      if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
      scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
      if (iComp == 1) then
        do k=0,2
          scr(iadC1+k) = scr(iadOpr+mInt+k)
          scr(iadC2+k) = scr(iadOpr+mInt+k)
        end do
      end if
      if (mInt == 0) cycle
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
      if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                      Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
    end do
    if (.not. NxtOpr) then
      if (.not. Short) call mma_deallocate(El_Work)
      exit
    end if
    iadTmt = iadOpr+nblock
    if (i <= 4) then
      iadTmp = iadTmt+nComp**2
      iscr = iadTmp+nComp
      if (iscr > maxscr) then
        if (.not. Short) call mma_deallocate(El_Work)
        call Error()
        return
      end if
    else
      iadTmp = iadTmt
    end if

    call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
    call prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),i,cScr,scr(iadTmt),scr(iadTmp), &
              ifallorb)
    nullify(cScr)
    if (.not. Short) call mma_deallocate(El_Work)
  end do

  ! Scan 'ONEINT' for magnetic hyperfine integrals

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! This check if the magnetic hypferfine integrals are available
  ! and only calculates hyperfine tensor matrix for UHF wavefunction
  ! in C1 symmetry. Since for UHF the spin density matrix is easily
  ! obtained.

  MAG_X2C = .false.
  irc = -1
  iopt = 1
  label = 'MAGXP  1'
  iComp = 1
  call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
  if (irc == 0) then
    MAG_X2C = .true.
    !if (Method == 'UHF-SCF ') then
    if (iUHF == 1) then
      if (UHF_HFC) then
        if (nIrrep == 1) then
          tNUC = 0
          nbast = nbas(0)
          call Get_iScalar('LP_nCenter',tNUC)
          call cmp_hfc(nbast,tNUC)
        else
          write(u6,'(/,/,6X,A)') 'Skipping Hyperfine tensor matrix for UHF with symmetry'
        end if
      else
        write(u6,'(/,/,6X,A)') 'Skipping Hyperfine tensor matrix'
      end if
    end if
  end if

  ! Scan 'ONEINT' for electric field integrals

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  !write(u6,*) ' Starting scan of ONEINT for various elec. field integrals'

  do iEF=0,2
    nComp = (iEF+1)*(iEF+2)/2

    iadEl = iadNuc+nComp
    iadLab = iadEl+nComp
    if (.not. Short) then
      call mma_allocate(El_Work,mDim,nComp,label='iadEl2')
      ip_Scr = ip_of_Work(Scr(1))
      iadEl = ip_of_Work(El_Work(1,1))-(ip_Scr-1)
    end if
    ! create vectors to store the sums of electronic and nuclear components over all centers
    call mma_allocate(ElSum,nComp,label='ElSum')
    call mma_allocate(NucSum,nComp,label='NucSum')
    ElSum(:) = Zero
    NucSum(:) = Zero

    ! loop over different operator origins (max. 99999)

    maxCen = 99999
    nCen = 0
    do i=1,maxCen
      call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
      call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl),1)
      call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
      write(label,'(a,i1,i5)') 'EF',iEF,i
      NxtOpr = .false.
      do iComp=1,nComp
        irc = -1
        iopt = 1
        call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
        if (irc == 0) mInt = idum(1)
        if (irc /= 0) cycle
        NxtOpr = .true.
        irc = -1
        iopt = 0
        iadOpr = iadLab+2*nComp
        call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
        if (irc /= 0) cycle
        if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
        scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
        if (iComp == 1) then
          do k=0,2
            scr(iadC1+k) = scr(iadOpr+mInt+k)
            scr(iadC2+k) = scr(iadOpr+mInt+k)
          end do
        end if
        if (mInt == 0) cycle
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
        if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                        Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
      end do
      if (.not. NxtOpr) exit

      call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
      call Prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),iEF,cScr,scr(iadTmt), &
                scr(iadTmp),ifallorb)
      nullify(cScr)
      ! add the components to the sums, and update the total number of centers
      do iComp=1,nComp
        do iOcc=1,mDim
          ElSum(iComp) = ElSum(iComp)+Scr(iadEl+(iComp-1)*mDim+iOcc-1)
        end do
      end do
      call DaXpY_(nComp,One,Scr(iadNuc),1,NucSum,1)
      nCen = i
    end do
    if (.not. Short) call mma_deallocate(El_Work)

    if (nCen > 0) then
      ! set the tolerance according to the total number of centers
      ! (assuming error scales with sqrt(ncen))
      iTol = 5
      iTol = iTol-nint(Half*log10(real(nCen,kind=wp)))
      ! set MAG_X2C to avoid tests of electric field properties when
      ! wavefunction is X2C transformed (there is no way to tell but we
      ! can kind of tell by reading MAG x2c integrals) This is a workaround.
      if (.not. MAG_X2C) then
        write(label,'(a,i1,a)') 'EF',iEF,'   el'
        call Add_Info(label,ElSum,nComp,iTol)
        write(label,'(a,i1,a)') 'EF',iEF,'  nuc'
        call Add_Info(label,NucSum,nComp,iTol)
      end if
    end if
    call mma_deallocate(ElSum)
    call mma_deallocate(NucSum)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  !write(u6,*) ' Starting scan of ONEINT for various contact term integrals'

  nComp = 1

  iadEl = iadNuc+nComp
  iadLab = iadEl+nComp
  if (.not. Short) then
    call mma_allocate(El_Work,mDim,nComp,label='iadEl2')
    ip_Scr = ip_of_Work(Scr(1))
    iadEl = ip_of_Work(El_Work(1,1))-(ip_Scr-1)
  end if
  ! create vectors to store the sums of electronic and nuclear components over all centers
  call mma_allocate(ElSum,nComp,label='ElSum')
  call mma_allocate(NucSum,nComp,label='NucSum')
  ElSum(:) = Zero
  NucSum(:) = Zero

  ! loop over different operator origins (max. 99999)

  maxCen = 99999
  nCen = 0
  do i=1,maxCen
    call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
    call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl),1)
    call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
    write(label,'(a,i5)') 'CNT',i
    NxtOpr = .false.

    iComp = 1
    irc = -1
    iopt = 1
    call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
    if (irc == 0) mInt = idum(1)
    if (irc == 0) then
      NxtOpr = .true.
      irc = -1
      iopt = 0
      iadOpr = iadLab+2*nComp
      call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
      if (irc == 0) then
        if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
        scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
        do k=0,2
          scr(iadC1+k) = scr(iadOpr+mInt+k)
          scr(iadC2+k) = scr(iadOpr+mInt+k)
        end do
        if (mInt /= 0) then
          call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
          if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab), &
                                                          nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
        end if
      end if
    end if
    if (.not. NxtOpr) exit

    call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
    call Prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),iEF,cScr,scr(iadTmt),scr(iadTmp), &
              ifallorb)
    nullify(cScr)
    ! add the components to the sums, and update the total number of centers
    do iComp=1,nComp
      do iOcc=1,mDim
        ElSum(iComp) = ElSum(iComp)+Scr(iadEl+(iComp-1)*mDim+iOcc-1)
      end do
    end do
    call DaXpY_(nComp,One,Scr(iadNuc),1,NucSum,1)
    nCen = i
  end do
  if (.not. Short) call mma_deallocate(El_Work)

  if (nCen > 0) then
    ! set the tolerance according to the total number of centers
    ! (assuming error scales with sqrt(ncen))
    iTol = 5
    iTol = iTol-nint(Half*log10(real(nCen,kind=wp)))
    if (.not. MAG_X2C) then
      write(label,'(a,a)') 'CNT','   el'
      call Add_Info(label,ElSum,nComp,iTol)
      write(label,'(a,a)') 'CNT','  nuc'
      call Add_Info(label,NucSum,nComp,iTol)
    end if
  end if
  call mma_deallocate(ElSum)
  call mma_deallocate(NucSum)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !write(u6,*) ' Starting scan of ONEINT diamagnetic shielding'

  nComp = 9
  lpole = 2

  iadEl = iadNuc+nComp
  iadLab = iadEl+nComp
  if (.not. Short) then
    call mma_allocate(El_Work,mDim,nComp,label='iadEl3')
    ip_Scr = ip_of_Work(Scr(1))
    iadEl = ip_of_Work(El_Work(1,1))-(ip_Scr-1)
  end if

  maxGG = 99
  maxCen = 99
  ! loop over different gauge origins (max.99)
  do j=1,maxGG
    call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
    call dcopy_(nComp*mDim,[Zero],0,Scr(iadEl),1)
    call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
    jRC = 0
    ! loop over different operator origins (max.99)
    do i=1,maxCen
      write(label,'(a,i2,i2)') 'DMS ',j,i
      NxtOpr = .false.
      do iComp=1,nComp
        irc = -1
        iopt = 1
        call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
        if (irc == 0) mInt = idum(1)
        if (irc /= 0) cycle
        NxtOpr = .true.
        irc = -1
        iopt = 0
        iadOpr = iadLab+2*nComp
        call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
        if (irc /= 0) cycle
        if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
        scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
        if (iComp == 1) then
          do k=0,2
            scr(iadC1+k) = scr(iadOpr+mInt+k)
          end do
        end if
        if (iComp == 2) then
          do k=0,2
            scr(iadC2+k) = scr(iadOpr+mInt+k)
          end do
        end if
        if (mInt == 0) cycle
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
        if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                        Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
      end do
      if (.not. NxtOpr) exit

      call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
      call prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),lpole,cScr,scr(iadTmt), &
                scr(iadTmp),ifallorb)
      nullify(cScr)
      jRC = 1
    end do
    if (jRC == 0) exit
  end do
  if (.not. Short) call mma_deallocate(El_Work)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iPL >= 2) then
    call CollapseOutput(0,'   Molecular properties:')
    write(u6,*)
  end if
  return

end subroutine Prpt_Internal

subroutine Error()
  write(u6,'(//1x,a/1x,a/1x,a,i8,a,i8)') ' Warning:',' Not enough scratch area to perform calculations',' Needed at least:',iscr, &
                                         '   available:',maxscr
  if (iPL >= 2) then
    call CollapseOutput(0,'   Molecular properties:')
    write(u6,*)
  end if
end subroutine Error

end subroutine Prpt_
