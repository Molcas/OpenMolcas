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
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas(0:nIrrep-1), nDim, n2Tot, MaxScr, iUHF
real(kind=wp), intent(in) :: Occ(nDim), Vec(n2Tot)
real(kind=wp), intent(out) :: Scr(MaxScr)
logical(kind=iwp), intent(in) :: var, Short, ifallorb
integer(kind=iwp) :: iadC1, iadC2, iadDen, iadDen_ab, iadEl, iadEl_Work, iadElSum, iadLab, iadNuc, iadNucSum, iadopr, iadtmp, &
                     iadtmt, idum(1), iInd1, iInd2, iOcc_ab, iopt, ip_Scr, ipD1ao, iPL, ipScr, ipVec, irc, iscr, iSmLbl, iTol, &
                     jRC, lpole, maxCen, maxGG, mBas(0:7), mDim, mInt, nbast, nblock, nCen, nComp, nfblock, nscr, tNUC
logical(kind=iwp) :: NxtOpr
character(len=8) :: label
real(kind=wp), parameter :: Thrs = 1.0e-6_wp
integer(kind=iwp), external :: ip_of_Work, iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
#include "hfc_logical.fh"
#include "WrkSpc.fh"

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

  if (iscr > maxscr) Go To 999
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Short) then
    iadDen = 1

    ! calculate the density matrix with all off-diagonal elements
    ! multipled by 2

    !if (iUHF ==  0) then
    call dcopy_(nblock,[Zero],0,Scr(iadDen),1)
    call GetMem('D1ao','Allo','Real',ipD1ao,nBlock)
    if (var) then
      call Get_D1ao_Var(Work(ipD1ao),nBlock)
    else
      call Get_D1ao(Work(ipD1ao),nBlock)
    end if
    do i=1,nBlock
      SCR(i) = Work(ipD1ao+i-1)
    end do
    call GetMem('Dens','Free','Real',ipD1ao,nBlock)
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
      call GetMem('iadEl1','Allo','Real',iadEl_Work,nComp*mDim)
      ip_Scr = ip_of_Work(Scr(1))
      iadEl = iadEl_Work-(ip_Scr-1)
    end if
    iscr = nscr+10+4*nComp
    if (iscr > maxscr) then
      if (.not. Short) call Free_Work(iadEl_Work)
      Go To 999
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
      if (irc /= 0) go to 101
      NxtOpr = .true.
      irc = -1
      iopt = 0
      iadOpr = iadLab+2*nComp
      call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
      if (irc /= 0) go to 101
      if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
      scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
      if (iComp == 1) then
        do k=0,2
          scr(iadC1+k) = scr(iadOpr+mInt+k)
          scr(iadC2+k) = scr(iadOpr+mInt+k)
        end do
      end if
      if (mInt == 0) Go To 101
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
      if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                      Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
  101 continue
    end do
    if (.not. NxtOpr) then
      if (.not. Short) call Free_Work(iadEl_Work)
      Go To 195
    end if
    iadTmt = iadOpr+nblock
    if (i <= 4) then
      iadTmp = iadTmt+nComp**2
      iscr = iadTmp+nComp
      if (iscr > maxscr) then
        if (.not. Short) call Free_Work(iadEl_Work)
        Go To 999
      end if
    else
      iadTmp = iadTmt
    end if

    call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
    call prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),i,cScr,scr(iadTmt),scr(iadTmp), &
              ifallorb)
    nullify(cScr)
    if (.not. Short) call Free_Work(iadEl_Work)
  end do

  ! Scan 'ONEINT' for magnetic hyperfine integrals

  195 continue
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
  else
    Go to 199
  end if

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

  ! Scan 'ONEINT' for electric field integrals

  199 continue
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  !write(u6,*) ' Starting scan of ONEINT for various elec. field integrals'

  do iEF=0,2
    nComp = (iEF+1)*(iEF+2)/2

    iadEl = iadNuc+nComp
    iadLab = iadEl+nComp
    if (.not. Short) then
      call GetMem('iadEl2','Allo','Real',iadEl_Work,nComp*mDim)
      ip_Scr = ip_of_Work(Scr(1))
      iadEl = iadEl_Work-(ip_Scr-1)
    end if
    ! create vectors to store the sums of electronic and nuclear components over all centers
    call GetMem('ElSum','Allo','Real',iadElSum,nComp)
    call GetMem('NucSum','Allo','Real',iadNucSum,nComp)
    call DZero(Work(iadElSum),nComp)
    call DZero(Work(iadNucSum),nComp)

    ! loop over different operator origins (max.99999)

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
        if (irc /= 0) go to 201
        NxtOpr = .true.
        irc = -1
        iopt = 0
        iadOpr = iadLab+2*nComp
        call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
        if (irc /= 0) Go To 201
        if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
        scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
        if (iComp == 1) then
          do k=0,2
            scr(iadC1+k) = scr(iadOpr+mInt+k)
            scr(iadC2+k) = scr(iadOpr+mInt+k)
          end do
        end if
        if (mInt == 0) Go To 201
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
        if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                        Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
  201   continue
      end do
      if (.not. NxtOpr) then
        if (.not. Short) call Free_Work(iadEl_Work)
        Go to 299
      end if

      call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
      call Prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),iEF,cScr,scr(iadTmt), &
                scr(iadTmp),ifallorb)
      nullify(cScr)
      ! add the components to the sums, and update the total number of centers
      do iComp=0,nComp-1
        iInd1 = iadElSum+iComp
        do iOcc=0,mDim-1
          Work(iInd1) = Work(iInd1)+Scr(iadEl+iComp*mDim+iOcc)
        end do
      end do
      call DaXpY_(nComp,One,Scr(iadNuc),1,Work(iadNucSum),1)
      nCen = i
    end do
    if (.not. Short) call Free_Work(iadEl_Work)

  299 continue
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
        call Add_Info(label,Work(iadElSum),nComp,iTol)
        write(label,'(a,i1,a)') 'EF',iEF,'  nuc'
        call Add_Info(label,Work(iadNucSum),nComp,iTol)
      end if
    end if
    call GetMem('ElSum','Free','Real',iadElSum,nComp)
    call GetMem('NucSum','Free','Real',iadNucSum,nComp)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  !write(u6,*) ' Starting scan of ONEINT for various contact term integrals'

  nComp = 1

  iadEl = iadNuc+nComp
  iadLab = iadEl+nComp
  if (.not. Short) then
    call GetMem('iadEl2','Allo','Real',iadEl_Work,nComp*mDim)
    ip_Scr = ip_of_Work(Scr(1))
    iadEl = iadEl_Work-(ip_Scr-1)
  end if
  ! create vectors to store the sums of electronic and nuclear components over all centers
  call GetMem('ElSum','Allo','Real',iadElSum,nComp)
  call GetMem('NucSum','Allo','Real',iadNucSum,nComp)
  call DZero(Work(iadElSum),nComp)
  call DZero(Work(iadNucSum),nComp)

  ! loop over different operator origins (max.99999)

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
    if (irc /= 0) go to 301
    NxtOpr = .true.
    irc = -1
    iopt = 0
    iadOpr = iadLab+2*nComp
    call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
    if (irc /= 0) Go To 301
    if (mInt /= 0) call CmpInt(scr(iadOpr),mInt,nBas,nIrrep,iSmLbl)
    scr(iadNuc+icomp-1) = scr(iadOpr+mInt+3)
    do k=0,2
      scr(iadC1+k) = scr(iadOpr+mInt+k)
      scr(iadC2+k) = scr(iadOpr+mInt+k)
    end do
    if (mInt == 0) Go To 301
    call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
    if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                    Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
  301 continue
    if (.not. NxtOpr) then
      if (.not. Short) call Free_Work(iadEl_Work)
      Go To 399
    end if

    call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
    call Prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),iEF,cScr,scr(iadTmt),scr(iadTmp), &
              ifallorb)
    nullify(cScr)
    ! add the components to the sums, and update the total number of centers
    do iComp=0,nComp-1
      iInd1 = iadElSum+iComp
      iInd2 = iadEl+iComp*mDim
      do iOcc=0,mDim-1
        Work(iInd1) = Work(iInd1)+Scr(iInd2+iOcc)
      end do
    end do
    call DaXpY_(nComp,One,Scr(iadNuc),1,Work(iadNucSum),1)
    nCen = i
  end do
  if (.not. Short) call Free_Work(iadEl_Work)

  399 continue
  if (nCen > 0) then
    ! set the tolerance according to the total number of centers
    ! (assuming error scales with sqrt(ncen))
    iTol = 5
    iTol = iTol-nint(Half*log10(real(nCen,kind=wp)))
    if (.not. MAG_X2C) then
      write(label,'(a,a)') 'CNT','   el'
      call Add_Info(label,Work(iadElSum),nComp,iTol)
      write(label,'(a,a)') 'CNT','  nuc'
      call Add_Info(label,Work(iadNucSum),nComp,iTol)
    end if
  end if
  call GetMem('ElSum','Free','Real',iadElSum,nComp)
  call GetMem('NucSum','Free','Real',iadNucSum,nComp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !write(u6,*) ' Starting scan of ONEINT diamagnetic shielding'

  nComp = 9
  lpole = 2

  iadEl = iadNuc+nComp
  iadLab = iadEl+nComp
  if (.not. Short) then
    call GetMem('iadEl3','Allo','Real',iadEl_Work,nComp*mDim)
    ip_Scr = ip_of_Work(Scr(1))
    iadEl = iadEl_Work-(ip_Scr-1)
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
        if (irc /= 0) go to 402
        NxtOpr = .true.
        irc = -1
        iopt = 0
        iadOpr = iadLab+2*nComp
        call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
        if (irc /= 0) go to 402
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
        if (mInt == 0) Go To 402
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim))
        if ((.not. Short) .and. (iUHF == 1)) call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen_ab),nDim,Occ(iOcc_ab),nblock, &
                                                        Scr(iadOpr),Scr(iadEl+(iComp-1)*mDim+nDim))
  402   continue
      end do
      if (.not. NxtOpr) Go To 4000

      call c_f_pointer(c_loc(scr(iadLab)),cScr,[1])
      call prop(short,label,scr(iadC1),scr(iadC2),nirrep,mBas,mDim,occ,Thrs,scr(iadEl),scr(iadNuc),lpole,cScr,scr(iadTmt), &
                scr(iadTmp),ifallorb)
      nullify(cScr)
      jRC = 1
    end do
  4000 continue
    if (jRC == 0) then
      if (.not. Short) call Free_Work(iadEl_Work)
      Go To 499
    end if
  end do
  if (.not. Short) call Free_Work(iadEl_Work)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  499 continue
  if (iPL >= 2) then
    call CollapseOutput(0,'   Molecular properties:')
    write(u6,*)
  end if
  return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  999 continue
  write(u6,'(//1x,a/1x,a/1x,a,i8,a,i8)') ' Warning:',' Not enough scratch area to perform calculations',' Needed at least:',iscr, &
                                        '   available:',maxscr
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (iPL >= 2) then
    call CollapseOutput(0,'   Molecular properties:')
    write(u6,*)
  end if

end subroutine Prpt_Internal

end subroutine Prpt_
