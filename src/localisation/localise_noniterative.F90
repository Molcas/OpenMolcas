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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine Localise_Noniterative(irc,Model,xNrm)
! Author: T.B. Pedersen
!
! Purpose: Non-iterative localisation of orbitals.
!          Models implemented:
!            Cholesky [MODEL='CHOL']
!            PAO      [MODEL='PAO ']

use Localisation_globals, only: AnaPAO, CMO, nBas, nFro, nOrb, nOrb2Loc, nSym, Occ, Thrs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
character(len=4), intent(in) :: Model
real(kind=wp), intent(out) :: xNrm
integer(kind=iwp) :: idum(1), iSym, kOff1, kOffC, kSav, l_Dens, l_Dv, l_R, Lu_, nOrPs
real(kind=wp) :: dum(1), yNrm
character(len=80) :: Txt
character(len=6) :: Namefile
character(len=4) :: myModel
logical(kind=iwp) :: Normalize
real(kind=wp), allocatable :: Dens(:), Dv(:), DvSav(:), R(:)
logical(kind=iwp), parameter :: Test_OrthoPAO = .false.
character(len=*), parameter :: SecNam = 'Localise_Noniterative'
integer(kind=iwp), external :: isFreeUnit

irc = 0
xNrm = Zero

myModel = Model
call UpCase(myModel)
if (myModel == 'CHOL') then
  !if (.not. Silent) then
  write(u6,'(/,1X,A)') 'Cholesky localisation'
  write(u6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (decomposition)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  l_Dens = nBas(1)**2
  do iSym=2,nSym
    l_Dens = max(l_Dens,nBas(iSym)**2)
  end do
  call mma_allocate(Dens,l_Dens,label='Dens')
  kOffC = 1
  do iSym=1,nSym
    if (nOrb2Loc(iSym) > 0) then
      kOff1 = kOffC+nBas(iSym)*nFro(iSym)
      call GetDens_Localisation(Dens,CMO(kOff1),nBas(iSym),nOrb2Loc(iSym))
      call ChoLoc(irc,Dens,CMO(kOff1),Thrs,yNrm,nBas(iSym),nOrb2Loc(iSym))
      xNrm = xNrm+yNrm*yNrm
      if (irc /= 0) then
        call mma_deallocate(Dens)
        irc = 1
        xNrm = -huge(xNrm)
        return
      end if
    end if
    kOffC = kOffC+nBas(iSym)**2
  end do
  xNrm = sqrt(xNrm)
  call mma_deallocate(Dens)
else if (myModel == 'PAO ') then
  !if (.not. Silent) then
  write(u6,'(/,1X,A)') 'PAO Cholesky localisation'
  write(u6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (decomposition)'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  l_Dv = nBas(1)**2
  l_R = nBas(1)**2
  do iSym=2,nSym
    l_R = l_R+nBas(iSym)**2
    l_Dv = max(l_Dv,nBas(iSym)**2)
  end do
  call mma_allocate(R,l_R,label='R')
  call mma_allocate(Dv,l_Dv,label='Dv')
  Normalize = .true.
  call GetRawPAOs(R,CMO,nBas,nOrb,nFro,nOrb2Loc,nSym,Normalize)
  if (AnaPAO) then
    call mma_allocate(DvSav,l_R,label='DvSav')
  end if
  kSav = 1
  kOffC = 1
  do iSym=1,nSym
    if (nOrb2Loc(iSym) > 0) then
      call GetDens_Localisation(Dv,R(kOffC),nBas(iSym),nBas(iSym))
      if (AnaPAO) then
        call dCopy_(nBas(iSym)**2,Dv,1,DvSav(kSav),1)
        kSav = kSav+nBas(iSym)**2
      end if
      kOff1 = kOffC+nBas(iSym)*nFro(iSym)
      call ChoLoc(irc,Dv,CMO(kOff1),Thrs,yNrm,nBas(iSym),nOrb2Loc(iSym))
      xNrm = xNrm+yNrm*yNrm
      if (irc /= 0) then
        if (AnaPAO) call mma_deallocate(DvSav)
        call mma_deallocate(R)
        call mma_deallocate(Dv)
        irc = 1
        xNrm = -huge(xNrm)
        return
      end if
    end if
    kOffC = kOffC+nBas(iSym)**2
  end do
  xNrm = sqrt(xNrm)
  if (AnaPAO) then
    call PAO_Analysis(DvSav,R,CMO)
    call mma_deallocate(DvSav)
  end if
  write(Namefile,'(A)') 'DPAORB'
  write(Txt,'(80X)')
  write(Txt,'(A)') 'Linearly dependent PAOs'
  Lu_ = isFreeUnit(11)
  call WrVec_Localisation(Namefile,Lu_,'CO',nSym,nBas,nBas,R,Occ,dum,idum,Txt)
  !if (.not. Silent) then
  write(u6,'(1X,A)') 'The DPAORB file has been written.'
  !end if
  write(Namefile,'(A)') 'IPAORB'
  write(Txt,'(80X)')
  write(Txt,'(A)') 'Linearly independent PAOs'
  Lu_ = isFreeUnit(11)
  call WrVec_Localisation(Namefile,Lu_,'CO',nSym,nBas,nBas,CMO,Occ,dum,idum,Txt)
  !if (.not. Silent) then
  write(u6,'(1X,A)') 'The IPAORB file has been written.'
  !end if
  call mma_deallocate(R)
  call mma_deallocate(Dv)
  nOrPs = 2 ! use 2 orthonorm. passes for num. accuracy
  call OrthoPAO_Localisation(CMO,nBas,nFro,nOrb2Loc,nSym,nOrPs,Test_OrthoPAO)
else
  write(Txt,'(A,A4)') 'Model = ',Model
  call SysAbendMsg(SecNam,'Unknown model',Txt)
end if

end subroutine Localise_Noniterative
