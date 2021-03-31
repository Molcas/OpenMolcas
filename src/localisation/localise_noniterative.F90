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

implicit real*8(a-h,o-z)
character*4 Model
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

character*21 SecNam
parameter(SecNam='Localise_Noniterative')

character*4 myModel
character*6 Namefile
character*80 Txt

logical Test_OrthoPAO, Normalize
parameter(Test_OrthoPAO=.false.)

dimension dum(1), idum(1)

irc = 0
xNrm = 0.0d0

myModel = Model
call UpCase(myModel)
if (myModel == 'CHOL') then
  !if (.not. Silent) then
  write(6,'(/,1X,A)') 'Cholesky localisation'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (decomposition)'
  write(6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  l_Dens = nBas(1)**2
  do iSym=2,nSym
    l_Dens = max(l_Dens,nBas(iSym)**2)
  end do
  call GetMem('Density','Allo','Real',ip_Dens,l_Dens)
  kOffC = ipCMO
  do iSym=1,nSym
    if (nOrb2Loc(iSym) > 0) then
      kOff1 = kOffC+nBas(iSym)*nFro(iSym)
      call GetDens_Localisation(Work(ip_Dens),Work(kOff1),nBas(iSym),nOrb2Loc(iSym))
      call ChoLoc(irc,Work(ip_Dens),Work(kOff1),Thrs,yNrm,nBas(iSym),nOrb2Loc(iSym))
      xNrm = xNrm+yNrm*yNrm
      if (irc /= 0) then
        call GetMem('Density','Free','Real',ip_Dens,l_Dens)
        irc = 1
        xNrm = -9.9d9
        return
      end if
    end if
    kOffC = kOffC+nBas(iSym)**2
  end do
  xNrm = sqrt(xNrm)
  call GetMem('Density','Free','Real',ip_Dens,l_Dens)
else if (myModel == 'PAO ') then
  !if (.not. Silent) then
  write(6,'(/,1X,A)') 'PAO Cholesky localisation'
  write(6,'(1X,A,1X,D12.4,A)') 'Convergence threshold:',Thrs,' (decomposition)'
  write(6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(6,'(1X,A,8(1X,I6))') 'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
  !end if
  l_Dv = nBas(1)**2
  l_R = nBas(1)**2
  do iSym=2,nSym
    l_R = l_R+nBas(iSym)**2
    l_Dv = max(l_Dv,nBas(iSym)**2)
  end do
  call GetMem('R','Allo','Real',ip_R,l_R)
  call GetMem('Dv','Allo','Real',ip_Dv,l_Dv)
  Normalize = .true.
  call GetRawPAOs(Work(ip_R),Work(ipCMO),nBas,nOrb,nFro,nOrb2Loc,nSym,Normalize)
  kSav = 0
  if (AnaPAO) then
    l_DvSav = l_R
    call GetMem('DvSav','Allo','Real',ip_DvSav,l_DvSav)
    kSav = ip_DvSav
  end if
  kOffR = ip_R
  kOffC = ipCMO
  do iSym=1,nSym
    if (nOrb2Loc(iSym) > 0) then
      call GetDens_Localisation(Work(ip_Dv),Work(kOffR),nBas(iSym),nBas(iSym))
      if (AnaPAO) then
        call dCopy_(nBas(iSym)**2,Work(ip_Dv),1,Work(kSav),1)
        kSav = kSav+nBas(iSym)**2
      end if
      kOff1 = kOffC+nBas(iSym)*nFro(iSym)
      call ChoLoc(irc,Work(ip_Dv),Work(kOff1),Thrs,yNrm,nBas(iSym),nOrb2Loc(iSym))
      xNrm = xNrm+yNrm*yNrm
      if (irc /= 0) then
        if (AnaPAO) then
          call GetMem('DvSav','Free','Real',ip_DvSav,l_DvSav)
        end if
        call GetMem('Dv','Free','Real',ip_Dv,l_Dv)
        call GetMem('R','Free','Real',ip_R,l_R)
        irc = 1
        xNrm = -9.9d9
        return
      end if
    end if
    kOffR = kOffR+nBas(iSym)**2
    kOffC = kOffC+nBas(iSym)**2
  end do
  xNrm = sqrt(xNrm)
  if (AnaPAO) then
    call PAO_Analysis(Work(ip_DvSav),Work(ip_R),Work(ipCMO))
    call GetMem('DvSav','Free','Real',ip_DvSav,l_DvSav)
  end if
  write(Namefile,'(A)') 'DPAORB'
  write(Txt,'(80X)')
  write(Txt,'(A)') 'Linearly dependent PAOs'
  Lu_ = isFreeUnit(11)
  call WrVec_Localisation(Namefile,Lu_,'CO',nSym,nBas,nBas,Work(ip_R),Work(ipOcc),dum,idum,Txt)
  !if (.not. Silent) then
  write(6,'(1X,A)') 'The DPAORB file has been written.'
  !end if
  write(Namefile,'(A)') 'IPAORB'
  write(Txt,'(80X)')
  write(Txt,'(A)') 'Linearly independent PAOs'
  Lu_ = isFreeUnit(11)
  call WrVec_Localisation(Namefile,Lu_,'CO',nSym,nBas,nBas,Work(ipCMO),Work(ipOcc),dum,idum,Txt)
  !if (.not. Silent) then
  write(6,'(1X,A)') 'The IPAORB file has been written.'
  !end if
  call GetMem('Dv','Free','Real',ip_Dv,l_Dv)
  call GetMem('R','Free','Real',ip_R,l_R)
  nOrPs = 2 ! use 2 orthonorm. passes for num. accuracy
  call OrthoPAO_Localisation(Work(ipCMO),nBas,nFro,nOrb2Loc,nSym,nOrPs,Test_OrthoPAO)
else
  write(Txt,'(A,A4)') 'Model = ',Model
  call SysAbendMsg(SecNam,'Unknown model',Txt)
end if

end subroutine Localise_Noniterative
