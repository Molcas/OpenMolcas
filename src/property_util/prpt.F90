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

subroutine Prpt()
!***********************************************************************
!                                                                      *
! Purpose: To set up all calling arguments for the subroutine          *
!          Prpt_ . For RASSCF to work MO coefficents and Occupation    *
!          numbers must be available on TMPORB.                        *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iDummy(1), iError, ipOcc, ipOcc_ab, ipScr, ipVec, ipVec_ab, iSym, iUHF, iWFType, Lu, MaxScr, n2Dim, n2Tot, &
                     nBas(8), nDim, nIrrep, nTriDim
real(kind=wp) :: Dummy(1)
logical(kind=iwp) :: ifallorb, Short, var
character(len=81) :: note
character(len=8) :: Method
character(len=4) :: PrpLst
character(len=2) :: lbl
integer(kind=iwp), external :: isFreeUnit
#include "WrkSpc.fh"

call GetEnvf('MOLCAS_PROPERTIES',PrpLst)
call UpCase(PrpLst)
if (PrPlst(1:3) == 'LON') then
  Short = .false.
  !ifallorb = .True.
else
  Short = .true.
  ifallorb = .false.
end if

! This variable is used so we know if the density we search for is labeled
! variational or not.
var = .false.

call Get_cArray('Relax Method',Method,8)

call Get_iScalar('nSym',nIrrep)

call Get_iArray('nBas',nBas,nIrrep)

nDim = 0
n2Dim = 0
nTriDim = 0
n2Tot = 0
do iSym=1,nIrrep
  nDim = nDim+nBas(iSym)
  nTriDim = nTriDim+nBas(iSym)*(nBas(iSym)+1)/2
  n2Tot = n2Tot+nBas(iSym)**2
end do

ipOcc = ip_Dummy ! dummy initialization
ipOcc_ab = ip_Dummy ! dummy initialization
ipVec = ip_Dummy ! dummy initialization
ipVec_ab = ip_Dummy ! dummy initialization

if ((Method == 'RHF-SCF ') .or. &
    (Method == 'IVO-SCF ') .or. &
    (Method == 'KS-DFT  ') .or. &
    (Method == 'UHF-SCF ')) then
  call Get_iScalar('SCF mode',iUHF)
else
  iUHF = 0
end if

if ((iUHF == 1) .or. (Method == 'RASSCFSA')) then
  call GetMem('Occ','Allo','Real',ipOcc,2*nDim)
  ipOcc_ab = ipOcc+nDim
else
  call GetMem('Occ','Allo','Real',ipOcc,nDim)
end if
if (Short) then
  ipVec = ip_Dummy
  lbl = 'O '
  n2Tot = 1
else
  if ((iUHF /= 1) .or. (Method == 'RASSCFSA')) then
    call GetMem('Vec','Allo','Real',ipVec,n2Tot)
  else
    call GetMem('Vec','Allo','Real',ipVec,2*n2Tot)
    ipVec_ab = ipVec+n2Tot
  end if
  lbl = 'CO'
end if

Lu = 10
Lu = IsFreeUnit(Lu)
if ((Method == 'RHF-SCF ') .or. &
    (Method == 'IVO-SCF ') .or. &
    (Method == 'KS-DFT  ') .or. &
    (Method == 'UHF-SCF ')) then
  if (iUHF /= 1) then
    call RdVec('SCFORB',Lu,Lbl,nIrrep,nBas,nBas,Work(ipVec),Work(ipOcc),Dummy,iDummy,'',0,iError)
  else
    call RdVec_('UHFORB',Lu,Lbl,iUHF,nIrrep,nBas,nBas,Work(ipVec),Work(ipVec_ab),Work(ipOcc),Work(ipOcc_ab),Dummy,Dummy,iDummy,'', &
                1,iError,iWFtype)
    if (Short) then
      do i=0,nDim-1
        Work(ipOcc+i) = Work(ipOcc+i)+Work(ipOcc_ab+i)
      end do
    end if
  end if
else if ((Method == 'RASSCF  ') .or. &
         (Method == 'CASSCF  ') .or. &
         (Method == 'CASDFT  ') .or. &
         (Method == 'CASSCFSA') .or. &
         (Method == 'CASPT2  ') .or. &
         (Method == 'RASSCFSA')) then
  if (Method == 'RASSCFSA') then
    call RdVec_('TMPORB',Lu,Lbl,iUHF,nIrrep,nBas,nBas,Work(ipVec),Work(ipVec_ab),Work(ipOcc),Work(ipOcc_ab),Dummy,Dummy,iDummy,'', &
                1,iError,iWFtype)
    if (Short) then
      do i=0,nDim-1
        Work(ipOcc+i) = Work(ipOcc+i)+Work(ipOcc_ab+i)
      end do
    end if
    var = .false.
  else
    call RdVec('TMPORB',Lu,Lbl,nIrrep,nBas,nBas,Work(ipVec),Work(ipOcc),Dummy,iDummy,note,0,iError)
    if (Note(2:4) == 'var') var = .true.
  end if
else if (Method == 'MBPT2   ') then
  ! MBPT2 has no occupation-numbers at the moment.
  call FZero(Work(ipOcc),nDim)
  var = .true.
else
  write(u6,*) 'Properties not supported for ',Method
end if

MaxScr = nTriDim+nDim*(nDim+1)/2+10+480+4*10
call GetMem('Scr','Allo','Real',ipScr,MaxScr)
call FZero(Work(ipScr),MaxScr)

call Prpt_(nIrrep,nBas,n2Dim,nDim,Work(ipOcc),n2Tot,Work(ipVec),MaxScr,Work(ipScr),var,Short,iUHF,ifallorb)

call GetMem('Scr','Free','Real',ipScr,MaxScr)
call GetMem('Occ','Free','Real',ipOcc,nDim)
if (.not. Short) call GetMem('Vec','Free','Real',ipVec,n2Tot)

return

end subroutine Prpt
