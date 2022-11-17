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

subroutine espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,ipIsMM,ipExt,iPL)
! Compute the expected values of the ESPF operators
! i.e. charges or charges+dipoles

use espf_global, only: MxExtPotComp
use Constants, only: Zero
use Definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp) :: natom, MltOrd, nMult, nGrdPt, ipTTT, ipMltp, ipGrid, ipIsMM, ipExt, iPL
#include "Molcas.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iAddPot, iAtom, iCur, iMult, iOff, ip_Charge, ipD2, ipEI, jMlt, kOrd, kPnt, ncmp
real(kind=wp) :: opnuc(1), SumOfChg, TotElecInt
character(len=LenIn) :: CName(MxAtom) ! IFG
character(len=3), parameter :: Axis(3) = [' x ',' y ',' z ']

if (iPL >= 5) then
  write(u6,*) ' In espf_mltp:',MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,ipIsMM
  call RecPrt('TTT',' ',Work(ipTTT),nGrdPt,nMult)
end if
call GetMem('Nuclear charge','Allo','Real',ip_Charge,natom)
call Get_Nuc_Charge_All(Work(ip_Charge),natom)
iMult = 0
do iAtom=0,natom-1
  if (iWork(ipIsMM+iAtom) == 0) then
    Work(ipMltp+iMult) = Work(ip_Charge+iAtom)
    if (MltOrd /= 1) then
      do iOff=1,MltOrd-1
        Work(ipMltp+iMult+iOff) = Zero
      end do
    end if
    iMult = iMult+MltOrd
  end if
end do
call GetMem('Nuclear charge','Free','Real',ip_Charge,natom)
opnuc = Zero
ncmp = 1
iAddPot = -2
call GetMem('dESPF2','Allo','Real',ipD2,nGrdPt)
call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipD2),nGrdPt,iAddPot)
if (iPL >= 5) call RecPrt('PV',' ',Work(ipD2),nGrdPt,1)
do jMlt=0,nMult-1
  do kPnt=0,nGrdPt-1
    Work(ipMltp+jMlt) = Work(ipMltp+jMlt)+Work(ipD2+kPnt)*Work(ipTTT+jMlt*nGrdPt+kPnt)
  end do
end do
call GetMem('dESPF2','Free','Real',ipD2,nGrdPt)

if (iPL >= 3) then
  write(u6,'(/,A,/)') '      Expectation values of the ESPF operators:'
  call GetMem('ElecInt','Allo','Real',ipEI,natom)
  call Get_CArray('Unique Atom Names',CName,LenIn*natom)
  SumOfChg = Zero
  jMlt = 0
  TotElecInt = Zero
  do iAtom=0,natom-1
    Work(ipEI+iAtom) = Zero
    iCur = ipExt+iAtom*MxExtPotComp
    if (iWork(ipIsMM+iAtom) == 1) cycle
    do kOrd=0,MltOrd-1
      if (kOrd == 0) then
        write(u6,1000) CName(iAtom+1),Work(ipMltp+jMlt)
        SumOfChg = SumOfChg+Work(ipMltp+jMlt)
      else
        write(u6,1001) Axis(kOrd),Work(ipMltp+jMlt+kOrd)
      end if
      Work(ipEI+iAtom) = Work(ipEI+iAtom)+Work(ipMltp+jMlt+kOrd)*Work(iCur+kOrd)
    end do
    jMlt = jMlt+MltOrd
    TotElecInt = TotElecInt+Work(ipEI+iAtom)
  end do
  write(u6,1002) SumOfChg
  write(u6,1003) TotElecInt
  do iAtom=0,natom-1
    if (iWork(ipIsMM+iAtom) == 0) write(u6,1004) CName(iAtom+1),Work(ipEI+iAtom)
  end do
  write(u6,'(A)')
  call GetMem('ElecInt','Free','Real',ipEI,natom)
end if

return

1000 format('        Charge on ',A,'      = ',F10.4)
1001 format('        + Dipole component ',A3,'= ',F10.4)
1002 format(/,'      Total ESPF charge     = ',F10.4,/)
1003 format(/,'      Total ESPF QM/MM interaction energy = ',F10.6,/)
1004 format('        ',A,' individual contribution =',F10.6)

end subroutine espf_mltp
