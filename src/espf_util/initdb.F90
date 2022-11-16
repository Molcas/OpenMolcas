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

subroutine InitDB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipDB,ipIsMM,iRMax,DeltaR,iGrdTyp,ipDGrd)
! Compute DB = d(ExtPot[TtT^-1]Tt)/dq

implicit real*8(A-H,O-Z)
#include "espf.fh"

! Various derivatives

call GetMem('DTTT','Allo','Real',ipDTTT,nMult*nGrdPt*nAtQM*3)
call GetMem('DT','Allo','Real',ipDT,nGrdPt*nMult*3*nAtQM)
call GetMem('DTT','Allo','Real',ipDTT,nMult*nMult*nAtQM*3)
call GetMem('DTTTT','Allo','Real',ipDTTTT,nMult*nMult*nAtQM*3)

call CalcDT(nMult,nGrdPt,natom,nAtQM,ipIsMM,iGrdTyp,Work(ipCord),Work(ipGrid),Work(ipDGrd),Work(ipT),Work(ipTT),Work(ipDT), &
            Work(ipDTT),Work(ipDTTTT),Work(ipDTTT))

call GetMem('DTTTT','Free','Real',ipDTTTT,nMult*nMult*nAtQM*3)
call GetMem('DTT','Free','Real',ipDTT,nMult*nMult*nAtQM*3)
call GetMem('DT','Free','Real',ipDT,nGrdPt*nMult*3*nAtQM)

! Finally dB = dTTT * V_ext + TTT * dV_ext

call CalcDB(nMult,nGrdPt,natom,nAtQM,ipIsMM,Work(ipTTT),Work(ipDTTT),Work(ipExt),Work(ipDB))
call GetMem('DTTT','Free','Real',ipDTTT,nMult*nGrdPt*nAtQM*3)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(iRMax)
  call Unused_real(DeltaR)
end if

end subroutine InitDB
