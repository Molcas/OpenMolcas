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
!====
subroutine CalcDB(nMult,nGrdPt,natom,nAtQM,ipIsMM,TTT,DTTT,ExtPot,DB)
! dB = dTTT * V_ext + TTT * dV_ext

implicit real*8(A-H,O-Z)
#include "espf.fh"
dimension TTT(nGrdPt,nMult), DTTT(nMult,nGrdPt,3,nAtQM),ExtPot(10,natom), DB(nGrdPt,3,nAtQM)

iPL = iPL_espf()

if (iPL >= 4) call RecPrt('TTT in calcdb',' ',TTT,nMult,nGrdPt)
nOrd = nMult/nAtQM
do iPnt=1,nGrdPt
  iQM = 0
  do iAt=1,natom
    if (iWork(ipIsMM+iAt-1) /= 0) goto 1
    iQM = iQM+1
    DB(iPnt,1,iQM) = TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(2,iAt)
    DB(iPnt,2,iQM) = TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(3,iAt)
    DB(iPnt,3,iQM) = TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(4,iAt)
    if (nOrd == 4) then
      DB(iPnt,1,iQM) = DB(iPnt,1,iQM)+TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot(5,iAt)+TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(8,iAt)+ &
                       TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(9,iAt)
      DB(iPnt,2,iQM) = DB(iPnt,2,iQM)+TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot(8,iAt)+TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(6,iAt)+ &
                       TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(10,iAt)
      DB(iPnt,3,iQM) = DB(iPnt,3,iQM)+TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot(9,iAt)+TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(10,iAt)+ &
                       TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(7,iAt)
    end if
    jQM = 0
    do jAt=1,natom
      if (iWork(ipIsMM+jAt-1) /= 0) goto 2
      jQM = jQM+1
      do iOrd=1,nOrd
        iMlt = nOrd*(jQM-1)+iOrd
        DB(iPnt,1,iQM) = DB(iPnt,1,iQM)+DTTT(iMlt,iPnt,1,iQM)*ExtPot(iOrd,jAt)
        DB(iPnt,2,iQM) = DB(iPnt,2,iQM)+DTTT(iMlt,iPnt,2,iQM)*ExtPot(iOrd,jAt)
        DB(iPnt,3,iQM) = DB(iPnt,3,iQM)+DTTT(iMlt,iPnt,3,iQM)*ExtPot(iOrd,jAt)
      end do
2     continue
    end do
1   continue
  end do
end do

! Some printing for debug

if (iPL >= 4) then
  do iQM=1,nAtQM
    write(6,*) 'dB/dq_i for i = ',iQM
    do jPnt=1,nGrdPt
      write(6,1234) jPnt,(DB(jPnt,iXYZ,iQM),iXYZ=1,3)
    end do
  end do
end if

return

1234 format(I6,3d13.6)

end
