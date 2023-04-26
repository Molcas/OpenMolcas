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

subroutine CalcDB(nMult,nGrdPt,natom,nAtQM,IsMM,TTT,DTTT,ExtPot,DB)
! dB = dTTT * V_ext + TTT * dV_ext

use espf_global, only: MxExtPotComp
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMult, nGrdPt, natom, nAtQM, IsMM(natom)
real(kind=wp), intent(in) :: TTT(nGrdPt,nMult), DTTT(nMult,nGrdPt,3,nAtQM), ExtPot(MxExtPotComp,natom)
real(kind=wp), intent(out) :: DB(nGrdPt,3,nAtQM)
integer(kind=iwp) :: iAt, iMlt, iOrd, iPL, iPnt, iQM, jAt, jPnt, jQM, nOrd
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()

if (iPL >= 4) call RecPrt('TTT in calcdb',' ',TTT,nMult,nGrdPt)
nOrd = nMult/nAtQM
do iPnt=1,nGrdPt
  iQM = 0
  do iAt=1,natom
    if (IsMM(iAt) /= 0) cycle
    iQM = iQM+1
    DB(iPnt,:,iQM) = TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(2:4,iAt)
    if (nOrd == 4) then
      DB(iPnt,1,iQM) = DB(iPnt,1,iQM)+TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot(5,iAt)+ &
                                      TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(8,iAt)+ &
                                      TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(9,iAt)
      DB(iPnt,2,iQM) = DB(iPnt,2,iQM)+TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot(8,iAt)+ &
                                      TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(6,iAt)+ &
                                      TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(10,iAt)
      DB(iPnt,3,iQM) = DB(iPnt,3,iQM)+TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot(9,iAt)+ &
                                      TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(10,iAt)+ &
                                      TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(7,iAt)
    end if
    jQM = 0
    do jAt=1,natom
      if (IsMM(jAt) /= 0) cycle
      jQM = jQM+1
      do iOrd=1,nOrd
        iMlt = nOrd*(jQM-1)+iOrd
        DB(iPnt,:,iQM) = DB(iPnt,:,iQM)+DTTT(iMlt,iPnt,:,iQM)*ExtPot(iOrd,jAt)
      end do
    end do
  end do
end do

! Some printing for debug

if (iPL >= 4) then
  do iQM=1,nAtQM
    write(u6,*) 'dB/dq_i for i = ',iQM
    do jPnt=1,nGrdPt
      write(u6,1234) jPnt,DB(jPnt,:,iQM)
    end do
  end do
end if

return

1234 format(I6,3d13.6)

end subroutine CalcDB
