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

subroutine InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipB,ipIsMM)
! Compute the electrostatic tensor matrix between QM atoms and grid points
! Then form the B = ExtPot[TtT^-1]Tt vector

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nMult, natom, nAtQM, nGrdPt, ipCord, ipGrid, ipT, ipTT, ipTTT, ipExt, ipB, ipIsMM
#include "WrkSpc.fh"
integer(kind=iwp) :: iCur, iMlt, iPL, iPnt, ipScr, iQM, J, jAt, jMlt, jPnt, kMlt, kPnt, nOrd
real(kind=wp) :: Det, R, R3, X, Y, Z
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()
nOrd = nMult/nAtQM

! T

do iPnt=1,nGrdPt
  iQM = 0
  do jAt=1,natom
    if (iWork(ipIsMM+jAt-1) == 1) cycle
    iQM = iQM+1
    X = Work(ipGrid+(IPnt-1)*3)-Work(ipCord+(jAt-1)*3)
    Y = Work(ipGrid+(IPnt-1)*3+1)-Work(ipCord+(jAt-1)*3+1)
    Z = Work(ipGrid+(IPnt-1)*3+2)-Work(ipCord+(jAt-1)*3+2)
    R = sqrt(X*X+Y*Y+Z*Z)
    iCur = ipT+(iPnt-1)*nMult+nOrd*(iQM-1)
    Work(iCur) = One/R
    if (nOrd > 1) then
      R3 = R*R*R
      Work(iCur+1) = X/R3
      Work(iCur+2) = Y/R3
      Work(iCur+3) = Z/R3
    end if
  end do
end do
if (iQM /= nAtQM) then
  write(u6,'(A,I4,A4,I4)') ' Error in espf/initb: iQM != nAtQM ',iQM,' != ',nAtQM
  call Abend()
end if

!call RecPrt('T',' ',Work(ipT),nGrdPt,nMult)
!nMax = Max(nGrdPt,nMult)
!call Allocate_Work(ipU,nMax*nMult)
!call Allocate_Work(ipW,nMult)
!call Allocate_Work(ipV,nMax*nMult)
!call Allocate_Work(ipScr,nMult)
!call SVD(nMax,nGrdPt,nMult,Work(ipT),Work(ipW),.true.,Work(ipU),.true.,Work(ipV),iErr,Work(ipScr))
!write(u6,*) 'iErr=',iErr
!call RecPrt('U',' ',Work(ipU),nMax,nMult)
!call RecPrt('w',' ',Work(ipW),nMult,1)
!call RecPrt('V',' ',Work(ipV),nMax,nMult)
!call Free_Work(ipW)
!call Free_Work(ipU)
!call Free_Work(ipV)
!call Free_Work(ipScr)

! TtT

do iMlt=1,nMult
  do jMlt=1,nMult
    iCur = ipTT+(iMlt-1)*nMult+(jMlt-1)
    Work(iCur) = Zero
    do kPnt=1,nGrdPt
      Work(iCur) = Work(iCur)+Work(ipT+(kPnt-1)*nMult+(iMlt-1))*Work(ipT+(kPnt-1)*nMult+(jMlt-1))
    end do
  end do
end do

! TtT^-1

call Allocate_Work(ipScr,nMult*nMult)
call minv(Work(ipTT),Work(ipScr),Det,nMult)
call dCopy_(nMult*nMult,Work(ipScr),1,Work(ipTT),1)
call Free_Work(ipScr)

! [TtT^-1]Tt

do iMlt=1,nMult
  do jPnt=1,nGrdPt
    iCur = ipTTT+(iMlt-1)*nGrdPt+(jPnt-1)
    Work(iCur) = Zero
    do kMlt=1,nMult
      Work(iCur) = Work(iCur)+Work(ipTT+(iMlt-1)*nMult+(kMlt-1))*Work(ipT+(jPnt-1)*nMult+(kMlt-1))
    end do
  end do
end do
if (iPL >= 4) call RecPrt('(TtT)^(-1)Tt matrix in InitB',' ',Work(ipTTT),nMult,nGrdPt)

! B = ExtPot[TtT^-1]Tt

do iPnt=1,nGrdPt
  iQM = 0
  iCur = ipB+(iPnt-1)
  Work(iCur) = Zero
  do jAt=1,natom
    if (iWork(ipIsMM+jAt-1) == 1) cycle
    iQM = iQM+1
    Work(iCur) = Work(iCur)+Work(ipExt+(jAt-1)*10)*Work(ipTTT+nOrd*(iQM-1)*nGrdPt+(iPnt-1))
    if (nOrd > 1) Work(iCur) = Work(iCur)+Work(ipExt+(jAt-1)*10+1)*Work(ipTTT+(nOrd*(iQM-1)+1)*nGrdPt+(iPnt-1))+ &
                               Work(ipExt+(jAt-1)*10+2)*Work(ipTTT+(nOrd*(iQM-1)+2)*nGrdPt+(iPnt-1))+ &
                               Work(ipExt+(jAt-1)*10+3)*Work(ipTTT+(nOrd*(iQM-1)+3)*nGrdPt+(iPnt-1))
  end do
end do
if (iPL >= 4) then
  write(u6,'(A)') ' In InitB (grid coordinates, B value)'
  do iPnt=1,nGrdPt
    write(u6,1234) iPnt,(Work(ipGrid+(iPnt-1)*3+J),J=0,2),Work(ipB+iPnt-1)
  end do
end if

return

1234 format(I4,4F12.6)

end subroutine InitB
