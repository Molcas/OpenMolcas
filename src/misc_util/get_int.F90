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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!**************************************************************
!
! Author :   F. Aquilante
!
!  Get_Int :  driver for the integral generator from Cholesky
!             vectors
!**************************************************************

subroutine Get_Int(rc,iOpt,iSymp,iSymq,iSymr,iSyms,Xint,lBuf,nMat)

implicit real*8(a-h,o-z)
integer rc, iOpt
integer iSymp, iSymq, iSymr, iSyms, Npq, Nrs
real*8 Xint(*)
integer lBuf, nMat
#include "RdOrd.fh"
#include "TwoDat.fh"
character*4 BaseNm
character*6 Fname
parameter(BaseNm='CHFV')
! Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1

! Check input parameters

rc = 0
if ((iOpt /= 1) .and. (iOpt /= 2)) then
  rc = rcRD06
  write(6,*) 'Get_Int: Invalid option'
  write(6,*) 'iOpt= ',iOpt
  call Abend()
end if
if ((iSymp < iSymq) .or. (iSymr < iSyms)) then
  rc = rcRD02
  write(6,*) 'Get_Int: invalid order of symmetry labels'
  call Abend()
end if
if (MulD2h(iSymp,iSymq) /= MulD2h(iSymr,iSyms)) then
  rc = rcRD01
  write(6,*) 'Get_Int: wrong symmetry labels, direct product is not total symmetric'
  call Abend()
end if
if (lBuf < 1) then
  rc = rcRD04
  write(6,*) 'Get_Int: invalid buffer size'
  write(6,*) 'lBuf=',lBuf
  call Abend()
end if

! Open files.
LuCVec(1) = 7
write(Fname,'(A4,I1,I1)') BaseNm,iSymp,iSymq
call DANAME_MF_WA(LuCVec(1),Fname)
if (iSymp /= iSymr) then
  LuCVec(2) = 7
  write(Fname,'(A4,I1,I1)') BaseNm,iSymr,iSyms
  call DANAME_MF_WA(LuCVec(2),Fname)
else
  LuCVec(2) = -1
end if

if (iSymp == iSymq) then
  Npq = nBas(iSymp)*(nBas(iSymp)+1)/2
else
  Npq = nBas(iSymp)*nBas(iSymq)
end if
if (iSymr == iSyms) then
  Nrs = nBas(iSymr)*(nBas(iSyms)+1)/2
else
  Nrs = nBas(iSymr)*nBas(iSyms)
end if

! For debug
!lBufs = lBuf
!lBuf = min(Nrs*min(Npq,2)+1,lBufs)

if (iOpt == 1) then
  pq1 = 1
  nMat = min(Npq,(lBuf-1)/Nrs)
else
  if ((pq1 < 1) .or. (pq1 > Npq)) then
    rc = 999999
    write(6,*) 'pq1 out of bounds: ',pq1
    call Abend()
    nMat = 99999999
  else
    nMat = min((Npq-pq1+1),(lBuf-1)/Nrs)
  end if
end if

if (nMat < 1) then
  return ! no more integrals to compute
else
  call GEN_INT(rc,iSymp,iSymq,iSymr,iSyms,pq1,nMat,Xint)
end if
pq1 = pq1+nMat

! Close files.
do i=1,2
  if (LuCVec(i) /= -1) then
    call DACLOS(LuCVec(i))
    LuCVec(i) = -1
  end if
end do

! LBuf = lBufs

return

end subroutine Get_Int
