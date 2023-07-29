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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

function Cho_VecBuf_Integrity_OK(Tol,Report)
!
! Thomas Bondo Pedersen, September 2012.
!
! Check Cholesky vector buffer integrity: compute norm and sum of
! vectors in the buffer and compare these values to the table
! generated at buffer initialization.

use Cholesky, only: CHVBFI, CHVBUF, InfVec, ip_CHVBFI_SYM, ip_CHVBUF_SYM, l_CHVBFI_SYM, LuPri, nDimRS, nSym, nVec_in_Buf
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: Cho_VecBuf_Integrity_OK
real(kind=wp) :: Tol
logical(kind=iwp) :: Report
integer(kind=iwp) :: ipV, iSym, jRed, jVec, n, nErr
logical(kind=iwp) :: OK
real(kind=wp) :: Nrm, Sm
real(kind=wp), external :: Cho_dSumElm, dDot_
! Statement functions
real(kind=wp) :: RefNrm, RefSm
integer(kind=iwp) :: i, j
RefNrm(i,j) = CHVBFI(ip_ChVBfI_Sym(j)+2*(i-1))
RefSm(i,j) = CHVBFI(ip_ChVBfI_Sym(j)+2*(i-1)+1)

nErr = 0
if (allocated(CHVBUF) .and. allocated(CHVBFI) .and. allocated(nDimRS)) then
  do iSym=1,nSym
    if ((nVec_in_Buf(iSym) > 0) .and. (l_ChVBfI_Sym(iSym) > 0)) then
      ipV = ip_ChVBuf_Sym(iSym)
      do jVec=1,nVec_in_Buf(iSym)
        jRed = InfVec(jVec,2,iSym)
        n = nDimRS(iSym,jRed)
        Nrm = sqrt(dDot_(n,CHVBUF(ipV),1,CHVBUF(ipV),1))
        Sm = Cho_dSumElm(CHVBUF(ipV),n)
        OK = (abs(Nrm-RefNrm(jVec,iSym)) < Tol) .and. (abs(Sm-RefSm(jVec,iSym)) < Tol)
        if (.not. OK) then
          nErr = nErr+1
          if (Report) then
            write(LuPri,'(A,I7,A,I2,A,I9)') 'Buffer corrupted: vector',jVec,' sym.',iSym,' dim.',n
            write(LuPri,'(3X,1P,3(A,D25.16))') 'Norm=',Nrm,' Reference=',RefNrm(jVec,iSym),' Diff=',Nrm-RefNrm(jVec,iSym)
            write(LuPri,'(3X,1P,3(A,D25.16))') 'Sum= ',Sm,' Reference=',RefSm(jVec,iSym),' Diff=',Sm-RefSm(jVec,iSym)
          end if
        end if
        ipV = ipV+n
      end do
    end if
  end do
end if
if (Report) then
  if (nErr /= 0) then
    write(LuPri,'(A,I7,A,1P,D25.16)') 'Buffer corrupted for ',nErr,' vectors. Tolerance=',Tol
  else
    write(LuPri,'(A,1P,D25.16)') 'Buffer integrity OK. Tolerance=',Tol
  end if
end if
Cho_VecBuf_Integrity_OK = nErr == 0

end function Cho_VecBuf_Integrity_OK
