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

subroutine Cho_ChkDia_A4(Diag,Dmax,iSym,nNeg,nNegT,nConv,xM,yM,zM)
!
! Purpose: check for negative diagonals.
!          Dmax is the max. diagonal (global), used for screening.
!          On exit,
!          nNeg = #negative zeroed diagonals.
!          nNegT = #negative diagonals
!          nConv = #diagonal elements < ThrCom
!          xM = max. element in Diag
!          yM = min. element in Diag
!          zM = max. abs. element in Diag

use Cholesky, only: Damp, iiBstR, IndRed, LuPri, nnBstR, SCDIAG, ThrCom, ThrNeg, TOONEG, WARNEG
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Diag(*)
real(kind=wp), intent(in) :: Dmax
integer(kind=iwp), intent(in) :: iSym
integer(kind=iwp), intent(out) :: nNeg, nNegT, nConv
real(kind=wp), intent(out) :: xM, yM, zM
integer(kind=iwp) :: i, j, j1, j2
real(kind=wp) :: Tst
character(len=*), parameter :: SecNam = 'Cho_ChkDia_A4'

nNeg = 0
nNegT = 0
nConv = 0
if (nnBstR(iSym,2) > 0) then
  xM = -9.9e9_wp
  yM = 9.9e9_wp
else
  xM = Zero
  yM = Zero
end if

j1 = iiBstR(iSym,2)+1
j2 = j1+nnBstR(iSym,2)-1
do j=j1,j2
  i = IndRed(j,2) ! addr in 1st reduced set
  xM = max(xM,Diag(i))
  yM = min(yM,Diag(i))
  if (Diag(i) < Zero) then
    nNegT = nNegT+1
    if (Diag(i) < ThrNeg) then
      nNeg = nNeg+1
      if (Diag(i) < TooNeg) then
        write(Lupri,'(A,A,I12,1X,1P,D16.8)') SecNam,': diagonal too negative: ',i,Diag(i)
        write(Lupri,'(A,A)') SecNam,': shutting down Cholesky decomposition!'
        call Cho_Quit('Diagonal too negative in '//SecNam,104)
      end if
      if (Diag(i) < WarNeg) write(Lupri,'(A,A,I12,1X,1P,D16.8,A)') SecNam,': Negative diagonal: ',i,Diag(i),' (zeroed)'
      Diag(i) = Zero
    end if
  end if
end do

zM = max(abs(xM),abs(yM))

do j=j1,j2
  i = IndRed(j,2)
  Tst = sqrt(abs(Diag(i))*Dmax)*Damp(2)
  if (Tst <= ThrCom) then
    nConv = nConv+1
    if (ScDiag) Diag(i) = Zero
  end if
end do

end subroutine Cho_ChkDia_A4
