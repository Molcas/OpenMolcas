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

subroutine CHO_REOVEC(IRS2F,N,LRDIM,WRK,LWRK)
!
! Purpose: reorder Cholesky vectors on disk to full storage.

use Index_Functions, only: iTri
use Cholesky, only: iBas, nBas, mmBstRT, nnBstRT
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, LRDIM, LWRK
integer(kind=iwp), intent(out) :: IRS2F(N,LRDIM)
real(kind=wp), intent(inout) :: WRK(LWRK)
integer(kind=iwp) :: IA, IB, IRS1, ISYMA, ISYMB, JA, JAB, JB
character(len=*), parameter :: SECNAM = 'CHO_REOVEC'
integer(kind=iwp), external :: CHO_ISAO

! Set up mapping from rs1 to full storage.
! ----------------------------------------

if (N < 3) call CHO_QUIT('Dimension error [1] in '//SECNAM,104)
if (LRDIM /= MMBSTRT) call CHO_QUIT('Dimension error [2] in '//SECNAM,104)
call CHO_RSTOF(IRS2F,N,MMBSTRT,1)
do IRS1=1,NNBSTRT(1)
  IA = IRS2F(1,IRS1)
  IB = IRS2F(2,IRS1)
  ISYMA = CHO_ISAO(IA)
  ISYMB = CHO_ISAO(IB)
  JA = IA-IBAS(ISYMA)
  JB = IB-IBAS(ISYMB)
  if (ISYMA == ISYMB) then
    JAB = ITRI(JA,JB)
  else
    JAB = NBAS(ISYMA)*(JB-1)+JA
  end if
  IRS2F(1,IRS1) = ISYMA
  IRS2F(2,IRS1) = ISYMB
  IRS2F(3,IRS1) = JAB
end do

! Set up index arrays and open files for storing full vectors.
! ------------------------------------------------------------

call CHO_REOINI()

! Reorder vectors on disk.
! ------------------------

call CHO_REOVC1(IRS2F,N,LRDIM,WRK,LWRK)

end subroutine CHO_REOVEC
