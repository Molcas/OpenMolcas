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

subroutine CHO_QUALIFY_2(DIAG,ISYM,ISHLAB,MEM,MEM0,LEFT)
!
! Purpose: qualify diagonals ("qualify until full, then largest").

use Cholesky, only: DiaMin, iiBstR, iiBstRSh, IndRed, iOffq, iQuAB, MaxQual, nnBstR, nnBstRSh, nQual
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: ISYM, ISHLAB, MEM
integer(kind=iwp), intent(inout) :: MEM0, LEFT
integer(kind=iwp) :: I, I1, I2, II, II1, IMAX, J, JJ, JJ1, K, K1, K2, KKMN, MAXQ, NDIM, NUMQ
real(kind=wp) :: XMAX, XMIN
character(len=*), parameter :: SECNAM = 'CHO_QUALIFY_2'

NDIM = NNBSTRSH(ISYM,ISHLAB,2)
if (NDIM > 0) then
  MAXQ = min(MAXQUAL-NQUAL(ISYM),LEFT/NNBSTR(ISYM,2))
  NUMQ = 0
  if (MAXQ > 0) then
    I1 = IIBSTR(ISYM,2)+IIBSTRSH(ISYM,ISHLAB,2)+1
    I2 = I1+NDIM-1
    if (MAXQ == 1) then ! qualify the largest > DIAMIN
      XMAX = DIAMIN(ISYM)
      IMAX = -1
      do I=I1,I2
        J = INDRED(I,2)
        if (DIAG(J) >= XMAX) then
          XMAX = DIAG(J)
          IMAX = I
        end if
      end do
      if (IMAX > 0) then
        NUMQ = NUMQ+1
        iQuAB(IOFFQ(ISYM)+NUMQ,ISYM) = IMAX
      end if
    else ! full search
      do I=I1,I2
        J = INDRED(I,2)
        if (DIAG(J) >= DIAMIN(ISYM)) then
          if (NUMQ < MAXQ) then
            NUMQ = NUMQ+1
            iQuAB(IOFFQ(ISYM)+NUMQ,ISYM) = I
          else if (NUMQ == MAXQ) then
            K1 = IOFFQ(ISYM)+1
            K2 = K1+NUMQ-1
            II1 = IQUAB(K1,ISYM)
            JJ1 = INDRED(II1,2)
            XMIN = DIAG(JJ1)
            KKMN = K1
            do K=K1+1,K2 ! find min. among qualified
              II = IQUAB(K,ISYM)
              JJ = INDRED(II,2)
              if (DIAG(JJ) < XMIN) then
                XMIN = DIAG(JJ)
                KKMN = K
              end if
            end do
            if (DIAG(J) > XMIN) iQuAB(KKMN,ISYM) = I ! replace
          else
            call CHO_QUIT('Logical error in '//SECNAM,104)
          end if
        end if
      end do
    end if
  end if
  NQUAL(ISYM) = NQUAL(ISYM)+NUMQ
  MEM0 = MEM0+NUMQ*NNBSTR(ISYM,2)
  LEFT = MEM-MEM0
end if

end subroutine CHO_QUALIFY_2
