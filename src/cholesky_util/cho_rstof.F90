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

subroutine CHO_RSTOF(IRS2F,N,LRDIM,IRED)
!
! Purpose: set up mapping between reduced set and SO indices
!          (i.e., full storage).
!
! IRS2F(1,irs) = alpha (SO index, not symmmetry reduced)
! IRS2F(2,irs) = beta  (SO index, not symmmetry reduced)

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Cholesky, only: iBas, iShlSO, iSOShl, mmBstRT, nBas, nBstSh, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: N, LRDIM, IRED
integer(kind=iwp), intent(out) :: IRS2F(N,LRDIM)
integer(kind=iwp) :: IA, IB, IRS, ISHLA, ISHLAB, ISHLB, ISYMA, ISYMAB, ISYMB, KA, KB, LA, LAB, LB
character(len=*), parameter :: SECNAM = 'CHO_RSTOF'
integer(kind=iwp), external :: CHO_F2SP, CHO_RS2F

if (N < 2) call CHO_QUIT('Dimension error [1] in '//SECNAM,104)
if (LRDIM /= MMBSTRT) call CHO_QUIT('Dimension error [2] in '//SECNAM,104)
IRS2F(:,1:MMBSTRT) = 0

do ISYMA=1,NSYM
  if (NBAS(ISYMA) > 0) then
    do ISYMB=1,ISYMA-1
      ISYMAB = MUL(ISYMA,ISYMB)
      do KB=1,NBAS(ISYMB)
        IB = IBAS(ISYMB)+KB
        LB = ISHLSO(IB)
        ISHLB = ISOSHL(IB)
        do KA=1,NBAS(ISYMA)
          IA = IBAS(ISYMA)+KA
          LA = ISHLSO(IA)
          ISHLA = ISOSHL(IA)
          if (ISHLA < ISHLB) then
            LAB = NBSTSH(ISHLB)*(LA-1)+LB
          else if (ISHLA == ISHLB) then
            LAB = ITRI(LA,LB)
          else
            LAB = NBSTSH(ISHLA)*(LB-1)+LA
          end if
          ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
          if (ISHLAB > 0) then
            IRS = CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
            if (IRS > 0) then
              IRS2F(1,IRS) = IA
              IRS2F(2,IRS) = IB
            end if
          end if
        end do
      end do
    end do
    ISYMB = ISYMA
    ISYMAB = 1
    do KA=1,NBAS(ISYMA)
      IA = IBAS(ISYMA)+KA
      LA = ISHLSO(IA)
      ISHLA = ISOSHL(IA)
      do KB=1,KA
        IB = IBAS(ISYMB)+KB
        LB = ISHLSO(IB)
        ISHLB = ISOSHL(IB)
        if (ISHLA < ISHLB) then
          LAB = NBSTSH(ISHLB)*(LA-1)+LB
        else if (ISHLA == ISHLB) then
          LAB = ITRI(LA,LB)
        else
          LAB = NBSTSH(ISHLA)*(LB-1)+LA
        end if
        ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
        if (ISHLAB > 0) then
          IRS = CHO_RS2F(LAB,ISHLAB,ISYMAB,IRED)
          if (IRS > 0) then
            IRS2F(1,IRS) = IA
            IRS2F(2,IRS) = IB
          end if
        end if
      end do
    end do
  end if
end do

end subroutine CHO_RSTOF
