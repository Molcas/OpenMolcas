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

subroutine CISX(IDX,D,DS,PS,PA,SCR)

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: NAC, NACPR2
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IDX(NAC)
real(kind=wp), intent(inout) :: D(*), DS(*), PS(*), PA(*)
real(kind=wp), intent(_OUT_) :: SCR(*)
integer(kind=iwp) :: I, ICASE, IJKLN, IJKLO, IJNEW, IJO, J, K, KLNEW, L, LLIM, NIJ, NIJKL
real(kind=wp) :: SGN, SGN0

! Convert from CI to SX ordering
! Note: A factor of 2 arises because matrices are folded

! one-body density
IJO = 0
NIJ = nTri_Elem(NAC)
NIJKL = nTri_Elem(NIJ)
do I=1,NAC
  do J=1,I
    IJNEW = iTri(IDX(I),IDX(J))
    IJO = IJO+1
    SCR(IJNEW) = D(IJO)
  end do
end do
D(1:NIJ) = SCR(1:NIJ)

! spin density
IJO = 0
NIJ = nTri_Elem(NAC)
NIJKL = nTri_Elem(NIJ)
do I=1,NAC
  do J=1,I
    IJNEW = iTri(IDX(I),IDX(J))
    IJO = IJO+1
    SCR(IJNEW) = DS(IJO)
  end do
end do
DS(1:NIJ) = SCR(1:NIJ)

! symmetrized (iCase=1) and antisymmetrized (iCase=2) two-body density
do ICASE=1,2
  IJKLO = 0
  SCR(1:NACPR2) = Zero
  do I=1,NAC
    do J=1,I
      if (IDX(J) > IDX(I)) then
        SGN0 = -One
      else
        SGN0 = One
      end if
      IJNEW = iTri(IDX(I),IDX(J))
      do K=1,I
        LLIM = K
        if (K == I) LLIM = J
        do L=1,LLIM
          IJKLO = IJKLO+1
          if (IDX(L) > IDX(K)) then
            SGN = -SGN0
          else
            SGN = SGN0
          end if
          KLNEW = iTri(IDX(K),IDX(L))
          IJKLN = iTri(IJNEW,KLNEW)
          if (ICASE == 1) then
            if (KLNEW > IJNEW) then
              if ((K == L) .and. (I /= J)) then
                SCR(IJKLN) = Two*PS(IJKLO)
              else if ((I == J) .and. (K /= L)) then
                SCR(IJKLN) = Half*PS(IJKLO)
              else
                SCR(IJKLN) = PS(IJKLO)
              end if
            else
              SCR(IJKLN) = PS(IJKLO)
            end if
          else
            SCR(IJKLN) = SGN*PA(IJKLO)
          end if
        end do
      end do
    end do
  end do
  if (ICASE == 1) then
    PS(1:NIJKL) = SCR(1:NIJKL)
  else
    PA(1:NIJKL) = SCR(1:NIJKL)
  end if

end do

end subroutine CISX
