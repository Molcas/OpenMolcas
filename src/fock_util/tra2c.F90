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

subroutine Tra2c(i,iSym,iBas,iAsh,j,jSym,jBas,jAsh,kl_Bas_pairs,ij_Orb_pairs,CMO_i,CMO_j,IJKL,C,TUKL)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: i, iSym, iBas, iAsh, j, jSym, jBas, jAsh, kl_Bas_pairs, ij_Orb_pairs
real(kind=wp), intent(in) :: CMO_i(iBas,iAsh), CMO_j(jBas,jAsh), IJKL(kl_Bas_pairs)
real(kind=wp), intent(out) :: C(ij_Orb_pairs)
real(kind=wp), intent(inout) :: TUKL(kl_Bas_pairs,ij_Orb_pairs)
integer(kind=iwp) :: iA, ijA, jA

!call RecPrt('IJKL',' ',IJKL,1,kl_Bas_Pairs)

if ((iSym == jSym) .and. (i /= j)) then
  ijA = 0
  do iA=1,iAsh
    do jA=1,iA
      ijA = ijA+1
      C(ijA) = CMO_i(i,iA)*CMO_i(j,jA)+CMO_i(j,iA)*CMO_i(i,jA)
    end do
  end do
  !call TriPrt('C',' ',C,iAsh)
else if (iSym == jSym) then
  ijA = 0
  do iA=1,iAsh
    do jA=1,iA
      ijA = ijA+1
      C(ijA) = CMO_i(i,iA)*CMO_i(i,jA)
    end do
  end do
  !call TriPrt('C',' ',C,iAsh)
else
  ijA = 0
  do iA=1,iAsh
    do jA=1,jAsh
      ijA = ijA+1
      C(ijA) = CMO_i(i,iA)*CMO_j(j,jA)
    end do
  end do
  !call RecPrt('C',' ',C,jAsh,iAsh)
end if
!write(u6,*) ' ij_Orb_Pairs=',ij_Orb_Pairs

call DNaXpY(ij_Orb_pairs,kl_Bas_pairs,C,1,IJKL,1,0,TUKL,1,kl_Bas_pairs)

return

end subroutine Tra2c
