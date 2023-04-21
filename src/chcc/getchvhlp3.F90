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

subroutine GetChVHlp3(L2,Tmp,cGrp,deGrp,LunAux)
! this routine does:
! read L2(m,c'de') into L2 from disk file
! structure of files: for each c',de' one file with
! name L2Name(cGrp,deGrp)
! @ citanie zatial odflaknute
!
! parameter description:
! L2     - Array for L2 (O)
! Tmp    - Temporary array of L2 size, used for mapping, if needed
! xGrp   - Groups of c, delta (I)
! LunAux - lun for auxiliary reading (I)

use chcc_global, only: DimGrpa, DimGrpbe, L2Name, nc
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: L2(*), Tmp(*)
integer(kind=iwp), intent(in) :: cGrp, deGrp, LunAux
integer(kind=iwp) :: length
character(len=6) :: LunName

! nacitanie (+expanzia, ak treba)
if (cGrp > deGrp) then
  LunName = L2Name(cGrp,deGrp)
  length = nc*DimGrpa(cGrp)*DimGrpbe(deGrp)
  call GetX(L2,length,LunAux,LunName,1,1)
else if (cGrp == deGrp) then
  LunName = L2Name(cGrp,deGrp)
  length = nc*DimGrpa(cGrp)*(DimGrpbe(deGrp)+1)/2
  call GetX(Tmp,length,LunAux,LunName,1,1)
  length = DimGrpa(cGrp)*(DimGrpbe(deGrp)+1)/2
  call Exp1(Tmp,L2,nc,length,DimGrpa(cGrp))
else
  LunName = L2Name(deGrp,cGrp)
  length = nc*DimGrpa(cGrp)*DimGrpbe(deGrp)
  call GetX(Tmp,length,LunAux,LunName,1,1)
  call Map3_132(Tmp,L2,nc,DimGrpbe(deGrp),DimGrpa(cGrp))
end if

return

end subroutine GetChVHlp3
