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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This program generates generally contracted basis sets from multiple *
! atomic densities. The densities are averaged with weights and the    *
! resulting eigenvectors, Natural Orbitals (NO), are used as basis     *
! functions where the corresponding eigenvalues, occupation numbers,   *
! are used as measure of importance.                                   *
! It is also possible to project out the contributions of vectors      *
! from the averaged density before generating the NO's.                *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine GenAno(ireturn)

use Genano_globals, only: isUHF, iProj, nSets
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i

call MkType()
call InpCtl_GenANO()
do i=1,nSets
   call RdCmo()
   call UpDens()
   if (isUHF == 1) then
     call RdCmo()
     call UpDens()
     isUHF = 0
   end if
end do
call SphAve()
if (iProj == 1) call Proj1()
if (iProj == 2) call Proj2()
call MkAno()
call Free_GenANO()

ireturn = 0

return

end subroutine GenAno
