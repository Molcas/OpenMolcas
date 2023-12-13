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

subroutine t3wresult(symi,symj,i,j,eaaa,eaab,eabb,ebbb)
! this routine writes:
!   0) value of SymiMin,Imin,SymJmin,Jmin from which accumulation started
!   1) value of SymI,SymJ,I,J
!   2) present stage of energies
! into T3tEne file and overwrites previous values

use CCT3_global, only: imin, jmin, symimin, symjmin
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: symi, symj, i, j
real(kind=wp), intent(in) :: eaaa, eaab, eabb, ebbb
integer(kind=iwp) :: lun

lun = 1
call Molcas_Open(lun,'T3tEne')
!open(unit=lun,file='T3tEne')

write(lun,97) symimin,imin,symjmin,jmin
write(lun,98) symi,symj
write(lun,98) i,j
write(lun,99) eaaa
write(lun,99) eaab
write(lun,99) eabb
write(lun,99) ebbb

close(lun)

return

97 format(2x,4(i4,2x))
98 format(2x,2(i4,2x))
99 format(2x,f22.16)

end subroutine t3wresult
