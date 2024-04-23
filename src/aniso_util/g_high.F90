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

subroutine g_high(esom,GRAD,s_som,dipsom,imltpl,d,Do_structure_abc,cryst,coord,gtens,maxes,iprint)
! this routine calculates the g-tensor and d-tensor in the basis of the any effective spin,
! (coming from 1 molecular term)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: imltpl, d, iprint
real(kind=wp), intent(in) :: esom(d), cryst(6), coord(3)
logical(kind=iwp), intent(in) :: GRAD, Do_structure_abc
real(kind=wp), intent(out) :: gtens(3), maxes(3,3)
complex(kind=wp), intent(in) :: s_som(3,d,d), dipsom(3,d,d)
integer(kind=iwp) :: i

! intializations
if (Iprint > 2) then
  call prMom('G_HIGH:  DIPSOM(l,i,j):',dipsom,d)
  call prMom('G_HIGH:   S_SOM(l,i,j):',s_som,d)
end if
!-----------------------------------------------------------------------
write(u6,'(/)')
write(u6,'(A)') repeat('%',95)
if (mod(d,2) == 0) then
  write(u6,'(5X,A,I2,A,I2,A)') 'CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET',iMLTPL, &
                               ' ( effective S = ',d-1,'/2)'
else
  write(u6,'(5X,A,I2,A,I1,A)') 'CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET',iMLTPL, &
                               ' ( effective S = ',(d-1)/2,')'
end if
write(u6,'(A)') repeat('%',95)
write(u6,'(A)') 'The pseudospin is defined in the basis of the following spin-orbit states:'
do i=1,d
  if (d > 9) then
    write(u6,'(a,i2,a,i2,a,f11.3,a)') 'spin-orbit state',i,'; energy(',i,') = ',ESOM(i),' cm-1.'
  else
    write(u6,'(a,i1,a,i1,a,f11.3,a)') 'spin-orbit state ',i,'; energy(',i,') = ',ESOM(i),' cm-1.'
  end if
end do
if (d == 2) write(u6,'(a,f17.10,a)') 'Tunnelling splitting:',ESOM(2)-ESOM(1),' cm-1.'

call G_HIGH_1(iMLTPL,d,ESOM,GRAD,S_SOM,dipsom,Do_structure_abc,cryst,coord,gtens,maxes,iprint)

return

end subroutine g_high
