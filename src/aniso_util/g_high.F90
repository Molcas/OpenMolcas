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

subroutine g_high(esom,GRAD,s_som,dipsom,imltpl,dim,Do_structure_abc,cryst,coord,gtens,maxes,iprint)
! this routine calculates the g-tensor and d-tensor in the basis of the any effective spin,
! (coming from 1 molecular term)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: imltpl, dim, iprint
logical, intent(in) :: Do_structure_abc, GRAD
real(kind=8), intent(in) :: esom(dim), cryst(6), coord(3)
real(kind=8), intent(out) :: gtens(3), maxes(3,3)
complex(kind=8), intent(in) :: dipsom(3,dim,dim), s_som(3,dim,dim)
! local variables:
integer :: i

! intializations
if (Iprint > 2) then
  call prMom('G_HIGH:  DIPSOM(l,i,j):',dipsom,dim)
  call prMom('G_HIGH:   S_SOM(l,i,j):',s_som,dim)
end if
!-----------------------------------------------------------------------
write(6,'(/)')
write(6,'(100A)') ('%',i=1,95)
if (mod(dim,2) == 0) then
  write(6,'(5X,A,I2,A,I2,A)') 'CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET',iMLTPL, &
                              ' ( effective S = ',dim-1,'/2)'
else
  write(6,'(5X,A,I2,A,I1,A)') 'CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET',iMLTPL, &
                              ' ( effective S = ',(dim-1)/2,')'
end if
write(6,'(100A)') ('%',i=1,95)
write(6,'(A)') 'The pseudospin is defined in the basis of the following spin-orbit states:'
do i=1,dim
  if (dim > 9) then
    write(6,'(a,i2,a,i2,a,f11.3,a)') 'spin-orbit state',i,'; energy(',i,') = ',ESOM(i),' cm-1.'
  else
    write(6,'(a,i1,a,i1,a,f11.3,a)') 'spin-orbit state ',i,'. energy(',i,') = ',ESOM(i),' cm-1.'
  end if
end do
if (dim == 2) write(6,'(a,f17.10,a)') 'Tunnelling splitting:',ESOM(2)-ESOM(1),' cm-1.'

call G_HIGH_1(iMLTPL,dim,ESOM,GRAD,S_SOM,dipsom,Do_structure_abc,cryst,coord,gtens,maxes,iprint)

return

end subroutine g_high
