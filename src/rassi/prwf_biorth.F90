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
subroutine prwf_biorth(nstate, njob, nconf, ndet, nasht, detocc, &
                       detcoeff,cithr)
! print CI expansion in determinant biorthonormal basis
use Definitions, only: iwp, wp, u6
implicit none
integer(kind=iwp), intent(in) :: nstate, njob, nconf, ndet, nasht
real(kind=wp), dimension(ndet), intent(in) :: detcoeff
character(len=(NASHT+1)), dimension(ndet), intent(in) :: detocc
integer(kind=iwp) :: i, ocsp
character(len=38) :: fomt
real*8 cithr
write(u6,*) ' ******* TRANSFORMED CI COEFFICIENTS *******'
write(u6,*) ' CI for state ', nstate
write(u6,*) ' This is on JobIph nr.', njob
write(u6,*) ' Its length NCI=', nconf
write(u6,*) ' Its length NDET=', ndet
if (nDet > 1) then
  ocsp = MAX(9,NASHT)
  write(fomt,'(A,I2,A)') '(I7,A16,A', ocsp, ',A5,G17.10,A5,G17.10)'
  write(u6,*) ' Occupation of active orbitals, and spin'
  write(u6,*) ' of open shells. (u,d: Spin up or down).'
  write(u6,'(A,A,A)') '    Det  ','                       ', &
                      '       Coef       Weight'
  do i=1,nDet
    if(ABS(detcoeff(i)).GT.cithr) then
      write(u6,fomt) i, '                 ', &
                    trim(detocc(i)), '     ', &
                    detcoeff(i), '     ', detcoeff(i)**2
    end if
  end do
  write(u6,*)('*',i=1,80)
end if
end subroutine prwf_biorth
