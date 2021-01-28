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

! Print the value of (unnormalized) vibrational wave functions.
!PAM97 For MOLCAS-4
subroutine prwf_vibrot(ndim,R)

use Vibrot_globals, only: J1A, J2A, nvib1, Vibwvs, iad12
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndim
real(kind=wp), intent(in) :: R(*)
integer(kind=iwp) :: ne1, ndim1, nwork, J1, Jad1, ist, iadr1, i, nv, nvsta, nvend
real(kind=wp), allocatable :: vib(:)

write(u6,*)
call CollapseOutput(1,'PRINTOUT OF VIBRATIONAL WAVE FUNCTIONS')
write(u6,*)

ne1 = nvib1+1
ndim1 = ndim+1
nwork = ne1*ndim1

call mma_allocate(Vib,nwork,label='Vib')

! Loop over rotational quantum numbers
do J1=J1A,J2A
  Jad1 = J1-J1A+1
  write(u6,1001) ' Rotational quantum number J=',J1

  ! Read vibrational functions for this J-value.
  ist = 1
  iadr1 = iad12(Jad1)
  do nv=1,ne1
    call DDafile(Vibwvs,2,vib(ist),ndim1,iadr1)
    ist = ist+ndim1
  end do

  ! Write out the wave functions:
  do nvsta=0,ne1-1,5
    nvend = min(nvsta+4,ne1-1)
    write(u6,*)
    write(u6,1002) 'Radial dist.',(nv,nv=nvsta,nvend)
    do i=1,ndim
      write(u6,1003) r(i),(vib(i+ndim1*nv),nv=nvsta,nvend)
    end do
  end do

! End of loop over J1.
end do

call mma_deallocate(vib)

call CollapseOutput(0,'PRINTOUT OF VIBRATIONAL WAVE FUNCTIONS')
write(u6,*)

return

1001 format(1x,a,i3)
1002 format(5x,a12,8x,'v=',5(i2,13x))
1003 format(1x,f12.6,5x,5f15.8)

end subroutine prwf_vibrot
