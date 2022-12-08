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
! Copyright (C) 2016, Sebastian Wouters                                *
!***********************************************************************
! Subroutine to write FCIDUMP file
! Written by Sebastian Wouters, Leuven, Aug 2016

subroutine FCIDUMP_OUTPUT(NACT,NELEC,TWOMS,ISYM,ORBSYM,ECONST,OEI,TEI,LINSIZE,NUM_TEI)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NACT, NELEC, TWOMS, ISYM, ORBSYM(NACT), LINSIZE, NUM_TEI
real(kind=wp), intent(in) :: ECONST, OEI(LINSIZE), TEI(NUM_TEI)
integer(kind=iwp) :: i, ij, ijkl, j, k, l, kl, writeout
real(kind=wp), parameter :: thr = 1.0e-16_wp
integer(kind=iwp), external :: isFreeUnit

writeout = isfreeunit(28)
!open(unit=writeout,file='FCIDUMP_CHEMPS2',action='write',status='replace')
call molcas_open(writeout,'FCIDUMP_CHEMPS2')
write(writeout,'(a11,i3,a7,i3,a5,i2,a1)') ' &FCI NORB=',NACT,',NELEC=',NELEC,',MS2=',TWOMS,','
write(writeout,'(a9)',advance='NO') '  ORBSYM='
do i=1,NACT
  write(writeout,'(i1,a1)',advance='NO') ORBSYM(i),','
end do
write(writeout,*)
write(writeout,'(a7,i1,a1)') '  ISYM=',ISYM,','
write(writeout,'(a2)') ' /'

ij = 0
ijkl = 0
do i=1,NACT
  do j=1,i
    ij = ij+1
    kl = 0
    kl_loop: do k=1,NACT
      do l=1,k
        kl = kl+1
        if (kl > ij) exit kl_loop
        ijkl = ijkl+1
        if (abs(TEI(ijkl)) >= thr) write(writeout,'(1x,es23.16e2,i4,i4,i4,i4)') TEI(ijkl),i,j,k,l
      end do
    end do kl_loop
  end do
end do

ij = 0
do i=1,NACT
  do j=1,i
    ij = ij+1
    if (abs(OEI(ij)) >= thr) write(writeout,'(1x,es23.16e2,i4,i4,i4,i4)') OEI(ij),i,j,0,0
  end do
end do

write(writeout,'(1x,es23.16e2,i4,i4,i4,i4)') ECONST,0,0,0,0

close(writeout)

return

end subroutine FCIDUMP_OUTPUT
