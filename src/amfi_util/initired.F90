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

subroutine initired()

!bs initialize all information for irreducible representations
!bs later on, it might be useful to have a switch for
!bs changing to other orders of IREDs like e.g. in TURBOMOLE
!
! HOW2ADD another symmetry:
!
! 1. add it in readbas.f to be accepted. Add the number of IRs
!
! 2. copy one of the symmetry-blocks in this subroutine and
!    edit the multiplication-table for the group
!
! 3. assign the right IRs to L_X, L_Y and L_Z
!
! that is  all. Good luck!!!

implicit real*8(a-h,o-z)
#include "para.fh"
#include "ired.fh"
character*3 symmetry

symmetry = 'D2H'  ! MOLCAS-Version
if (symmetry == 'D2H') then
  mult(2,1) = 2
  mult(3,1) = 3
  mult(4,1) = 4
  mult(5,1) = 5
  mult(6,1) = 6
  mult(7,1) = 7
  mult(8,1) = 8

  mult(3,2) = 4
  mult(4,2) = 3
  mult(5,2) = 6
  mult(6,2) = 5
  mult(7,2) = 8
  mult(8,2) = 7

  mult(4,3) = 2
  mult(5,3) = 7
  mult(6,3) = 8
  mult(7,3) = 5
  mult(8,3) = 6

  mult(5,4) = 8
  mult(6,4) = 7
  mult(7,4) = 6
  mult(8,4) = 5

  mult(6,5) = 2
  mult(7,5) = 3
  mult(8,5) = 4

  mult(7,6) = 4
  mult(8,6) = 3

  mult(8,7) = 2

  do ired=1,8
    mult(ired,ired) = 1
  end do
  do irun=2,8
    do jrun=1,irun-1
      mult(jrun,irun) = mult(irun,jrun)
    end do
  end do
  !BS write(6,*)
  !BS write(6,*) 'multiplicitation table (Atkins,Child and Phillips)'
  !BS write(6,*)
  !BS do ired=1,8
  !BS   write(6,'(8I5)') (mult(jred,ired),jred=1,8)
  !BS   write(6,*)
  !BS end do

  IRLX = 4
  IRLY = 3
  IRLZ = 2
  !bs assume same order of ireds as Atkins Child and Phillips use..
  !bs would lead to an order with 1 to 1, 2 to 2 ...
  !bs however, this is the molecule/ seward order.
  iredorder(1) = 1
  iredorder(2) = 4
  iredorder(3) = 6
  iredorder(4) = 7
  iredorder(5) = 8
  iredorder(6) = 5
  iredorder(7) = 3
  iredorder(8) = 2
  do ired=1,8
    iredorderinv(iredorder(ired)) = ired
  end do
  ipow2ired(0,0,0) = iredorder(1)
  ipow2ired(1,1,0) = iredorder(2)
  ipow2ired(1,0,1) = iredorder(3)
  ipow2ired(0,1,1) = iredorder(4)
  ipow2ired(1,1,1) = iredorder(5)
  ipow2ired(0,0,1) = iredorder(6)
  ipow2ired(0,1,0) = iredorder(7)
  ipow2ired(1,0,0) = iredorder(8)
  !write(6,*) 'interacting IRs '
  do ired=1,8
    IRwithLX(ired) = iredorder(mult(IRLX,iredorderinv(ired)))
    IRwithLY(ired) = iredorder(mult(IRLY,iredorderinv(ired)))
    IRwithLZ(ired) = iredorder(mult(IRLZ,iredorderinv(ired)))
    !write(6,*) IRwithLX(ired),IRwithLY(ired),IRwithLZ(ired)
  end do
else if (symmetry == 'C2V') then
  !bs 1. A1 2. A2 3. B1 4. B2
  mult(2,1) = 2
  mult(3,1) = 3
  mult(4,1) = 4

  mult(3,2) = 4
  mult(4,2) = 3

  mult(4,3) = 2

  do ired=1,4
    mult(ired,ired) = 1
  end do
  do irun=2,4
    do jrun=1,irun-1
      mult(jrun,irun) = mult(irun,jrun)
    end do
  end do
  write(6,*)
  write(6,*) 'multiplicitation table '
  write(6,*)
  do ired=1,4
    write(6,'(4I5)') (mult(jred,ired),jred=1,4)
    write(6,*)
  end do

  IRLX = 4
  IRLY = 3
  IRLZ = 2
  !bs this is the molecule/seward order.
  iredorder(1) = 1
  iredorder(2) = 4
  iredorder(3) = 2
  iredorder(4) = 3
  do ired=1,4
    iredorderinv(iredorder(ired)) = ired
  end do
  ipow2ired(0,0,0) = iredorder(1)
  ipow2ired(1,1,0) = iredorder(2)
  ipow2ired(1,0,1) = iredorder(3)
  ipow2ired(0,1,1) = iredorder(4)
  ipow2ired(1,1,1) = iredorder(2)
  ipow2ired(0,0,1) = iredorder(1)
  ipow2ired(0,1,0) = iredorder(4)
  ipow2ired(1,0,0) = iredorder(3)
  !write(6,*) 'interacting IRs '
  do ired=1,4
    IRwithLX(ired) = iredorder(mult(IRLX,iredorderinv(ired)))
    IRwithLY(ired) = iredorder(mult(IRLY,iredorderinv(ired)))
    IRwithLZ(ired) = iredorder(mult(IRLZ,iredorderinv(ired)))
    !write(6,*) IRwithLX(ired),IRwithLY(ired),IRwithLZ(ired)
  end do
else if (symmetry == 'D2 ') then
  !bs 1. A1 2. B1 3. B2 4. B3
  mult(2,1) = 2
  mult(3,1) = 3
  mult(4,1) = 4

  mult(3,2) = 4
  mult(4,2) = 3
  mult(4,3) = 2

  do ired=1,4
    mult(ired,ired) = 1
  end do
  do irun=2,4
    do jrun=1,irun-1
      mult(jrun,irun) = mult(irun,jrun)
    end do
  end do
  write(6,*)
  write(6,*) 'multiplicitation table '
  write(6,*)
  do ired=1,4
    write(6,'(4I5)') (mult(jred,ired),jred=1,4)
    write(6,*)
  end do

  IRLX = 4
  IRLY = 3
  IRLZ = 2
  iredorder(1) = 1
  iredorder(2) = 2
  iredorder(3) = 3
  iredorder(4) = 4
  do ired=1,4
    iredorderinv(iredorder(ired)) = ired
  end do
  ipow2ired(0,0,0) = iredorder(1)
  ipow2ired(1,1,0) = iredorder(2)
  ipow2ired(1,0,1) = iredorder(3)
  ipow2ired(0,1,1) = iredorder(4)
  ipow2ired(1,1,1) = iredorder(1)
  ipow2ired(0,0,1) = iredorder(2)
  ipow2ired(0,1,0) = iredorder(3)
  ipow2ired(1,0,0) = iredorder(4)
  !write(6,*) 'interacting IRs '
  do ired=1,4
    IRwithLX(ired) = iredorder(mult(IRLX,iredorderinv(ired)))
    IRwithLY(ired) = iredorder(mult(IRLY,iredorderinv(ired)))
    IRwithLZ(ired) = iredorder(mult(IRLZ,iredorderinv(ired)))
    !write(6,*) IRwithLX(ired),IRwithLY(ired),IRwithLZ(ired)
  end do
else if (symmetry == 'C2H') then
  !bs assume 1.Ag 2.Au 3.Bg 4.Bu
  mult(2,1) = 2
  mult(3,1) = 3
  mult(4,1) = 4

  mult(3,2) = 4
  mult(4,2) = 3

  mult(4,3) = 2

  do ired=1,4
    mult(ired,ired) = 1
  end do
  do irun=2,4
    do jrun=1,irun-1
      mult(jrun,irun) = mult(irun,jrun)
    end do
  end do
  write(6,*)
  write(6,*) 'multiplicitation table '
  write(6,*)
  do ired=1,4
    write(6,'(4I5)') (mult(jred,ired),jred=1,4)
    write(6,*)
  end do

  IRLX = 3
  IRLY = 3
  IRLZ = 1
  iredorder(1) = 1
  iredorder(2) = 2
  iredorder(3) = 3
  iredorder(4) = 4
  do ired=1,4
    iredorderinv(iredorder(ired)) = ired
  end do
  ipow2ired(0,0,0) = iredorder(1)
  ipow2ired(1,1,0) = iredorder(1)
  ipow2ired(1,0,1) = iredorder(3)
  ipow2ired(0,1,1) = iredorder(3)
  ipow2ired(1,1,1) = iredorder(2)
  ipow2ired(0,0,1) = iredorder(2)
  ipow2ired(0,1,0) = iredorder(4)
  ipow2ired(1,0,0) = iredorder(4)
  !write(6,*) 'interacting IRs '
  do ired=1,4
    IRwithLX(ired) = iredorder(mult(IRLX,iredorderinv(ired)))
    IRwithLY(ired) = iredorder(mult(IRLY,iredorderinv(ired)))
    IRwithLZ(ired) = iredorder(mult(IRLZ,iredorderinv(ired)))
    !write(6,*) IRwithLX(ired),IRwithLY(ired),IRwithLZ(ired)
  end do
else if (symmetry == 'CS ') then
  write(6,*) 'CS in initired '
  !bs assume 1.A' 2.A'
  mult(2,1) = 2

  do ired=1,2
    mult(ired,ired) = 1
  end do
  do irun=2,2
    do jrun=1,irun-1
      mult(jrun,irun) = mult(irun,jrun)
    end do
  end do
  write(6,*)
  write(6,*) 'multiplicitation table '
  write(6,*)
  do ired=1,2
    write(6,'(2I5)') (mult(jred,ired),jred=1,2)
    write(6,*)
  end do

  IRLX = 2
  IRLY = 2
  IRLZ = 1
  iredorder(1) = 1
  iredorder(2) = 2
  do ired=1,2
    iredorderinv(iredorder(ired)) = ired
  end do
  ipow2ired(0,0,0) = iredorder(1)
  ipow2ired(1,1,0) = iredorder(1)
  ipow2ired(1,0,1) = iredorder(2)
  ipow2ired(0,1,1) = iredorder(2)
  ipow2ired(1,1,1) = iredorder(2)
  ipow2ired(0,0,1) = iredorder(2)
  ipow2ired(0,1,0) = iredorder(1)
  ipow2ired(1,0,0) = iredorder(1)
  !write(6,*) 'interacting IRs '
  do ired=1,2
    IRwithLX(ired) = iredorder(mult(IRLX,iredorderinv(ired)))
    IRwithLY(ired) = iredorder(mult(IRLY,iredorderinv(ired)))
    IRwithLZ(ired) = iredorder(mult(IRLZ,iredorderinv(ired)))
    !write(6,*) IRwithLX(ired),IRwithLY(ired),IRwithLZ(ired)
  end do
end if

return

end subroutine initired
