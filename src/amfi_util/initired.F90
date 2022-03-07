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

subroutine initired(symmetry)
!bs initialize all information for irreducible representations
!bs later on, it might be useful to have a switch for
!bs changing to other orders of IREDs like e.g. in TURBOMOLE
!
! HOW2ADD another symmetry:
!
! 1. add it in readbas to be accepted. Add the number of IRs
!
! 2. copy one of the symmetry-blocks in this subroutine and
!    edit the nired for the group
!
! 3. assign the right IRs to L_X, L_Y and L_Z
!
! that is  all. Good luck!!!

use AMFI_global, only: ipow2ired
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

implicit none
character(len=3), intent(in) :: symmetry
integer(kind=iwp) :: ired, iredorder(8), jred, nired
character(len=5) :: frmt

nired = 0

select case (symmetry)

  case ('D2H') ! MOLCAS-Version
    nired = 8
    !BS write(u6,*)
    !BS write(u6,*) 'multiplication table (Atkins,Child and Phillips)'
    !BS write(u6,*)
    !BS write(frmt,'("(",I1,"I5)")') nired
    !BS do ired=1,nired
    !BS   write(u6,frmt) (Mul(jred,ired),jred=1,nired)
    !BS   write(u6,*)
    !BS end do

    !IRLX = 4
    !IRLY = 3
    !IRLZ = 2
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
    !bs irreducible representation of the cartesian functions in D2H
    !bs order taken from Tables for Group theory:
    !bs Atkins, Child and Phillips   Oxford University Press 1970
    !bs 1. AG: only even powers (0,0,0)
    !bs 2. B1G: (1,1,0)    L_z
    !bs 3. B2G: (1,0,1)    L_y
    !bs 4. B3G: (0,1,1)    L_x
    !bs 5. AU:  (1,1,1)
    !bs 6. B1U: (0,0,1)
    !bs 7. B2U: (0,1,0)
    !bs 8. B3U: (1,0,0)
    ipow2ired(0,0,0) = iredorder(1)
    ipow2ired(1,1,0) = iredorder(2)
    ipow2ired(1,0,1) = iredorder(3)
    ipow2ired(0,1,1) = iredorder(4)
    ipow2ired(1,1,1) = iredorder(5)
    ipow2ired(0,0,1) = iredorder(6)
    ipow2ired(0,1,0) = iredorder(7)
    ipow2ired(1,0,0) = iredorder(8)

  case ('C2V')
    nired = 4
    !bs 1. A1 2. A2 3. B1 4. B2
    write(u6,*)
    write(u6,*) 'multiplication table '
    write(u6,*)
    write(frmt,'("(",I1,"I5)")') nired
    do ired=1,nired
      write(u6,frmt) (Mul(jred,ired),jred=1,nired)
      write(u6,*)
    end do

    !IRLX = 4
    !IRLY = 3
    !IRLZ = 2
    !bs this is the molecule/seward order.
    iredorder(1) = 1
    iredorder(2) = 4
    iredorder(3) = 2
    iredorder(4) = 3
    ipow2ired(0,0,0) = iredorder(1)
    ipow2ired(1,1,0) = iredorder(2)
    ipow2ired(1,0,1) = iredorder(3)
    ipow2ired(0,1,1) = iredorder(4)
    ipow2ired(1,1,1) = iredorder(2)
    ipow2ired(0,0,1) = iredorder(1)
    ipow2ired(0,1,0) = iredorder(4)
    ipow2ired(1,0,0) = iredorder(3)

  case ('D2 ')
    nired = 4
    !bs 1. A1 2. B1 3. B2 4. B3
    write(u6,*)
    write(u6,*) 'multiplication table '
    write(u6,*)
    write(frmt,'("(",I1,"I5)")') nired
    do ired=1,nired
      write(u6,frmt) (Mul(jred,ired),jred=1,nired)
      write(u6,*)
    end do

    !IRLX = 4
    !IRLY = 3
    !IRLZ = 2
    iredorder(1) = 1
    iredorder(2) = 2
    iredorder(3) = 3
    iredorder(4) = 4
    ipow2ired(0,0,0) = iredorder(1)
    ipow2ired(1,1,0) = iredorder(2)
    ipow2ired(1,0,1) = iredorder(3)
    ipow2ired(0,1,1) = iredorder(4)
    ipow2ired(1,1,1) = iredorder(1)
    ipow2ired(0,0,1) = iredorder(2)
    ipow2ired(0,1,0) = iredorder(3)
    ipow2ired(1,0,0) = iredorder(4)

  case ('C2H')
    nired = 4
    !bs assume 1.Ag 2.Au 3.Bg 4.Bu
    write(u6,*)
    write(u6,*) 'multiplication table '
    write(u6,*)
    write(frmt,'("(",I1,"I5)")') nired
    do ired=1,nired
      write(u6,frmt) (Mul(jred,ired),jred=1,nired)
      write(u6,*)
    end do

    !IRLX = 3
    !IRLY = 3
    !IRLZ = 1
    iredorder(1) = 1
    iredorder(2) = 2
    iredorder(3) = 3
    iredorder(4) = 4
    ipow2ired(0,0,0) = iredorder(1)
    ipow2ired(1,1,0) = iredorder(1)
    ipow2ired(1,0,1) = iredorder(3)
    ipow2ired(0,1,1) = iredorder(3)
    ipow2ired(1,1,1) = iredorder(2)
    ipow2ired(0,0,1) = iredorder(2)
    ipow2ired(0,1,0) = iredorder(4)
    ipow2ired(1,0,0) = iredorder(4)

  case ('CS ')
    nired = 2
    !bs assume 1.A' 2.A'
    write(u6,*)
    write(u6,*) 'multiplication table '
    write(u6,*)
    write(frmt,'("(",I1,"I5)")') nired
    do ired=1,nired
      write(u6,frmt) (Mul(jred,ired),jred=1,nired)
      write(u6,*)
    end do

    !IRLX = 2
    !IRLY = 2
    !IRLZ = 1
    iredorder(1) = 1
    iredorder(2) = 2
    ipow2ired(0,0,0) = iredorder(1)
    ipow2ired(1,1,0) = iredorder(1)
    ipow2ired(1,0,1) = iredorder(2)
    ipow2ired(0,1,1) = iredorder(2)
    ipow2ired(1,1,1) = iredorder(2)
    ipow2ired(0,0,1) = iredorder(2)
    ipow2ired(0,1,0) = iredorder(1)
    ipow2ired(1,0,0) = iredorder(1)
end select

!write(u6,*) 'interacting IRs '
!do ired=1,nired
!  iredorderinv(iredorder(ired)) = ired
!end do
!do ired=1,nired
!  IRwithLX(ired) = iredorder(Mul(IRLX,iredorderinv(ired)))
!  IRwithLY(ired) = iredorder(Mul(IRLY,iredorderinv(ired)))
!  IRwithLZ(ired) = iredorder(Mul(IRLZ,iredorderinv(ired)))
!  !write(u6,*) IRwithLX(ired),IRwithLY(ired),IRwithLZ(ired)
!end do

return

end subroutine initired
