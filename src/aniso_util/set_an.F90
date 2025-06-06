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

subroutine Set_an(ln,a)
!----------------------------------------!
! anm( ln, rank )
! are tabulated from the book:
!
! A. Abragam & B. Bleaney
! EPR of Transition Ions
! Oxford University Press, 1970, Table 20
!----------------------------------------!

use Constants, only: Zero, One, Two, Four, Five, Seven, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ln
real(kind=wp), intent(out) :: a(6)
real(kind=wp) :: an(28,6)

an(:,:) = Zero
a(:) = Zero

!---------  Ce  -------------------------------!
an(1,2) = -Two/35.0_wp               !       -2 / 5*7
an(1,4) = Two/315.0_wp               !        2 / 9*5*7
an(1,6) = Zero
!---------  Pr  -------------------------------!
an(2,2) = -52.0_wp/2475.0_wp         !    -4*13 / 9*25*11
an(2,4) = -Four/5445.0_wp            !       -4 / 9*5*121
an(2,6) = 272.0_wp/4459455.0_wp      !    16*17 / 81*5*7*121*13
!---------  Nd  -------------------------------!
an(3,2) = -Seven/1089.0_wp           !       -7 / 9*121
an(3,4) = -136.0_wp/467181.0_wp      !    -8*17 / 27*11*121*13
an(3,6) = -1615.0_wp/42513471.0_wp   ! -5*17*19 / 27*7*1331*169
!---------  Pm  -------------------------------!
an(4,2) = 14.0_wp/1815.0_wp          !       14 / 3*5*121
an(4,4) = 952.0_wp/2335905.0_wp      !   8*7*17 / 27*5*1331*13
an(4,6) = 2584.0_wp/42513471.0_wp    !  8*17*19 / 27*7*1331*169
!---------  Sm  -------------------------------!
an(5,2) = 13.0_wp/315.0_wp           !       13 / 9*5*7
an(5,4) = 26.0_wp/10395.0_wp         !     2*13 / 27*5*7*11
an(5,6) = Zero
!---------  Eu  -------------------------------!
an(6,2) = Zero
an(6,4) = Zero
an(6,6) = Zero
!---------  Gd  -------------------------------!
an(7,2) = Zero
an(7,4) = Zero
an(7,6) = Zero
!---------  Tb  -------------------------------!
an(8,2) = -One/99.0_wp               !  -1 / 9*11
an(8,4) = Two/16335.0_wp             !   2 / 27*5*121
an(8,6) = -One/891891.0_wp           !  -1 / 81*7*121*13
!---------  Dy  -------------------------------!
an(9,2) = -Two/315.0_wp              !  -2 / 9*5*7
an(9,4) = -Eight/135135.0_wp         !  -8 / 27*5*7*11*13
an(9,6) = Four/3864861.0_wp          !   4 / 27*7*121*169
!---------  Ho  -------------------------------!
an(10,2) = -One/450.0_wp             !  -1 / 2*9*25
an(10,4) = -One/30030.0_wp           !  -1 / 2*3*5*7*11*13
an(10,6) = -Five/3864861.0_wp        !  -5 / 27*7*121*169
!---------  Er  -------------------------------!
an(11,2) = Four/1575.0_wp            !   4 / 9*25*7
an(11,4) = Two/45045.0_wp            !   2 / 9*5*7*11*13
an(11,6) = Eight/3864861.0_wp        !   8 / 27*7*121*169
!---------  Tm  -------------------------------!
an(12,2) = One/99.0_wp               !   1 / 9*11
an(12,4) = Eight/49005.0_wp          !   8 / 81*5*121
an(12,6) = -Five/891891.0_wp         !  -5 / 81*7*121*13
!---------  Yb  -------------------------------!
an(13,2) = Two/63.0_wp               !   2 / 9*7
an(13,4) = -Two/1155.0_wp            !  -2 / 3*5*7*11
an(13,6) = Four/27027.0_wp           !   4 / 27*7*11*13
!---------  Lu  -------------------------------!
an(14,2) = Zero
an(14,4) = Zero
an(14,6) = Zero
!---------  Lu  -------------------------------!

!   Actinides (III)

!---------  Th  -------------------------------!
an(15,2) = -Two/35.0_wp              ! -2 / 5*7
an(15,4) = Two/315.0_wp              !  2 / 9*5*7
an(15,6) = Zero
!---------  Pa  -------------------------------!
an(16,2) = -52.0_wp/2475.0_wp        ! -4*13 / 9*25*11
an(16,4) = -Four/5445.0_wp           !    -4 / 9*5*121
an(16,6) = 272.0_wp/4459455.0_wp     ! 16*17 / 81*5*7*121*13
!---------  U   -------------------------------!
an(17,2) = -Seven/1089.0_wp          !       -7 / 9*121
an(17,4) = -136.0_wp/467181.0_wp     !    -8*17 / 27*11*121*13
an(17,6) = -1615.0_wp/42513471.0_wp  ! -5*17*19 / 27*7*1331*169
!---------  Np  -------------------------------!
an(18,2) = 14.0_wp/1815.0_wp         !       14 / 3*5*121
an(18,4) = 952.0_wp/2335905.0_wp     !   8*7*17 / 27*5*1331*13
an(18,6) = 2584.0_wp/42513471.0_wp   !  8*17*19 / 27*7*1331*169
!---------  Pu  -------------------------------!
an(19,2) = 13.0_wp/315.0_wp          !       13 / 9*5*7
an(19,4) = 26.0_wp/10395.0_wp        !     2*13 / 27*5*7*11
an(19,6) = Zero
!---------  Am  -------------------------------!
an(20,2) = Zero
an(20,4) = Zero
an(20,6) = Zero
!---------  Cm  -------------------------------!
an(21,2) = Zero
an(21,4) = Zero
an(21,6) = Zero
!---------  Bk  -------------------------------!
an(22,2) = -One/99.0_wp              !  -1 / 9*11
an(22,4) = Two/16335.0_wp            !   2 / 27*5*121
an(22,6) = -One/891891.0_wp          !  -1 / 81*7*121*13
!---------  Cf  -------------------------------!
an(23,2) = -Two/315.0_wp             !  -2 / 9*5*7
an(23,4) = -Eight/135135.0_wp        !  -8 / 27*5*7*11*13
an(23,6) = Four/3864861.0_wp         !   4 / 27*7*121*169
!---------  Es  -------------------------------!
an(24,2) = -One/450.0_wp             !  -1 / 2*9*25
an(24,4) = -One/30030.0_wp           !  -1 / 2*3*5*7*11*13
an(24,6) = -Five/3864861.0_wp        !  -5 / 27*7*121*169
!---------  Fm  -------------------------------!
an(25,2) = Four/1575.0_wp            !   4 / 9*25*7
an(25,4) = Two/45045.0_wp            !   2 / 9*5*7*11*13
an(25,6) = Eight/3864861.0_wp        !   8 / 27*7*121*169
!---------  Md  -------------------------------!
an(26,2) = One/99.0_wp               !   1 / 9*11
an(26,4) = Eight/49005.0_wp          !   8 / 81*5*121
an(26,6) = -Five/891891.0_wp         !  -5 / 81*7*121*13
!---------  No  -------------------------------!
an(27,2) = Two/63.0_wp               !   2 / 9*7
an(27,4) = -Two/1155.0_wp            !  -2 / 3*5*7*11
an(27,6) = Four/27027.0_wp           !   4 / 27*7*11*13
!---------  Lr  -------------------------------!
an(28,2) = Zero
an(28,4) = Zero
an(28,6) = Zero
!---------  End -------------------------------!

a(:) = an(ln,:)

return

end subroutine Set_an
