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

subroutine gugadrt_paras_calculate()

use gugadrt_global, only: ja_sys, jb_sys, jc_sys, jroute_sys, n_electron, norb_all, norb_dz, spin
use Constants, only: Half

implicit none
!data inlptb_new/
!         1  2  3  4  5  6  7  8  9  10  11  12 13
!! a^r=1
!     *    -1, 0, 0, 0, 0, 0,-7, 0,-9,-10,  0,-12, 0,
!! a^l=2
!     *     0, 0, 0, 0, 0,-6, 0,-8, 0,  0,-11,  0, 0,
!! b_r=3
!     *    4, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
!! b_l=4
!     *    5, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
!! b^r=5
!     *    0, 7, 8,10,11, 0, 0, 0, 0,  0,  0,  0, 0,
!! b^l=6
!     *    0, 0, 9, 0,12, 0, 0, 0, 0,  0,  0,  0, 9,
!! c^'=7
!     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3,
!! c^"=8
!     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,
!! d_l^r=9
!     *    6, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0,
!! d^r^r=10
!     *    0,-2, 0,-4, 0, 0, 0, 0, 0,  0,  0,  0, 0,
!! d^r^l=11
!     *    0, 0,-3, 0,-5, 0, 0, 0, 0,  0,  0,  0, 0,
!! c^'" =12
!     *    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0/

!****************************************************************
!   ar      =1 (+a^r)     drr     =2 (+d^rr)   drl     =3 (+d^rl)
!   arbr    =4 (+d^rr)    arbl    =5 (+d^rl)   ard_l^r =6 (+a^l)
!   drrb^r  =7 (+a^r)     drlb^r  =8 (+a^l)    drlb^l  =9 (+a^r)
!   arbrb^r =10 (+a^r)    arblb^r =11 (+a^l)   arblb^l =12 (+a^r)
!   drl     =13 (*)
!****************************************************************
!data inlptb_new
!    */ -1,  0,  4,  5,  0,  0,  1,  1,  6,  0,  0,  1,
!    *   0,  0,  0,  0,  7,  0,  2,  2,  0, -2,  0,  2,
!    *   0,  0,  0,  0,  8,  9,  3,  3,  0,  0, -3,  3,
!    *   0,  0,  0,  0, 10,  0,  4,  4,  0, -4,  0,  4,
!    *   0,  0,  0,  0, 11, 12,  5,  5,  0,  0, -5,  5,
!    *   0, -6,  0,  0,  0,  0,  6,  6,  0,  0,  0,  6,
!    *  -7,  0,  0,  0,  0,  0,  7,  7,  0,  0,  0,  7,
!    *   0, -8,  0,  0,  0,  0,  8,  8,  0,  0,  0,  8,
!    *  -9,  0,  0,  0,  0,  0,  9,  9,  0,  0,  0,  9,
!    * -10,  0,  0,  0,  0,  0, 10, 10,  0,  0,  0, 10,
!    *   0,-11,  0,  0,  0,  0, 11, 11,  0,  0,  0, 11,
!    * -12,  0,  0,  0,  0,  0, 12, 12,  0,  0,  0, 12,
!    *   0,  0,  0,  0,  0,  9,  3, 13,  0,  0,  0,  0/
!data indord_plptype/0,0,0,1,1,3,3,1,1,2,2,2/   !severe_new_error_1

ja_sys = int(n_electron*Half-spin)-norb_dz
jb_sys = int(spin+spin)
jc_sys = norb_all-ja_sys-jb_sys

if (jb_sys == 0) jroute_sys = 1
if (jb_sys == 1) jroute_sys = 2
if (jb_sys > 1) jroute_sys = 3

return

end subroutine gugadrt_paras_calculate
