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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This procedure classifies the atomic orbitals of atoms into 6 classes*
! and returns the count. Which shells that are returned are specified  *
! by the the option switch opt:                                        *
!                                                                      *
!  1 DeepCore                                                          *
!  2 Core                                                              *
!  4 SoftCore                                                          *
!  8 DeepValence                                                       *
! 16 Valence                                                           *
! 32 ExtraValence                                                      *
!                                                                      *
! The numbers are added up to get more than one shell reported.        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: May 2003                                                    *
!                                                                      *
!***********************************************************************

subroutine OrbType(Z,List,opt)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Z, opt
integer(kind=iwp), intent(out) :: List(0:3)
integer(kind=iwp) :: DeepCore(0:3), Core(0:3), SoftCore(0:3), DeepValence(0:3), Valence(0:3), ExtraValence(0:3)
integer(kind=iwp), parameter :: s = 0, p = 1, d = 2, f = 3

!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
DeepCore(:) = 0
Core(:) = 0
SoftCore(:) = 0
DeepValence(:) = 0
Valence(:) = 0
ExtraValence(:) = 0
List(:) = 0

!----------------------------------------------------------------------*
! How many shells are there for this atom                              *
!----------------------------------------------------------------------*
select case (Z)
  case (0) ! Dummy
  case (1:2) ! H-He
    Valence(s) = 1
  case (3:4) ! Li-Be
    Core(s) = 1
    Valence(s) = 1
    ExtraValence(p) = 1
  case (5:10) ! B-Ne
    Core(s) = 1
    Valence(s) = 1
    Valence(p) = 1
  case (11:12) ! Na-Mg
    DeepCore(s) = 1
    SoftCore(s) = 1
    SoftCore(p) = 1
    Valence(s) = 1
    ExtraValence(p) = 1
  case (13:18) ! Al-Ar
    DeepCore(s) = 1
    Core(s) = 1
    Core(p) = 1
    Valence(s) = 1
    Valence(p) = 1
  case (19:20) ! K-Ca
    DeepCore(s) = 2
    DeepCore(p) = 1
    SoftCore(s) = 1
    SoftCore(p) = 1
    Valence(s) = 1
    ExtraValence(p) = 1
  case (21:30) ! Sc-Zn
    DeepCore(s) = 2
    DeepCore(p) = 1
    Core(s) = 1
    Core(p) = 1
    Valence(s) = 1
    Valence(d) = 1
    ExtraValence(p) = 1
  case (31:36) ! Ga-Kr
    DeepCore(s) = 2
    DeepCore(p) = 1
    Core(s) = 1
    Core(p) = 1
    Core(d) = 1
    Valence(s) = 1
    Valence(p) = 1
  case (37:38) ! Rb-Sr
    DeepCore(s) = 3
    DeepCore(p) = 2
    DeepCore(d) = 1
    SoftCore(s) = 1
    SoftCore(p) = 1
    Valence(s) = 1
    ExtraValence(p) = 1
  case (39:48) ! Y-Cd
    DeepCore(s) = 3
    DeepCore(p) = 2
    DeepCore(d) = 1
    Core(s) = 1
    Core(p) = 1
    Valence(s) = 1
    Valence(d) = 1
    ExtraValence(p) = 1
  case (49:54) ! In-Xe
    DeepCore(s) = 3
    DeepCore(p) = 2
    DeepCore(d) = 1
    Core(s) = 1
    Core(p) = 1
    Core(d) = 1
    Valence(s) = 1
    Valence(p) = 1
  case (55:56) ! Cs-Ba
    DeepCore(s) = 4
    DeepCore(p) = 3
    DeepCore(d) = 2
    SoftCore(s) = 1
    SoftCore(p) = 1
    Valence(s) = 1
    ExtraValence(p) = 1
  case (57:70) ! La-Yb
    DeepCore(s) = 4
    DeepCore(p) = 3
    DeepCore(d) = 2
    Core(s) = 1
    Core(p) = 1
    Valence(s) = 1
    Valence(f) = 1
    ExtraValence(p) = 1
  case (71:80) ! Lu-Hg
    DeepCore(s) = 4
    DeepCore(p) = 3
    DeepCore(d) = 2
    Core(s) = 1
    Core(p) = 1
    SoftCore(f) = 1
    Valence(s) = 1
    Valence(d) = 1
    ExtraValence(p) = 1
  case (81:86) ! Tl-Rn
    DeepCore(s) = 4
    DeepCore(p) = 3
    DeepCore(d) = 2
    Core(s) = 1
    Core(p) = 1
    Core(f) = 1
    SoftCore(d) = 1
    Valence(s) = 1
    Valence(p) = 1
  case (87:88) ! Fr-Ra
    DeepCore(s) = 5
    DeepCore(p) = 4
    DeepCore(d) = 3
    DeepCore(f) = 1
    SoftCore(s) = 1
    SoftCore(p) = 1
    Valence(s) = 1
    ExtraValence(p) = 1
  case (89:102) ! Ac-No
    DeepCore(s) = 5
    DeepCore(p) = 4
    DeepCore(d) = 3
    DeepCore(f) = 1
    Core(s) = 1
    Core(p) = 1
    Valence(s) = 1
    Valence(f) = 1
    ExtraValence(p) = 1
  case (103:112) ! Lr-Cn
    DeepCore(s) = 5
    DeepCore(p) = 4
    DeepCore(d) = 3
    DeepCore(f) = 1
    Core(s) = 1
    Core(p) = 1
    SoftCore(f) = 1
    Valence(s) = 1
    Valence(d) = 1
    ExtraValence(p) = 1
  case default
    write(u6,*) 'orbtype: element ',Z,' not yet implemented'
    call Abend()
end select
!----------------------------------------------------------------------*
! Fill up the list to be returned                                      *
!----------------------------------------------------------------------*
if (btest(opt,0)) List(:) = List(:)+DeepCore(:)
if (btest(opt,1)) List(:) = List(:)+Core(:)
if (btest(opt,2)) List(:) = List(:)+SoftCore(:)
if (btest(opt,3)) List(:) = List(:)+DeepValence(:)
if (btest(opt,4)) List(:) = List(:)+Valence(:)
if (btest(opt,5)) List(:) = List(:)+ExtraValence(:)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine OrbType
