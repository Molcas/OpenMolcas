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
integer(kind=iwp), intent(out) :: List(4)
integer(kind=iwp) :: DeepCore(4), Core(4), SoftCore(4), DeepValence(4), Valence(4), ExtraValence(4), i

!----------------------------------------------------------------------*
! Is this a legal element?                                             *
!----------------------------------------------------------------------*
!NIKO if ((Z < 1) .or. (Z > 112)) then
if ((Z < 0) .or. (Z > 112)) then
  write(u6,*) 'orbtype: do only know elements 1-112'
  call Abend()
end if
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

! Dummy
if (Z == 0) then
  Valence(1) = 0 ! this is a redundant operation.
! H-He
else if (Z <= 2) then
  Valence(1) = 1
! Li-Be
else if (Z <= 4) then
  Core(1) = 1
  Valence(1) = 1
  ExtraValence(2) = 1
! B-Ne
else if (Z <= 10) then
  Core(1) = 1
  Valence(1) = 1
  Valence(2) = 1
! Na-Mg
else if (Z <= 12) then
  DeepCore(1) = 1
  SoftCore(1) = 1
  SoftCore(2) = 1
  Valence(1) = 1
  ExtraValence(2) = 1
! Al-Ar
else if (Z <= 18) then
  DeepCore(1) = 1
  Core(1) = 1
  Core(2) = 1
  Valence(1) = 1
  Valence(2) = 1
! K-Ca
else if (Z <= 20) then
  DeepCore(1) = 2
  DeepCore(2) = 1
  SoftCore(1) = 1
  SoftCore(2) = 1
  Valence(1) = 1
  ExtraValence(2) = 1
! Sc-Zn
else if (Z <= 30) then
  DeepCore(1) = 2
  DeepCore(2) = 1
  Core(1) = 1
  Core(2) = 1
  Valence(1) = 1
  Valence(3) = 1
  ExtraValence(2) = 1
! Ga-Kr
else if (Z <= 36) then
  DeepCore(1) = 2
  DeepCore(2) = 1
  Core(1) = 1
  Core(2) = 1
  Core(3) = 1
  Valence(1) = 1
  Valence(2) = 1
! Rb-Sr
else if (Z <= 38) then
  DeepCore(1) = 3
  DeepCore(2) = 2
  DeepCore(3) = 1
  SoftCore(1) = 1
  SoftCore(2) = 1
  Valence(1) = 1
  ExtraValence(2) = 1
! Y-Cd
else if (Z <= 48) then
  DeepCore(1) = 3
  DeepCore(2) = 2
  DeepCore(3) = 1
  Core(1) = 1
  Core(2) = 1
  Valence(1) = 1
  Valence(3) = 1
  ExtraValence(2) = 1
! In-Xe
else if (Z <= 54) then
  DeepCore(1) = 3
  DeepCore(2) = 2
  DeepCore(3) = 1
  Core(1) = 1
  Core(2) = 1
  Core(3) = 1
  Valence(1) = 1
  Valence(2) = 1
! Cs-Ba
else if (Z <= 56) then
  DeepCore(1) = 4
  DeepCore(2) = 3
  DeepCore(3) = 2
  SoftCore(1) = 1
  SoftCore(2) = 1
  Valence(1) = 1
  ExtraValence(1) = 1
! La-Yb
else if (Z <= 70) then
  DeepCore(1) = 4
  DeepCore(2) = 3
  DeepCore(3) = 2
  Core(1) = 1
  Core(2) = 1
  Valence(1) = 1
  Valence(4) = 1
  ExtraValence(2) = 1
! Lu-Hg
else if (Z <= 80) then
  DeepCore(1) = 4
  DeepCore(2) = 3
  DeepCore(3) = 2
  Core(1) = 1
  Core(2) = 1
  SoftCore(4) = 1
  Valence(1) = 1
  Valence(3) = 1
  ExtraValence(2) = 1
! Tl-Rn
else if (Z <= 86) then
  DeepCore(1) = 4
  DeepCore(2) = 3
  DeepCore(3) = 2
  Core(1) = 1
  Core(2) = 1
  Core(4) = 1
  SoftCore(3) = 1
  Valence(1) = 1
  Valence(2) = 1
! Fr-Ra
else if (Z <= 88) then
  DeepCore(1) = 5
  DeepCore(2) = 4
  DeepCore(3) = 3
  DeepCore(4) = 1
  SoftCore(1) = 1
  SoftCore(2) = 1
  Valence(1) = 1
  ExtraValence(2) = 1
! Ac-No
else if (Z <= 102) then
  DeepCore(1) = 5
  DeepCore(2) = 4
  DeepCore(3) = 3
  DeepCore(4) = 1
  Core(1) = 1
  Core(2) = 1
  Valence(1) = 1
  Valence(4) = 1
  ExtraValence(2) = 1
! Lr-Cn
else if (Z <= 112) then
  DeepCore(1) = 5
  DeepCore(2) = 4
  DeepCore(3) = 3
  DeepCore(4) = 1
  Core(1) = 1
  Core(2) = 1
  SoftCore(4) = 1
  Valence(1) = 1
  Valence(3) = 1
  ExtraValence(2) = 1
else
  write(u6,*) 'orbtype: element',Z,' not yet implemented'
  call Abend()
end if
!----------------------------------------------------------------------*
! Fill up the list to be returned                                      *
!----------------------------------------------------------------------*
if (btest(opt,0)) then
  do i=1,4
    List(i) = List(i)+DeepCore(i)
  end do
end if
if (btest(opt,1)) then
  do i=1,4
    List(i) = List(i)+Core(i)
  end do
end if
if (btest(opt,2)) then
  do i=1,4
    List(i) = List(i)+SoftCore(i)
  end do
end if
if (btest(opt,3)) then
  do i=1,4
    List(i) = List(i)+DeepValence(i)
  end do
end if
if (btest(opt,4)) then
  do i=1,4
    List(i) = List(i)+Valence(i)
  end do
end if
if (btest(opt,5)) then
  do i=1,4
    List(i) = List(i)+ExtraValence(i)
  end do
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine OrbType
