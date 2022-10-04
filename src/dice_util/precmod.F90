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
! Copyright (C) 2022, George Benthien                                  *
!***********************************************************************
! String utilites from George Benthien
! Approved for using by George Benthien
! http://gbenthien.net/strings/str-index.html

module precision

! Real kinds

integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

! Integer kinds

integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

!Complex kinds

integer, parameter :: kc4 = kr4                            ! single precision complex
integer, parameter :: kc8 = kr8                            ! double precision complex

end module precision
