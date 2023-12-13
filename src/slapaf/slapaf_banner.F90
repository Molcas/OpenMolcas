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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

subroutine Slapaf_banner()

use Definitions, only: u6

implicit none

! figlet -f larry3d.flf Slap af
write(u6,100) " ____     ___                                      ___"
write(u6,100) "/\  _`\  /\_ \                                   /'___\"
write(u6,100) "\ \,\L\_\\//\ \      __     _____          __   /\ \__/"
write(u6,100) " \/_\__ \  \ \ \   /'__`\  /\ '__`\      /'__`\ \ \ ,__\"
write(u6,100) "   /\ \L\ \ \_\ \_/\ \L\.\_\ \ \L\ \    /\ \L\.\_\ \ \_/"
write(u6,100) "   \ `\____\/\____\ \__/.\_\\ \ ,__/    \ \__/.\_\\ \_\"
write(u6,100) "    \/_____/\/____/\/__/\/_/ \ \ \/      \/__/\/_/ \/_/"
write(u6,100) "                              \ \_\"
write(u6,100) "                               \/_/"

100 format(16x,a)

end subroutine Slapaf_banner
