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

subroutine Gateway_banner()

use Definitions, only: u6

implicit none

! figlet -f larry3d.flf Gateway
write(u6,100) " ____              __"
write(u6,100) "/\  _`\           /\ \__"
write(u6,100) "\ \ \L\_\     __  \ \ ,_\    __   __  __  __     __     __  __"
write(u6,100) " \ \ \L_L   /'__`\ \ \ \/  /'__`\/\ \/\ \/\ \  /'__`\  /\ \/\ \"
write(u6,100) "  \ \ \/, \/\ \L\.\_\ \ \_/\  __/\ \ \_/ \_/ \/\ \L\.\_\ \ \_\ \"
write(u6,100) "   \ \____/\ \__/.\_\\ \__\ \____\\ \___x___/'\ \__/.\_\\/`____ \"
write(u6,100) "    \/___/  \/__/\/_/ \/__/\/____/ \/__//__/   \/__/\/_/ `/___/> \"
write(u6,100) "                                                            /\___/"
write(u6,100) "                                                            \/__/"

100 format(16x,a)

end subroutine Gateway_banner
