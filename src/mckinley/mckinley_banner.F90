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

subroutine McKinley_banner()

use Definitions, only: u6

implicit none

! figlet -f larry3d.flf McKinley
write(u6,100) "                 __  __               ___"
write(u6,100) " /'\_/`\        /\ \/\ \  __         /\_ \"
write(u6,100) "/\      \    ___\ \ \/'/'/\_\    ___ \//\ \      __   __  __"
write(u6,100) "\ \ \__\ \  /'___\ \ , < \/\ \ /' _ `\ \ \ \   /'__`\/\ \/\ \"
write(u6,100) " \ \ \_/\ \/\ \__/\ \ \\`\\ \ \/\ \/\ \ \_\ \_/\  __/\ \ \_\ \"
write(u6,100) "  \ \_\\ \_\ \____\\ \_\ \_\ \_\ \_\ \_\/\____\ \____\\/`____ \"
write(u6,100) "   \/_/ \/_/\/____/ \/_/\/_/\/_/\/_/\/_/\/____/\/____/ `/___/> \"
write(u6,100) "                                                          /\___/"
write(u6,100) "                                                          \/__/"

100 format(25x,a)

end subroutine McKinley_banner
