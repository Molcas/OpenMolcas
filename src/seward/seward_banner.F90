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

subroutine Seward_Banner()

use Definitions, only: u6

implicit none

! figlet -f larry3d.flf Seward
write(u6,100) " ____                                           __"
write(u6,100) "/\  _`\                                        /\ \"
write(u6,100) "\ \,\L\_\     __   __  __  __     __     _ __  \_\ \"
write(u6,100) " \/_\__ \   /'__`\/\ \/\ \/\ \  /'__`\  /\`'__\/'_` \"
write(u6,100) "   /\ \L\ \/\  __/\ \ \_/ \_/ \/\ \L\.\_\ \ \//\ \L\ \"
write(u6,100) "   \ `\____\ \____\\ \___x___/'\ \__/.\_\\ \_\\ \___,_\"
write(u6,100) "    \/_____/\/____/ \/__//__/   \/__/\/_/ \/_/ \/__,_ /"
write(u6,100) ""
write(u6,100) ""

100 format(21x,a)

end subroutine Seward_Banner
