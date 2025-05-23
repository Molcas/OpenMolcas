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

subroutine MCLR_banner()

use Definitions, only: u6

implicit none

! figlet -f larry3d.flf MCLR
write(u6,100) "        ____     __       ____"
write(u6,100) " /'\_/`\/\  _`\  /\ \     /\  _`\"
write(u6,100) "/\      \ \ \/\_\\ \ \    \ \ \L\ \"
write(u6,100) "\ \ \__\ \ \ \/_/_\ \ \  __\ \ ,  /"
write(u6,100) " \ \ \_/\ \ \ \L\ \\ \ \L\ \\ \ \\ \"
write(u6,100) "  \ \_\\ \_\ \____/ \ \____/ \ \_\ \_\"
write(u6,100) "   \/_/ \/_/\/___/   \/___/   \/_/\/ /"
write(u6,100) ""
write(u6,100) ""

100 format(16x,a)

end subroutine MCLR_banner
