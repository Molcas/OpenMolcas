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

module segtab

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: IBVPT(26) = [0,0,0,0,1,1,2,2,1,1,2,1,1,2,2,1,2,2,3,3,3,3,3,3,3,3], &
                                IC1(26) = [0,1,2,3,0,2,0,1,0,1,1,2,3,0,1,2,2,3,1,3,2,3,0,1,2,3], &
                                IC2(26) = [0,1,2,3,1,3,2,3,0,1,2,2,3,0,1,1,2,3,0,2,0,1,0,1,2,3], &
                                ISVC(26) = [1,1,1,1,1,6,1,5,1,2,4,7,2,1,7,3,2,2,1,5,1,6,1,1,1,1], &
                                ITVPT(26) = [0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,1,1,2,2,3,3,3,3]

public :: IBVPT, IC1, IC2, ISVC, ITVPT

end module segtab
