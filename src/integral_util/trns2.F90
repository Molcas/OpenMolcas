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
! Copyright (C) 1990, Roland Lindh                                     *
!***********************************************************************

subroutine Trns2(Win,Wout,nvec,nc)
!***********************************************************************
!                                                                      *
! Object: utility routine to transform a AO batch in case of redun-    *
!         dancy of type aA=bB or cC=dD.                                *
!                                                                      *
! Called from: Trns1                                                   *
!                                                                      *
! Calling    : DCopy  (ESSL)                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             June '90                                                 *
!***********************************************************************

implicit none
integer nVec, nc
real*8 Win(nvec,nc,nc), Wout(nvec,nc,nc)
integer ic, id

do ic=1,nc
  do id=1,nc
    call dcopy_(nvec,Win(1,ic,id),1,Wout(1,id,ic),1)
  end do
end do

return

end subroutine Trns2
