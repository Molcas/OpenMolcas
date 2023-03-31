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
       subroutine mv0zero (dd,length,mat)
!
!      mat = 0
!
#include "ccsd1.fh"
       integer dd
       integer length
       real*8  mat(1:dd)
!
!      help variables
!
       integer init
       real*8  zero
       data    zero/0.0d0/
!
!
       if (mhkey.eq.1) then
!      ESSL

       call dcopy_(length,[zero],0,mat,1)
!
       else
!      Fortran matrix handling
!
       do 10 init=1,length
       mat(init) = zero
 10    continue
!
       end if
!
       return
       end
