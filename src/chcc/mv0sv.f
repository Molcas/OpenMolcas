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
       subroutine mv0sv (dd,length,mat,f)
!
!      A = A .f
!
#include "chcc1.fh"
       integer dd
       integer length
       real*8  mat(1:dd)
       real*8 f
!
!      help variables
!
       integer init
!
!
       if (mhkey.eq.1) then
!      ESSL
!
       do  5 init=1,length
       mat(init) = mat(init)*f
  5    continue
!
       else
!      Fortran matrix handling
!
       do 10 init=1,length
       mat(init) = mat(init)*f
 10    continue
!
       end if
!
!
        return
        end
