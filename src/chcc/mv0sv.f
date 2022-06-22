************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine mv0sv (dd,length,mat,f)
C
c      A = A .f
c
#include "chcc1.fh"
       integer dd
       integer length
       real*8  mat(1:dd)
       real*8 f
C
C      help variables
C
       integer init
C
C
       if (mhkey.eq.1) then
c      ESSL
c
       do  5 init=1,length
       mat(init) = mat(init)*f
  5    continue
c
       else
c      Fortran matrix handling
C
       do 10 init=1,length
       mat(init) = mat(init)*f
 10    continue
c
       end if
c
c
        return
        end
