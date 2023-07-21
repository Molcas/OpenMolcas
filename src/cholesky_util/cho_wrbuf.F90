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

subroutine CHO_WRBUF(LENGTH,BUF,IBUF,LENBUF,IUNIT)
!
! Purpose: write buffer to disk.

implicit real*8(a-h,o-z)
real*8 BUF(LENBUF)
integer IBUF(4,LENBUF)

write(IUNIT) LENGTH,BUF,IBUF

end subroutine CHO_WRBUF
