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
! Copyright (C) 2004, Giovanni Ghigo                                   *
!***********************************************************************

subroutine put_tra_comm(IBD2M,NSYMX,NORBX,NOSHX,LUINTMX)

implicit real*8(a-h,o-z)
integer IBD2M(3,36*36), NSYMX, NORBX(8), NOSHX(8), LUINTMX
#include "intgrl.fh"

do j=1,36*36
  do i=1,3
    IAD2M(i,j) = IBD2M(i,j)
  end do
end do
NSYMZ = NSYMX
do i=1,8
  NORBZ(i) = NORBX(i)
  NOSHZ(i) = NOSHX(i)
end do
LUINTMZ = LUINTMX

return

end subroutine put_tra_comm
