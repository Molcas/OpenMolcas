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

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IBD2M(3,36*36), NSYMX, NORBX(8), NOSHX(8), LUINTMX
#include "intgrl.fh"

IAD2M(:,:) = IBD2M(:,:)
NSYMZ = NSYMX
NORBZ(:) = NORBX(:)
NOSHZ(:) = NOSHX(:)
LUINTMZ = LUINTMX

return

end subroutine put_tra_comm
