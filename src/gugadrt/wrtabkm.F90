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

subroutine wrtabkm(iabkm,labkm,nabcbit,iabcbit,jatmp,jbtmp,jmtmp,kktmp)
! this subroutine unpack ja jb jm kt from compress arrays

implicit real*8(a-h,o-z)
dimension iabkm(1:labkm)

call packnod(iabkm,1,jatmp,nabcbit,iabcbit,labkm)
call packnod(iabkm,2,jbtmp,nabcbit,iabcbit,labkm)
call packnod(iabkm,3,jmtmp,nabcbit,iabcbit,labkm)
call packnod(iabkm,4,kktmp,nabcbit,iabcbit,labkm)

return

end subroutine wrtabkm
