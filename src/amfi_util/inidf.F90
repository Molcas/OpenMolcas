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

subroutine inidf()
!bs initializes the df on module with double factorials

use AMFI_global, only: df, dffrac
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: irun, jbm

df(0) = One
df(1) = One
do irun=2,ubound(df,1)
  df(irun) = real(irun,kind=wp)*df(irun-2)
end do
do jbm=0,ubound(df,1)
  dffrac(:,jbm) = df(:)/df(jbm)
end do

return

end subroutine inidf
