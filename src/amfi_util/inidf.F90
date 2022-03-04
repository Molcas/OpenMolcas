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
!bs initializes the df on common block  with double facultatives

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "para.fh"
#include "dofuc.fh"
integer(kind=iwp) :: ibm, irun, jbm

df(0) = One
df(1) = One
do irun=2,ndfmx
  df(irun) = real(irun,kind=wp)*df(irun-2)
end do
do jbm=0,ndfmx-1
  do ibm=jbm,ndfmx
    dffrac(ibm,jbm) = df(ibm)/df(jbm)
  end do
end do
do jbm=1,ndfmx
  do ibm=0,jbm-1
    dffrac(ibm,jbm) = One/dffrac(jbm,ibm)
  end do
end do

return

end subroutine inidf
