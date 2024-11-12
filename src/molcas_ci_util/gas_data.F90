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

module gas_data
  use definitions, only: iwp
  implicit none
  Private
# include "Molcas.fh"
  Logical(kind=iwp), Public:: iDoGas=.FALSE.
  INTEGER(kind=iwp), Public:: NGAS=0
  INTEGER(kind=iwp), Public:: NGSSH(mxGAS,mxSym)=0
  INTEGER(kind=iwp), Public:: IGSOCCX(mxGAS,2)=0

end module gas_data
