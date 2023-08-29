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

function CHO_ISAO(IAO)
!
! Purpose: return symmetry of AO number IAO (in global list).

use Cholesky, only: IBAS, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_ISAO
integer(kind=iwp), intent(in) :: IAO
integer(kind=iwp), external :: CHO_IRANGE
#ifdef _DEBUGPRINT_
character(len=8), parameter :: SECNAM = 'CHO_ISAO'

if ((IAO > NBAST) .or. (IAO < 1)) then
  write(LUPRI,'(//,1X,A,A,I10)') SECNAM,': AO index out of bounds: ',IAO
  write(LUPRI,'(A,I10,A,/)') 'Maximum possible: NBAST = ',NBAST,'(from common block)'
  if (NBAST < 1) then
    call CHO_QUIT('Initialization error detected in '//SECNAM,102)
  else
    call CHO_QUIT('Internal error detected in '//SECNAM,103)
  end if
end if
#endif

CHO_ISAO = CHO_IRANGE(IAO,IBAS,NSYM,.false.)

end function CHO_ISAO
