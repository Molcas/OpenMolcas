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

function CHO_ISAOSH(IAO,ISHL)
!
! Purpose: return symmetry of AO number IAO in shell ISHL.

use Cholesky, only: iBasSh, nSym
#ifdef _DEBUGPRINT_
use Cholesky, only: nBstSh
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: CHO_ISAOSH
integer(kind=iwp), intent(in) :: IAO, ISHL
integer(kind=iwp), external :: CHO_IRANGE
#ifdef _DEBUGPRINT_
character(len=*), parameter :: SECNAM = 'CHO_ISAOSH'

if ((ISHL > NSHELL) .or. (ISHL < 1)) then
  write(LUPRI,'(//,1X,A,A,I10)') SECNAM,': shell index out of bounds: ',ISHL
  write(LUPRI,'(A,I10,A,/)') 'Maximum possible: NSHELL = ',NSHELL,'(from common block)'
  if (NSHELL < 1) then
    call CHO_QUIT('Initialization error detected in '//SECNAM,102)
  else
    call CHO_QUIT('Internal error detected in '//SECNAM,103)
  end if
else if ((IAO > NBSTSH(ISHL)) .or. (IAO < 1)) then
  write(LUPRI,'(//,1X,A,A,I10)') SECNAM,': AO index out of bounds: ',IAO,' shell: ',ISHL
  write(LUPRI,'(A,I10,A,/)') 'Maximum possible: NBSTSH(ISHL) = ',NBSTSH(ISHL),'(from common block)'
  if (NBSTSH(ISHL) < 1) then
    call CHO_QUIT('Initialization error detected in '//SECNAM,102)
  else
    call CHO_QUIT('Internal error detected in '//SECNAM,103)
  end if
end if
#endif

CHO_ISAOSH = CHO_IRANGE(IAO,IBASSH(:,ISHL),NSYM,.false.)

end function CHO_ISAOSH
