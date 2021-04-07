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

subroutine ORTHO_MOTRA(nSym,nBas,nDel,Ovlp,CMO)
! Objective: Orthonormalize input vectors
!           (Gram-Schmidt orthogonaliztion)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nDel(nSym)
real(kind=wp), intent(in) :: Ovlp(*)
real(kind=wp), intent(inout) :: CMO(*)
#include "trafo_motra.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: IJ, IM, ISYM, LW1, LW2, LW3, NORBI

! Allocate work space

call GETMEM('SCR1','ALLO','REAL',LW1,N2MAX)
call GETMEM('SCR2','ALLO','REAL',LW2,N2MAX)
call GETMEM('SCR3','ALLO','REAL',LW3,N2MAX)

! Loop over symmetries and orthogonalize

IJ = 1
IM = 1
do ISYM=1,NSYM
  NORBI = NBAS(ISYM)-NDEL(ISYM)
  if (NORBI > 0) then
    call SQUARE(OVLP(IJ),WORK(LW3),1,NBAS(ISYM),NBAS(ISYM))
    !call MXMA(WORK(LW3),1,NBAS(ISYM),CMO(IM),1,NBAS(ISYM),WORK(LW2),1,NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),NORBI)
    call DGEMM_('N','N',NBAS(ISYM),NORBI,NBAS(ISYM),One,WORK(LW3),NBAS(ISYM),CMO(IM),NBAS(ISYM),Zero,WORK(LW2),NBAS(ISYM))
    !call MXMA(CMO(IM),NBAS(ISYM),1,WORK(LW2),1,NBAS(ISYM),WORK(LW1),1,NORBI,NORBI,NBAS(ISYM),NORBI)
    call DGEMM_('T','N',NORBI,NORBI,NBAS(ISYM),One,CMO(IM),NBAS(ISYM),WORK(LW2),NBAS(ISYM),Zero,WORK(LW1),NORBI)
    call ORTHOX_MOTRA(WORK(LW1),CMO(IM),NORBI,NBAS(ISYM))
  end if
  IJ = IJ+NBAS(ISYM)*(NBAS(ISYM)+1)/2
  IM = IM+NBAS(ISYM)*NBAS(ISYM)
end do

! Deallocate work space and exit

call GETMEM('SCR3','FREE','REAL',LW3,N2MAX)
call GETMEM('SCR2','FREE','REAL',LW2,N2MAX)
call GETMEM('SCR1','FREE','REAL',LW1,N2MAX)

return

end subroutine ORTHO_MOTRA
