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

subroutine Desymmetrize(SOInt,nSOInt,Scr,nScr,AOInt,nBas,nBas1,SymInv,nSym,iSyLbl)

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nSOInt, nScr, nSym, nBas(0:nSym-1), nBas1, iSyLbl
real(kind=wp), intent(in) :: SOInt(nSOInt), SymInv(nBas1**2)
real(kind=wp), intent(out) :: Scr(nScr), AOInt(nBas1,nBas1)
integer(kind=iwp) :: ijSym, iOffPi, iOffPj, iOffSO, iSym, jSym

AOInt(:,:) = Zero

iOffSO = 1
iOffPi = 1
do iSym=0,nSym-1
  iOffPj = 1
  do jSym=0,iSym
    ijSym = Mul(iSym+1,jSym+1)-1
    if (btest(iSyLbl,ijSym)) then
      if (nBas(iSym)*nBas(jSym) == 0) cycle
      if (iSym == jSym) then
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Diagonal Block'
        call RecPrt('SOInt',' ',SOInt(iOffSO),nBas(iSym),nBas(iSym))
#       endif

        call DGEMM_('N','T',nBas(iSym),nBas1,nBas(iSym),One,SOInt(iOffSO),nBas(iSym),SymInv(iOffPi),nBas1,Zero,Scr,nBas(iSym))

        call DGEMM_('N','N',nBas1,nBas1,nBas(iSym),One,SymInv(iOffPi),nBas1,Scr,nBas(iSym),One,AOInt,nBas1)

      else
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Off-Diagonal Block',iSym,jSym
        call RecPrt('SOInt',' ',SOInt(iOffSO),nBas(iSym),nBas(jSym))
        call RecPrt('Pi',' ',SymInv(iOffPi),nbas1,nBas(iSym))
        call RecPrt('Pj',' ',SymInv(iOffPj),nbas1,nBas(jSym))
#       endif

        call DGEMM_('N','T',nBas(iSym),nBas1,nBas(jSym),One,SOInt(iOffSO),nBas(iSym),SymInv(iOffPj),nBas1,Zero,Scr,nBas(iSym))

        call DGEMM_('N','N',nBas1,nBas1,nBas(iSym),One,SymInv(iOffPi),nBas1,Scr,nBas(iSym),One,AOInt,nBas1)

        call DGEMM_('T','T',nBas1,nBas1,nBas(iSym),One,Scr,nBas(iSym),SymInv(iOffPi),nBas1,One,AOInt,nBas1)

      end if

      iOffSO = iOffSO+nBas(iSym)*nBas(jSym)
    end if
    iOffPj = iOffPj+nBas(jSym)*nBas1
  end do
  iOffPi = iOffPi+nBas(iSym)*nBas1
end do
#ifdef _DEBUGPRINT_
call RecPrt('AOInt',' ',AOInt,nBas1,nBas1)
#endif

return

end subroutine Desymmetrize
