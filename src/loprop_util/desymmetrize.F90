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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 SOInt(nSOInt), Scr(nScr), AOInt(nBas1,nBas1),SymInv(nBas1**2)
integer nBas(0:nSym-1)

call FZero(AOInt,nBas1**2)

iOffSO = 1
iOffPi = 1
do iSym=0,nSym-1
  iOffPj = 1
  do jSym=0,iSym
    ijSym = ieor(iSym,jSym)
    if (iand(iSyLbl,2**ijSym) == 0) Go To 20
    if (nBas(iSym)*nBas(jSym) == 0) Go To 30
    if (iSym == jSym) then
#     ifdef _DEBUGPRINT_
      write(6,*) 'Diagonal Block'
      call RecPrt('SOInt',' ',SOInt(iOffSO),nBas(iSym),nBas(iSym))
#     endif

      call DGEMM_('N','T',nBas(iSym),nBas1,nBas(iSym),1.0d0,SOInt(iOffSO),nBas(iSym),SymInv(iOffPi),nBas1,0.0d0,Scr,nBas(iSym))

      call DGEMM_('N','N',nBas1,nBas1,nBas(iSym),One,SymInv(iOffPi),nBas1,Scr,nBas(iSym),One,AOInt,nBas1)

    else
#     ifdef _DEBUGPRINT_
      write(6,*) 'Off-Diagonal Block',iSym,jSym
      call RecPrt('SOInt',' ',SOInt(iOffSO),nBas(iSym),nBas(jSym))
      call RecPrt('Pi',' ',SymInv(iOffPi),nbas1,nBas(iSym))
      call RecPrt('Pj',' ',SymInv(iOffPj),nbas1,nBas(jSym))
#     endif

      call DGEMM_('N','T',nBas(iSym),nBas1,nBas(jSym),1.0d0,SOInt(iOffSO),nBas(iSym),SymInv(iOffPj),nBas1,0.0d0,Scr,nBas(iSym))

      call DGEMM_('N','N',nBas1,nBas1,nBas(iSym),One,SymInv(iOffPi),nBas1,Scr,nBas(iSym),One,AOInt,nBas1)

      call DGEMM_('T','T',nBas1,nBas1,nBas(iSym),One,Scr,nBas(iSym),SymInv(iOffPi),nBas1,One,AOInt,nBas1)

    end if

30  continue
    iOffSO = iOffSO+nBas(iSym)*nBas(jSym)
20  continue
    iOffPj = iOffPj+nBas(jSym)*nBas1
  end do
  iOffPi = iOffPi+nBas(iSym)*nBas1
end do
#ifdef _DEBUGPRINT_
call RecPrt('AOInt',' ',AOInt,nBas1,nBas1)
#endif

return

end subroutine Desymmetrize
