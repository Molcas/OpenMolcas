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

subroutine Mp2Diag()
!***********************************************************************
!                                                                      *
! Called from: WfCtl_MP2                                               *
!                                                                      *
!***********************************************************************

use MBPT2_Global, only: ip_DiaA, ipInt1, ipScr1, mAdDel, mAdFro, mAdOcc, mAdVir
use Constants, only: One, Four
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: iA, iB, iI, iJ, ipIntC, iSym, iSym1, iSym2, lint, nMaxOrb
real(kind=wp) :: E_a, E_i, Ediff, xibja, xijab, xijba
#include "WrkSpc.fh"
#include "corbinf.fh"
! statement function
integer(kind=iwp) :: i, j, k, iDiaA
iDiaA(i,j,k) = ip_DiaA(k)+j-1+(nOcc(k)+nFro(k))*(i-1)

!                                                                      *
!***********************************************************************
!                                                                      *
nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do

lint = nMaxOrb
call GetMem('Int1','Allo','Real',ipInt1,lInt)
call GetMem('IntC','Allo','Real',ipIntC,lInt)
call GetMem('Scr1','Allo','Real',ipScr1,lInt)

do iSym=1,nSym

  do iA=1,nExt(iSym)+nDel(iSym)
    iB = iA
    call Exch(iSym,iSym,iSym,iSym,iA+nOcc(iSym)+nFro(iSym),iB+nOcc(iSym)+nFro(iSym),Work(ipInt1),Work(ipScr1))
    call Coul(iSym,iSym,iSym,iSym,iA+nOcc(iSym)+nFro(iSym),iB+nOcc(iSym)+nFro(iSym),Work(ipIntC),Work(ipScr1))
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' *  A,B = ',iA,iB
    call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
    call RecPrt('IntC:','(8F10.6)',Work(ipIntC),nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
#   endif

    ! A temporary storage for the coulomb integrals
    ! in triangular form, they will later be put in
    ! IntC as a square matrix.
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Calculate the diagonal

    do iI=1,nFro(iSym)+nOcc(iSym)
      iJ = iI
      xijab = Work(ipInt1+iJ-1+(iI-1)*(nOrb(iSym)+nDel(iSym)))
      xijba = xijab
      xibja = Work(ipIntC+iI-1+(iJ-1)*(nOrb(iSym)+nDel(iSym)))
      if (iA <= nExt(iSym)) then
        E_a = Work(mAdVir(iSym)+iA-1)
      else
        E_a = Work(mAdDel(iSym)+iA-nExt(iSym)-1)
      end if
      if (iI > nFro(iSym)) then
        E_i = Work(mAdOcc(iSym)+iI-nFro(iSym)-1)
      else
        E_i = Work(mAdFro(iSym)+iI-1)
      end if
      Ediff = E_a-E_i
      !-----------------------------------------------------------------
      !write(u6,*) 'xijab',xijab
      !write(u6,*) 'xijba',xijba
      !write(u6,*) 'xibja',xibja
      !-----------------------------------------------------------------
      Work(iDiaA(iA,iI,iSym)) = Work(iDiaA(iA,iI,iSym))+One/(Ediff+Four*xijab-xijba-xibja)

    end do
  end do
end do

call GetMem('Int1','Free','Real',ipInt1,lInt)
call GetMem('IntC','Free','Real',ipIntC,lInt)
call GetMem('Scr1','Free','Real',ipScr1,lInt)
#ifdef _DEBUGPRINT_
do iSym=1,nSym
  write(u6,*) 'Symmetry nr',iSym
  call RecPrt('Diag(ia|ia)','',work(iDiaA(1,1,iSym)),nOcc(iSym),nExt(iSym))
end do
#endif

return

end subroutine Mp2Diag
