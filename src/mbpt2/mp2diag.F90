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

use MBPT2_Global, only: DiaA, EOcc, EVir, mAdDel, mAdFro, mAdOcc, mAdVir
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Four
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: iA, iB, iI, iJ, iSym, iSym1, iSym2, nMaxOrb
real(kind=wp) :: E_a, E_i, Ediff, xibja, xijab, xijba
real(kind=wp), allocatable :: Int1(:), IntC(:), Scr1(:)
#include "corbinf.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do

call mma_allocate(Int1,nMaxOrb,label='Int1')
call mma_allocate(IntC,nMaxOrb,label='IntC')
call mma_allocate(Scr1,nMaxOrb,label='Scr1')

do iSym=1,nSym

  do iA=1,nExt(iSym)+nDel(iSym)
    iB = iA
    call Exch(iSym,iSym,iSym,iSym,iA+nOcc(iSym)+nFro(iSym),iB+nOcc(iSym)+nFro(iSym),Int1,Scr1)
    call Coul(iSym,iSym,iSym,iSym,iA+nOcc(iSym)+nFro(iSym),iB+nOcc(iSym)+nFro(iSym),IntC,Scr1)
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' *  A,B = ',iA,iB
    call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
    call RecPrt('IntC:','(8F10.6)',IntC,nOrb(iSym)+nDel(iSym),nOrb(iSym)+nDel(iSym))
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
      xijab = Int1(iJ+(iI-1)*(nOrb(iSym)+nDel(iSym)))
      xijba = xijab
      xibja = IntC(iI+(iJ-1)*(nOrb(iSym)+nDel(iSym)))
      if (iA <= nExt(iSym)) then
        E_a = EVir(mAdVir(iSym)+iA-1)
      else
        E_a = EVir(mAdDel(iSym)+iA-nExt(iSym)-1)
      end if
      if (iI > nFro(iSym)) then
        E_i = EOcc(mAdOcc(iSym)+iI-nFro(iSym)-1)
      else
        E_i = EOcc(mAdFro(iSym)+iI-1)
      end if
      Ediff = E_a-E_i
      !-----------------------------------------------------------------
      !write(u6,*) 'xijab',xijab
      !write(u6,*) 'xijba',xijba
      !write(u6,*) 'xibja',xibja
      !-----------------------------------------------------------------
      DiaA%SB(iSym)%A2(iI,iA) = DiaA%SB(iSym)%A2(iI,iA)+One/(Ediff+Four*xijab-xijba-xibja)

    end do
  end do
end do

call mma_deallocate(Int1)
call mma_deallocate(IntC)
call mma_deallocate(Scr1)
#ifdef _DEBUGPRINT_
do iSym=1,nSym
  write(u6,*) 'Symmetry nr',iSym
  call RecPrt('Diag(ia|ia)','',DiaA%SB(iSym)%A2,nOcc(iSym),nExt(iSym))
end do
#endif

return

end subroutine Mp2Diag
