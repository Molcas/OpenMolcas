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

subroutine MP2Ap(iSymIA,iSymJB,AP,P)
! A subroutine that calculates A*p_k in the PCG-algorithm

use MBPT2_Global, only: EOcc, EVir, iPoVec, mAdDel, mAdFro, mAdOcc, mAdVir
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSymIA, iSymJB
real(kind=wp), intent(inout) :: AP(*)
real(kind=wp), intent(in) :: P(*)
integer(kind=iwp) :: iA, iB, iI, iJ, index1, index2, iSym1, iSym2, nB, nJ, nMaxOrb
real(kind=wp) :: E_a, E_i, Ediff, Fac, xiajb, xibja, xijab
real(kind=wp), allocatable :: Int1(:), Int2(:), IntC(:), Scr1(:)
#include "corbinf.fh"

nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do
call mma_allocate(Int1,nMaxOrb,label='Int1')
call mma_allocate(Int2,nMaxOrb,label='Int2')
call mma_allocate(IntC,nMaxOrb,label='IntC')
call mma_allocate(Scr1,nMaxOrb,label='Scr1')

nB = nExt(iSymJB)+nDel(iSymJB)
do iA=1,nExt(iSymIA)+nDel(iSymIA)
  if (iSymIA == iSymJB) nB = iA
  do iB=1,nB
    call Exch(iSymIA,iSymIA,iSymJB,iSymJB,iA+nFro(iSymIA)+nOcc(iSymIA),iB+nFro(iSymJB)+nOcc(iSymJB),Int1,Scr1)
    if (iSymIA /= iSymJB) then
      call Exch(iSymJB,iSymIA,iSymIA,iSymJB,iA+nFro(iSymIA)+nOcc(iSymIA),iB+nFro(iSymJB)+nOcc(iSymJB),Int2,Scr1)
    end if
    call Coul(iSymIA,iSymJB,iSymIA,iSymJB,iA+nFro(iSymIA)+nOcc(iSymIA),iB+nFro(iSymJB)+nOcc(iSymJB),IntC,Scr1)
    !write(u6,*)
    !write(u6,*) ' *  A,B = ',iA,iB
    !call RecPrt('Int1:','(8F10.6)',Int1,nOrb(iSymIA)+nDel(iSymIA),nOrb(iSymJB)+nDel(iSymJB))
    !if (iSymIA /= iSymJB) then
    !  call RecPrt('Int2:','(8F10.6)',Int2,nOrb(iSymJB)+nDel(iSymJB),nOrb(iSymIA)+nDel(iSymIA))
    !end if
    !call RecPrt('IntC:','(8F10.6)',IntC,nOrb(iSymIA)+nDel(iSymIA),nOrb(iSymJB)+nDel(iSymJB))

    ! This loop will traverse an triangular AIxBJ-matrix
    do iI=1,nFro(iSymIA)+nOcc(iSymIA)
      nJ = nFro(iSymJB)+nOcc(iSymJB)
      if ((iA == iB) .and. (iSymIA == iSymJB)) nJ = iI
      do iJ=1,nJ
        Fac = One
        ! We will count the diagonal twice so the factor half will get it right
        if ((iA == iB) .and. (iI == iJ) .and. (iSymIA == iSymJB)) then
          Fac = Half
        end if
        ! For multiplication with the BJ element of P
        index1 = iPoVec(iSymJB)+iJ+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1)
        index2 = iPoVec(iSymIA)+iI+(nFro(iSymIA)+nOcc(iSymIA))*(iA-1)
        xiajb = Int1(iI+(iJ-1)*(nOrb(iSymIA)+nDel(iSymIA)))
        xijab = IntC(iI+(iJ-1)*(nOrb(iSymIA)+nDel(iSymIA)))
        if (iSymIA == iSymJB) then
          xibja = Int1(iJ+(iI-1)*(nOrb(iSymJB)+nDel(iSymJB)))

        else
          xibja = Int2(iJ+(iI-1)*(nOrb(iSymJB)+nDel(iSymJB)))
        end if

        AP(index1) = AP(index1)+Fac*(Four*xiajb-xibja-xijab)*P(iPoVec(iSymIA)+iI+(nFro(iSymIA)+nOcc(iSymIA))*(iA-1))
        !---- Debug comments -------------------------------------------
        !write(u6,*) 'Symm A B',iSymIA,iSymJB
        !write(u6,*) 'IAJB',iI,iA,iJ,iB
        !write(u6,*) 'jiba',xijab
        !write(u6,*) 'adress AP',iPoVec(iSymJB)+iJ-1+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1)
        !write(u6,*) 'Contr coul',Fac*(xijab)*P(iPoVec(iSymIA)+iI+(nFro(iSymIA)+nOcc(iSymIA))*(iA-1))
        !---------------------------------------------------------------
        !write(u6,*) 'IA JB',iOrbI,iVirA,iOrbJ,iVirB
        !write(u6,*) 'Index1:',index1
        !write(u6,*) 'A-elementet',Fac*(Four*xijab-xijba-xibja)
        !write(u6,*) 'P-vector',P(iPoVec(JB)+iOrbJ+nOcc(JB)*(iVirB-1))
        !write(u6,*) 'Index2:',index2
        !write(u6,*) 'A-elementet',Fac*(Four*xijab-xijba-xibja)
        !write(u6,*) 'P-vector',P(iPoVec(IA)+iOrbI+nOcc(IA)*(iVirA-1))

        AP(index2) = AP(index2)+Fac*(Four*xiajb-xibja-xijab)*P(iPoVec(iSymJB)+iJ+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1))
        if ((iA == iB) .and. (iI == iJ) .and. (iSymIA == iSymJB)) then
          if (iA <= nExt(iSymIA)) then
            E_a = EVir(mAdVir(iSymIA)+iA-1)
          else
            E_a = EVir(mAdDel(iSymIA)+iA-nExt(iSymIA)-1)
          end if
          if (iI > nFro(iSymIA)) then
            E_i = EOcc(mAdOcc(iSymIA)+iI-nFro(iSymIA)-1)
          else
            E_i = EOcc(mAdFro(iSymIA)+iI-1)
          end if
          Ediff = E_a-E_i
          AP(index1) = AP(index1)+EDiff*P(iPoVec(iSymJB)+iJ+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1))

        end if
      end do !iOrbJ
    end do   !iOrbI

  end do     !iOrbB
end do       !iOrbA

call mma_deallocate(Int1)
call mma_deallocate(Int2)
call mma_deallocate(IntC)
call mma_deallocate(Scr1)

return

end subroutine MP2Ap
