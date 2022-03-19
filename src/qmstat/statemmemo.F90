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

! MO-basis route.
subroutine StateMMEmo(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME,iCent,ipAvRed,Cha,Dip,Qua)

use Index_Functions, only: nTri3_Elem, nTri_Elem
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: nAObas, nMObas, nState, nTyp, iCi, iBigT, iMME(nTri3_Elem(MxMltp)), iCent(nTri_Elem(nAObas)), ipAvRed
real(kind=wp) :: Cha(MxStOT,MxQCen), Dip(MxStOT,3,MxQCen), Qua(MxStOT,6,MxQCen)
integer(kind=iwp) :: i, iB1, iB2, ipAOG, ipAOG_s, ipMOG, ipMOG_s, ipO, iS1, iS2, iTEMP, iTyp, j, kaunta, kaunter, kk, nSizeA, nSizeM
real(kind=wp) :: PerAake

kaunter = 0
nSizeA = nTri_Elem(nAObas)
nSizeM = nTri_Elem(nMObas)
call GetMem('Transition','Allo','Real',ipMOG,nSizeM)
call GetMem('SqMO','Allo','Real',ipMOG_s,nMObas**2)
call GetMem('TEMP','Allo','Real',iTEMP,nAObas*nMObas)
call GetMem('SqAO','Allo','Real',ipAOG_s,nAObas**2)
call GetMem('TransitionA','Allo','Real',ipAOG,nSizeA)
call GetMem('OnTheWay','Allo','Real',ipO,nTyp)

! Loop over state pairs.

do iS1=1,nState
  do iS2=1,iS1
    kaunter = kaunter+1

    ! Collect the proper piece of the TDM in MO-basis.

    call dCopy_(nSizeM,Work(iBigT+nSizeM*(kaunter-1)),1,Work(ipMOG),1)

    ! Additional transformation step from MO to AO.

    call Square(Work(ipMOG),Work(ipMOG_s),1,nMObas,nMObas)
    kk = 0
    do i=1,nMObas
      do j=1,nMObas
        if (i /= j) Work(ipMOG_s+kk) = Half*Work(ipMOG_s+kk)
        kk = kk+1
      end do
    end do
    call Dgemm_('N','N',nAObas,nMObas,nMObas,One,Work(ipAvRed),nAObas,Work(ipMOG_s),nMObas,Zero,Work(iTEMP),nAObas)
    call Dgemm_('N','T',nAObas,nAObas,nMObas,One,Work(iTEMP),nAObas,Work(ipAvRed),nAObas,Zero,Work(ipAOG_s),nAObas)
    kk = 0
    do i=1,nAObas
      do j=1,nAObas
        if (i /= j) Work(ipAOG_s+kk) = Two*Work(ipAOG_s+kk)
        kk = kk+1
      end do
    end do
    call SqToTri_Q(Work(ipAOG_s),Work(ipAOG),nAObas)

    ! Loop over AO-basis pairs.

    kaunta = 0
    do iB1=1,nAObas
      do iB2=1,iB1
        PerAake = Work(ipAOG+kaunta)
        do iTyp=1,nTyp
          Work(ipO+iTyp-1) = Work(iMME(iTyp)+kaunta)*PerAake
        end do
        kaunta = kaunta+1
        Cha(kaunter,iCent(kaunta)) = Cha(kaunter,iCent(kaunta))+Work(ipO)
        Dip(kaunter,1,iCent(kaunta)) = Dip(kaunter,1,iCent(kaunta))+Work(ipO+1)
        Dip(kaunter,2,iCent(kaunta)) = Dip(kaunter,2,iCent(kaunta))+Work(ipO+2)
        Dip(kaunter,3,iCent(kaunta)) = Dip(kaunter,3,iCent(kaunta))+Work(ipO+3)
        Qua(kaunter,1,iCent(kaunta)) = Qua(kaunter,1,iCent(kaunta))+Work(ipO+4)
        Qua(kaunter,2,iCent(kaunta)) = Qua(kaunter,2,iCent(kaunta))+Work(ipO+5)
        ! The reason why 7 and 6 are interchanged is that
        ! QMSTAT uses the ordering xx,xy,yy,xz,yz,zz while
        ! Seward uses the ordering xx,xy,xz,yy,yz,zz.
        Qua(kaunter,3,iCent(kaunta)) = Qua(kaunter,3,iCent(kaunta))+Work(ipO+7)
        Qua(kaunter,4,iCent(kaunta)) = Qua(kaunter,4,iCent(kaunta))+Work(ipO+6)
        Qua(kaunter,5,iCent(kaunta)) = Qua(kaunter,5,iCent(kaunta))+Work(ipO+8)
        Qua(kaunter,6,iCent(kaunta)) = Qua(kaunter,6,iCent(kaunta))+Work(ipO+9)
      end do
    end do
  end do
end do
call GetMem('Transition','Free','Real',ipMOG,nSizeM)
call GetMem('SqMO','Free','Real',ipMOG_s,nMObas**2)
call GetMem('TEMP','Free','Real',iTEMP,nAObas*nMObas)
call GetMem('SqAO','Free','Real',ipAOG_s,nAObas**2)
call GetMem('TransitionA','Free','Real',ipAOG,nSizeA)
call GetMem('OnTheWay','Free','Real',ipO,nTyp)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iCi)

end subroutine StateMMEmo
