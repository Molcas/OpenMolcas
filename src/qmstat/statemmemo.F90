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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6), iCent(MxBas**2)
dimension Cha(MxStOT,MxQCen), Dip(MxStOT,3,MxQCen)
dimension Qua(MxStOT,6,MxQCen)

kaunter = 0
nSizeA = nAObas*(nAObas+1)/2
nSizeM = nMObas*(nMObas+1)/2
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

    call dCopy_(nSizeM,Work(iBigT+nSizeM*(kaunter-1)),iONE,Work(ipMOG),iONE)

    ! Additional transformation step from MO to AO.

    call Square(Work(ipMOG),Work(ipMOG_s),iONE,nMObas,nMObas)
    kk = 0
    do i=1,nMObas
      do j=1,nMObas
        if (i /= j) Work(ipMOG_s+kk) = 0.5d0*Work(ipMOG_s+kk)
        kk = kk+1
      end do
    end do
    call Dgemm_('N','N',nAObas,nMObas,nMObas,ONE,Work(ipAvRed),nAObas,Work(ipMOG_s),nMObas,ZERO,Work(iTEMP),nAObas)
    call Dgemm_('N','T',nAObas,nAObas,nMObas,ONE,Work(iTEMP),nAObas,Work(ipAvRed),nAObas,ZERO,Work(ipAOG_s),nAObas)
    kk = 0
    do i=1,nAObas
      do j=1,nAObas
        if (i /= j) Work(ipAOG_s+kk) = 2.0d0*Work(ipAOG_s+kk)
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
