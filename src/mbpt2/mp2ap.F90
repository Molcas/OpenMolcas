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

subroutine MP2Ap(iSymIA,iSymJB,ip_AP,ip_P)
! A subroutine that calculates A*p_k in the PCG-algorithm

implicit real*8(a-h,o-z)
!defining One etc.
#include "real.fh"
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "trafo.fh"
#include "files_mbpt2.fh"
#include "corbinf.fh"

nMaxOrb = 0
do iSym1=1,nSym
  do iSym2=1,nSym
    nMaxOrb = max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))*(nOrb(iSym2)+nDel(iSym2)))
  end do
end do
lint = nMaxOrb
call GetMem('Int1','Allo','Real',ipInt1,lInt)
call GetMem('Int2','Allo','Real',ipInt2,lInt)
call GetMem('IntC','Allo','Real',ipIntC,lInt)
call GetMem('Scr1','Allo','Real',ipScr1,lInt)

nB = nExt(iSymJB)+nDel(iSymJB)
do iA=1,nExt(iSymIA)+nDel(iSymIA)
  if (iSymIA == iSymJB) nB = iA
  do iB=1,nB
    call Exch(iSymIA,iSymIA,iSymJB,iSymJB,iA+nFro(iSymIA)+nOcc(iSymIA),iB+nFro(iSymJB)+nOcc(iSymJB),Work(ipInt1),Work(ipScr1))
    if (iSymIA /= iSymJB) then
      call Exch(iSymJB,iSymIA,iSymIA,iSymJB,iA+nFro(iSymIA)+nOcc(iSymIA),iB+nFro(iSymJB)+nOcc(iSymJB),Work(ipInt2),Work(ipScr1))
    end if
    call Coul(iSymIA,iSymJB,iSymIA,iSymJB,iA+nFro(iSymIA)+nOcc(iSymIA),iB+nFro(iSymJB)+nOcc(iSymJB),Work(ipIntC),Work(ipScr1))
    !write(6,*)
    !write(6,*) ' *  A,B = ',iA,iB
    !call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymIA)+nDel(iSymIA),nOrb(iSymJB)+nDel(iSymJB))
    !if (iSymIA.Ne.iSymJB) Then
    !  call RecPrt('Int2:','(8F10.6)',Work(ipInt2),nOrb(iSymJB)+nDel(iSymJB),nOrb(iSymIA)+nDel(iSymIA))
    !end if
    !call RecPrt('IntC:','(8F10.6)',Work(ipIntC),nOrb(iSymIA)+nDel(iSymIA),nOrb(iSymJB)+nDel(iSymJB))

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
        index1 = ip_Ap+iPoVec(iSymJB)+iJ-1+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1)
        index2 = ip_Ap+iPoVec(iSymIA)+iI-1+(nFro(iSymIA)+nOcc(iSymIA))*(iA-1)
        xiajb = Work(ipInt1+iI-1+(iJ-1)*(nOrb(iSymIA)+nDel(iSymIA)))
        xijab = Work(ipIntC+iI-1+(iJ-1)*(nOrb(iSymIA)+nDel(iSymIA)))
        if (iSymIA == iSymJB) then
          xibja = Work(ipInt1+iJ-1+(iI-1)*(nOrb(iSymJB)+nDel(iSymJB)))

        else
          xibja = Work(ipInt2+iJ-1+(iI-1)*(nOrb(iSymJB)+nDel(iSymJB)))
        end if

        Work(index1) = Work(index1)+Fac*(Four*xiajb-xibja-xijab)*Work(ip_P+iPoVec(iSymIA)+iI-1+(nFro(iSymIA)+nOcc(iSymIA))*(iA-1))
        !---- Debug comments -------------------------------------------
        !write(6,*) 'Symm A B',iSymIA,iSymJB
        !write(6,*) 'IAJB',iI,iA,iJ,iB
        !write(6,*) 'jiba',xijab
        !write(6,*) 'adress AP',iPoVec(iSymJB)+iJ-1+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1)
        !write(6,*) 'Contr coul',Fac*(xijab)*Work(ip_P+iPoVec(iSymIA)+iI-1+(nFro(iSymIA)+nOcc(iSymIA))*(iA-1))
        !---------------------------------------------------------------
        !write(6,*) 'IA JB',iOrbI,iVirA,iOrbJ,iVirB
        !write(6,*) 'Index1:',index1-iAp
        !write(6,*) 'A-elementet',Fac*(Four*xijab-xijba-xibja)
        !write(6,*) 'P-vector',Work(iP+iPoVec(JB)+iOrbJ-1+nOcc(JB)*(iVirB-1))
        !write(6,*) 'Index2:',index2-iAp
        !write(6,*) 'A-elementet',Fac*(Four*xijab-xijba-xibja)
        !write(6,*) 'P-vector',Work(iP+iPoVec(IA)+iOrbI-1+nOcc(IA)*(iVirA-1))

        Work(index2) = Work(index2)+Fac*(Four*xiajb-xibja-xijab)*Work(ip_P+iPoVec(iSymJB)+iJ-1+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1))
        if ((iA == iB) .and. (iI == iJ) .and. (iSymIA == iSymJB)) then
          if (iA <= nExt(iSymIA)) then
            E_a = Work(mAdVir(iSymIA)+iA-1)
          else
            E_a = Work(mAdDel(iSymIA)+iA-nExt(iSymIA)-1)
          end if
          if (iI > nFro(iSymIA)) then
            E_i = Work(mAdOcc(iSymIA)+iI-nFro(iSymIA)-1)
          else
            E_i = Work(mAdFro(iSymIA)+iI-1)
          end if
          Ediff = E_a-E_i
          Work(index1) = Work(index1)+EDiff*Work(ip_P+iPoVec(iSymJB)+iJ-1+(nFro(iSymJB)+nOcc(iSymJB))*(iB-1))

        end if
      end do !iOrbJ
    end do   !iOrbI

  end do     !iOrbB
end do       !iOrbA

call GetMem('Int1','Free','Real',ipInt1,lInt)
call GetMem('Int2','Free','Real',ipInt2,lInt)
call GetMem('IntC','Free','Real',ipIntC,lInt)
call GetMem('Scr1','Free','Real',ipScr1,lInt)

return

end subroutine MP2Ap
