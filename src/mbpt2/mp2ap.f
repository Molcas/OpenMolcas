************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*  A subroutine that calculates A*p_k in the PCG-algorithm
*
      Subroutine MP2Ap(iSymIA,iSymJB,ip_AP,ip_P)
*
      Implicit real*8 (a-h,o-z)
*defining One etc.
#include "real.fh"
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "trafo.fh"
#include "files_mbpt2.fh"
#include "corbinf.fh"

*
      nMaxOrb=0
      Do iSym1 = 1,nSym
         Do iSym2 = 1,nSym
            nMaxOrb = Max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))
     &                               *(nOrb(iSym2)+nDel(iSym2)))
         End Do
      End Do
      lint = nMaxOrb
      Call GetMem('Int1','Allo','Real',ipInt1,lInt)
      Call GetMem('Int2','Allo','Real',ipInt2,lInt)
      Call GetMem('IntC','Allo','Real',ipIntC,lInt)
      Call GetMem('Scr1','Allo','Real',ipScr1,lInt)

      nB = nExt(iSymJB) + nDel(iSymJB)
      Do iA = 1, nExt(iSymIA) + nDel(iSymIA)
         if(iSymIA.eq.iSymJB) nB = iA
         Do iB = 1, nB
            Call Exch(iSymIA,iSymIA,iSymJB,iSymJB,
     &                iA + nFro(iSymIA) + nOcc(iSymIA),
     &                iB + nFro(iSymJB) + nOcc(iSymJB),
     &                Work(ipInt1), Work(ipScr1))
            If(iSymIA.Ne.iSymJB) Then
               Call Exch(iSymJB, iSymIA, iSymIA, iSymJB,
     &                   iA + nFro(iSymIA) + nOcc(iSymIA),
     &                   iB + nFro(iSymJB) + nOcc(iSymJB),
     &                   Work(ipInt2), Work(ipScr1))
            End If
            Call Coul(iSymIA,iSymJB,iSymIA,iSymJB,
     &                iA + nFro(iSymIA) + nOcc(iSymIA),
     &                iB + nFro(iSymJB) + nOcc(iSymJB),
     &                Work(ipIntC), Work(ipScr1))
*            Write(6,*)
*            Write(6,*) ' *  A,B = ',iA,iB
*            Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
*     &                  nOrb(iSymIA)+nDel(iSymIA),
*     &                  nOrb(iSymJB)+nDel(iSymJB))
*            If (iSymIA.Ne.iSymJB) Then
*               Call RecPrt('Int2:','(8F10.6)',Work(ipInt2),
*     &                     nOrb(iSymJB)+nDel(iSymJB),
*     &                     nOrb(iSymIA)+nDel(iSymIA))
*            End If
*            Call RecPrt('IntC:','(8F10.6)',Work(ipIntC),
*     &                  nOrb(iSymIA)+nDel(iSymIA),
*     &                  nOrb(iSymJB)+nDel(iSymJB))
**********************************************************************
*
*
*           This loop will traverse an triangular AIxBJ-matrix
            Do iI= 1,nFro(iSymIA) + nOcc(iSymIA)
               nJ=nFro(iSymJB) + nOcc(iSymJB)
               If((iA.eq.iB).and.(iSymIA.eq.iSymJB)) nJ=iI
               Do iJ=1,nJ
                  Fac = One
*                 We will count the diagonal twice so the factor
*                 half will get it right
                  If((iA.eq.iB).and.(iI.eq.iJ)
     &                .and.(iSymIA.eq.iSymJB)) Then
                     Fac = Half
                  Endif
*                  For multiplication with the BJ element of P
                  index1 = ip_Ap + iPoVec(iSymJB) +iJ-1 +
     &                      (nFro(iSymJB) + nOcc(iSymJB))*(iB-1)
                  index2 = ip_Ap + iPoVec(iSymIA) + iI-1 +
     &                      (nFro(iSymIA) +nOcc(iSymIA))*(iA-1)
                  xiajb = Work(ipInt1  + iI-1 + (iJ-1)*
     &                         (nOrb(iSymIA)+nDel(iSymIA)))
                  xijab = Work(ipIntC  + iI-1 +
     &                        (iJ-1)*(nOrb(iSymIA)+nDel(iSymIA)))
                  If (iSymIA.eq.iSymJB) then
                  xibja = Work(ipInt1 + iJ-1 +
     &                        (iI-1)*(nOrb(iSymJB)+nDel(iSymJB)))


                  Else
                     xibja = Work(ipInt2 + iJ-1 +
     &                     (iI-1)*(nOrb(iSymJB)+nDel(iSymJB)))
                  EndIf
*
                  Work(index1)=Work(index1) +
     &                          Fac*(Four*xiajb-xibja-xijab)*
     &                          work(ip_P + iPoVec(iSymIA) +
     &                           iI-1 +
     &                          (nFro(iSymIA) + nOcc(iSymIA))*(iA-1))
*------------ Debug comments -------------------------------------------
*                        Write(6,*) 'Symm A B', iSymIA, iSymJB
*                        Write(6,*) 'IAJB', iI,iA,iJ,iB
*                        Write(6,*) 'jiba', xijab
*                        Write(6,*) 'adress AP',
*     &                             iPoVec(iSymJB) +iJ-1 +
*     &                       (nFro(iSymJB) + nOcc(iSymJB))*(iB-1)
*                        Write(6,*) 'Contr coul',
*     &                             Fac*(xijab)*
*     &                          work(ip_P + iPoVec(iSymIA) +
*     &                           iI-1 +
*     &                          (nFro(iSymIA) + nOcc(iSymIA))*(iA-1))
*------------------------------------------------------------------------
*                  write(6,*) 'IA JB',iOrbI,iVirA,iOrbJ,iVirB
*                  write(6,*) 'Index1:',index1-iAp
*                  write(6,*) 'A-elementet', Fac*(Four*xijab-xijba-xibja)
*                  write(6,*) 'P-vector',work(iP+ iPoVec(JB) + iOrbJ-1+
*     &                          nOcc(JB)*(iVirB-1))
*                  write(6,*) 'Index2:',index2-iAp
*                  write(6,*) 'A-elementet', Fac*(Four*xijab-xijba-xibja)
*                  write(6,*) 'P-vector',work(iP+ iPoVec(IA) + iOrbI-1+
*     &                          nOcc(IA)*(iVirA-1))

                  Work(index2)= Work(index2) +
     &                           Fac*(Four*xiajb-xibja-xijab)*
     &                           work(ip_P + iPoVec(iSymJB) +
     &                            iJ-1 +
     &                           (nFro(iSymJB) + nOcc(iSymJB))*(iB-1))
                  If((iA.eq.iB).and.(iI.eq.iJ).and.
     &                (iSymIA.eq.iSymJB)) then
                     If(iA .le. nExt(iSymIA)) Then
                        E_a = Work(mAdVir(iSymIA) + iA -1)
                     Else
                        E_a = Work(mAdDel(iSymIA) + iA - nExt(iSymIA)-1)
                     End If
                     If(iI.gt.nFro(iSymIA)) Then
                        E_i = Work(mAdOcc(iSymIA) + iI - nFro(iSymIA)-1)
                     Else
                        E_i = Work(mAdFro(iSymIA) + iI -1)
                     End If
                     Ediff = E_a - E_i
                     Work(index1) = Work(index1) + EDiff *
     &                          work(ip_P + iPoVec(iSymJB) +
     &                           iJ-1 +
     &                          (nFro(iSymJB) + nOcc(iSymJB))*(iB-1))

                  Endif
               EndDo            !iOrbJ
            EndDo               !iOrbI
*

         EndDo                  !iOrbB
      EndDo                     !iOrbA
*
      Call GetMem('Int1','Free','Real',ipInt1,lInt)
      Call GetMem('Int2','Free','Real',ipInt2,lInt)
      Call GetMem('IntC','Free','Real',ipIntC,lInt)
      Call GetMem('Scr1','Free','Real',ipScr1,lInt)
      Return
      End
