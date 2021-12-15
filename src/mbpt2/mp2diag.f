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
      SubRoutine Mp2Diag()
*                                                                      *
************************************************************************
*                                                                      *
*          Called from: WfCtl_MP2
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
#include "mp2grad.fh"
#include "trafo.fh"
#include "files_mbpt2.fh"
#include "corbinf.fh"
*#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iDiaA(i,j,k) = ip_DiaA(k) + j-1 +
     &                           (nOcc(k)+nFro(k))*(i-1)
*                                                                      *
************************************************************************
*                                                                      *
      nMaxOrb=0
      Do iSym1 = 1,nSym
         Do iSym2 = 1,nSym
            nMaxOrb = Max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))
     &                               *(nOrb(iSym2)+nDel(iSym2)))
         End Do
      End Do
*
      lint = nMaxOrb
      Call GetMem('Int1','Allo','Real',ipInt1,lInt)
      Call GetMem('IntC','Allo','Real',ipIntC,lInt)
      Call GetMem('Scr1','Allo','Real',ipScr1,lInt)
*
      Do iSym = 1, nSym
*
         Do iA = 1, nExt(iSym) + nDel(iSym)
            iB = iA
            Call Exch(iSym,iSym,iSym,iSym,
     &                iA + nOcc(iSym) + nFro(iSym),
     &                iB + nOcc(iSym) + nFro(iSym),
     &                Work(ipInt1),Work(ipScr1))
            Call Coul(iSym,iSym,iSym,iSym,
     &                iA + nOcc(iSym) + nFro(iSym),
     &                iB + nOcc(iSym) + nFro(iSym),
     &                Work(ipIntC),Work(ipScr1))
#ifdef _DEBUGPRINT_
            Write(6,*)
            Write(6,*) ' *  A,B = ',iA,iB
            Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
     &                  nOrb(iSym)+nDel(iSym),
     &                  nOrb(iSym)+nDel(iSym))
            Call RecPrt('IntC:','(8F10.6)',Work(ipIntC),
     &                  nOrb(iSym)+nDel(iSym),
     &                  nOrb(iSym)+nDel(iSym))
#endif
*
*           A temporary storage for the coulomb integrals
*           in triangular form, they will later be put in
*           IntC as a square matrix.
*                                                                      *
************************************************************************
*                                                                      *
*           Calculate the diagonal
*
            Do iI = 1, nFro(iSym) + nOcc(iSym)
               iJ = iI
               xijab = Work(ipInt1 + iJ-1 +
     &                         (iI-1)*(nOrb(iSym)+nDel(iSym)))
               xijba = xijab
               xibja = Work(ipIntC + iI-1 +
     &                   (iJ-1)*(nOrb(iSym)+nDel(iSym)))
               If(iA .le. nExt(iSym)) Then
                  E_a = Work(mAdVir(iSym) + iA -1)
               Else
                  E_a = Work(mAdDel(iSym) + iA - nExt(iSym) -1)
               End If
               If(iI.gt.nFro(iSym)) Then
                  E_i = Work(mAdOcc(iSym) + iI - nFro(iSym) -1)
               Else
                  E_i = Work(mAdFro(iSym) + iI -1)
               End If
               Ediff = E_a - E_i
*------------------------------------------------------------
*               write(6,*) 'xijab', xijab
*               write(6,*) 'xijba', xijba
*               write(6,*) 'xibja', xibja
*------------------------------------------------------------
               Work(iDiaA(iA,iI,iSym)) =
     &           Work(iDiaA(iA,iI,iSym)) + 1.0d0/( Ediff +
     &            Four*xijab-xijba-xibja)
     &
            EndDo
         EndDo
      EndDo

      Call GetMem('Int1','Free','Real',ipInt1,lInt)
      Call GetMem('IntC','Free','Real',ipIntC,lInt)
      Call GetMem('Scr1','Free','Real',ipScr1,lInt)
#ifdef _DEBUGPRINT_
      Do iSym = 1,nSym
          write (6,*) 'Symmetry nr', iSym
          Call RecPrt('Diag(ia|ia)','',
     &                work(iDiaA(1,1,iSym)),
     &                nOcc(iSym),nExt(iSym))
      EndDo
#endif
      return
      End
