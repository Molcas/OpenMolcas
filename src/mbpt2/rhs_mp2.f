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
*                                                                      *
************************************************************************
*                                                                      *
*     The RHS for the MP2-gradients
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine RHS_MP2 ()
*
      Implicit real*8 (a-h,o-z)
#include "files_mbpt2.fh"
#include "trafo.fh"
#include "corbinf.fh"
*defining One etc.
#include "real.fh"
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "chomp2_cfg.fh"
*                                                                      *
************************************************************************
*                                                                      *
      IAD13=0
      LIADOUT=3*36*36
*
*     Read transformed integrals from a logical unit to a
*     buffer
*
      Call iDAFILE(LuIntM,2,IADOUT,LIADOUT,IAD13)
*
*     Sort startadresses for Virtual and occupied orbital
*     energies in nice vectors
*
      nVirTot=0
      nDelTot=0
      Do i = 1, nSym
         nVirTot = nVirTot + nExt(i)
         nDelTot = nDelTot + nDel(i)
      End Do
*
*     The frozen energies are stashed at the end!
*
      mAdOcc(1) = ipEOcc
      do iSym = 2, nSym
         mAdOcc(iSym) = mAdOcc(iSym-1) + nOcc(iSym-1)
      enddo
*
*     The deleted "energies" are stached at the end!
*
      mAdVir(1) = ipEVir
      do iSym=2, nSym
         mAdVir(iSym) = mAdVir(iSym-1) + nExt(iSym-1)
      enddo
*
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
*----- Print occupied and virtual energies
      Do iSym = 1, nSym
         If (nOcc(iSym).gt.0) Then
         Write(6,*)
         Write(6,*) 'Occupied energies, Sym', iSym
         Write(6,*)
         Do iOcc = 0, nOcc(iSym)-1
            Write(6,*) Work(mAdOcc(iSym)+iOcc)
         End Do
         End If
         If (nExt(iSym).gt.0) Then
         Write(6,*)
         Write(6,*) 'Virtual energies, Sym', iSym
         Write(6,*)
         Do iExt = 0, nExt(iSym)-1
            Write(6,*) Work(mAdVir(iSym)+iExt)
         End Do
         End If
         If (nFro(iSym).gt.0) Then
         Write(6,*)
         Write(6,*) 'Frozen energies, Sym', iSym
         Write(6,*)
         Do iFro = 0, nFro(iSym)-1
            Write(6,*) Work(mAdFro(iSym)+iFro)
         End Do
         End If
         If (nDel(iSym).gt.0) Then
         Write(6,*)
         Write(6,*) 'Deleted energies, Sym', iSym
         Write(6,*)
         Do iDel = 0, nDel(iSym)-1
            Write(6,*) Work(mAdDel(iSym)+iDel)
            Write(6,*)
         End Do
         End If
      End Do
#endif
*
      EMP2 = Zero
      VECL2= One
*
*     Find the largest block needed to store exchange integrals
*     either using IJ or AB as fix indeces.
*
      nMaxOrb=0
      Do iSym1 = 1,nSym
         Do iSym2 = 1,nSym
            nMaxOrb = Max(nMaxOrb,(nOrb(iSym1)+nDel(iSym1))
     &                               *(nOrb(iSym2)+nDel(iSym2)))
         End Do
      End Do
*
      lInt = nMaxOrb
*
*     Allocate space for the block of exchange integrals
*     type 1 and type 2 (ia|jb) and (ib|ja)
*
      Call GetMem('Int1','Allo','Real',ipInt1,lInt)
      Call GetMem('Int2','Allo','Real',ipInt2,lInt)
*
*     Allocate a scratchvector for Coul and Exch
*     subroutines.
*
      Call GetMem('Scr1','Allo','Real',ipScr1,lInt)
*
*     Refer to the three-index storing of the symmetryblocks
*     by this routine
*
#ifdef _DEBUGPRINT_
      Write (6,*) 'Before RHS_Mp2_help1'
      Do iSym = 1, nSym
         nI = nOcc(iSym) + nFro(iSym)
         nA = nExt(iSym) + nDel(iSym)
         nB = nOrb(iSym) + nDel(iSym)
         Call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),
     &               nI,nA)
         Call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
      End Do
#endif
      Do iSymI = 1, nSym
         Do iSymJ = 1, iSymI
            Do iSymA = 1,nSym
               iSymB = IEOR(iSymA-1,(IEOR(iSymI-1,iSymJ-1))) + 1
               If (iSymB.le.iSymA) then
*                 Check if there is any doubly excited configurations in
*                 the current symmetry.
                  if((nOrb(iSymI)+nDel(iSymI))
     &              *(nOrb(iSymJ)+nDel(iSymJ))
     &              *(nOrb(iSymA)+nDel(iSymA))
     &              *(nOrb(iSymB)+nDel(iSymB)).ne.0) then
*
                     Call RHS_Mp2_help1(iSymA, iSymB, iSymI, iSymJ)
*
                  endif         !The constructed symmetry was empty
               Endif            !No <occ occ | vir vir> in this symmetry
            enddo               !ASym
         enddo                  !JSym
      enddo                     !ISym
#ifdef _DEBUGPRINT_
      Write (6,*) 'After RHS_Mp2_help1'
      Do iSym = 1, nSym
         nI = nOcc(iSym) + nFro(iSym)
         nA = nExt(iSym) + nDel(iSym)
         nB = nOrb(iSym) + nDel(iSym)
         Call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),
     &               nI,nA)
         Call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Do iSym =1, nSym
         Do i = 1, nOcc(iSym)+nExt(iSym)
            Do j = 1, nFro(iSym)
               Work(ip_Density(iSym) + i+nFro(iSym)-1 +
     &              (j-1)*(nOrb(iSym) + nDel(iSym))) =
     &         Work(ip_Density(iSym) + j-1 +
     &              (i+nFro(iSym)-1)*(nOrb(iSym) + nDel(iSym)))
            End Do
            Do j= 1, nDel(iSym)
              Work(ip_Density(iSym) + j+nOrb(iSym)-1 +
     &              (i+nFro(iSym)-1)*(nOrb(iSym) + nDel(iSym))) =
     &           Work(ip_Density(iSym) + i+nFro(iSym)-1 +
     &              (j+nOrb(iSym)-1)*(nOrb(iSym) + nDel(iSym)))
            End Do
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*) 'After RHS_Mp2_help1 xx'
      Do iSym = 1, nSym
         nI = nOcc(iSym) + nFro(iSym)
         nA = nExt(iSym) + nDel(iSym)
         nB = nOrb(iSym) + nDel(iSym)
         Call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),
     &               nI,nA)
         Call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
      End Do
#endif
*
*     Since some terms in the Lagrangian is dependent of Pab and Pij we
*     have to loop through all the symmetries again since they are not done
*     until the loop is over.
*
      Do iSymI = 1, nSym
         Do iSymJ = 1, iSymI
            Do iSymA = 1,nSym
               iSymB = IEOR(iSymA-1,(IEOR(iSymI-1,iSymJ-1))) + 1
               if (iSymB.le.iSymA) then
*                 Check if there is any doubly excited configurations in
*                 the current symmetry.
                  if((nOrb(iSymI)+nDel(iSymI))
     &              *(nOrb(iSymJ)+nDel(iSymJ))
     &              *(nOrb(iSymA)+nDel(iSymA))
     &              *(nOrb(iSymB)+nDel(iSymB)).ne.0) then
                     Call RHS_Mp2_help2(iSymA,iSymB,iSymI,iSymJ)
                  endif         !The constructed symmetry was empty
               endif            !No <occ occ | vir vir> in this symmetry
            enddo               !ASym
         enddo                  !JSym
      enddo                     !ISym
#ifdef _DEBUGPRINT_
      Write (6,*) 'After RHS_Mp2_help1 yy'
      Do iSym = 1, nSym
         nI = nOcc(iSym) + nFro(iSym)
         nA = nExt(iSym) + nDel(iSym)
         nB = nOrb(iSym) + nDel(iSym)
         Call RecPrt('MP2Lagr',' ',Work(ip_Mp2Lagr(iSym)),
     &               nI,nA)
         Call RecPrt('Dens',' ',Work(ip_Density(iSym)),nB,nB)
      End Do
#endif
*
*     Need two extra integral fields due to symmetrization
*
      If(.not. NoGamma) Then
         Call GetMem('Int1_2','Allo','Real',ipInt1_2,lInt)
         Call GetMem('Int2_2','Allo','Real',ipInt2_2,lInt)
*
*        Construct and backtransform nonsep 2pdm.
*
         Call Gamma_new()
*
         Call GetMem('Int1_2','Free','Real',ipInt1_2,lInt)
         Call GetMem('Int2_2','Free','Real',ipInt2_2,lInt)
      End If
*
      Call GetMem('Int1','Free','Real',ipInt1,lInt)
      Call GetMem('Int2','Free','Real',ipInt2,lInt)
      Call GetMem('Scr1','Free','Real',ipScr1,lInt)
*
#ifdef _DEBUGPRINT_
      write (6,*) 'EMP2 is ', EMP2
      write (6,*) ' '
      Do iSym = 1, nSym
         Write(6,*) 'Density matrix for Symm:', iSym
         Call RecPrt('MP2Density','',Work(ip_Density(iSym)),
     &               nOrb(iSym) + nDel(iSym), nOrb(iSym) + nDel(iSym))
      End Do
      Do iSym = 1, nSym
         Write(6,*) 'WDensity matrix for Symm:', iSym
         Call RecPrt('MP2WDensity','',Work(ip_WDensity(iSym)),
     &               nOrb(iSym) + nDel(iSym), nOrb(iSym) + nDel(iSym))
      End Do
      Do iSym = 1, nSym
         Write(6,*) 'Lagrangian matrix for symm', iSym
         Call RecPrt('Lagr2','',Work(ip_Mp2Lagr(iSym)),
     &               nFro(iSym) + nOcc(iSym),
     &               nExt(iSym) + nDel(iSym))
      End Do
#endif
*
      VECL2=Sqrt(One/VECL2)
*
      Return
      end
