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
      Subroutine Build_Mp2Dens(TriDens,nTriDens,ip_Density,CMO,mSym,
     &                         nOrbAll,nOccAll,Diagonalize)

      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "corbinf.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
*
      Integer   nTriDens
      Real*8    TriDens(nTriDens)
      Integer   ip_Density(8)
      Real*8    CMO(*)
      Integer   mSym
      Integer   nOrbAll(8), nOccAll(8)
      Logical   Diagonalize

      Integer   ipSymRec(8)
      Integer   ipSymTri(8)
      Integer   ipSymLin(8)
      Character(LEN=30) Note

      Real*8, Allocatable:: AORecBlock(:), TmpRecBlock(:)
      Real*8, Allocatable:: AOTriBlock(:)
      Real*8, Allocatable:: MOTriBlock(:)
      Real*8, Allocatable:: EigenVecBlock(:), EigenValBlock(:)
      Real*8, Allocatable:: EigenVecTot(:), EigenValTot(:)
      Real*8, Allocatable:: Energies(:)
      Integer, Allocatable:: IndT(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
*                                                                      *
************************************************************************
*                                                                      *
      nOrbAllTot = nOrbAll(1)
      nOrbAllMax = nOrbAll(1)
      lRecTot = nOrbAll(1)*nOrbAll(1)
      Do iSym = 2, mSym
         nOrbAllTot = nOrbAllTot + nOrbAll(iSym)
         nOrbAllMax = Max(nOrbAllMax, nOrbAll(iSym))
         lRecTot = lRecTot + nOrbAll(iSym)*nOrbAll(iSym)
      End Do
*

*     A blockmatrix of the size Orb(iSym) X Orb(iSym) is
*     allocated and set to zero
      Call mma_allocate(AORecBlock,nOrbAllMax**2,Label='AORecBlock')
      Call mma_allocate(TmpRecBlock,nOrbAllMax**2,Label='TmpRecBlock')
      Call mma_allocate(AOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2,
     &                  Label='AOTriBlock')

      If(Diagonalize) Then
         Call mma_allocate(MOTriBlock,nOrbAllMax*(nOrbAllMax+1) / 2,
     &                     Label='MOTriBlock')
         Call mma_allocate(EigenVecBlock,nOrbAllMax*nOrbAllMax,
     &                     Label='EigenVecBlock')
         Call mma_allocate(EigenValBlock,nOrbAllMax,
     &                     Label='EigenValBlock')
         Call mma_allocate(EigenVecTot,lRecTot,Label='EigenVecTot')
         Call mma_allocate(EigenValTot,nOrbAllTot,Label='EigenValTot')
         Call mma_allocate(Energies,nOrbAllTot,Label='Energies')
         Call mma_allocate(IndT,7,mSym,Label='IndT')
         Energies(:)=Zero
      End If

      AORecBlock(:)=Zero
      TmpRecBlock(:)=Zero
      AOTriBlock(:)=Zero

*
*     Setup a pointer to a symmetryblock in rect or tri representation.
      ipSymRec(1) = 0
      ipSymTri(1) = 0
      ipSymLin(1) = 0
      Do iSym = 2, 8
         ipSymTri(iSym) = ipSymTri(iSym-1)
     &                  + (nOrbAll(iSym-1)) * (nOrbAll(iSym-1)+1)/2
         ipSymRec(iSym) = ipSymRec(iSym-1)
     &                  + (nOrbAll(iSym-1))**2
         ipSymLin(iSym) = ipSymLin(iSym-1)
     &                  + (nOrbAll(iSym-1))
      End Do
*
      Do iSym = 1, mSym

         If(nOrbAll(iSym).ne.0) Then
*
*           Set initial values for the eigenvectors to the CMO-matrix.
*           If we transform the CMO-matrix in the same way as the MO-matrix
*           needs to be treated to be diagonal we get a transformation from
*           AO to canonical basis. Also sets the energies to zero since they
*           have no physical relevance in this basis.
            If(Diagonalize) Then
               Do i = 1, nOrbAll(iSym)**2
                  EigenVecBlock(i) = CMO(ipSymRec(iSym)+i)
               End Do
            End If
*     Transform the symmetryblock to AO-density
            Call DGEMM_('N','N',nOrbAll(iSym),nOrbAll(iSym),
     &                 nOrbAll(iSym),1.0d0 , CMO(ipSymRec(iSym)+1),
     &                 nOrbAll(iSym),Work(ip_Density(iSym)),
     &                 nOrbAll(iSym), 0.0d0, TmpRecBlock,
     &                 nOrbAll(iSym))
            Call DGEMM_('N','T',nOrbAll(iSym),nOrbAll(iSym),
     &                 nOrbAll(iSym),
     &                 1.0d0,TmpRecBlock,nOrbAll(iSym),
     &                       CMO(ipSymRec(iSym)+1),nOrbAll(iSym),
     &                 0.0d0,AORecBlock,nOrbAll(iSym))
*            Call RecPrt('AODens:','(20F8.5)',AORecBlock,
*     &                  nOrb(iSym),nOrb(iSym))
*            Call RecPrt('MODens:','(20F8.5)',Work(ip_MORecBlock),
*     &                  nOrb(iSym), nOrb(iSym))
            Call Fold_Mat(1,nOrbAll(iSym),AORecBlock,AOTriBlock)
            call dcopy_(nOrbAll(iSym)*(nOrbAll(iSym)+1)/2,
     &                 AOTriBlock,1,
     &                 TriDens(1+ipSymTri(iSym)),1)

            If(Diagonalize) Then
*     Make a normal folded matrix
*
               index = 0
               Do i = 1, nOrbAll(iSym)
                  Do j = 1,i
                     MOTriBlock(1+index) =  Work(ip_Density(iSym)+
     &                    j-1+(i-1)*(nOrbAll(iSym)))
                     index = index + 1
                  End Do
               End Do
*
               Call NIDiag(MOTriBlock,EigenVecBlock,
     &                     nOrbAll(iSym),nOrbAll(iSym),0)



               Do i = 1, nOrbAll(iSym)
                  EigenValBlock(i) =MOTriBlock(iTri(i,i))
               End Do
*
               Call JacOrd3(EigenValBlock,EigenVecBlock,
     &                      nOrbAll(iSym),nOrbAll(iSym))

#ifdef _DEBUGPRINT_
               Write(6,*) 'The sorted eigenvalues are'
               Do i = 1, nOrbAll(iSym)
                  Write(6,*) EigenValBlock(i)
               End Do
               Write(6,*) 'Eigenvectors sorted'
               Do i=1, nOrbAll(iSym)**2
                  Write(6,*) EigenVecBlock(i)
               End Do
#endif
*
               call dcopy_(nOrbAll(iSym)**2,EigenVecBlock,1,
     &                    EigenVecTot(1+ipSymRec(iSym)),1)
               call dcopy_(nOrbAll(iSym),EigenValBlock,1,
     &                    EigenValTot(1+ipSymLin(iSym)),1)

            End If

*
         End If
      End Do

*     Put the MP2 natural orbitals we just produced on disk to the file
*     MP2ORB.

      If(Diagonalize) Then
         LuMP2=50
         LuMP2=IsFreeUnit(LuMP2)
*        Build the TypeIndex array
         Do iSym=1,mSym
            IndT(1,iSym)=nFro(iSym)
            IndT(2,iSym)=nOcc(iSym)
            IndT(3,iSym)=0
            IndT(4,iSym)=0
            IndT(5,iSym)=0
            IndT(6,iSym)=nOrb(iSym)-nFro(iSym)-nOcc(iSym)-nDel(iSym)
            IndT(7,iSym)=nDel(iSym)
         End Do
         Note='*  Natural MP2 orbitals'
         Call WrVec('MP2ORB',LuMP2,'COEI',mSym,nOrbAll,nOrbAll,
     &              EigenVecTot,EigenValTot,Energies,IndT,Note)
*        Create a molden-file
         iUHF = 0
         Call Molden_Interface(iUHF,'MP2ORB','MD_MP2')
*
      End If



      Call mma_deallocate(AORecBlock)
      Call mma_deallocate(TmpRecBlock)
      Call mma_deallocate(AOTriBlock)
*
      If(Diagonalize) Then
         Call mma_deallocate(MOTriBlock)
         Call mma_deallocate(EigenVecBlock)
         Call mma_deallocate(EigenValBlock)
         Call mma_deallocate(EigenVecTot)
         Call mma_deallocate(EigenValTot)
         Call mma_deallocate(Energies)
         Call mma_deallocate(IndT)
      End If

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nOccAll)
      End
