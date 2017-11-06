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
      Subroutine Build_Mp2Dens(ip_TriDens,ip_Density,CMO,
     &                          mSym,nOrbAll,nOccAll,Diagonalize)

      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "corbinf.fh"
*
      Integer   ip_AOTriBlock
      Real*8    CMO(*)
      Integer   ipSymRec(8)
      Integer   ipSymTri(8)
      Integer   ipSymLin(8)
      Integer   nOrbAll(8), nOccAll(8), ip_Density(8)
      Logical   Diagonalize,AddFragments
      Character*30 Note
*
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
*
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
      Call GetMem('AORecBlock','Allo','Real',ip_AORecBlock,
     &            nOrbAllMax**2)
      Call GetMem('TmpRecBlock','Allo','Real',ip_TmpRecBlock,
     &            nOrbAllMax**2)
      Call GetMem('AOTriBlock','Allo','Real',ip_AOTriBlock,
     &            nOrbAllMax*(nOrbAllMax+1) / 2)
      If(Diagonalize) Then
         Call GetMem('MOTriBlock','Allo','Real',ip_MOTriBlock,
     &            nOrbAllMax*(nOrbAllMax+1) / 2)
         Call GetMem('EigenVecBlock','Allo','Real',ip_EigenVecBlock,
     &               nOrbAllMax*nOrbAllMax)
         Call GetMem('EigenValBlock','Allo','Real',ip_EigenValBlock,
     &               nOrbAllMax)
         Call GetMem('EigenVectors','Allo','Real',ip_EigenVecTot,
     &               lRecTot)
         Call GetMem('EigenValues','Allo','Real',ip_EigenValTot,
     &               nOrbAllTot)
         Call GetMem('Energies','Allo','Real',ip_Energies,nOrbAllTot)
         Call GetMem('IndT','Allo','Inte',ip_IndT,7*mSym)
         Call FZero(Work(ip_Energies),nOrbAllTot)
      End If


      Call FZero(Work(ip_AORecBlock),nOrbAllMax**2)
      Call FZero(Work(ip_TmpRecBlock),nOrbAllMax**2)
      Call FZero(Work(ip_AOTriBlock),nOrbAllMax
     &                            * (nOrbAllMax+1)/2)

*
*     Setup a pointer to a symmetryblock in rect or tri representation.
      ipSymRec(1) = 0
      ipSymTri(1) = 0
      ipSymLin(1) = 0
      Do iSym = 2, 8
         ipSymTri(iSym) = ipSymTri(iSym-1)
     &                  + (nOrbAll(iSym-1))
     &                  * (nOrbAll(iSym-1)+1)/2
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
                  Work(ip_EigenVecBlock + i-1)
     &               = CMO(ipSymRec(iSym)+i)
               End Do
            End If
*     Transform the symmetryblock to AO-density
            Call DGEMM_('N','N',nOrbAll(iSym),nOrbAll(iSym),
     &                 nOrbAll(iSym),1.0d0 , CMO(ipSymRec(iSym)+1),
     &                 nOrbAll(iSym),Work(ip_Density(iSym)),
     &                 nOrbAll(iSym), 0.0d0, Work(ip_TmpRecBlock),
     &                 nOrbAll(iSym))
            Call DGEMM_('N','T',nOrbAll(iSym),nOrbAll(iSym),
     &                 nOrbAll(iSym),1.0d0,Work(ip_TmpRecBlock),
     &                 nOrbAll(iSym),CMO(ipSymRec(iSym)+1),
     &                 nOrbAll(iSym),0.0d0, Work(ip_AORecBlock),
     &                 nOrbAll(iSym))
*            Call RecPrt('AODens:','(20F8.5)',Work(ip_AORecBlock),
*     &                  nOrb(iSym),nOrb(iSym))
*            Call RecPrt('MODens:','(20F8.5)',Work(ip_MORecBlock),
*     &                  nOrb(iSym), nOrb(iSym))
            Call Fold_Mat(1,nOrbAll(iSym),Work(ip_AORecBlock),
     &                    Work(ip_AOTriBlock))
            call dcopy_(nOrbAll(iSym)*(nOrbAll(iSym)+1)/2,
     &                 Work(ip_AOTriBlock),1,
     &                 Work(ip_TriDens+ipSymTri(iSym)),1)

            If(Diagonalize) Then
*     Make a normal folded matrix
*
               index = 0
               Do i = 1, nOrbAll(iSym)
                  Do j = 1,i
                     Work(ip_MOTriBlock+index) =  Work(ip_Density(iSym)+
     &                    j-1+(i-1)*(nOrbAll(iSym)))
                     index = index + 1
                  End Do
               End Do
*
               Call NIDiag(Work(ip_MOTriBlock),Work(ip_EigenVecBlock),
     &                     nOrbAll(iSym),nOrbAll(iSym),0)



               Do i = 1, nOrbAll(iSym)
                  Work(ip_EigenValBlock+i-1) =
     &                    Work(ip_MOTriBlock+iTri(i,i)-1)
               End Do
*
               Call JacOrd3(Work(ip_EigenValBlock),
     &                      Work(ip_EigenVecBlock),
     &                      nOrbAll(iSym),nOrbAll(iSym))

#ifdef _DEBUG_
               Write(6,*) 'The sorted eigenvalues are'
               Do i = 1, nOrbAll(iSym)
                  Write(6,*) Work(ip_EigenValBlock+i-1)
               End Do
               Write(6,*) 'Eigenvectors sorted'
               Do i=1, nOrbAll(iSym)**2
                  Write(6,*) Work(ip_EigenVecBlock+i-1)
               End Do
#endif
*
               call dcopy_(nOrbAll(iSym)**2,
     &                    Work(ip_EigenVecBlock),1,
     &                    Work(ip_EigenVecTot+ipSymRec(iSym)),1)
               call dcopy_(nOrbAll(iSym),
     &                    Work(ip_EigenValBlock),1,
     &                    Work(ip_EigenValTot+ipSymLin(iSym)),1)

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
         iOff=ip_IndT
         Do iSym=1,mSym
            iWork(iOff+0)=nFro(iSym)
            iWork(iOff+1)=nOcc(iSym)
            iWork(iOff+2)=0
            iWork(iOff+3)=0
            iWork(iOff+4)=0
            iWork(iOff+5)=nOrb(iSym)-nFro(iSym)-nOcc(iSym)-nDel(iSym)
            iWork(iOff+6)=nDel(iSym)
            iOff=iOff+7
         End Do
         Note='*  Natural MP2 orbitals'
         Call WrVec('MP2ORB',LuMP2,'COEI',mSym,nOrbAll,nOrbAll,
     &              Work(ip_EigenVecTot),Work(ip_EigenValTot),
     &              Work(ip_Energies),iWork(ip_IndT),Note)
*        Create a molden-file
         AddFragments = .True.
         iUHF = 0
         Call Molden_Interface(iUHF,'MP2ORB','MD_MP2',AddFragments)
*
      End If



      Call GetMem('AORecBlock','Free','Real',ip_AORecBlock,
     &            nOrbAllMax**2)
      Call GetMem('TmpRecBlock','Free','Real',ip_TmpRecBlock,
     &            nOrbAllMax**2)
      Call GetMem('AOTriBlock','Free','Real',ip_AOTriBlock,
     &            nOrbAllMax*(nOrbAllMax+1) / 2)
*
      If(Diagonalize) Then
         Call GetMem('MOTriBlock','Free','Real',ip_MOTriBlock,
     &            nOrbAllMax*(nOrbAllMax+1) / 2)
         Call GetMem('EigenVecBlock','Free','Real',ip_EigenVecBlock,
     &               nOrbAllMax*nOrbAllMax)
         Call GetMem('EigenValBlock','Free','Real',ip_EigenValBlock,
     &               nOrbAllMax)
         Call GetMem('EigenVectors','Free','Real',ip_EigenVecTot,
     &               lRecTot)
         Call GetMem('EigenValues','Free','Real',ip_EigenValTot,
     &               nOrbAllTot)
         Call GetMem('Energies','Free','Real',ip_Energies,nOrbAllTot)
         Call GetMem('IndT','Free','Inte',ip_IndT,7*mSym)
      End If

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nOccAll)
      End
