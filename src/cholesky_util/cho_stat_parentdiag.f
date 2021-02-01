************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_Stat_ParentDiag()
C
C     Thomas Bondo Pedersen, February 2006.
C
C     Purpose: print statistics about the diagonals from which the
C              Cholesky vectors are computed. I.e., whether they are
C              one-center or two-center diagonals. Does not work with
C              symmetry!
C
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Character*19 SecNam
      Parameter (SecNam = 'Cho_Stat_ParentDiag')

      Character*(LENIN8) AtomLabel(MxBas)

      Logical Debug

      Integer  Cho_iSumElm
      External Cho_iSumElm

      Parameter (numAt = 6)
      Real*8  Ratio(numAt)
      Real*8  RClass(3)
      Integer iClass(4)
      integer nPseudo


      nBas_Start(i)=iWork(ip_nBas_Start-1+i)
      nBas_per_Atom(i)=iWork(ip_nBas_per_Atom-1+i)
      iBF2Atom(i)=iWork(ip_iBF2Atom-1+i)
      Coord(i,j)=Work(ip_Coord-1+3*(j-1)+i)
      mapRS2F(i,j)=iWork(ip_mapRS2F-1+2*(j-1)+i)
      nPC1(i)=iWork(ip_PC1-1+i)
      RC2(i)=Work(ip_RC2-1+i)

C     Set debug flag.
C     ---------------

#if defined (_DEBUGPRINT_)
      Debug = .True.
#else
      Debug = .False.
#endif

C     Symmetry not allowed.
C     ---------------------

      If (nSym .ne. 1) Return

C     Check if anything to do.
C     ------------------------

      If (NumCho(1) .ne. NumChT) Then
         Call SysAbendMsg(SecNam,'NumChT .ne. NumCho(1)',' ')
      End If
      If (NumChT .lt. 1) Return

C     Get number of atoms.
C     --------------------

C     Call Get_nAtoms_All(nAtom)
      Call Get_iScalar('Bfn Atoms',nAtom)
      Call Get_iScalar('Pseudo atoms',nPseudo)

C     Get atomic labels and basis function labels.
C     --------------------------------------------

      Call Get_cArray('Unique Basis Names',AtomLabel,LENIN8*nBasT)

C     Allocate and get index arrays for indexation of basis functions on
C     each atom.
C     ------------------------------------------------------------------

      l_nBas_per_Atom = nAtom
      l_nBas_Start    = nAtom
      Call GetMem('nB_per_Atom','Allo','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Allo','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                 AtomLabel,nBasT,nAtom,Debug)

C     Allocate and get nuclear coordinates.
C     -------------------------------------

      l_Coord = 3*nAtom
      Call GetMem('Coord','Allo','Real',ip_Coord,l_Coord)
C     Call Get_dArray('Unique Coordinates',Work(ip_Coord),l_Coord)
      Call Get_dArray('Bfn Coordinates',Work(ip_Coord),l_Coord)

C     Allocate and set mapping from basis function to atom.
C     -----------------------------------------------------

      l_iBF2Atom = nBasT
      Call GetMem('iBF2Atom','Allo','Inte',ip_iBF2Atom,l_iBF2Atom)
      ip0 = ip_iBF2Atom - 1
      Do iAtom = 1,nAtom
         i1 = nBas_Start(iAtom)
         i2 = i1 + nBas_per_Atom(iAtom) - 1
         Do i = i1,i2
            iWork(ip0+i) = iAtom
         End Do
      End Do

      If (Debug) Then
         Write(Lupri,*)
         Write(Lupri,*) SecNam,': mapping from basis function to atom:'
         Do i = 1,nBasT
            Write(Lupri,*) 'Basis function ',i,' is centered on atom ',
     &                     iBF2Atom(i),' labeled ',AtomLabel(i)(1:LENIN)
         End Do
      End If

C     Allocate and set mapping from 1st reduced set to full storage.
C     --------------------------------------------------------------

      l_mapRS2F = 2*nnBstRT(1)
      Call GetMem('mapRS2F','Allo','Inte',ip_mapRS2F,l_mapRS2F)
      Call Cho_RStoF(iWork(ip_mapRS2F),2,nnBstRT(1),1)

C     Allocate counters.
C     ------------------

      l_PC1 = nAtom
      l_RC2 = NumChT
      Call GetMem('PC1','Allo','Inte',ip_PC1,l_PC1)
      Call GetMem('RC2','Allo','Real',ip_RC2,l_RC2)

C     Compute statistics.
C     -------------------

      Rmin = 1.0d15
      Rmax = -1.0d15
      Rave = 0.0d0
      Call iCopy(l_PC1,[0],0,iWork(ip_PC1),1)
      Call dCopy_(l_RC2,[0.0d0],0,Work(ip_RC2),1)
      Do iVec = 1,NumChT
         iRS1 = InfVec(iVec,1,1)
         iA = mapRS2F(1,iRS1)
         iB = mapRS2F(2,IRS1)
         iAtA = iBF2Atom(iA)
         iAtB = iBF2Atom(iB)
         If (iAtA .eq. iAtB) Then
            iWork(ip_PC1-1+iAtA) = iWork(ip_PC1-1+iAtA) + 1
         Else
            R = sqrt((Coord(1,iAtA)-Coord(1,iAtB))**2
     &              +(Coord(2,iAtA)-Coord(2,iAtB))**2
     &              +(Coord(3,iAtA)-Coord(3,iAtB))**2)
            Work(ip_RC2-1+iVec) = R
            Rmin = min(Rmin,R)
            Rmax = max(Rmax,R)
            Rave = Rave + R
         End If
      End Do
      nTot1 = Cho_iSumElm(iWork(ip_PC1),nAtom)
      nTot2 = NumChT - nTot1
      If (nTot2 .gt. 0) Then
         Rave = Rave/dble(nTot2)
      Else If (nTot2 .lt. 0) Then
         Call SysAbendMsg(SecNam,'Setup error!','nTot2 < 0')
      End If

C     Print overall statisctics.
C     --------------------------

      Write(Lupri,'(//,2X,A,/,2X,A)')
     & 'Parent Diagonals','----------------'
      Write(Lupri,'(/,A,I9,A,F7.2,A)')
     & 'Number of vectors from 1-center diagonals:',nTot1,
     & ' (',1.0d2*dble(nTot1)/dble(NumChT),'%)'
      Write(Lupri,'(A,I9,A,F7.2,A)')
     & 'Number of vectors from 2-center diagonals:',nTot2,
     & ' (',1.0d2*dble(nTot2)/dble(NumChT),'%)'

C     Print statistics for 1-center vectors.
C     --------------------------------------

      Write(Lupri,'(/,1X,A)')
     & 'Vectors from 1-center diagonals:'
      nBatch = (nAtom-nPseudo-1)/numAt + 1 ! exclude pseudo-atoms
      Do iBatch = 1,nBatch
         If (iBatch .eq. nBatch) Then
!           exclude pseudo-atoms
            nAt = nAtom - nPseudo - numAt*(nBatch-1)
         Else
            nAt = numAt
         End If
         iAt0 = numAt*(iBatch-1)
         iAt1 = iAt0 + 1
         iAt2 = iAt0 + nAt
         Do i = 1,nAt
            iAt = iAt0 + i
            If (nBas_per_Atom(iAt) .lt. 1) Then
               If (nPC1(iAt) .gt. 0) Then
                  Call SysAbendMsg(SecNam,
     &                 'No basis functions, but >0 vectors !?!?',' ')
               End If
               Ratio(i) = 0.0d0
            Else
               Ratio(i) = dble(nPC1(iAt))/dble(nBas_per_Atom(iAt))
            End If
         End Do
         Write(Lupri,'(/,A,6(6X,A))')
     &   'Label              ',
     &    (AtomLabel(nBas_Start(i))(1:LENIN),i=iAt1,iAt2)
         Write(Lupri,'(A,6(1X,I9))')
     &   'Center no.         ',(i,i=iAt1,iAt2)
         Write(Lupri,'(A,6(1X,I9))')
     &   'Vectors (M)        ',(nPC1(i),i=iAt1,iAt2)
         Write(Lupri,'(A,6(1X,I9))')
     &   'Basis functions (N)',(nBas_per_Atom(i),i=iAt1,iAt2)
         Write(Lupri,'(A,6(1X,F9.2))')
     &   'Ratio (M/N)        ',(Ratio(i),i=1,nAt)
         If (iBatch .ne. nBatch) Write(Lupri,*)
      End Do

C     Print statistics for 2-center vectors.
C     --------------------------------------

      If (nTot2 .gt. 0) Then
         Write(Lupri,'(/,1X,A)')
     &   'Vectors from 2-center diagonals:'
         RClass(1) = Rave - (Rave-Rmin)/2.0d0
         RClass(2) = Rave
         RClass(3) = Rave + (Rmax-Rave)/2.0d0
         Call iCopy(4,[0],0,iClass,1)
         nChk = 0
         Do iVec = 1,NumChT
            If (RC2(iVec) .gt. 0.0d0) Then
               nChk = nChk + 1
               If (RC2(iVec) .le. RClass(1)) Then
                  iClass(1) = iClass(1) + 1
               Else If (RC2(iVec) .le. RClass(2)) Then
                  iClass(2) = iClass(2) + 1
               Else If (RC2(iVec) .le. RClass(3)) Then
                  iClass(3) = iClass(3) + 1
               Else
                  iClass(4) = iClass(4) + 1
               End If
            End If
         End Do
         If (nChk .ne. nTot2) Then
C FAQ
C --- Replace the call to SysAbendMsg with just a warning.
C
C            Call SysAbendMsg(SecNam,'nChk .ne. nTot2',' ')
C
C --- TODO/FIX  figure out if with ghost atoms is just a mistmatch
C               or there is really a bug
C
            write(6,*) SecNam//': Warning! (nChk .ne. nTot2); could be',
     &' due to the presence of ghost atoms.'
         End If
         Write(Lupri,'(/,A,1P,3D15.5)')
     &   'Min, average, and max center distance: ',Rmin,Rave,Rmax
         Write(Lupri,'(A,D12.2,A,I9)')
     &   '#vectors with center distance                R <= ',
     &   RClass(1),': ',iClass(1)
         Write(Lupri,'(A,D12.2,A,D12.2,A,I9)')
     &   '#vectors with center distance ',RClass(1),' < R <= ',
     &   RClass(2),': ',iClass(2)
         Write(Lupri,'(A,D12.2,A,D12.2,A,I9)')
     &   '#vectors with center distance ',RClass(2),' < R <= ',
     &   RClass(3),': ',iClass(3)
         Write(Lupri,'(A,D12.2,A,12X,A,I9)')
     &   '#vectors with center distance ',RClass(3),' < R    ',
     &   ': ',iClass(4)
      End If

C     Deallocations.
C     --------------

      Call GetMem('RC2','Free','Real',ip_RC2,l_RC2)
      Call GetMem('PC1','Free','Inte',ip_PC1,l_PC1)
      Call GetMem('mapRS2F','Free','Inte',ip_mapRS2F,l_mapRS2F)
      Call GetMem('iBF2Atom','Free','Inte',ip_iBF2Atom,l_iBF2Atom)
      Call GetMem('Coord','Free','Real',ip_Coord,l_Coord)
      Call GetMem('nB_Start','Free','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call GetMem('nB_per_Atom','Free','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)

      End
