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
      SubRoutine Cho_SetAtomShl(irc,iAtomShl,n)
      use ChoArr, only: iSOShl
C
C     Purpose: set mapping from shell to atom (i.e., center).
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Integer iAtomShl(n)
#include "cholesky.fh"
#include "choprint.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Character(LEN=14), Parameter:: SecNam = 'Cho_SetAtomShl'

      Character(LEN=LENIN8) AtomLabel(MxBas)

      Integer, Parameter:: Info_Debug = 4

#if defined (_DEBUGPRINT_)
      Logical, Parameter:: Debug = .True.
#else
      Logical, Parameter:: Debug = .False.
#endif

      Integer, Allocatable:: nBas_per_Atom(:)
      Integer, Allocatable:: nBas_Start(:)

      If (Debug) Then
         Write(Lupri,*) '>>> Enter ',SecNam
      End If

C     Check.
C     ------

      irc = 0
      If (nSym .ne. 1) Then ! does not work with symmetry
         irc = 1
         If (Debug) Then
            Write(Lupri,*) '>>> Exit ',SecNam,
     &      ' (error exit: symmetry not allowed!)'
         End If
         Return
      End If
      If (n .lt. nShell) Then
         Call Cho_Quit(SecNam//': iAtomShl not allocated correctly!',
     &                 104)
      End If

C     Get number of atoms.
C     --------------------

C     Call Get_nAtoms_All(nAtom)
      Call Get_iScalar('Bfn Atoms',nAtom)

C     Get atomic labels and basis function labels.
C     --------------------------------------------

      Call Get_cArray('Unique Basis Names',AtomLabel,LENIN8*nBasT)

C     Allocate and get index arrays for indexation of basis functions on
C     each atom.
C     ------------------------------------------------------------------

      Call mma_allocate(nBas_per_Atom,nAtom,Label='nBas_per_Atom')
      Call mma_allocate(nBas_Start,nAtom,Label='nBas_Start')
      Call BasFun_Atom(nBas_per_Atom,nBas_Start,
     &                 AtomLabel,nBasT,nAtom,Debug)

C     Set shell-to-atom mapping.
C     --------------------------

      Do iAtom = 1,nAtom
         i1 = nBas_Start(iAtom)
         i2 = i1 + nBas_per_Atom(iAtom) - 1
         Do i = i1,i2
            iAtomShl(iSOShl(i)) = iAtom
         End Do
      End Do

      If (iPrint.ge.Info_Debug .or. Debug) Then
         nErr = 0
         Write(Lupri,*)
         Write(Lupri,*) SecNam,': shell-to-atom mapping:'
         numSh  = 7
         nBatch = (nShell-1)/numSh + 1
         Do iBatch = 1,nBatch
            If (iBatch .eq. nBatch) Then
               nSh = nShell - numSh*(nBatch-1)
            Else
               nSh = numSh
            End If
            iSh0 = numSh*(iBatch-1)
            iSh1 = iSh0 + 1
            iSh2 = iSh0 + nSh
            Write(Lupri,'(/,A,7(1X,I9))')
     &      'Shell:',(iSh,iSh=iSh1,iSh2)
            Write(Lupri,'(A,7(1X,I9))')
     &      'Atom :',(iAtomShl(iSh),iSh=iSh1,iSh2)
            Do iSh = iSh1,iSh2
               If (iAtomShl(iSh).lt.1 .or. iAtomShl(iSh).gt.nAtom) Then
                  nErr = nErr + 1
               End If
            End Do
         End Do
         If (nErr .ne. 0) Then
            Call Cho_Quit(SecNam//': shell-to-atom init failed!',104)
         End If
      End If

C     Deallocations.
C     --------------

      Call mma_deallocate(nBas_Start)
      Call mma_deallocate(nBas_per_Atom)

      If (Debug) Then
         Write(Lupri,*) '>>> Exit ',SecNam
      End If

      End
