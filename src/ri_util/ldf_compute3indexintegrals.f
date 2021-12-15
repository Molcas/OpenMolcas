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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_Compute3IndexIntegrals_1(AB,C,tau,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: compute 3-index integrals of the form (u_A v_B|K_C)
C              where |u_A v_B) are AO products on atom pair AB and K_C
C              are one-center auxiliary functions on atom C. Note that
C              *all* auxiliary functions on C are included (including
C              linearly dependent ones).
C
C     tau is the prescreening threshold. Note that the integral
C     prescreening arrays must have been set up before calling this
C     routine.
C
      Implicit None
      Integer AB
      Integer C
      Real*8  tau
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_int3.fh"
#include "localdf_bas.fh"
#include "ldf_integral_prescreening_info.fh"

      Character*28 SecNam
      Parameter (SecNam='LDF_Compute3IndexIntegrals_1')

      External Int_LDF_3Indx_1

      Integer  LDF_nBas_Atom, LDF_nBasAux_Atom
      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBas_Atom, LDF_nBasAux_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
#if defined (_DEBUGPRINT_)
      Integer  LDF_nAtom
      External LDF_nAtom
#endif

      Real*8 tau2

      Integer ip_iOffRow, l_iOffRow
      Integer ip_iOffCol, l_iOffCol
      Integer A, B
      Integer nu, nv
      Integer M
      Integer l_xInt
      Integer nShellA, nShellB, nAuxShellC
      Integer ipA, ipB, ipC
      Integer n, ip0
      Integer iS, jS, kS, lS
      Integer iShell, jShell, kShell, lShell
      Integer K, ijK, jiK
      Integer ip_SewWrk, l_SewWrk

      Integer i, j
      Integer nBasSh
      Integer iOffCol
      Integer iOffRow
      Integer AP_Atoms
      Real*8  Imax
      Real*8  Gmax
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iOffCol(i)=iWork(ip_iOffCol-1+i)
      iOffRow(i,j)=iWork(ip_iOffRow-1+nShellA*(j-1)+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      Imax(i,j)=Work(iWork(ip_IDiag+2*(AB-1)+1)-1+nShellA*(j-1)+i)
      Gmax(i)=Work(iWork(ip_GDiag_1C+2*(C-1)+1)-1+i)

#if defined (_USE_APD_INTEGRALS_)
      If (.True.) Then
         Call WarningMessage(0,SecNam//': Using APD integrals!')
         Call LDF_APD3IndexIntegrals_1(AB,C,l_xInt_,xInt)
         Return
      End If
#endif

#if defined (_DEBUGPRINT_)
      If (AB.lt.1 .or. AB.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': AB out of bounds!')
         Call LDF_Quit(1)
      End If
      M=LDF_nShell_Atom(AP_Atoms(1,AB))*LDF_nShell_Atom(AP_Atoms(2,AB))
      If (iWork(ip_IDiag+2*(AB-1)).ne.M) Then
         Call WarningMessage(2,SecNam//': IDiag not properly set!')
         Call LDF_Quit(1)
      End If
      If (C.lt.1 .or. C.gt.LDF_nAtom()) Then
         Call WarningMessage(2,SecNam//': C out of bounds!')
         Call LDF_Quit(1)
      End If
      M=LDF_nAuxShell_Atom(C)
      If (iWork(ip_GDiag_1C+2*(C-1)).ne.M) Then
         Call WarningMessage(2,SecNam//': GDiag_1C not properly set!')
         Call LDF_Quit(1)
      End If
#endif

      ! Square screening threshold
      tau2=tau**2

      ! Get atoms of atom pair AB
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)

      ! Get number of basis functions on atoms A and B
      nu=LDF_nBas_Atom(A)
      nv=LDF_nBas_Atom(B)

      ! Get number of auxiliary functions on C
      M=LDF_nBasAux_Atom(C)

      ! Check integral array dimension
      nRow=nu*nv
      l_xInt=nRow*M
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! Get number of shells on each atom A,B
      nShellA=LDF_nShell_Atom(A)
      nShellB=LDF_nShell_Atom(B)

      ! Get pointers to lists of shells on A,B
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1

      ! Allocate row offset array
      l_iOffRow=nShellA*nShellB
      Call GetMem('3IiOffR','Allo','Inte',ip_iOffRow,l_iOffRow)

      ! Set row offset array
      n=0
      Do jS=1,nShellB
         jShell=iWork(ipB+jS)
         ip0=ip_iOffRow-1+nShellA*(jS-1)
         Do iS=1,nShellA
            iWork(ip0+iS)=n
            iShell=iWork(ipA+iS)
            n=n+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do
#if defined (_DEBUGPRINT_)
      If (n.ne.nRow) Then
         Call WarningMessage(2,SecNam//': row dimension problem!')
         Call LDF_Quit(1)
      End If
#endif

      ! Get number of auxiliary shells on C
      nAuxShellC=LDF_nAuxShell_Atom(C)

      ! Get pointer to auxiliary shells on C
      ipC=LDF_lAuxShell_Atom(C)-1

      ! Allocate columns offset array
      l_iOffCol=nAuxShellC
      Call GetMem('3IiOffC','Allo','Inte',ip_iOffCol,l_iOffCol)

      ! Set column offset array
      n=0
      ip0=ip_iOffCol-1
      Do kS=1,nAuxShellC
         iWork(ip0+kS)=n
         kShell=iWork(ipC+kS)
         n=n+nBasSh(kShell)
      End Do
#if defined (_DEBUGPRINT_)
      If (n.ne.M) Then
         Call WarningMessage(2,SecNam//': column dimension problem!')
         Call LDF_Quit(1)
      End If
#endif

      ! Allocate memory for Seward
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      kShell=nShell_Valence+nShell_Auxiliary+1
      SHC=kShell
      If (A.eq.B) Then
         ! compute (u_A v_B|K_C) for u_A >= v_B (shell blocks)
         Do lS=1,nAuxShellC
            lShell=iWork(ipC+lS)
            SHD=lShell
            iCol0=iOffCol(lS)
            Do jS=1,nShellB
               jShell=iWork(ipB+jS)
               SHB=jShell
               Do iS=jS,nShellA
                  If (Imax(iS,jS)*Gmax(lS).ge.tau2) Then
                     iShell=iWork(ipA+iS)
                     SHA=iShell
                     iRow0=iOffRow(iS,jS)
                     Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                              xInt,l_xInt,Int_LDF_3Indx_1)
                  End If
               End Do
            End Do
         End Do
         ! Now copy to get lower triangle u_A < v_B (shell blocks)
         Do K=0,M-1
            ip0=nRow*K
            Do jS=2,nShellB
               jShell=iWork(ipB+jS)
               Do iS=1,jS-1
                  iShell=iWork(ipA+iS)
                  Do j=1,nBasSh(jShell)
                     Do i=1,nBasSh(iShell)
                        ijK=ip0+iOffRow(iS,jS)+nBasSh(iShell)*(j-1)+i
                        jiK=ip0+iOffRow(jS,iS)+nBasSh(jShell)*(i-1)+j
                        xInt(ijK)=xInt(jiK)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      Else If (A.gt.B) Then
         ! compute (u_A v_B|K_C)
         Do lS=1,nAuxShellC
            lShell=iWork(ipC+lS)
            SHD=lShell
            iCol0=iOffCol(lS)
            Do jS=1,nShellB
               jShell=iWork(ipB+jS)
               SHB=jShell
               Do iS=1,nShellA
                  If (Imax(iS,jS)*Gmax(lS).ge.tau2) Then
                     iShell=iWork(ipA+iS)
                     SHA=iShell
                     iRow0=iOffRow(iS,jS)
                     Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                              xInt,l_xInt,Int_LDF_3Indx_1)
                  End If
               End Do
            End Do
         End Do
      Else
         Call WarningMessage(2,SecNam//': A<B')
         Call LDF_Quit(1)
      End If

      ! Deallocate
      Call xRlsMem_Ints()
      Call GetMem('3IiOffC','Free','Inte',ip_iOffCol,l_iOffCol)
      Call GetMem('3IiOffR','Free','Inte',ip_iOffRow,l_iOffRow)

      nRow=0
      iRow0=0
      iCol0=0
      SHA=0
      SHB=0
      SHC=0
      SHD=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Compute3IndexIntegrals_2(AB,CD,tau,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: compute 3-index integrals of the form (u_A v_B|K_CD)
C              where |u_A v_B) are AO products on atom pair AB and K_CD
C              are two-center auxiliary functions on atom pair CD.
C
C     tau is the prescreening threshold. Note that the integral
C     prescreening arrays must have been set up before calling this
C     routine.
C
      Implicit None
      Integer AB
      Integer CD
      Real*8  tau
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_int.fh"
#include "localdf_bas.fh"
#include "ldf_integral_prescreening_info.fh"

      Character*28 SecNam
      Parameter (SecNam='LDF_Compute3IndexIntegrals_2')

      Integer  LDF_nBasAux_Pair
      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Real*8 tau2

      Integer M
      Integer A, B
      Integer l_xInt
      Integer nShellA, nShellB
      Integer n
      Integer ipA, ipB
      Integer iS, jS, klS
      Integer iShell, jShell, kShell, lShell
      Integer ip0, ij0, ji0, ij, ji
      Integer K
      Integer ip_SewWrk, l_SewWrk

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer i2CList
      Real*8  Gmax
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      i2CList(i,j)=iWork(ip_2CList-1+l_2CList_1*(j-1)+i)
      Gmax(i)=Work(iWork(ip_GDiag_2C+2*(CD-1)+1)-1+i)

#if defined (_USE_APD_INTEGRALS_)
      If (.True.) Then
         Call WarningMessage(0,SecNam//': Using APD integrals!')
         Call LDF_APD3IndexIntegrals_2(AB,CD,l_xInt_,xInt)
         Return
      End If
#endif

#if defined (_DEBUGPRINT_)
      If (AB.lt.1 .or. AB.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': AB out of bounds!')
         Call LDF_Quit(1)
      End If
      If (CD.lt.1 .or. CD.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': CD out of bounds!')
         Call LDF_Quit(1)
      End If
      If (iWork(ip_GDiag_2C+2*(CD-1)).lt.1) Then
         Call WarningMessage(2,SecNam//': GDiag_2C not properly set!')
         Call LDF_Quit(1)
      End If
#endif

      ! Get number of two-center auxiliary functions on CD
      M=AP_2CFunctions(1,CD)

      ! Return if nothing to do
      If (M.lt.1) Return

      ! Square prescreening threshold
      tau2=tau**2

      ! Get atoms of atom pair AB
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)

      ! Get row dimension
      nRow_uvJ=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      If (nRow_uvJ.lt.1) Return

      ! Check integral array dimension
      l_xInt=nRow_uvJ*M
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! Set column index array including one-center functions.
      ! Then shift indices to exclude one-center functions.
      Call LDF_SetIndxG(CD)
      Call LDF_ColMod(LDF_nBasAux_Pair(CD)-M)

      ! Get number of shells on A and B
      nShellA=LDF_nShell_Atom(A)
      nShellB=LDF_nShell_Atom(B)

      ! Allocate row offset array
      l_iOff=nShellA*nShellB
      Call GetMem('3I2iOff','Allo','Inte',ip_iOff,l_iOff)

      ! Get pointers to shell lists
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1

      ! Set row offset array
      n=0
      Do jS=1,nShellB
         jShell=iWork(ipB+jS)
         ip0=ip_iOff-1+nShellA*(jS-1)
         Do iS=1,nShellA
            iWork(ip0+iS)=n
            iShell=iWork(ipA+iS)
            n=n+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do
#if defined (_DEBUGPRINT_)
      If (n.ne.nRow_uvJ) Then
         Call WarningMessage(2,SecNam//': row dimension problem!')
         Call LDF_Quit(1)
      End If
#endif

      ! Allocate memory for Seward
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      Do klS=1,l_2CList_2
         kShell=i2CList(1,klS)
         lShell=i2CList(2,klS)
         SPAB=i2CList(3,klS)
         Call LDF_CI_uvJ_PS(AB,kShell,lShell,Gmax(klS),tau2,l_xInt,xInt)
      End Do

      ! Free Seward memory
      Call xRlsMem_Ints()

      ! Postprocessing for A=B:
      ! Only lower (block) triangle computed => transpose to get upper
      If (A.eq.B) Then
         Do K=0,M-1
            ip0=nRow_uvJ*K
            Do jS=2,nShellB
               jShell=iWork(ipB+jS)
               Do iS=1,jS-1
                  iShell=iWork(ipA+iS)
                  ij0=ip0+iWork(ip_iOff-1+nShellA*(jS-1)+iS)
                  ji0=ip0+iWork(ip_iOff-1+nShellB*(iS-1)+jS)
                  Do j=1,nBasSh(jShell)
                     Do i=1,nBasSh(iShell)
                        ij=ij0+nBasSh(iShell)*(j-1)+i
                        ji=ji0+nBasSh(jShell)*(i-1)+j
                        xInt(ij)=xInt(ji)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If

      ! Deallocate row offset array
      Call GetMem('3I2iOff','Free','Inte',ip_iOff,l_iOff)
      ip_iOff=0
      l_iOff=0
      nRow_uvJ=0
      iOffuv=0

      ! Unset indices
      Call LDF_UnsetIndxG()
      SHA=0
      SHB=0
      SHC=0
      SHD=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ColMod(M)
      Implicit None
      Integer M
#include "WrkSpc.fh"
#include "localdf_int.fh"

      Integer i, j, ip0, ij0, ij

      ip0=ip_IndxG2-1
      Do j=1,l_IndxG2_2
         ij0=ip0+l_IndxG2_1*(j-1)
         Do i=1,l_IndxG2_1
            ij=ij0+i
            iWork(ij)=max(iWork(ij)-M,0)
         End Do
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CI_uvJ_PS(iAtomPair,iShell,jShell,Gmax,tau2,
     &                         l_xInt,xInt)
C
C     Thomas Bondo Pedersen, November 2010.
C     - modified version of LDF_CI_uvJ
C
C     Compute integrals (uv|J) where J is an auxiliary function in
C     shell pair iShell,jShell [iShell=dummy-shell if 1-center].
C
      Implicit None
      Integer iAtomPair
      Integer iShell, jShell
      Real*8  Gmax, tau2
      Integer l_xInt
      Real*8  xInt(l_xInt)
#include "WrkSpc.fh"
#include "localdf_int.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_integral_prescreening_info.fh"

      Character*13 SecNam
      Parameter (SecNam='LDF_CI_uvJ_PS')

      External Int_LDF_uvJ

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer kAtom, lAtom
      Integer nShell_kAtom, nShell_lAtom
      Integer kS, lS
      Integer kShell, lShell
      Integer ipk, ipl
#if defined (_DEBUGPRINT_)
      Integer M
#endif

      Integer i, j
      Integer AP_Atoms
      Real*8  Imax
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      Imax(i,j)=Work(iWork(ip_IDiag+2*(iAtomPair-1)+1)-1
     &                             +nShell_kAtom*(j-1)+i)

#if defined (_DEBUGPRINT_)
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': iAtomPair out of bounds!')
         Call LDF_Quit(1)
      End If
      M=LDF_nShell_Atom(AP_Atoms(1,iAtomPair))
     & *LDF_nShell_Atom(AP_Atoms(2,iAtomPair))
      If (iWork(ip_IDiag+2*(iAtomPair-1)).ne.M) Then
         Call WarningMessage(2,SecNam//': IDiag not properly set!')
         Call LDF_Quit(1)
      End If
#endif

      kAtom=AP_Atoms(1,iAtomPair)
      lAtom=AP_Atoms(2,iAtomPair)
      nShell_kAtom=LDF_nShell_Atom(kAtom)
      nShell_lAtom=LDF_nShell_Atom(lAtom)
      ipk=LDF_lShell_Atom(kAtom)-1
      ipl=LDF_lShell_Atom(lAtom)-1

      SHA=iShell
      SHB=jShell

      If (kAtom.eq.lAtom) Then
         ! NOTE: only lower (block) triangular part is computed!!!
         Do lS=1,nShell_lAtom
            lShell=iWork(ipl+lS)
            SHD=lShell
            Do kS=lS,nShell_kAtom
               If (Imax(kS,lS)*Gmax.ge.tau2) Then
                  kShell=iWork(ipk+kS)
                  SHC=kShell
                  ! iOffuv = row offset
                  iOffuv=iWork(ip_iOff-1+nShell_kAtom*(lS-1)+kS)
                  Call Eval_IJKL(iShell,jShell,kShell,lShell,xInt,
     &                           l_xInt,Int_LDF_uvJ)
               End If
            End Do
         End Do
      Else If (kAtom.gt.lAtom) Then
         Do lS=1,nShell_lAtom
            lShell=iWork(ipl+lS)
            SHD=lShell
            Do kS=1,nShell_kAtom
               If (Imax(kS,lS)*Gmax.ge.tau2) Then
                  kShell=iWork(ipk+kS)
                  SHC=kShell
                  ! iOffuv = row offset
                  iOffuv=iWork(ip_iOff-1+nShell_kAtom*(lS-1)+kS)
                  Call Eval_IJKL(iShell,jShell,kShell,lShell,xInt,
     &                           l_xInt,Int_LDF_uvJ)
               End If
            End Do
         End Do
      Else
         Call WarningMessage(2,SecNam//': kAtom<lAtom')
         Call LDF_Quit(1)
      End If

      End
