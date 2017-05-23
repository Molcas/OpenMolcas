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
      Subroutine LDF_Compute2IndexIntegrals_11(A,B,tau,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: compute 2-index integrals (J|K) where |J) are auxiliary
C              one-center functions on A and |K) on B.
C
C     tau is the prescreening threshold. Note that the integral
C     prescreening arrays must have been set up before calling this
C     routine.
C
      Implicit None
      Integer A
      Integer B
      Real*8  tau
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int3.fh"
#include "ldf_integral_prescreening_info.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_Compute2IndexIntegrals_11')

      External Int_LDF_2Indx_11

      Integer  LDF_nBasAux_Atom, LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBasAux_Atom, LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
#if defined (_DEBUG_)
      Integer  LDF_nAtom
      External LDF_nAtom
#endif

      Real*8 tau2

      Integer nJ, nK
      Integer l_xInt
      Integer nAuxShellA, nAuxShellB
      Integer ipA, ipB
      Integer dShell
      Integer ip_SewWrk, l_SewWrk
      Integer jS, kS
      Integer jShell, kShell
      Integer JK, KJ
      Integer J, K
#if defined (_DEBUG_)
      Integer M
#endif

      Integer i
      Integer nBasSh
      Real*8  Gmax_A, Gmax_B
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      Gmax_A(i)=Work(iWork(ip_GDiag_1C+2*(A-1)+1)-1+i)
      Gmax_B(i)=Work(iWork(ip_GDiag_1C+2*(B-1)+1)-1+i)

#if defined (_USE_APD_INTEGRALS_)
      If (.True.) Then
         Call WarningMessage(0,SecNam//': Using APD integrals!')
         Call LDF_APD2IndexIntegrals_11(A,B,l_xInt_,xInt)
         Return
      End If
#endif

#if defined (_DEBUG_)
      If (A.lt.1 .or. A.gt.LDF_nAtom()) Then
         Call WarningMessage(2,SecNam//': A out of bounds!')
         Call LDF_Quit(1)
      End If
      M=LDF_nAuxShell_Atom(A)
      If (iWork(ip_GDiag_1C+2*(A-1)).ne.M) Then
         Call WarningMessage(2,
     &                       SecNam//': GDiag_1C not properly set [A]!')
         Call LDF_Quit(1)
      End If
      If (B.lt.1 .or. B.gt.LDF_nAtom()) Then
         Call WarningMessage(2,SecNam//': B out of bounds!')
         Call LDF_Quit(1)
      End If
      M=LDF_nAuxShell_Atom(B)
      If (iWork(ip_GDiag_1C+2*(B-1)).ne.M) Then
         Call WarningMessage(2,
     &                       SecNam//': GDiag_1C not properly set [B]!')
         Call LDF_Quit(1)
      End If
#endif

      ! Square threshold
      tau2=tau**2

      ! Get number of auxiliary functions on A and B
      nJ=LDF_nBasAux_Atom(A)
      nK=LDF_nBasAux_Atom(B)

      ! Check integral dimension
      l_xInt=nJ*nK
      If (l_xInt.lt.1) Return
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      End If

      ! Set row dimension (used for extracting integrals)
      ! Stored in localdf_int3.fh
      nRow=nJ

      ! Get number of aux shells on A,B
      nAuxShellA=LDF_nAuxShell_Atom(A)
      nAuxShellB=LDF_nAuxShell_Atom(B)

      ! Get pointers to aux shell lists on A,B
      ipA=LDF_lAuxShell_Atom(A)-1
      ipB=LDF_lAuxShell_Atom(B)-1

      ! Set dummy shell
      dShell=nShell_Valence+nShell_Auxiliary+1

      ! Allocate Seward memory
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      SHA=dShell
      SHC=dShell
      If (A.eq.B) Then
         iCol0=0
         Do kS=1,nAuxShellB
            kShell=iWork(ipB+kS)
            SHD=kShell
            iRow0=iCol0
            Do jS=kS,nAuxShellA
               jShell=iWork(ipA+jS)
               If (Gmax_A(jS)*Gmax_B(kS).ge.tau2) Then
                  SHB=jShell
                  Call Eval_IJKL(dShell,jShell,dShell,kShell,xInt,
     &                           l_xInt,Int_LDF_2Indx_11)
               End If
               iRow0=iRow0+nBasSh(jShell)
            End Do
            iCol0=iCol0+nBasSh(kShell)
         End Do
         iCol0=nBasSh(iWork(ipB+1))
         Do kS=2,nAuxShellB
            kShell=iWork(ipB+kS)
            iRow0=0
            Do jS=1,kS-1
               jShell=iWork(ipA+jS)
               Do K=iCol0+1,iCol0+nBasSh(kShell)
                  Do J=iRow0+1,iRow0+nBasSh(jShell)
                     JK=nRow*(K-1)+J
                     KJ=nRow*(J-1)+K
                     xInt(JK)=xInt(KJ)
                  End Do
               End Do
               iRow0=iRow0+nBasSh(jShell)
            End Do
            iCol0=iCol0+nBasSh(kShell)
         End Do
      Else
         iCol0=0
         Do kS=1,nAuxShellB
            kShell=iWork(ipB+kS)
            SHD=kShell
            iRow0=0
            Do jS=1,nAuxShellA
               jShell=iWork(ipA+jS)
               If (Gmax_A(jS)*Gmax_B(kS).ge.tau2) Then
                  SHB=jShell
                  Call Eval_IJKL(dShell,jShell,dShell,kShell,xInt,
     &                           l_xInt,Int_LDF_2Indx_11)
               End If
               iRow0=iRow0+nBasSh(jShell)
            End Do
            iCol0=iCol0+nBasSh(kShell)
         End Do
      End If

      ! Release Seward memory
      Call xRlsMem_Ints()

      SHA=0
      SHB=0
      SHC=0
      SHD=0
      nRow=0
      iRow0=0
      iCol0=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Compute2IndexIntegrals_12(A,CD,tau,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: compute 2-index integrals (J|K) where |J) are auxiliary
C              one-center functions on atom A and |K) are two-center
C              auxiliary functions on atom pair CD.
C
C     tau is the prescreening threshold. Note that the integral
C     prescreening arrays must have been set up before calling this
C     routine.
C
      Implicit None
      Integer A
      Integer CD
      Real*8  tau
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_integral_prescreening_info.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_Compute2IndexIntegrals_12')

      External Int_LDF_2Indx_12

      Integer  LDF_nBasAux_Pair
      Integer  LDF_nBasAux_Atom, LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBasAux_Pair
      External LDF_nBasAux_Atom, LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
#if defined (_DEBUG_)
      Integer  LDF_nAtom
      External LDF_nAtom
#endif

      Real*8 tau2

      Integer nJ, nK
      Integer l_xInt
      Integer nAuxShellA, ipA
      Integer dShell
      Integer ip_SewWrk, l_SewWrk
      Integer klS, kShell, lShell
      Integer jS, jShell
#if defined (_DEBUG_)
      Integer M
#endif

      Integer i, j
      Integer nBasSh
      Integer AP_2CFunctions
      Integer i2CList
      Real*8  Gmax_A
      Real*8  Gmax_CD
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      i2CList(i,j)=iWork(ip_2CList-1+l_2CList_1*(j-1)+i)
      Gmax_A(i)=Work(iWork(ip_GDiag_1C+2*(A-1)+1)-1+i)
      Gmax_CD(i)=Work(iWork(ip_GDiag_2C+2*(CD-1)+1)-1+i)

#if defined (_USE_APD_INTEGRALS_)
      If (.True.) Then
         Call WarningMessage(0,SecNam//': Using APD integrals!')
         Call LDF_APD2IndexIntegrals_12(A,CD,l_xInt_,xInt)
         Return
      End If
#endif

#if defined (_DEBUG_)
      If (A.lt.1 .or. A.gt.LDF_nAtom()) Then
         Call WarningMessage(2,SecNam//': A out of bounds!')
         Call LDF_Quit(1)
      End If
      M=LDF_nAuxShell_Atom(A)
      If (iWork(ip_GDiag_1C+2*(A-1)).ne.M) Then
         Call WarningMessage(2,SecNam//': GDiag_1C not properly set!')
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

      ! Square threshold
      tau2=tau**2

      ! Get number of auxiliary two-center functions on CD
      nK=AP_2CFunctions(1,CD)
      If (nK.lt.1) Return

      ! Get number of auxiliary one-center functions on A
      nJ=LDF_nBasAux_Atom(A)
      If (nJ.lt.1) Return

      ! Check integral array dimension
      l_xInt=nJ*nK
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      End If

      ! Get number of auxiliary shells on A
      nAuxShellA=LDF_nAuxShell_Atom(A)

      ! Get pointer to list of aux shells on A
      ipA=LDF_lAuxShell_Atom(A)-1

      ! Set dummy shell
      dShell=nShell_Valence+nShell_Auxiliary+1

      ! Set column index array including one-center functions.
      ! Set row dimension for integral extraction
      ! Then shift indices to exclude one-center functions.
      Call LDF_SetIndxG(CD)
      nRow_G=nJ
      Call LDF_ColMod(LDF_nBasAux_Pair(CD)-nK)

      ! Allocate Seward memory
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      SHA=dShell
      Do klS=1,l_2CList_2
         kShell=i2CList(1,klS)
         lShell=i2CList(2,klS)
         SPAB=i2CList(3,klS)
         SHC=kShell
         SHD=lShell
         iOffuv=0
         Do jS=1,nAuxShellA
            jShell=iWork(ipA+jS)
            If (Gmax_A(jS)*Gmax_CD(klS).ge.tau2) Then
               SHB=jShell
               Call Eval_IJKL(dShell,jShell,kShell,lShell,xInt,l_xInt,
     &                        Int_LDF_2Indx_12)
            End If
            iOffuv=iOffuv+nBasSh(jShell)
         End Do
      End Do

      ! Release Seward memory
      Call xRlsMem_Ints()

      ! Unset index arrays
      Call LDF_UnsetIndxG()

      SHA=0
      SHB=0
      SHC=0
      SHD=0
      SPAB=0
      SPCD=0
      nRow_G=0
      nRow_uvJ=0
      iOffuv=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Compute2IndexIntegrals_22(AB,CD,tau,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: compute 2-index integrals (J|K) where |J) are auxiliary
C              two-center functions on atom pair AB and |K) are two-
C              center auxiliary functions on atom pair CD.
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
#include "localdf_bas.fh"
#include "localdf_int2.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_integral_prescreening_info.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_Compute2IndexIntegrals_22')

      External Int_LDF_JK_2P

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Real*8 tau2

      Integer nJ, nK
      Integer l_xInt
      Integer ip_SewWrk, l_SewWrk
      Integer ijS, klS
      Integer iShell, jShell, kShell, lShell
      Integer ii, jj, kk, ll
      Integer ij0, kl0
      Integer ij, kl
      Integer ijkl, klij

      Integer i, j
      Integer nBasSh
      Integer AP_2CFunctions
      Integer iRow
      Integer iCol
      Real*8  Gmax_AB
      Real*8  Gmax_CD
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      iRow(i,j)=iWork(ip_AB_IndxG2-1+l_AB_IndxG2_1*(j-1)+i)
      iCol(i,j)=iWork(ip_CD_IndxG2-1+l_CD_IndxG2_1*(j-1)+i)
      Gmax_AB(i)=Work(iWork(ip_GDiag_2C+2*(AB-1)+1)-1+i)
      Gmax_CD(i)=Work(iWork(ip_GDiag_2C+2*(CD-1)+1)-1+i)

#if defined (_USE_APD_INTEGRALS_)
      If (.True.) Then
         Call WarningMessage(0,SecNam//': Using APD integrals!')
         Call LDF_APD2IndexIntegrals_22(AB,CD,l_xInt_,xInt)
         Return
      End If
#endif

#if defined (_DEBUG_)
      If (AB.lt.1 .or. AB.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': AB out of bounds!')
         Call LDF_Quit(1)
      End If
      If (iWork(ip_GDiag_2C+2*(AB-1)).lt.1) Then
         Call WarningMessage(2,
     &                      SecNam//': GDiag_2C not properly set [AB]!')
         Call LDF_Quit(1)
      End If
      If (CD.lt.1 .or. CD.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': CD out of bounds!')
         Call LDF_Quit(1)
      End If
      If (iWork(ip_GDiag_2C+2*(CD-1)).lt.1) Then
         Call WarningMessage(2,
     &                      SecNam//': GDiag_2C not properly set [CD]!')
         Call LDF_Quit(1)
      End If
#endif

      ! Square threshold
      tau2=tau**2

      ! Get row and column dimensions
      nJ=AP_2CFunctions(1,AB)
      nK=AP_2CFunctions(1,CD)

      ! Check integral dimension
      l_xInt=nJ*nK
      If (l_xInt.lt.1) Return
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      End If

      ! Set indices, redefine row/col dim, and remove one-center
      ! functions
      Call LDF_SetIndx_JK_2P(AB,CD)
      nAB=nJ
      nCD=nK
      Call LDF_ColMod2(LDF_nBasAux_Pair(AB)-nJ,LDF_nBasAux_Pair(CD)-nK)

      ! Allocate memory for Seward
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      If (AB.eq.CD) Then
         Do klS=0,l_CD_2CList_2-1
            kShell=iWork(ip_CD_2CList+3*klS)
            lShell=iWork(ip_CD_2CList+3*klS+1)
            SPCD=iWork(ip_CD_2CList+3*klS+2)
            SHC=kShell
            SHD=lShell
            Do ijS=klS,l_AB_2CList_2-1
               If (Gmax_AB(ijS+1)*Gmax_CD(klS+1).ge.tau2) Then
                  iShell=iWork(ip_AB_2CList+3*ijS)
                  jShell=iWork(ip_AB_2CList+3*ijS+1)
                  SPAB=iWork(ip_AB_2CList+3*ijS+2)
                  SHA=iShell
                  SHB=jShell
                  Call Eval_IJKL(iShell,jShell,kShell,lShell,xInt,
     &                           l_xInt,Int_LDF_JK_2P)
               End If
            End Do
         End Do
         Do klS=1,l_CD_2CList_2-1
            kShell=iWork(ip_CD_2CList+3*klS)
            lShell=iWork(ip_CD_2CList+3*klS+1)
            SPCD=iWork(ip_CD_2CList+3*klS+2)
            Do ijS=0,klS-1
               iShell=iWork(ip_AB_2CList+3*ijS)
               jShell=iWork(ip_AB_2CList+3*ijS+1)
               SPAB=iWork(ip_AB_2CList+3*ijS+2)
               Do ll=1,nBasSh(lShell)
                  kl0=nBasSh(kShell)*(ll-1)
                  Do kk=1,nBasSh(kShell)
                     kl=kl0+kk
                     If (iCol(kl,SPCD).gt.0) Then
                        Do jj=1,nBasSh(jShell)
                           ij0=nBasSh(iShell)*(jj-1)
                           Do ii=1,nBasSh(iShell)
                              ij=ij0+ii
                              If (iRow(ij,SPAB).gt.0) Then
                                 ijkl=nAB*(iCol(kl,SPCD)-1)
     &                               +iRow(ij,SPAB)
                                 klij=nAB*(iRow(ij,SPAB)-1)
     &                               +iCol(kl,SPCD)
                                 xInt(ijkl)=xInt(klij)
                              End If
                           End Do
                        End Do
                     End If
                  End Do
               End Do
            End Do
         End Do
      Else
         Do klS=0,l_CD_2CList_2-1
            kShell=iWork(ip_CD_2CList+3*klS)
            lShell=iWork(ip_CD_2CList+3*klS+1)
            SPCD=iWork(ip_CD_2CList+3*klS+2)
            SHC=kShell
            SHD=lShell
            Do ijS=0,l_AB_2CList_2-1
               If (Gmax_AB(ijS+1)*Gmax_CD(klS+1).ge.tau2) Then
                  iShell=iWork(ip_AB_2CList+3*ijS)
                  jShell=iWork(ip_AB_2CList+3*ijS+1)
                  SPAB=iWork(ip_AB_2CList+3*ijS+2)
                  SHA=iShell
                  SHB=jShell
                  Call Eval_IJKL(iShell,jShell,kShell,lShell,xInt,
     &                           l_xInt,Int_LDF_JK_2P)
               End If
            End Do
         End Do
      End If

      ! Release Seward memory
      Call xRlsMem_Ints()

      ! Unset indices
      Call LDF_UnsetIndx_JK_2P()

      SHA=0
      SHB=0
      SHC=0
      SHD=0
      SPAB=0
      SPCD=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ColMod2(MAB,MCD)
      Implicit None
      Integer MAB, MCD
#include "WrkSpc.fh"
#include "localdf_int2.fh"

      Integer i, j, ip0, ij0, ij

      ip0=ip_AB_IndxG2-1
      Do j=1,l_AB_IndxG2_2
         ij0=ip0+l_AB_IndxG2_1*(j-1)
         Do i=1,l_AB_IndxG2_1
            ij=ij0+i
            iWork(ij)=max(iWork(ij)-MAB,0)
         End Do
      End Do

      ip0=ip_CD_IndxG2-1
      Do j=1,l_CD_IndxG2_2
         ij0=ip0+l_CD_IndxG2_1*(j-1)
         Do i=1,l_CD_IndxG2_1
            ij=ij0+i
            iWork(ij)=max(iWork(ij)-MCD,0)
         End Do
      End Do

      End
