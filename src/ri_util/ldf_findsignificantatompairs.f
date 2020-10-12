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
      Subroutine LDF_FindSignificantAtomPairs(irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: find significant atom pairs based on Cauchy-Schwarz
C              screening. The integral diagonal is computed as a
C              biproduct. Results are stored in ldf_atom_pair_info.fh
C
C     Return codes: 0 for success, 1 otherwise.
C
      Implicit None
      Integer irc
#include "localdf.fh"
#include "WrkSpc.fh"

      Character*28 SecNam
      Parameter (SecNam='LDF_FindSignificantAtomPairs')

      Integer nRSAP
      Integer ip_RSAP, l_RSAP

      Real*8 tau
      Real*8 CutInt_, CutInt

      ! Init return code
      irc=0

      ! Set CutInt so low that all diagonals are computed
      ! Save old value (in order to reset it after diag calc)
      Call LDF_GetCutInt(CutInt_)
      CutInt=1.0d-99
      Call LDF_SetCutInt(CutInt)

      ! Get rough list of significant atom pairs based on Cauchy-Schwarz
      ! and estimated integral diagonals.
      tau=Thr_Prescreen**2
      nRSAP=0
      ip_RSAP=0
      Call LDF_RoughSAP(tau,nRSAP,ip_RSAP,irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_RoughSAP returned code',irc
         irc=1
         Return
      End If

      ! Set up list of significant atom pairs, computing integral
      ! diagonals as a biproduct.
      tau=Thr_Prescreen**2
      Call LDF_SAP(tau,nRSAP,iWork(ip_RSAP),irc)
      If (irc.ne.0) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': LDF_SAP returned code',irc
         irc=1
         Return
      End If

      ! Deallocate rough list of significant atom pairs
      l_RSAP=2*nRSAP
      Call GetMem('LDF_AP','Free','Inte',ip_RSAP,l_RSAP)

      ! Reset CutInt
      Call LDF_SetCutInt(CutInt_)

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_RoughSAP(tau,nAtomPair,ip_AtomPair,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Get rough estimate of significant atom pairs.
C
      Implicit None
      Real*8  tau ! screening threshold
      Integer nAtomPair
      Integer ip_AtomPair
      Integer irc
#include "WrkSpc.fh"

      Integer nShell, nAtom
      Integer ip_Dmax, l_Dmax
      Integer ip_Tmax, l_Tmax
      Integer nShell_i, nShell_j
      Integer ip_i, ip_j
      Integer iAtom, jAtom
      Integer iS, jS
      Integer ijAtom, jiAtom
      Integer l_AtomPair

      Real*8 Dmax_All

      Integer  LDF_nShell, LDF_nAtom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell, LDF_nAtom, LDF_nShell_Atom, LDF_lShell_Atom

      Integer i, j
      Integer iShell_i, iShell_j
      Real*8  Tmax, Dmax
      iShell_i(i)=iWork(ip_i-1+i)
      iShell_j(i)=iWork(ip_j-1+i)
      Tmax(i,j)=Work(ip_Tmax-1+nShell*(j-1)+i)
      Dmax(i,j)=Work(ip_Dmax-1+nAtom*(j-1)+i)

      ! Init return code
      irc=0

      ! Get total number of valence shells
      nShell=LDF_nShell()

      ! Get total number of atoms
      nAtom=LDF_nAtom()

      ! Allocate atom pair diagonal
      l_Dmax=nAtom*nAtom
      Call GetMem('LDF_Dmax','Allo','Real',ip_Dmax,l_Dmax)

      ! Allocate and get estimated shell pair diagonals
      l_Tmax=nShell*nShell
      Call GetMem('LDF_Tmax','Allo','Real',ip_Tmax,l_Tmax)
      Call Shell_MxSchwz(nShell,Work(ip_Tmax))

      ! Find max for each atom pair
      Call Cho_dZero(Work(ip_Dmax),l_Dmax)
      Do jAtom=1,nAtom
         nShell_j=LDF_nShell_Atom(jAtom)
         ip_j=LDF_lShell_Atom(jAtom)
         ijAtom=ip_Dmax-1+nAtom*(jAtom-1)+jAtom
         Do jS=1,nShell_j
            Do iS=jS,nShell_j
               Work(ijAtom)=max(Work(ijAtom),
     &                          Tmax(iShell_j(iS),iShell_j(jS)))
            End Do
         End Do
         Do iAtom=jAtom+1,nAtom
            nShell_i=LDF_nShell_Atom(iAtom)
            ip_i=LDF_lShell_Atom(iAtom)
            ijAtom=ip_Dmax-1+nAtom*(jAtom-1)+iAtom
            Do jS=1,nShell_j
               Do iS=1,nShell_i
                  Work(ijAtom)=max(Work(ijAtom),
     &                             Tmax(iShell_i(iS),iShell_j(jS)))
               End Do
            End Do
            jiAtom=ip_Dmax-1+nAtom*(iAtom-1)+jAtom
            Work(jiAtom)=Work(ijAtom)
         End Do
      End Do

      ! Deallocate shell pair diagonal
      Call GetMem('LDF_Tmax','Free','Real',ip_Tmax,l_Tmax)

      ! Find global max among atom pairs
      Dmax_All=Dmax(1,1)
      Do iAtom=2,nAtom
         Do jAtom=1,iAtom
            Dmax_All=max(Dmax_All,Dmax(iAtom,jAtom))
         End Do
      End Do

      ! Compute nAtomPair
      nAtomPair=0
      Do iAtom=1,nAtom
         Do jAtom=1,iAtom
            If (Dmax_All*Dmax(iAtom,jAtom) .gt. tau) Then
               nAtomPair=nAtomPair+1
            End If
         End Do
      End Do

      ! Allocate atom pair list
      l_AtomPair=2*nAtomPair
      Call GetMem('LDF_AP','Allo','Inte',ip_AtomPair,l_AtomPair)

      ! Set atom pair list
      ijAtom=0
      Do iAtom=1,nAtom
         Do jAtom=1,iAtom
            If (Dmax_All*Dmax(iAtom,jAtom) .gt. tau) Then
               iWork(ip_AtomPair+2*ijAtom)=iAtom
               iWork(ip_AtomPair+2*ijAtom+1)=jAtom
               ijAtom=ijAtom+1
            End If
         End Do
      End Do

      ! Deallocate atom pair diagonal
      Call GetMem('LDF_Dmax','Free','Real',ip_Dmax,l_Dmax)

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_SAP(tau,nRSAP,ID_RSAP,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Get significant atom pairs.
C
      Implicit None
      Real*8  tau
      Integer nRSAP
      Integer ID_RSAP(2,nRSAP)
      Integer irc
#include "WrkSpc.fh"

      Integer ip_TmpDiag, l_TmpDiag
      Integer iRSAP
      Integer iAtom, jAtom
      Integer ni, nj

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      ! Init return code
      irc=0

      ! Return if nothing to do
      If (nRSAP.lt.1) Return

      ! Allocate tmp memory for diagonal
      l_TmpDiag=0
      Do iRSAP=1,nRSAP
         iAtom=ID_RSAP(1,iRSAP)
         jAtom=ID_RSAP(2,iRSAP)
         ni=LDF_nBas_Atom(iAtom)
         If (iAtom.eq.jAtom) Then
            l_TmpDiag=l_TmpDiag+ni*(ni+1)/2
         Else If (iAtom.gt.jAtom) Then
            nj=LDF_nBas_Atom(jAtom)
            l_TmpDiag=l_TmpDiag+ni*nj
         Else
            Call WarningMessage(2,'LDF_SAP: iAtom<jAtom')
            Call LDF_Quit(1)
         End If
      End Do
      Call GetMem('TmpDiag','Allo','Real',ip_TmpDiag,l_TmpDiag)

      ! Compute diagonal
      Call LDF_ComputeRSAPDiagonal(nRSAP,ID_RSAP,l_TmpDiag,
     &                             Work(ip_TmpDiag))

      ! Set significant atom pairs in ldf_atom_pair_info.fh
      ! Store pointers to diagonal blocks
      Call LDF_SetAPI(nRSAP,ID_RSAP,l_TmpDiag,Work(ip_TmpDiag))

      ! Deallocate TmpDiag
      Call GetMem('TmpDiag','Free','Real',ip_TmpDiag,l_TmpDiag)

c Avoid unused argument warnings
      If (.False.) Call Unused_real(tau)
      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_ComputeRSAPDiagonal(nRSAP,ID_RSAP,l_Diag,Diag)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Compute diagonal for atom pairs listed in ID_RSAP.
C
      Implicit None
      Integer nRSAP
      Integer ID_RSAP(2,nRSAP)
      Integer l_Diag
      Real*8  Diag(l_Diag)
#include "WrkSpc.fh"

      Integer ip_iOff, l_iOff
      Integer ip_SewWrk, l_SewWrk
      Integer l, ID
      Integer iRSAP
      Integer iAtom, jAtom
      Integer ni, nj

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer i
      Integer iOff
      iOff(i)=iWork(ip_iOff-1+i)

      ! Compute diagonal offset array
      l_iOff=nRSAP+1
      Call GetMem('iOff','Allo','Inte',ip_iOff,l_iOff)
      l=1
      Do iRSAP=1,nRSAP
         iWork(ip_iOff-1+iRSAP)=l
         iAtom=ID_RSAP(1,iRSAP)
         jAtom=ID_RSAP(2,iRSAP)
         ni=LDF_nBas_Atom(iAtom)
         If (iAtom.eq.jAtom) Then
            l=l+ni*(ni+1)/2
         Else
            nj=LDF_nBas_Atom(jAtom)
            l=l+ni*nj
         End If
      End Do
      iWork(ip_iOff+nRSAP)=l
#if defined (_DEBUGPRINT_)
      If (l_Diag.lt.(iOff(nRSAP+1)-1)) Then
         Call WarningMessage(2,
     &                   'LDF_ComputeRSAPDiagonal: dimension mismatch!')
         Call LDF_Quit(1)
      End If
#endif

      ! Compute diagonal integrals (parallelzation over atom pairs)
      Call Init_Tsk(ID,nRSAP)
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)
      Call Cho_dZero(Diag,iOff(nRSAP+1)-1)
      Do While (Rsv_Tsk(ID,iRSAP))
         l=iOff(iRSAP+1)-iOff(iRSAP)
         Call LDF_ComputeAPDiagonal(ID_RSAP(1,iRSAP),
     &                              ID_RSAP(2,iRSAP),
     &                              l,Diag(iOff(iRSAP)))
      End Do
      Call GAdGOP(Diag,iOff(nRSAP+1)-1,'+')
      Call xRlsMem_Ints()
      Call Free_Tsk(ID)

      ! Deallocation
      Call GetMem('iOff','Free','Inte',ip_iOff,l_iOff)

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_ComputeAPDiagonal(iAtom,jAtom,l_Diag,Diag)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Compute diagonal for atom pair iAtom,jAtom
C
      Implicit None
      Integer iAtom, jAtom
      Integer l_Diag
      Real*8  Diag(l_Diag)
#include "WrkSpc.fh"

      Character*21 SecNam
      Parameter (SecNam='LDF_ComputeAPDiagonal')

      Integer iOff, l
      Integer iSi, iSj
      Integer ni, nj
      Integer iS, jS
      Integer ip_i, ip_j

      Integer  LDF_lShell_Atom, LDF_nShell_Atom, LDF_nBasSh_Atom
      External LDF_lShell_Atom, LDF_nShell_Atom, LDF_nBasSh_Atom

      Integer i
      Integer iShell_i, iShell_j
      iShell_i(i)=iWork(ip_i-1+i)
      iShell_j(i)=iWork(ip_j-1+i)

      ip_i=LDF_lShell_Atom(iAtom)
      ip_j=LDF_lShell_Atom(jAtom)

      iOff=1
      If (iAtom.eq.jAtom) Then
         Do iSi=1,LDF_nShell_Atom(iAtom)
            ni=LDF_nBasSh_Atom(iSi,iAtom)
            iS=iShell_i(iSi)
            Do iSj=1,iSi-1
               nj=LDF_nBasSh_Atom(iSj,iAtom)
               jS=iShell_i(iSj)
               l=ni*nj
               Call LDF_CAPD(iS,jS,l,Diag(iOff))
               iOff=iOff+l
            End Do
            l=ni*(ni+1)/2
            Call LDF_CAPD(iS,iS,l,Diag(iOff))
            iOff=iOff+l
         End Do
      Else If (iAtom.gt.jAtom) Then
         Do iSj=1,LDF_nShell_Atom(jAtom)
            nj=LDF_nBasSh_Atom(iSj,jAtom)
            jS=iShell_j(iSj)
            Do iSi=1,LDF_nShell_Atom(iAtom)
               ni=LDF_nBasSh_Atom(iSi,iAtom)
               iS=iShell_i(iSi)
               l=ni*nj
               Call LDF_CAPD(iS,jS,l,Diag(iOff))
               iOff=iOff+l
            End Do
         End Do
      Else
         Call WarningMessage(2,SecNam//': iAtom<jAtom')
         Call LDF_Quit(1)
      End If

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_CAPD(iS,jS,l_Diag,Diag)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: compute diagonal (iS jS | iS jS)
C
      Implicit None
      Integer iS, jS
      Integer l_Diag
      Real*8  Diag(l_Diag)
#include "localdf_int.fh"

      External Integral_WrOut_LDF_Diag

      SHA=iS
      SHB=jS
      SHC=iS
      SHD=jS

      Call Eval_IJKL(iS,jS,iS,jS,Diag,l_Diag,Integral_WrOut_LDF_Diag)

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_SetAPI(nRSAP,ID_RSAP,l_TmpDiag,TmpDiag)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: set atom pair info
C
C     - tbp, March 2011: use block matrix for diagonals (i.e. diagonal
C                        blocks no longer stored in triangular format)
C
      Implicit None
      Integer nRSAP
      Integer ID_RSAP(2,nRSAP)
      Integer l_TmpDiag
      Real*8  TmpDiag(l_TmpDiag)
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

#if defined (_DEBUGPRINT_)
      Integer iAtomPair
      Integer uv
      Integer ip_LT, l_LT
      Integer ip_Q, l_Q
      Real*8  Mx
#endif
      Integer i0, ii, l
      Integer iRSAP
      Integer iAtom, jAtom
      Integer ip_APDmax, l_APDmax
      Integer ni, nj
      Integer ip

      Real*8 Dmax_All, tau

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Integer i

C     Count significant atom pairs.
C     =============================

      ! Set screening threshold
      tau=Thr_Prescreen**2

      ! Allocate atom pair max diag
      l_APDmax=nRSAP
      Call GetMem('APDmax','Allo','Real',ip_APDmax,l_APDmax)

      ! Find max. diagonal, globally and per atom pair
      Dmax_All=-9.9d9
      i0=1
      Do iRSAP=1,nRSAP
         iAtom=ID_RSAP(1,iRSAP)
         jAtom=ID_RSAP(2,iRSAP)
         ni=LDF_nBas_Atom(iAtom)
         If (iAtom.eq.jAtom) Then
            l=ni*(ni+1)/2
         Else If (iAtom.gt.jAtom) Then
            nj=LDF_nBas_Atom(jAtom)
            l=ni*nj
         Else
            Call WarningMessage(2,'LDF_SetAPI: iAtom<jAtom [1]')
            Call LDF_Quit(1)
            l=0
         End If
         ii=ip_APDmax-1+iRSAP
         Work(ii)=TmpDiag(i0)
         Do i=i0+1,i0+l-1
            Work(ii)=max(Work(ii),TmpDiag(i))
         End Do
         Dmax_All=max(Dmax_All,Work(ii))
         i0=i0+l
      End Do

      ! Count atom pairs
      NumberOfAtomPairs=0
      i0=ip_APDmax-1
      Do iRSAP=1,nRSAP
         If ((Dmax_All*Work(i0+iRSAP)).gt.tau) Then
            NumberOfAtomPairs=NumberOfAtomPairs+1
         End If
      End Do

C     Set atom pair info and diagonal blocks.
C     =======================================

      ! Allocate and set AP_Atoms
      l_AP_Atoms=2*NumberOfAtomPairs
      Call GetMem('LDFAPA','Allo','Inte',ip_AP_Atoms,l_AP_Atoms)
      i0=ip_APDmax-1
      i=0
      Do iRSAP=1,nRSAP
         iAtom=ID_RSAP(1,iRSAP)
         jAtom=ID_RSAP(2,iRSAP)
         If ((Dmax_All*Work(i0+iRSAP)).gt.tau) Then
            iWork(ip_AP_Atoms+2*i)=iAtom
            iWork(ip_AP_Atoms+2*i+1)=jAtom
            i=i+1
         End If
      End Do
      If (i.ne.NumberOfAtomPairs) Then
         Call WarningMessage(2,'LDF_SetAPI: i != NumberOfAtomPairs [1]')
         Call LDF_Quit(1)
      End If

      ! Allocate diagonals (as block matrices)
      Call LDF_AllocateBlockMatrix('APD',ip_AP_Diag)
      l_AP_Diag=NumberOfAtomPairs
      Call LDF_AllocateBlockMatrix('APB',ip_AP_DiagBak)
      l_AP_DiagBak=NumberOfAtomPairs

      ! Set diagonals
      i0=ip_APDmax-1
      ii=1
      i=0
      Do iRSAP=1,nRSAP
         iAtom=ID_RSAP(1,iRSAP)
         jAtom=ID_RSAP(2,iRSAP)
         If (iAtom.eq.jAtom) Then
            ni=LDF_nBas_Atom(iAtom)
            l=ni*(ni+1)/2
            If ((Dmax_All*Work(i0+iRSAP)).gt.tau) Then
               ip=iWork(ip_AP_Diag+i)
               Call LDF_LT2Q(iAtom,TmpDiag(ii),Work(ip))
               Call dCopy_(ni**2,Work(ip),1,
     &                          Work(iWork(ip_AP_DiagBak+i)),1)
               i=i+1
            End If
         Else If (iAtom.gt.jAtom) Then
            ni=LDF_nBas_Atom(iAtom)
            nj=LDF_nBas_Atom(jAtom)
            l=ni*nj
            If ((Dmax_All*Work(i0+iRSAP)).gt.tau) Then
               ip=iWork(ip_AP_Diag+i)
               Call dCopy_(l,TmpDiag(ii),1,Work(ip),1)
               ip=iWork(ip_AP_DiagBak+i)
               Call dCopy_(l,TmpDiag(ii),1,Work(ip),1)
               i=i+1
            End If
         Else
            Call WarningMessage(2,'LDF_SetAPI: iAtom<jAtom [2]')
            Call LDF_Quit(1)
            l=0 ! avoid compiler warning
         End If
         ii=ii+l
      End Do
      If (i.ne.NumberOfAtomPairs) Then
         Call WarningMessage(2,'LDF_SetAPI: i != NumberOfAtomPairs [2]')
         Call LDF_Quit(1)
      End If

      ! Deallocate atom pair max diag
      l_APDmax=nRSAP
      Call GetMem('APDmax','Free','Real',ip_APDmax,l_APDmax)

C     Allocate and initialize number of linearly dependent one-center
C     functions and number of two-center functions.
C     ===============================================================

      l_AP_1CLinDep=2*NumberOfAtomPairs
      l_AP_2CFunctions=2*NumberOfAtomPairs
      Call GetMem('AP1CLD','Allo','Inte',ip_AP_1CLinDep,l_AP_1CLinDep)
      Call GetMem('AP2CFN','Allo','Inte',ip_AP_2CFunctions,
     &                                    l_AP_2CFunctions)
      Call iZero(iWork(ip_AP_1CLinDep),l_AP_1CLinDep)
      Call iZero(iWork(ip_AP_2CFunctions),l_AP_2CFunctions)

#if defined (_DEBUGPRINT_)
C     Debug: "unit test" of LDF_Q2LT and LDF_LT2Q
C     ===========================================
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=iWork(ip_AP_Atoms+2*(iAtomPair-1))
         jAtom=iWork(ip_AP_Atoms+2*(iAtomPair-1)+1)
         If (iAtom.eq.jAtom) Then
            If (LDF_nBas_Atom(iAtom).gt.0) Then
               l_LT=LDF_nBas_Atom(iAtom)*(LDF_nBas_Atom(iAtom)+1)/2
               Call GetMem('LT','Allo','Real',ip_LT,l_LT)
               l_Q=LDF_nBas_Atom(iAtom)**2
               Call GetMem('Q','Allo','Real',ip_Q,l_Q)
               ip=iWork(ip_AP_Diag-1+iAtomPair)
               Call LDF_Q2LT(iAtom,Work(ip),Work(ip_LT))
               Call LDF_LT2Q(iAtom,Work(ip_LT),Work(ip_Q))
               Call dAXPY_(l_Q,-1.0d0,Work(ip),1,Work(ip_Q),1)
               Mx=0.0d0
               Do uv=0,l_Q-1
                  Mx=max(Mx,abs(Work(ip_Q+uv)))
               End Do
               If (Mx.gt.1.0d-14) Then
                  Call WarningMessage(2,
     &                          'LDF_SetAPI: LDF_Q2LT/LDF_LT2Q failure')
                  Write(6,'(A,1P,D20.10)') 'Max asymmetry:',Mx
                  Call LDF_Quit(1)
               End If
               Call GetMem('Q','Free','Real',ip_Q,l_Q)
               Call GetMem('LT','Free','Real',ip_LT,l_LT)
            End If
         End If
      End Do
#endif

      End
