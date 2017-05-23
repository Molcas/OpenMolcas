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
* Copyright (C) 2010-2012, Thomas Bondo Pedersen                       *
************************************************************************
      Subroutine LDF_CheckPairIntegrals(Mode,iAtomPair,l_C,C,irc)
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: check fitting coefficients by comparing exact and
C              approximate integrals for atom pair iAtomPair.
C
      Integer Mode
      Integer iAtomPair
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "localdf.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*22 SecNam
      Parameter (SecNam='LDF_CheckPairIntegrals')

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_nBasAux_Pair

      Integer iAtom, jAtom
      Integer nBas_iAtom, nBas_jAtom
      Integer nShell_iAtom, nShell_jAtom
      Integer M

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)
      nBas_iAtom=LDF_nBas_Atom(iAtom)
      nBas_jAtom=LDF_nBas_Atom(jAtom)
      nShell_iAtom=LDF_nShell_Atom(iAtom)
      nShell_jAtom=LDF_nShell_Atom(jAtom)
      M=LDF_nBasAux_Pair(iAtomPair)
      Call Cho_Head(SecNam//': Integral Check','-',80,6)
      Write(6,'(A,I9)')
     & 'Atom Pair..............................',iAtomPair
      Write(6,'(A,2I9)')
     & 'Atoms..................................',iAtom,jAtom
      Write(6,'(A,2I9)')
     & 'Number of basis functions..............',nBas_iAtom,nBas_jAtom
      Write(6,'(A,2I9)')
     & 'Number of shells.......................',
     & nShell_iAtom,nShell_jAtom
      Write(6,'(A,I9)')
     & 'Number of auxiliary functions..........',M
      Write(6,'(A,1P,D15.6)')
     & 'Target accuracy........................',Thr_Accuracy
      Call xFlush(6)
      If (Mode.eq.1) Then
         Call LDF_CheckPairIntegrals_Robust(iAtomPair,l_C,C,irc)
         If (irc.ne.0) Then
            If (irc.eq.1) Then
               Write(6,'(A)')
     &         '(Delta(AB)|Delta(AB)) matrix not symmetric'
            Else If (irc.eq.2) Then
               Write(6,'(A)')
     &         '(Delta(AB)|Delta(AB)) matrix not positive semidefinite'
            Else If (irc.eq.3) Then
               Write(6,'(A)')
     &         '(Delta(AB)|Delta(AB)) matrix diagonal not consistent'
            Else
               Write(6,'(A,A,I10,A)')
     &         'Non-zero return code from ',
     &         'LDF_CheckPairIntegrals_Robust:',irc,' (unkown)'
            End If
         End If
      Else If (Mode.eq.2) Then
         Call LDF_CheckPairIntegrals_Nonrobust(iAtomPair,l_C,C,irc)
      Else If (Mode.eq.3) Then
         Call LDF_CheckPairIntegrals_HlfNHlf(iAtomPair,l_C,C,irc)
      Else
         Call WarningMessage(2,SecNam//': illegal Mode')
         Call LDF_Quit(1)
      End If
      If (irc.eq.0) Then
         Write(6,'(A,A,I10)')
     &   SecNam,': pair integrals all right for atom pair',iAtomPair
         Call xFlush(6)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPairIntegrals_Robust(AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: check fitting coefficients by comparing exact and
C     approximate integrals for atom pair iAtomPair.
C
C     Checks:
C     1) Symmetry: irc=1 is returned if the diff ints are not symmetric.
C     2) Diagonal: irc=2 is returned if the diagonal diff ints are not
C                  equal to the updated diagonal.
C     3) Accuracy: if 2C functions are included and the max diagonal is
C                  greater than the target accuracy, irc=3 is returned.
C                  Exception: if this is constrained LDF, irc=0 is
C                  returned and a warning is printed.
C     4) PSD: irc=4 is returned if the diff ints do not constitute a
C             positive semidefinite matrix (which implies that the
C             Cauchy-Schwarz inequality is fulfilled).
C
C     Note: this is a debug code and has not been written with
C     efficiency in mind!
C
C     Note: this routine uses the robust integral representation.
C
      Implicit None
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "localdf.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_CheckPairIntegrals_Robust')

      Logical  isSymmetric
      External isSymmetric

      Integer  LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      External LDF_AtomPair_DiagDim, LDF_nBasAux_Pair

      Real*8 Tol
      Parameter (Tol=1.0d-12)

      Integer nAB
      Integer M, Mm
      Integer ip_Int, l_Int
      Integer ip_3IndxInt, l_3IndxInt
      Integer ip_G, l_G
      Integer ip_D
      Integer kl
      Integer n
      Integer ip_Vec, l_Vec
      Integer nAcc

      Real*8 D_max

      ! Init return code
      irc=0

      ! Get and check dimensions
      nAB=LDF_AtomPair_DiagDim(AB)
      M=LDF_nBasAux_Pair(AB)
      Mm=max(M,1)
      If (l_C.lt.(nAB*M)) Then
         Call WarningMessage(2,SecNam//': insufficient array dimension')
         Call LDF_Quit(1)
      End If

      ! Exit with proper return code if nothing to do
      If (M.lt.1) Then
         If (nAB.lt.1) Then
            irc=0
            Return
         End If
      Else
         If (nAB.lt.1) Then
            irc=-1
            Return
         End If
      End If

      ! Allocate and compute exact integrals (uv|kl)
      l_Int=nAB**2
      Call GetMem('CPII','Allo','Real',ip_Int,l_Int)
      Call LDF_ComputeValenceIntegrals(AB,AB,l_Int,Work(ip_Int))

      ! Check symmetry
      If (.not.isSymmetric(Work(ip_Int),nAB,Tol)) Then
         Call WarningMessage(2,
     &        SecNam//': (AB|AB) integrals not symmetric')
         Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
         Call LDF_Quit(1)
      End If

      ! Set indices for G matrix calculation
      Call LDF_SetIndxG(AB)

      ! Allocate and compute 3 index integrals (uv|J)
      l_3IndxInt=nAB*Mm
      Call GetMem('CPI3I','Allo','Real',ip_3IndxInt,l_3IndxInt)
      Call LDF_ComputeIntegrals_uvJ(AB,l_3IndxInt,Work(ip_3IndxInt))

      ! Allocate and compute G matrix (J|K)
      l_G=M**2
      Call GetMem('CPIG','Allo','Real',ip_G,l_G)
      Call LDF_ComputeGMat(AB,M,Work(ip_G))

      ! (uv|kl) := - sum_K (uv|K)*C(kl,K) + (uv|kl)
      Call dGeMM_('N','T',nAB,nAB,M,
     &            -1.0d0,Work(ip_3IndxInt),nAB,C,nAB,
     &            1.0d0,Work(ip_Int),nAB)

      ! (kl|J) := sum_K C(kl,K)*(K|J) - (kl|J)
      Call dGeMM_('N','N',nAB,M,M,
     &            1.0d0,C,nAB,Work(ip_G),Mm,
     &            -1.0d0,Work(ip_3IndxInt),nAB)

      ! (uv|kl) = sum_J C(uv,J)*(kl|J) + (uv|kl)
      Call dGeMM_('N','T',nAB,nAB,M,
     &            1.0d0,C,nAB,Work(ip_3IndxInt),nAB,
     &            1.0d0,Work(ip_Int),nAB)

      ! Deallocation
      Call GetMem('CPIG','Free','Real',ip_G,l_G)
      Call GetMem('CPI3I','Free','Real',ip_3IndxInt,l_3IndxInt)
      Call LDF_UnsetIndxG()

      ! Check 1) Symmetry
      If (irc.eq.0) Then
         If (.not.isSymmetric(Work(ip_Int),nAB,Tol)) Then
            Call WarningMessage(2,
     &        SecNam//': (Delta(AB)|Delta(AB)) integrals not symmetric')
            Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
            irc=1
         End If
      End If

      ! Check 2) Diagonal
      If (irc.eq.0) Then
         ip_D=iWork(ip_AP_Diag-1+AB)
         n=nAB+1
         kl=0
         Do While (kl.lt.nAB .and. irc.eq.0)
            If (abs(Work(ip_D+kl)-Work(ip_Int+n*kl)).gt.Tol) Then
               Call WarningMessage(2,
     &          SecNam//': (Delta(AB)|Delta(AB)) diagonal inconsistent')
               Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
               irc=2
            End If
            kl=kl+1
         End Do
      End If

      ! Check 3) Accuracy
      If (irc.eq.0) Then
         If (LDF2) Then
            D_max=0.0d0
            nAcc=0
            n=nAB+1
            Do kl=0,nAB-1
               If (Work(ip_Int+n*kl).gt.Thr_Accuracy) Then
                  nAcc=nAcc+1
                  D_max=max(D_max,Work(ip_Int+n*kl))
               End If
            End Do
            If (nAcc.ne.0) Then
               Call WarningMessage(2,
     &                   SecNam//': error greater than target accuracy')
               Write(6,'(A,1P,D20.10)') 'Max diagonal:',D_max
               If (LDF_Constraint.eq.-1) Then ! unconstrained LDF
                  irc=3
               End If
            End If
         End If
      End If

      ! Check 4) PSD
      If (irc.eq.0) Then
         l_Vec=nAB**2
         Call GetMem('CPIV','Allo','Real',ip_Vec,l_Vec)
         Call CD_InCore(Work(ip_Int),nAB,Work(ip_Vec),nAB,n,Tol,irc)
         Call GetMem('CPIV','Free','Real',ip_Vec,l_Vec)
         If (irc.ne.0) Then
            Call WarningMessage(2,
     &              SecNam//': (Delta(AB)|Delta(AB)) integrals not PSD')
            Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
            irc=4
         End If
      End If

      ! Deallocation
      Call GetMem('CPII','Free','Real',ip_Int,l_Int)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPairIntegrals_Nonrobust(iAtomPair,l_C,C,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Check fitting coefficients by comparing exact and
C     approximate integrals for atom pair iAtomPair. A report is
C     printed on Cauchy-Schwarz (C-S) violations of the residual matrix.
C     If all integrals are represented with an accuracy better than (or
C     equal to) the max diagonal, irc=0 on exit; else, irc=1.
C     If 2-center functions are included, all errors must be below or
C     equal to the target accuracy. If not, irc=1 is returned.
C
C     Note: this is a debug code and has not been written with
C     efficiency in mind!
C
C     Note: this routine uses the non-robust integral representation.
C
      Implicit None
      Integer iAtomPair
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "WrkSpc.fh"
#include "localdf.fh"
#include "localdf_int.fh"
#include "localdf_bas.fh"
#include "localdf_print.fh"
#include "ldf_atom_pair_info.fh"

      Real*8 Small
      Parameter (Small=1.0d-3)

      Integer ThrPrint
      Parameter (ThrPrint=5)

      Character*32 SecNam
      Parameter (SecNam='LDF_CheckPairIntegrals_Nonrobust')

      External Int_LDF_SQ

      Integer  LDF_nBas_Atom
      Integer  LDF_nShell_Atom, LDF_lShell_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom, LDF_nBasAux_Pair

      Real*8 tC0, tC1, tW0, tW1
      Real*8 Err_max, Err_min, Err_amx, Err_amn, Err_ave, Err_aav
      Real*8 Dij, Dkl, Dijkl, Dmax
      Real*8 CSDevMax, CSDev

      Integer M, Mm
      Integer ip_G, l_G
      Integer iAtom, jAtom
      Integer nBas_iAtom, nBas_jAtom, nB2
      Integer ip_Int, l_Int
      Integer ip_V, l_V
      Integer ip_C_, l_C_
      Integer nShell_iAtom, nShell_jAtom
      Integer ipi, ipj
      Integer iS, jS, iShell, jShell
      Integer kS, lS, kShell, lShell
      Integer iS1, kS1, ijS, klS
      Integer maxSP
      Integer ip_SPB, l_SPB, n
      Integer ip_SPD, l_SPD
      Integer k, l
      Integer nij, nkl
      Integer nErr_CS, nAcc_CS
      Integer ip_D
      Integer ipInt
      Integer ip_SewWrk, l_SewWrk

      Logical WriteFirstLine

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer SPB
      Integer SPD
      Integer iTri
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      SPB(i,j)=iWork(ip_SPB-1+nShell_iAtom*(j-1)+i)
      SPD(i,j)=iWork(ip_SPD-1+nShell_iAtom*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      Call CWTime(tC0,tW0)
      WriteFirstLine=iPrint.ge.ThrPrint

      ! Init return code
      irc=0

      ! Get number of auxiliary functions
      M=LDF_nBasAux_Pair(iAtomPair)
      Mm=max(M,1)

      ! Get atoms
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Get number of basis functions on each atom
      nBas_iAtom=LDF_nBas_Atom(iAtom)
      nBas_jAtom=LDF_nBas_Atom(jAtom)
      nB2=nBas_iAtom*nBas_jAtom

      If (l_C.lt.nB2*M) Then
         Call WarningMessage(2,SecNam//': l_C<nB2*M')
         Call LDF_Quit(1)
      End If

      ! Exit with proper return code if nothing to do
      If (M.lt.1) Then
         If (nB2.lt.1) Then
            irc=0
            Return
         End If
      Else
         If (nB2.lt.1) Then
            irc=1
            Return
         End If
      End If

      ! Get number of shells on each atom
      nShell_iAtom=LDF_nShell_Atom(iAtom)
      nShell_jAtom=LDF_nShell_Atom(jAtom)

      ! Get pointer to diagonal for this atom pair
      ip_D=iWork(ip_AP_Diag-1+iAtomPair)-1

      ! Allocate and compute G matrix
      ! G[JK]=(J|K)
      l_G=Mm**2
      Call GetMem('ChkG','Allo','Real',ip_G,l_G)
      Call LDF_SetIndxG(iAtomPair)
      Call LDF_ComputeGMat(iAtomPair,M,Work(ip_G))
      Call LDF_UnsetIndxG()

      ! Get pointers to atomic shell lists
      ipi=LDF_lShell_Atom(iAtom)-1
      ipj=LDF_lShell_Atom(jAtom)-1

      ! Find largest shell pair
      maxSP=0
      Do jS=1,nShell_jAtom
         jShell=iWork(ipj+jS)
         Do iS=1,nShell_iAtom
            iShell=iWork(ipi+iS)
            maxSP=max(maxSP,nBasSh(iShell)*nBasSh(jShell))
         End Do
      End Do

      ! Allocate memory for integrals corresponding to largest shell
      ! pair
      l_Int=maxSP**2
      Call GetMem('ChkInt','Allo','Real',ip_Int,l_Int)

      ! Allocate tmp arrays to store shell pair blocks of C
      l_C_=maxSP*Mm
      Call GetMem('ChkC','Allo','Real',ip_C_,l_C_)

      ! Allocate V array to store intermediate
      ! V[J,kl]=sum_K (J|K)*C[kl,K]
      l_V=Mm*maxSP
      Call GetMem('ChkV','Allo','Real',ip_V,l_V)

      ! Allocate and set up offset array to shell pair blocks in C
      l_SPB=nShell_iAtom*nShell_jAtom
      Call GetMem('ChkSPB','Allo','Inte',ip_SPB,l_SPB)
      n=0
      Do jS=1,nShell_jAtom
         jShell=iWork(ipj+jS)
         Do iS=1,nShell_iAtom
            iWork(ip_SPB-1+nShell_iAtom*(jS-1)+iS)=n
            iShell=iWork(ipi+iS)
            n=n+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do

      ! Allocate and set up offset array to shell pair blocks in
      ! diagonal
      l_SPD=nShell_iAtom*nShell_jAtom
      Call GetMem('ChkSPD','Allo','Inte',ip_SPD,l_SPD)
      n=0
      If (iAtom.eq.jAtom) Then
         Do iS=1,nShell_iAtom
            iShell=iWork(ipi+iS)
            Do jS=1,iS-1
               iWork(ip_SPD-1+nShell_iAtom*(jS-1)+iS)=n
               iWork(ip_SPD-1+nShell_iAtom*(iS-1)+jS)=n
               jShell=iWork(ipi+jS)
               n=n+nBasSh(iShell)*nBasSh(jShell)
            End Do
            iWork(ip_SPD-1+nShell_iAtom*(iS-1)+iS)=n
            n=n+nBasSh(iShell)*(nBasSh(iShell)+1)/2
         End Do
      Else If (iAtom.gt.jAtom) Then
         Call iCopy(nShell_iAtom*nShell_jAtom,iWork(ip_SPB),1,
     &                                        iWork(ip_SPD),1)
      Else
         Call WarningMessage(2,SecNam//': iAtom<jAtom')
         Call LDF_Quit(1)
      End If

      ! Allocate memory for seward
      Call GetMem('MaxMem','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compare integrals one shell quadruple at a time
      nErr_CS=0
      nAcc_CS=0
      Err_max=-9.9d9
      Err_min=9.9d9
      Err_amx=-9.9d9
      Err_amn=9.9d9
      Err_ave=0.0d0
      Err_aav=0.0d0
      Dmax=-9.9d9
      CSDevMax=0.0d0
      n=0
      Do lS=1,nShell_jAtom
         lShell=iWork(ipj+lS)
         If (iAtom.eq.jAtom) Then
            kS1=lS
         Else
            kS1=1
         End If
         Do kS=kS1,nShell_iAtom
            kShell=iWork(ipi+kS)
            klS=iTri(kShell,lShell)
            nkl=nBasSh(kShell)*nBasSh(lShell)
            Do K=0,M-1
               Call dCopy_(nkl,C(nB2*K+SPB(kS,lS)+1),1,
     &                        Work(ip_C_+nkl*K),1)
            End Do
            Call dGeMM_('N','T',M,nkl,M,
     &                  1.0d0,Work(ip_G),Mm,Work(ip_C_),nkl,
     &                  0.0d0,Work(ip_V),Mm)
            Do jS=1,nShell_jAtom
               jShell=iWork(ipj+jS)
               If (iAtom.eq.jAtom) Then
                  iS1=jS
               Else
                  iS1=1
               End If
               Do iS=iS1,nShell_iAtom
                  iShell=iWork(ipi+iS)
                  ijS=iTri(iShell,kShell)
                  If (ijS.ge.klS) Then
                     nij=nBasSh(iShell)*nBasSh(jShell)
                     Do K=0,M-1
                        Call dCopy_(nij,C(nB2*K+SPB(iS,jS)+1),1,
     &                                 Work(ip_C_+nij*K),1)
                     End Do
                     SHA=iShell
                     SHB=jShell
                     SHC=kShell
                     SHD=lShell
                     Call Cho_dZero(Work(ip_Int),nij*nkl)
                     Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                              Work(ip_Int),l_Int,
     &                              Int_LDF_SQ)
                     Call dGeMM_('N','N',nij,nkl,M,
     &                           -1.0d0,Work(ip_C_),nij,Work(ip_V),Mm,
     &                           1.0d0,Work(ip_Int),nij)
                     ipInt=ip_Int-1
                     Do l=1,nBasSh(lShell)
                        Do k=1,nBasSh(kShell)
                           If (iAtom.eq.jAtom) Then
                              Dkl=Work(ip_D+SPD(kS,lS)+iTri(k,l))
                           Else
                              Dkl=Work(ip_D+SPD(kS,lS)
     &                                 +nBasSh(kShell)*(l-1)+k)
                           End If
                           Dmax=max(Dmax,Dkl)
                           Do j=1,nBasSh(jShell)
                              Do i=1,nBasSh(iShell)
                                 If (iAtom.eq.jAtom) Then
                                    Dij=Work(ip_D+SPD(iS,jS)+iTri(i,j))
                                 Else
                                    Dij=Work(ip_D+SPD(iS,jS)
     &                                       +nBasSh(iShell)*(j-1)+i)
                                 End If
                                 Dijkl=max(sqrt(Dij*Dkl),1.0d-14)
                                 ipInt=ipInt+1
                                 CSDev=abs(Work(ipInt))-Dijkl
                                 CSDevMax=max(CSDevMax,CSDev)
                                 If (abs(Work(ipInt)).gt.Dijkl) Then
                                    If (WriteFirstLine) Then
                                       Write(6,'(A,A)')
     &                                 'Cauchy-Schwarz check [shells,',
     &                                 'indices,error,CS upper bound]:'
                                       WriteFirstLine=.False.
                                    End If
                                    nErr_CS=nErr_CS+1
                                    If (abs(CSDev).gt.Small) Then
                                       If (iPrint.ge.ThrPrint) Then
                                          Write(6,'(8I8,1P,2D21.12)')
     &                                    iS,jS,kS,lS,i,j,k,l,
     &                                    Work(ipInt),Dijkl
                                       End If
                                    Else
                                       nAcc_CS=nAcc_CS+1
                                       If (iPrint.ge.ThrPrint) Then
                                          Write(6,'(8I8,1P,2D21.12,A)')
     &                                    iS,jS,kS,lS,i,j,k,l,
     &                                   Work(ipInt),Dijkl,' (accepted)'
                                       End If
                                    End If
                                 End If
                                 Err_max=max(Err_max,Work(ipInt))
                                 Err_min=min(Err_min,Work(ipInt))
                                 Err_amx=max(Err_amx,abs(Work(ipInt)))
                                 Err_amn=min(Err_amn,abs(Work(ipInt)))
                                 Err_ave=Err_ave+Work(ipInt)
                                 Err_aav=Err_aav+abs(Work(ipInt))
                              End Do
                           End Do
                        End Do
                     End Do
                     n=n+nij*nkl
                  End If
               End Do
            End Do
         End Do
      End Do
      If (n.eq.0) Then
         Call WarningMessage(2,SecNam//': n=0 (division by zero)')
         Call LDF_Quit(1)
      End If
      Err_ave=Err_ave/dble(n)
      Err_aav=Err_aav/dble(n)
      Write(6,'(A,1P,2D15.6)')
     & 'Max and Min Error......................',Err_max,Err_min
      Write(6,'(A,1P,2D15.6)')
     & 'Max and Min Abs Error..................',Err_amx,Err_amn
      Write(6,'(A,1P,2D15.6)')
     & 'Average and Average Abs Error..........',Err_ave,Err_aav
      Write(6,'(A,1P,D15.6)')
     & 'Max Cauchy-Schwarz deviation...........',CSDevMax
      Write(6,'(A,I9)')
     & 'Number of integrals checked............',n
      Write(6,'(A,I9)')
     & 'Number of Cauchy-Schwarz violations....',nErr_CS
      Write(6,'(A,I9)')
     & 'Number of acceptable C-S violations....',nAcc_CS
      nErr_CS=nErr_CS-nAcc_CS ! unacceptable C-S violations
      Write(6,'(A,I9,A,1P,D15.6,A)')
     & 'Number of unacceptable C-S violations..',nErr_CS,
     & ' (deviation >',Small,')'
      If (Err_amx.gt.Dmax .and. abs(Err_amx-Dmax).gt.Small) Then
         Write(6,'(A,1P,D15.6,A,A,D15.6,A)')
     &   'Max abs error deviates from max diagonal (',Dmax,') by',
     &   ' more than ',Small,' => not accepted!'
         irc=1
      Else
         If (nErr_CS.gt.0) Then
            Write(6,'(A,1P,D15.6,A,A,D15.6,A)')
     &      'Max abs error deviates from max diagonal (',Dmax,') by',
     &      ' less than ',Small,' => accepted!'
         End If
         irc=0
      End If
      If (LDF2) Then
         If (Err_amx.gt.Thr_Accuracy .and.
     &       abs(Err_amx-Thr_Accuracy).gt.Small) Then
            Write(6,'(A,1P,D15.6,A,A,D15.6,A)')
     &      'Max abs error deviates from target accuracy (',
     &      Thr_Accuracy,') by',
     &      ' more than ',Small,' => not accepted!'
            irc=1
         End If
      End If
      Call xFlush(6)

      ! Deallocate seward memory
      Call xRlsMem_Ints()

      ! Deallocate offset array to shell pair blocks in
      ! diagonal
      Call GetMem('ChkSPD','Free','Inte',ip_SPD,l_SPD)

      ! Deallocate offset array to shell pair blocks in C
      Call GetMem('ChkSPB','Free','Inte',ip_SPB,l_SPB)

      ! Deallocate V array
      Call GetMem('ChkV','Free','Real',ip_V,l_V)

      ! Deallocate C array
      Call GetMem('ChkC','Free','Real',ip_C_,l_C_)

      ! Deallocate integral array
      Call GetMem('ChkInt','Free','Real',ip_Int,l_Int)

      ! Deallocate G matrix
      Call GetMem('ChkG','Free','Real',ip_G,l_G)

      ! Print timing
      Call CWTime(tC1,tW1)
      Write(6,'(A,2(1X,F12.2),A)')
     & 'CPU and wall time for integral check:',
     & tC1-tC0,tW1-tW0,' seconds'

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CheckPairIntegrals_HlfNHlf(AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, January 2012.
C     - based on LDF_CheckPairIntegrals_Robust
C
C     Purpose: check fitting coefficients by comparing exact and
C     approximate integrals for atom pair iAtomPair.
C
C     Checks:
C     1) Symmetry: irc=1 is returned if the diff ints are not symmetric.
C     2) Diagonal: irc=2 is returned if the diagonal diff ints are not
C                  equal to the updated diagonal.
C     3) Accuracy: if 2C functions are included and the max diagonal is
C                  greater than the target accuracy, irc=3 is returned.
C                  Exception: if this is constrained LDF, irc=0 is
C                  returned and a warning is printed.
C     4) PSD: irc=4 is returned if the diff ints do not constitute a
C             positive semidefinite matrix (which implies that the
C             Cauchy-Schwarz inequality is fulfilled).
C
C     Note: this is a debug code and has not been written with
C     efficiency in mind!
C
C     Note: this routine uses the half-and-half integral representation.
C
      Implicit None
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "localdf.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*30 SecNam
      Parameter (SecNam='LDF_CheckPairIntegrals_HlfNHlf')

      Logical  isSymmetric
      External isSymmetric

      Integer  LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      External LDF_AtomPair_DiagDim, LDF_nBasAux_Pair

      Real*8 Tol
      Parameter (Tol=1.0d-12)

      Integer nAB
      Integer M, Mm
      Integer ip_Int, l_Int
      Integer ip_3IndxInt, l_3IndxInt
      Integer ip_D
      Integer kl
      Integer n
      Integer ip_Vec, l_Vec
      Integer nAcc

      Real*8 D_max

      ! Init return code
      irc=0

      ! Get and check dimensions
      nAB=LDF_AtomPair_DiagDim(AB)
      M=LDF_nBasAux_Pair(AB)
      Mm=max(M,1)
      If (l_C.lt.(nAB*M)) Then
         Call WarningMessage(2,SecNam//': insufficient array dimension')
         Call LDF_Quit(1)
      End If

      ! Exit with proper return code if nothing to do
      If (M.lt.1) Then
         If (nAB.lt.1) Then
            irc=0
            Return
         End If
      Else
         If (nAB.lt.1) Then
            irc=-1
            Return
         End If
      End If

      ! Allocate and compute exact integrals (uv|kl)
      l_Int=nAB**2
      Call GetMem('CPII','Allo','Real',ip_Int,l_Int)
      Call LDF_ComputeValenceIntegrals(AB,AB,l_Int,Work(ip_Int))

      ! Check symmetry
      If (.not.isSymmetric(Work(ip_Int),nAB,Tol)) Then
         Call WarningMessage(2,
     &        SecNam//': (AB|AB) integrals not symmetric')
         Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
         Call LDF_Quit(1)
      End If

      ! Set indices for G matrix calculation
      Call LDF_SetIndxG(AB)

      ! Allocate and compute 3 index integrals (uv|J)
      l_3IndxInt=nAB*Mm
      Call GetMem('CPI3I','Allo','Real',ip_3IndxInt,l_3IndxInt)
      Call LDF_ComputeIntegrals_uvJ(AB,l_3IndxInt,Work(ip_3IndxInt))

      ! (uv|kl) := -1/2*sum_K (uv|K)*C(kl,K) + (uv|kl)
      Call dGeMM_('N','T',nAB,nAB,M,
     &            -0.5d0,Work(ip_3IndxInt),nAB,C,nAB,
     &            1.0d0,Work(ip_Int),nAB)

      ! (uv|kl) := -1/2*sum_J C(uv,J)*(kl|J) + (uv|kl)
      Call dGeMM_('N','T',nAB,nAB,M,
     &            -0.5d0,C,nAB,Work(ip_3IndxInt),nAB,
     &            1.0d0,Work(ip_Int),nAB)

      ! Deallocation
      Call GetMem('CPI3I','Free','Real',ip_3IndxInt,l_3IndxInt)
      Call LDF_UnsetIndxG()

      ! Check 1) Symmetry
      If (irc.eq.0) Then
         If (.not.isSymmetric(Work(ip_Int),nAB,Tol)) Then
            Call WarningMessage(2,
     &        SecNam//': (Delta(AB)|Delta(AB)) integrals not symmetric')
            Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
            irc=1
         End If
      End If

      ! Check 2) Diagonal
      If (irc.eq.0) Then
         ip_D=iWork(ip_AP_Diag-1+AB)
         n=nAB+1
         kl=0
         Do While (kl.lt.nAB .and. irc.eq.0)
            If (abs(Work(ip_D+kl)-Work(ip_Int+n*kl)).gt.Tol) Then
               Call WarningMessage(2,
     &          SecNam//': (Delta(AB)|Delta(AB)) diagonal inconsistent')
               Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
               irc=2
            End If
            kl=kl+1
         End Do
      End If

      ! Check 3) Accuracy
      If (irc.eq.0) Then
         If (LDF2) Then
            D_max=0.0d0
            nAcc=0
            n=nAB+1
            Do kl=0,nAB-1
               If (Work(ip_Int+n*kl).gt.Thr_Accuracy) Then
                  nAcc=nAcc+1
                  D_max=max(D_max,Work(ip_Int+n*kl))
               End If
            End Do
            If (nAcc.ne.0) Then
               Call WarningMessage(2,
     &                   SecNam//': error greater than target accuracy')
               Write(6,'(A,1P,D20.10)') 'Max diagonal:',D_max
               If (LDF_Constraint.eq.-1) Then ! unconstrained LDF
                  irc=3
               End If
            End If
         End If
      End If

      ! Check 4) PSD
      If (irc.eq.0) Then
         l_Vec=nAB**2
         Call GetMem('CPIV','Allo','Real',ip_Vec,l_Vec)
         Call CD_InCore(Work(ip_Int),nAB,Work(ip_Vec),nAB,n,Tol,irc)
         Call GetMem('CPIV','Free','Real',ip_Vec,l_Vec)
         If (irc.ne.0) Then
            Call WarningMessage(2,
     &              SecNam//': (Delta(AB)|Delta(AB)) integrals not PSD')
            Write(6,'(A,1P,D20.10)') 'Tolerance=',Tol
            irc=4
         End If
      End If

      ! Deallocation
      Call GetMem('CPII','Free','Real',ip_Int,l_Int)

      End
