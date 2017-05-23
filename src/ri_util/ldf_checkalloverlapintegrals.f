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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_CheckAllOverlapIntegrals(doPrint,Tol2C,
     &                                        MAE,AB_MAE,
     &                                        MRNrm,AB_MRNrm)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compare exact and fitted overlap for all atom pairs, i.e.
C              sum_J C(uv,J)*\int J(r)dr - S(uv)
C              The routine does NOT judge whether the MAE (or another
C              measure) is too large. However, it does check that
C              overlaps corresponding to 2C auxiliary functions are
C              represented exactly (to within the input tolerance
C              Tol2C).
C
C     Input:
C        doPrint  -- if true, print info
C        Tol2C    -- error tolerance for overlaps corresponding to
C                    2C auxiliary functions
C     Output:
C        MAE      -- max abs error
C        AB_MAE   -- atom pair corresponding to MAE
C        MRNrm    -- max relative error norm
C        AB_MRNrm -- atom pair corresponding to MRNRm
C
      Implicit None
      Logical doPrint
      Real*8  Tol2C
      Real*8  MAE
      Integer AB_MAE
      Real*8  MRNrm
      Integer AB_MRNrm
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"

      Character*28 SecNam
      Parameter (SecNam='LDF_CheckAllOverlapIntegrals')

      Real*8   dDot_, Cho_dSumElm, LDF_AtomicDistance
      external ddot_, Cho_dSumElm, LDF_AtomicDistance

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair_wLD, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair_wLD, LDF_nBasAux_Pair

      Integer ip_Stat, l_Stat
      Integer ip_SBlocks
      Integer ip_AuxIntVec
      Integer AB
      Integer A, B
      Integer nAB
      Integer ip_C, l_C, l_Cmax
      Integer ipS
      Integer n2CErr, irc
      Integer AB_MAE_2CT

      Real*8  Nrm, Sm
      Real*8  SNrm, SSm
      Real*8  RMS, RNrm
      Real*8  MAE_2C, MAE_2CT

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      MAE=-9.9d9
      AB_MAE=-1
      MRNrm=-9.9d9
      AB_MRNrm=-1
      MAE_2C=0.0d0
      MAE_2CT=0.0d0
      AB_MAE_2CT=-1
      l_Cmax=0
      Do AB=1,NumberOfAtomPairs
         l_Cmax=max(l_Cmax,LDF_nBas_Atom(AP_Atoms(1,AB))
     &                    *LDF_nBas_Atom(AP_Atoms(2,AB))
     &                    *LDF_nBasAux_Pair_wLD(AB))
      End Do
      If (l_Cmax.lt.1) Return
      Call GetMem('Coeff','Allo','Real',ip_C,l_Cmax)
      Call LDF_AllocateBlockMatrix('Ovl',ip_SBlocks)
      Call LDF_GetBlockedOverlapMatrix(0,ip_SBlocks)
      Call LDF_AllocateAuxBasVector('Int',ip_AuxIntVec)
      Call LDF_ComputeAuxInt(ip_AuxIntVec)
      l_Stat=6
      Call GetMem('S_Stat','Allo','Real',ip_Stat,l_Stat)
      If (doPrint) Call Cho_Head('LDF Overlap Check','-',80,6)
      n2CErr=0
      Do AB=1,NumberOfAtomPairs
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         l_C=nAB*LDF_nBasAux_Pair_wLD(AB)
         If (l_C.gt.0) Then
            Call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l_C)
            ipS=iWork(ip_SBlocks-1+AB)
            If (doPrint) Then
               SNrm=sqrt(dDot_(nAB,Work(ipS),1,Work(ipS),1))
               SSm=Cho_dSumElm(Work(ipS),nAB)
            End If
            Call LDF_ComputeOverlapFromAuxInt(AB,
     &                                1.0d0,l_C,Work(ip_C),ip_AuxIntVec,
     &                                   -1.0d0,nAB,Work(ipS))
            If (LDF2) Then
               Call LDF_Check2COverlap(.False.,AB,nAB,Work(ipS),Tol2C,
     &                                 MAE_2C,irc)
               n2CErr=n2CErr+irc
               If (MAE_2C.gt.MAE_2CT) Then
                  MAE_2CT=MAE_2C
                  AB_MAE_2CT=AB
               End If
            End If
            If (nAB.lt.1) Then
               Call Cho_dZero(Work(ip_Stat),l_Stat)
               Nrm=0.0d0
               RMS=0.0d0
               Sm=0.0d0
            Else
               ipS=iWork(ip_SBlocks-1+AB)
               Call Statistics(Work(ipS),nAB,Work(ip_Stat),
     &                         1,2,3,4,5,6,0)
               Nrm=dDot_(nAB,Work(ipS),1,Work(ipS),1)
               RMS=sqrt(Nrm/dble(nAB))
               Nrm=sqrt(Nrm)
               Sm=Cho_dSumElm(Work(ipS),nAB)
            End If
            If (SNrm.gt.0.0d0) Then
               RNrm=Nrm/SNrm
            Else
               RNrm=0.0d0
            End If
            If (doPrint) Then
               Write(6,'(/,2X,A,10X,I10,2X,A,2I10,2X,A,1P,D20.10)')
     &         'Atom pair.........',AB,
     &         'Atoms.............',A,B,
     &         'Atomic distance...',LDF_AtomicDistance(A,B)
               Write(6,'(2X,A,5X,I15,2X,A,5X,I15,A)')
     &         'Dimension.........',nAB,
     &         'Auxiliary basis...',
     &         LDF_nBasAux_Pair(AB),' (w/o LinDep)'
               Write(6,'(1P,2(2X,A,D20.10))')
     &         'Norm of S.........',SNrm,
     &         'Sum of S..........',SSm
               Write(6,'(1P,2(2X,A,D20.10))')
     &         'Norm of error.....',Nrm,
     &         'Sum of error......',Sm
               Write(6,'(1P,2(2X,A,D20.10))')
     &         'Average error.....',Work(ip_Stat),
     &         'Standard deviation',Work(ip_Stat+5)
               Write(6,'(1P,2(2X,A,D20.10))')
     &         'Minimum error.....',Work(ip_Stat+2),
     &         'Maximum error.....',Work(ip_Stat+3)
               Write(6,'(1P,2(2X,A,D20.10))')
     &         'RMS error.........',RMS,
     &         'Abs average error.',Work(ip_Stat+1)
               Write(6,'(1P,2(2X,A,D20.10))')
     &         'Relative error nrm',RNrm,
     &         'Maximum abs error.',Work(ip_Stat+4)
               Write(6,'(1P,2X,A,D20.10)')
     &         'Max abs 2C error..',MAE_2C
               Call xFlush(6)
            End If
            If (Work(ip_Stat+4).gt.MAE) Then
               MAE=Work(ip_Stat+4)
               AB_MAE=AB
            End If
            If (RNrm.gt.MRNrm) Then
               MRNrm=RNrm
               AB_MRNrm=AB
            End If
         End If
      End Do
      Call GetMem('S_Stat','Free','Real',ip_Stat,l_Stat)
      Call LDF_DeallocateAuxBasVector('Int',ip_AuxIntVec)
      Call LDF_DeallocateBlockMatrix('Ovl',ip_SBlocks)
      Call GetMem('Coeff','Free','Real',ip_C,l_Cmax)

      If (doPrint) Then
         Write(6,*)
         If (AB_MAE.gt.0) Then
            Write(6,'(2X,A,1P,D20.10,1X,A,I10,2X,A,D20.10)')
     &      'Max abs error..........',MAE,'@AB=',AB_MAE,'Distance=',
     &      LDF_AtomicDistance(AP_Atoms(1,AB_MAE),AP_Atoms(2,AB_MAE))
         End If
         If (AB_MRNrm.gt.0) Then
            Write(6,'(2X,A,1P,D20.10,1X,A,I10,2X,A,D20.10)')
     &     'Max relative error norm',MRNrm,'@AB=',AB_MRNrm,'Distance=',
     &     LDF_AtomicDistance(AP_Atoms(1,AB_MRNrm),AP_Atoms(2,AB_MRNrm))
         End If
         If (AB_MAE_2CT.gt.0) Then
            Write(6,'(2X,A,1P,D20.10,1X,A,I10,2X,A,D20.10)')
     &      'Max abs 2C error.......',MAE_2CT,'@AB=',AB_MAE_2CT,
     &      'Distance=',
     &      LDF_AtomicDistance(AP_Atoms(1,AB_MAE_2CT),
     &                         AP_Atoms(2,AB_MAE_2CT))
         End If
         Write(6,'(2X,A,10X,I10,16X,A,1P,D20.10)')
     &   'Number of 2C errors....',n2CErr,'Tolerance=',Tol2C
         Call xFlush(6)
      End If

      If (n2CErr.gt.0) Then
         Call WarningMessage(2,SecNam//': too large 2C overlap errors')
         Write(6,'(A,1P,D20.10,2X,A,I10)')
     &   'Tolerance=',Tol2C,
     &   'Number of errors=',n2CErr
         Call LDF_Quit(1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ComputeOverlapFromAuxInt(AB,
     &                                        alpha,l_C,C,ip_AuxIntVec,
     &                                        beta,l_S,S)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute fitted overlap for atom pair AB, i.e.
C              S(uA vB) := alpha * sum_J C(uA vB,J) * \int J(r)dr
C                         + beta * S(uA vB)
C
      Implicit None
      Integer AB
      Real*8  alpha
      Integer l_C
      Real*8  C(l_C)
      Integer ip_AuxIntVec
      Real*8  beta
      Integer l_S
      Real*8  S(l_S)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

#if defined (_DEBUG_)
      Character*28 SecNam
      Parameter (SecNam='LDF_ComputeOverlapFromAuxInt')
#endif

      Integer  LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom

      Integer A, B
      Integer nAB
      Integer ipC
      Integer ipJ
      Integer M

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      If (nAB.lt.1) Return
#if defined (_DEBUG_)
      If (nAB.gt.l_S) Then
         Call WarningMessage(2,SecNam//': dimension mismatch [1]')
         Call LDF_Quit(1)
      End If
#endif

      If (beta.ne.1.0d0) Then
         If (abs(beta).lt.1.0d-24) Then
            Call fZero(S,nAB)
         Else
            Call dScal_(nAB,beta,S,1)
         End If
      End If
      ipC=1
      ipJ=iWork(ip_AuxIntVec-1+A)
      M=LDF_nBasAux_Atom(A)
      Call dGeMV_('N',nAB,M,
     &           alpha,C(ipC),nAB,Work(ipJ),1,
     &           1.0d0,S,1)
      If (B.ne.A) Then
         ipC=ipC+nAB*M
         ipJ=iWork(ip_AuxIntVec-1+B)
         M=LDF_nBasAux_Atom(B)
         Call dGeMV_('N',nAB,M,
     &              alpha,C(ipC),nAB,Work(ipJ),1,
     &              1.0d0,S,1)
      End If
      If (AP_2CFunctions(1,AB).gt.0) Then
         ipC=ipC+nAB*M
         ipJ=iWork(ip_AuxIntVec-1+LDF_nAtom()+AB)
         M=AP_2CFunctions(1,AB)
         Call dGeMV_('N',nAB,M,
     &              alpha,C(ipC),nAB,Work(ipJ),1,
     &              1.0d0,S,1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Check2COverlap(doPrint,AB,l_S,S,Tol2C,MAE,irc)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: check that overlaps corresponding to 2C fitting functions
C              are represented to within the tolerance Tol2C. On exit,
C              irc is the number of too large errors and MAE is the max
C              abs 2C error.
C
      Implicit None
      Logical doPrint
      Integer AB
      Integer l_S
      Real*8  S(l_S)
      Real*8  Tol2C
      Real*8  MAE
      Integer irc
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_info.fh"
#include "ldf_atom_pair_info.fh"

      Character*18 SecNam
      Parameter (SecNam='LDF_Check2COverlap')

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Integer A, B
      Integer nAB
      Integer nSA, nSB
      Integer ip_Offset, l_Offset
      Integer iSA, iSB
      Integer iShellA
      Integer i2C
      Integer iA, iB
      Integer ipA, ip

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer List2C
      Integer iOff
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      List2C(i,j)=iWork(AP_2CFunctions(2,AB)-1+4*(j-1)+i)
      iOff(i,j)=iWork(ip_Offset-1+nSA*(j-1)+i)

      irc=0
      MAE=0.0d0
      If (AP_2CFunctions(1,AB).lt.1) Return

      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      If (l_S.lt.nAB) Then
         Call WarningMessage(2,SecNam//': l_S < nAB')
         Call LDF_Quit(1)
      End If

      nSA=LDF_nShell_Atom(A)
      nSB=LDF_nShell_Atom(B)
      l_Offset=nSA*nSB
      Call GetMem('Offset','Allo','Inte',ip_Offset,l_Offset)
      Call LDF_uvOffset(AB,nSA,nSB,iWork(ip_Offset))
      ipA=LDF_lShell_Atom(A)-1
      Do i2C=1,AP_2CFunctions(1,AB)
         iSA=List2C(1,i2C)
         iA=List2C(2,i2C)
         iSB=List2C(3,i2C)
         iB=List2C(4,i2C)
         iShellA=iWork(ipA+iSA)
         ip=iOff(iSA,iSB)+nBasSh(iShellA)*(iB-1)+iA
         If (abs(S(ip)).gt.Tol2C) Then
            irc=irc+1
         End If
         MAE=max(MAE,abs(S(ip)))
      End Do
      Call GetMem('Offset','Free','Inte',ip_Offset,l_Offset)

      If (doPrint) Then
         Write(6,'(2X,A,I10,2X,A,1P,D20.10)')
     &   'AB=',AB,'Max abs 2C overlap error=',MAE
         Call xFlush(6)
      End If

      End
