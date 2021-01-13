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
      Subroutine LDF_AddChargeConstraintCorrection(AB,l_C,C)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: adadd the correction from the charge constraint,
C
C              C(uv,J) += lambda(uv)*C1(J)
C
C              where uv and J belong to atom pair AB,
C
C              lambda(uv) = eta * [ S(uv) - sum_J C(uv,J)*n(J) ]
C
C              (here C(uv,J) are the coefficients on input)
C
C              eta = 1/[sum_J C1(J)*n(J)]
C
C              n(J) = \int J(r)dr
C
C              and C1(J) is the solution of
C
C              sum_K C1(K)*(K|J) = n(J)
C
C              The n(J) vector must be computed before entering this
C              routine (available through
C              ldf_charge_constraint_info.fh). The overlap matrix S(uv)
C              is computed directly in this routine.
C
C              Overlap and lambda matrices are stored in work at
C              locations defined in ldf_charge_constraint_info.fh
C              and are available on exit from this routine.
C
C     Implementation outline:
C
C     1. Compute G matrix, G(JK)=(J|K)
C     2. Compute overlap matrix for atom pair AB
C     3. Set up rhs of C1 equation, i.e. n(J) without linearly dependent
C        1C functions and including 2C functions (if any)
C     4. Solve equations for C1 using Lapack routine dPoSV
C     5. compute eta
C     6. compute lambda(uv)
C     7. update coefficients C+=lambda*C1
C     8. clean C (i.e. make sure that C(uv,J)=delta(uv,J) for J
C        corresponding to J)
C
      Implicit None
      Integer AB
      Integer l_C
      Real*8  C(l_C)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_charge_constraint_info.fh"

      Character*33 SecNam
      Parameter (SecNam='LDF_AddChargeConstraintCorrection')

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      real*8 ddot_
      external ddot_

      Character*8 Label
      Character*1 UpLo

      Integer A, B
      Integer nA, nB, nAB, M
      Integer ip_G, l_G
      Integer ip_n, l_n
      Integer ip_C1, l_C1
      Integer nRHS
      Integer irc
      Integer ip_S
      Integer ip_lambda

      Real*8  eta

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Check that charge constraint info is set
      If (.not.ChargeConstraintSet) Then
         Call WarningMessage(2,
     &                       SecNam//': charge constraint info not set')
         Call LDF_Quit(1)
      End If

      ! Set dimensions
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nA=LDF_nBas_Atom(A)
      nB=LDF_nBas_Atom(B)
      nAB=nA*nB
      M=LDF_nBasAux_Pair(AB)
      If (nAB.lt.1 .or. M.lt.1) Then
         ip_lambda=0
         Return
      End If
      If ((nAB*M).gt.l_C) Then
         Call WarningMessage(2,SecNam//': array dimension error')
         Call LDF_Quit(1)
      End If

      ! Allocate and compute G matrix
      l_G=M**2
      Call GetMem('CLDFG','Allo','Real',ip_G,l_G)
      Call LDF_SetIndxG(AB)
      Call LDF_ComputeGMat(AB,M,Work(ip_G))
      Call LDF_UnsetIndxG()

      ! Compute overlap matrix S(uv)
      ip_S=ip_CC_Overlap
      Label='Mltpl  0'
      Call LDF_SetOneEl(Label)
      Call LDF_ComputeOverlapBlock(AB,nAB,Work(ip_S))
      Call LDF_UnsetOneEl(Label)

      ! Allocate and set up n vector
      l_n=M
      Call GetMem('CLDFn','Allo','Real',ip_n,l_n)
      Call LDF_CC_GetRHS(AB,nAB,Work(ip_S),M,Work(ip_n))

      ! RHS vector (n) is overwritten by solution in dPOSV, so
      ! allocate C1 and copy RHS to this location
      l_C1=M
      Call GetMem('CLDFC1','Allo','Real',ip_C1,l_C1)
      Call dCopy_(M,Work(ip_n),1,Work(ip_C1),1)

      ! Solve positive definite linear equations to obtain C1:
      ! sum_K (J|K)*C1(K) = n(J)
      UpLo='L'
      nRHS=1
      irc=0
      call dposv_(UpLo,M,nRHS,Work(ip_G),M,Work(ip_C1),M,irc)
      If (irc.ne.0) Then
         Call WarningMessage(2,
     &                      SecNam//': non-zero return code from dPOSV')
         Write(6,'(A,I10)') 'Return code:',irc
         If (irc.gt.0) Then
            Write(6,'(A)')
     &      '   => G matrix not positive definite'
         Else
            Write(6,'(A,I2,A)')
     &      '   => argument no.',-irc,' has an illegal value'
         End If
         Call LDF_Quit(1)
      End If

      ! Compute eta = 1/sum_J C1(J)*n(J)
      eta=dDot_(M,Work(ip_C1),1,Work(ip_n),1)
      If (abs(eta).lt.1.0d-14) Then
         Call WarningMessage(2,SecNam//': division by zero (eta)')
         Call LDF_Quit(1)
      End If
      eta=1.0d0/eta

      ! Compute lambda(uv)=eta*(S(uv)-sum_K C(uv,K)*n(K))
      ip_lambda=ip_CC_lambda
      Call dCopy_(nAB,Work(ip_S),1,Work(ip_lambda),1)
      Call dGeMV_('N',nAB,M,
     &           -eta,C,nAB,Work(ip_n),1,
     &           eta,Work(ip_lambda),1)

      ! Explicitly zero lambda elements corresponding to 2C functions
      ! Should be zero automatically, but this removes small numerical
      ! inaccuracies.
      Call LDF_CleanLambda(AB,nAB,Work(ip_lambda))

      ! Add charge constraint to C
      ! C(uv,J) = C(uv,J) + lambda(uv)*C1(J)
      Call dGeR(nAB,M,1.0d0,Work(ip_lambda),1,Work(ip_C1),1,C,nAB)

      ! Clean rows corresponding to 2C functions, i.e.
      ! C(uv,J)=delta(uv,J) where J is a 2C function
      Call LDF_CleanC(AB,C,nAB,M)

      ! Deallocation
      Call GetMem('CLDFC1','Free','Real',ip_C1,l_C1)
      Call GetMem('CLDFn','Free','Real',ip_n,l_n)
      Call GetMem('CLDFG','Free','Real',ip_G,l_G)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CleanLambda(AB,nAB,lambda)
      Implicit None
      Integer AB
      Integer nAB
      Real*8  lambda(nAB)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*15 SecNam
      Parameter (SecNam='LDF_CleanLambda')

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer A, B
      Integer nRow_Map, nCol_Map
      Integer ip_Map, l_Map
      Integer i2CF

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer Map
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      Map(i,j)=iWork(ip_Map-1+nRow_Map*(j-1)+i)

      If (AP_2CFunctions(1,AB).gt.0) Then
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         If (nAB.lt.LDF_nBas_Atom(A)*LDF_nBas_Atom(B)) Then
            Call WarningMessage(2,
     &                         SecNam//': insufficient array dimension')
            Call LDF_Quit(1)
         End If
         nRow_Map=AP_2CFunctions(1,AB)
         If (A.eq.B) Then
            nCol_Map=2
         Else
            nCol_Map=1
         End If
         l_Map=nRow_Map*nCol_Map
         Call GetMem('CLDFMap','Allo','Inte',ip_Map,l_Map)
         Call LDF_Map2CF(AB,nRow_Map,nCol_Map,iWork(ip_Map))
         Do i2CF=1,AP_2CFunctions(1,AB)
            lambda(Map(i2CF,1))=0.0d0
         End Do
         If (A.eq.B) Then
            Do i2CF=1,AP_2CFunctions(1,AB)
               lambda(Map(i2CF,2))=0.0d0
            End Do
         End If
         Call GetMem('CLDFMap','Free','Inte',ip_Map,l_Map)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CC_GetRHS(AB,l_S,S,M,n)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: get rhs vector for the charge constraint equation,
C                 n(J) = \int J(r)dr
C              excluding linearly dependent functions. For 2C functions,
C              the integral is copied from the overlap matrix S.
C              The linearly independent 1C functions are copied from the
C              blocked aux bas vector stored in
C              ldf_charge_constraint_info.fh (which must be set up
C              before calling this routine).
C
      Implicit None
      Integer AB
      Integer l_S
      Real*8  S(l_S)
      Integer M
      Real*8  n(M)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_charge_constraint_info.fh"

      Character*13 SecNam
      Parameter (SecNam='LDF_CC_GetRHS')

      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom, LDF_nBasAux_Atom
      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom, LDF_nBasAux_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom
      Logical  LDF_isLinDep
      External LDF_isLinDep

#if defined (_DEBUGPRINT_)
      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair
#endif

      Integer A, B
      Integer ip, ipn
      Integer nSA, nSB
      Integer ipA, ipB
      Integer iS, iShell, ii
      Integer MA, MB
      Integer ip_iOff, l_iOff
      Integer i2C
      Integer iSA, iSB, iA, iB, iShellA

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      Integer List2C
      Integer iOff
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      List2C(i,j)=iWork(AP_2CFunctions(2,AB)-1+4*(j-1)+i)
      iOff(i,j)=iWork(ip_iOff-1+nSA*(j-1)+i)

      ! Check that charge constraint info is set
      If (.not.ChargeConstraintSet) Then
         Call WarningMessage(2,
     &                       SecNam//': charge constraint info not set')
         Call LDF_Quit(1)
      End If

      ! Get atoms
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)

      ! Copy 1C part
      If (AP_1CLinDep(1,AB).lt.1) Then
         ip=iWork(ip_CC_AuxIntVec-1+A)
         ipn=1
         MA=LDF_nBasAux_Atom(A)
         Call dCopy_(MA,Work(ip),1,n(ipn),1)
         ipn=ipn+MA
         If (B.ne.A) Then
            ip=iWork(ip_CC_AuxIntVec-1+B)
            MB=LDF_nBasAux_Atom(B)
            Call dCopy_(MB,Work(ip),1,n(ipn),1)
            ipn=ipn+MB
         End If
      Else
         ip=iWork(ip_CC_AuxIntVec-1+A)
         ipn=1
         nSA=LDF_nAuxShell_Atom(A)
         ipA=LDF_lAuxShell_Atom(A)-1
         Do iS=1,nSA
            iShell=iWork(ipA+iS)
            Do ii=1,nBasSh(iShell)
               If (.not.LDF_isLinDep(ii,iS,A,AB)) Then
                  n(ipn)=Work(ip)
                  ipn=ipn+1
               End If
               ip=ip+1
            End Do
         End Do
         If (B.ne.A) Then
            ip=iWork(ip_CC_AuxIntVec-1+B)
            nSB=LDF_nAuxShell_Atom(B)
            ipB=LDF_lAuxShell_Atom(B)-1
            Do iS=1,nSB
               iShell=iWork(ipB+iS)
               Do ii=1,nBasSh(iShell)
                  If (.not.LDF_isLinDep(ii,iS,B,AB)) Then
                     n(ipn)=Work(ip)
                     ipn=ipn+1
                  End If
                  ip=ip+1
               End Do
            End Do
         End If
      End If

      ! Copy 2C part
      If (AP_2CFunctions(1,AB).gt.0) Then
         nSA=LDF_nShell_Atom(A)
         nSB=LDF_nShell_Atom(B)
         l_iOff=nSA*nSB
         Call GetMem('iOff','Allo','Inte',ip_iOff,l_iOff)
         Call LDF_uvOffset(AB,nSA,nSB,iWork(ip_iOff))
         ipA=LDF_lShell_Atom(A)-1
         Do i2C=1,AP_2CFunctions(1,AB)
            iSA=List2C(1,i2C)
            iA=List2C(2,i2C)
            iSB=List2C(3,i2C)
            iB=List2C(4,i2C)
            iShellA=iWork(ipA+iSA)
            ip=iOff(iSA,iSB)+nBasSh(iShellA)*(iB-1)+iA
            n(ipn)=S(ip)
            ipn=ipn+1
         End Do
         Call GetMem('iOff','Free','Inte',ip_iOff,l_iOff)
      End If

#if defined (_DEBUGPRINT_)
      ! Check dimension
      ipn=ipn-1
      If (ipn.ne.LDF_nBasAux_Pair(AB) .or. ipn.ne.M) Then
         Call WarningMessage(2,SecNam//': dimension error')
         Call LDF_Quit(1)
      End If
#endif

      End
