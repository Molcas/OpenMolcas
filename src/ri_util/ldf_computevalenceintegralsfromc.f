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
      Subroutine LDF_ComputeValenceIntegralsFromC(Mode,tau,
     &                                            AB,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: Compute valence integrals (u_A v_B | k_C l_D) for atom
C              pairs AB and CD using LDF coefficients and either robust,
C              nonrobust, or half-and-half integral representation.
C              Mode=1: robust, ({AB}|CD)+(AB|{CD})-({AB}|{CD})
C              Mode=2: robust, ({AB}|{CD})
C              Mode=3: half-and-half, [({AB]|CD)+(AB|{CD})]/2
C
C     Note: this routine is mainly for debugging purposes.
C     Integral prescreening info should be set up before calling this
C     routine (although it does not have to be).
C
      Implicit None
      Integer Mode
      Real*8  tau
      Integer AB
      Integer CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*32 SecNam
      Parameter (SecNam='LDF_ComputeValenceIntegralsFromC')

      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet
#if defined (_DEBUG_)
      Logical  isSymmetric
      External isSymmetric
#endif

      Integer  LDF_nBas_Atom
      Integer  LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nBas_Atom
      External LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD

      Logical IPI_set_here

      Integer A, B, C, D
      Integer nBas_A, nBas_B, nBas_C, nBas_D
      Integer nAB, nCD
      Integer MAB, MCD
      Integer ip_CAB, l_CAB
      Integer ip_CCD, l_CCD
      Integer ip, l
      Integer M
      Integer ipC
      Integer l_xInt

      Real*8 Factor

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Check integral mode
      If (Mode.lt.1 .or. Mode.gt.3) Then
         Call WarningMessage(2,SecNam//': integral mode out of bounds')
         Call LDF_Quit(1)
      End If

      ! Get atoms
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)

      ! Get number of basis functions on each atom
      nBas_A=LDF_nBas_Atom(A)
      nBas_B=LDF_nBas_Atom(B)
      nBas_C=LDF_nBas_Atom(C)
      nBas_D=LDF_nBas_Atom(D)

      ! Get product dimensions
      nAB=nBas_A*nBas_B
      nCD=nBas_C*nBas_D

      ! Quick return?
      If (nAB.lt.1 .or. nCD.lt.1) Return

      ! Check integral array dimension
      l_xInt=nAB*nCD
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! Init integral array
      Call Cho_dZero(xInt,l_xInt)

      ! Get number of auxiliary basis functions on each pair
      MAB=LDF_nBasAux_Pair_wLD(AB)
      MCD=LDF_nBasAux_Pair_wLD(CD)

      ! Quick return?
      If (Mode.eq.1 .or. Mode.eq.3) Then
         If (MAB.lt.1 .and. MCD.lt.1) Return
      Else If (Mode.eq.2) Then
         If (MAB.lt.1 .or. MCD.lt.1) Return
      Else
         Write(6,'(A,A)')
     &   SecNam,': I should never end up at this place!'
         Call LDF_Quit(1)
      End If

      ! Set prescreening info (if not already done)
      If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
         Call LDF_SetIntegralPrescreeningInfo()
         IPI_set_here=.True.
      Else
         IPI_set_here=.False.
      End If

      ! Allocate and read coefficients
      l_CAB=nAB*MAB
      Call GetMem('IntCAB','Allo','Real',ip_CAB,l_CAB)
      Call LDF_CIO_ReadC_wLD(AB,Work(ip_CAB),l_CAB)
      If (AB.eq.CD) Then
         l_CCD=l_CAB
         ip_CCD=ip_CAB
      Else
         l_CCD=nCD*MCD
         Call GetMem('IntCCD','Allo','Real',ip_CCD,l_CCD)
         Call LDF_CIO_ReadC_wLD(CD,Work(ip_CCD),l_CCD)
      End If

      ! 3-index contributions
      If (Mode.eq.1 .or. Mode.eq.3) Then
C-tbp: commented following lines:
C symmetry is broken with constrained fitting, unless both
C terms are done explicitly for AB=CD. Not clear why at the
C moment (Feb. 4, 2012)
C-tbp    If (AB.eq.CD) Then
C-tbp       Factor=2.0d0
C-tbp    Else
            Factor=1.0d0
C-tbp    End If
         ! Compute sum_K_C (u_A v_B|K_C)*C(k_C l_D, K_C)
         ipC=ip_CCD
         M=LDF_nBasAux_Atom(C)
         l=nAB*M
         Call GetMem('(AB|C)','Allo','Real',ip,l)
         Call LDF_Compute3IndexIntegrals_1(AB,C,tau,l,Work(ip))
         Call dGeMM_('N','T',nAB,nCD,M,
     &               Factor,Work(ip),nAB,Work(ipC),nCD,
     &               1.0d0,xInt,nAB)
         Call GetMem('(AB|C)','Free','Real',ip,l)
         ipC=ipC+nCD*M
         If (D.ne.C) Then
            ! Compute sum_K_D (u_A v_B|K_D)*C(k_C l_D, K_D)
            M=LDF_nBasAux_Atom(D)
            l=nAB*M
            Call GetMem('(AB|D)','Allo','Real',ip,l)
            Call LDF_Compute3IndexIntegrals_1(AB,D,tau,l,Work(ip))
            Call dGeMM_('N','T',nAB,nCD,M,
     &                  Factor,Work(ip),nAB,Work(ipC),nCD,
     &                  1.0d0,xInt,nAB)
            Call GetMem('(AB|D)','Free','Real',ip,l)
            ipC=ipC+nCD*M
         End If
         If (AP_2CFunctions(1,CD).gt.0) Then
            ! Compute sum_K_CD (u_A v_B|K_CD)*C(k_C l_D, K_CD)
            M=AP_2CFunctions(1,CD)
            l=nAB*M
            Call GetMem('(AB|[CD])','Allo','Real',ip,l)
            Call LDF_Compute3IndexIntegrals_2(AB,CD,tau,l,Work(ip))
            Call dGeMM_('N','T',nAB,nCD,M,
     &                  Factor,Work(ip),nAB,Work(ipC),nCD,
     &                  1.0d0,xInt,nAB)
            Call GetMem('(AB|[CD])','Free','Real',ip,l)
         End If
C-tbp: commented following line
C symmetry is broken with constrained fitting, unless both
C terms are done explicitly for AB=CD. Not clear why at the
C moment (Feb. 4, 2012)
C-tbp    If (AB.ne.CD) Then
            ! Compute sum_J_A C(u_A v_B, J_A)*(k_C l_D|J_A)
            ipC=ip_CAB
            M=LDF_nBasAux_Atom(A)
            l=nCD*M
            Call GetMem('(CD|A)','Allo','Real',ip,l)
            Call LDF_Compute3IndexIntegrals_1(CD,A,tau,l,Work(ip))
            Call dGeMM_('N','T',nAB,nCD,M,
     &                  1.0d0,Work(ipC),nAB,Work(ip),nCD,
     &                  1.0d0,xInt,nAB)
            Call GetMem('(CD|A)','Free','Real',ip,l)
            ipC=ipC+nAB*M
            If (B.ne.A) Then
               ! Compute sum_J_B C(u_A v_B, J_B)*(k_C l_D|J_B)
               M=LDF_nBasAux_Atom(B)
               l=nCD*M
               Call GetMem('(CD|B)','Allo','Real',ip,l)
               Call LDF_Compute3IndexIntegrals_1(CD,B,tau,l,Work(ip))
               Call dGeMM_('N','T',nAB,nCD,M,
     &                     1.0d0,Work(ipC),nAB,Work(ip),nCD,
     &                     1.0d0,xInt,nAB)
               Call GetMem('(CD|B)','Free','Real',ip,l)
               ipC=ipC+nAB*M
            End If
            If (AP_2CFunctions(1,AB).gt.0) Then
               ! Compute sum_J_AB C(u_A v_B, J_AB)*(k_C l_D|J_AB)
               M=AP_2CFunctions(1,AB)
               l=nCD*M
               Call GetMem('(CD|[AB])','Allo','Real',ip,l)
               Call LDF_Compute3IndexIntegrals_2(CD,AB,tau,l,Work(ip))
               Call dGeMM_('N','T',nAB,nCD,M,
     &                     1.0d0,Work(ipC),nAB,Work(ip),nCD,
     &                     1.0d0,xInt,nAB)
               Call GetMem('(CD|[AB])','Free','Real',ip,l)
            End If
C-tbp    End If
         If (Mode.eq.3) Call dScal_(l_xInt,0.5d0,xInt,1)
      End If

#if defined (_DEBUG_)
      If (AB.eq.CD) Then
         If (.not.isSymmetric(xInt,nAB,1.0d-12)) Then
            Call WarningMessage(0,
     &                          SecNam//': integrals not symmetric [1]')
            Write(6,'(A,I10)') 'AB=',AB
         End If
      End If
#endif

      ! 2-center contributions
      If (Mode.eq.1 .or. Mode.eq.2) Then
         ! Set factor for 2-center contributions
         If (Mode.eq.1) Then
            Factor=-1.0d0
         Else
            Factor=1.0d0
         End If
         l=max(LDF_nBasAux_Atom(A),
     &         LDF_nBasAux_Atom(B),
     &         AP_2CFunctions(1,AB))*nCD
         Call GetMem('Intermediate','Allo','Real',ip,l)
         ipC=ip_CAB
         M=LDF_nBasAux_Atom(A)
         Call LDF_CVIFC_1(A,CD,tau,l_CCD,Work(ip_CCD),l,Work(ip))
         Call dGeMM_('N','N',nAB,nCD,M,
     &               Factor,Work(ipC),nAB,Work(ip),max(M,1),
     &               1.0d0,xInt,nAB)
         ipC=ipC+nAB*M
         If (B.ne.A) Then
            M=LDF_nBasAux_Atom(B)
            Call LDF_CVIFC_1(B,CD,tau,l_CCD,Work(ip_CCD),l,Work(ip))
            Call dGeMM_('N','N',nAB,nCD,M,
     &                  Factor,Work(ipC),nAB,Work(ip),max(M,1),
     &                  1.0d0,xInt,nAB)
            ipC=ipC+nAB*M
         End If
         If (AP_2CFunctions(1,AB).gt.0) Then
            M=AP_2CFunctions(1,AB)
            Call LDF_CVIFC_2(AB,CD,tau,l_CCD,Work(ip_CCD),l,Work(ip))
            Call dGeMM_('N','N',nAB,nCD,M,
     &                  Factor,Work(ipC),nAB,Work(ip),max(M,1),
     &                  1.0d0,xInt,nAB)
         End If
         Call GetMem('Intermediate','Free','Real',ip,l)
      End If

#if defined (_DEBUG_)
      If (AB.eq.CD) Then
         If (.not.isSymmetric(xInt,nAB,1.0d-12)) Then
            Call WarningMessage(0,
     &                          SecNam//': integrals not symmetric [2]')
            Write(6,'(A,I10)') 'AB=',AB
         End If
      End If
#endif

      ! If integral prescreening info was set here, unset info.
      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      ! Deallocation
      If (CD.ne.AB) Then
         Call GetMem('IntCCD','Free','Real',ip_CCD,l_CCD)
      End If
      Call GetMem('IntCAB','Free','Real',ip_CAB,l_CAB)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CVIFC_1(A,CD,tau,l_CCD,CCD,l_X_,X)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute X intermediate.
C
C        X(J_A,k_C l_D) = sum_K_C  (J_A | K_C)* CCD(k_C l_D,K_C)
C                       + sum_K_D  (J_A | K_D)* CCD(k_C l_D,K_D)
C                       + sum_K_CD (J_A | K_CD)*CCD(k_C l_D,K_CD)
C
      Implicit None
      Integer A
      Integer CD
      Real*8  tau
      Integer l_CCD
      Real*8  CCD(l_CCD)
      Integer l_X_
      Real*8  X(l_X_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*11 SecNam
      Parameter (SecNam='LDF_CVIFC_1')

      Integer  LDF_nBas_Atom
      Integer  LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nBas_Atom
      External LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD

      Integer C, D
      Integer nCD
      Integer MA, MC, MD, MCD
      Integer ip, l
      Integer ipC
      Integer l_X

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      MA=LDF_nBasAux_Atom(A)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)
      nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      MC=LDF_nBasAux_Atom(C)
      MD=LDF_nBasAux_Atom(D)
      MCD=AP_2CFunctions(1,CD)

      If (MA.lt.1) Return
      If (nCD.lt.1) Return

      If (C.eq.D) Then
         If (l_CCD.lt.nCD*(MC+MCD)) Then
            Call WarningMessage(2,SecNam//': Illegal CCD dimension')
            Call LDF_Quit(1)
         End If
         If (LDF_nBasAux_Pair_wLD(CD).ne.(MC+MCD)) Then
            Call WarningMessage(2,
     &                          SecNam//': Pair auxbas dimension error')
            Call LDF_Quit(1)
         End If
      Else
         If (l_CCD.lt.nCD*(MC+MD+MCD)) Then
            Call WarningMessage(2,SecNam//': Illegal CCD dimension')
            Call LDF_Quit(1)
         End If
         If (LDF_nBasAux_Pair_wLD(CD).ne.(MC+MD+MCD)) Then
            Call WarningMessage(2,
     &                          SecNam//': Pair auxbas dimension error')
            Call LDF_Quit(1)
         End If
      End If

      l_X=MA*nCD
      If (l_X.gt.l_X_) Then
         Call WarningMessage(2,SecNam//': Insufficient X dimension')
         Call LDF_Quit(1)
      End If
      Call Cho_dZero(X,l_X)

      If (LDF_nBasAux_Pair_wLD(CD).lt.1) Return

      l=MA*max(MC,MD,MCD)
      Call GetMem('CVIFC1','Allo','Real',ip,l)
      ipC=1
      Call LDF_Compute2IndexIntegrals_11(A,C,tau,l,Work(ip))
      Call dGeMM_('N','T',MA,nCD,MC,
     &            1.0d0,Work(ip),MA,CCD(ipC),nCD,
     &            1.0d0,X,MA)
      ipC=ipC+nCD*MC
      If (D.ne.C) Then
         Call LDF_Compute2IndexIntegrals_11(A,D,tau,l,Work(ip))
         Call dGeMM_('N','T',MA,nCD,MD,
     &               1.0d0,Work(ip),MA,CCD(ipC),nCD,
     &               1.0d0,X,MA)
         ipC=ipC+nCD*MD
      End If
      If (MCD.gt.0) Then
         Call LDF_Compute2IndexIntegrals_12(A,CD,tau,l,Work(ip))
         Call dGeMM_('N','T',MA,nCD,MCD,
     &               1.0d0,Work(ip),MA,CCD(ipC),nCD,
     &               1.0d0,X,MA)
      End If
      Call GetMem('CVIFC1','Free','Real',ip,l)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CVIFC_2(AB,CD,tau,l_CCD,CCD,l_X_,X)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute X intermediate.
C
C        X(J_AB,k_C l_D) = sum_K_C  (J_AB | K_C)* CCD(k_C l_D,K_C)
C                        + sum_K_D  (J_AB | K_D)* CCD(k_C l_D,K_D)
C                        + sum_K_CD (J_AB | K_CD)*CCD(k_C l_D,K_CD)
C
      Implicit None
      Integer AB
      Integer CD
      Real*8  tau
      Integer l_CCD
      Real*8  CCD(l_CCD)
      Integer l_X_
      Real*8  X(l_X_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*11 SecNam
      Parameter (SecNam='LDF_CVIFC_2')

      Integer  LDF_nBas_Atom
      Integer  LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nBas_Atom
      External LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD

      Integer C, D
      Integer nCD
      Integer MAB, MC, MD, MCD
      Integer ip, l
      Integer ipC
      Integer l_X

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      MAB=AP_2CFunctions(1,AB)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)
      nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      MC=LDF_nBasAux_Atom(C)
      MD=LDF_nBasAux_Atom(D)
      MCD=AP_2CFunctions(1,CD)

      If (MAB.lt.1) Return
      If (nCD.lt.1) Return

      If (C.eq.D) Then
         If (l_CCD.lt.nCD*(MC+MCD)) Then
            Call WarningMessage(2,SecNam//': Illegal CCD dimension')
            Call LDF_Quit(1)
         End If
         If (LDF_nBasAux_Pair_wLD(CD).ne.(MC+MCD)) Then
            Call WarningMessage(2,
     &                          SecNam//': Pair auxbas dimension error')
            Call LDF_Quit(1)
         End If
      Else
         If (l_CCD.lt.nCD*(MC+MD+MCD)) Then
            Call WarningMessage(2,SecNam//': Illegal CCD dimension')
            Call LDF_Quit(1)
         End If
         If (LDF_nBasAux_Pair_wLD(CD).ne.(MC+MD+MCD)) Then
            Call WarningMessage(2,
     &                          SecNam//': Pair auxbas dimension error')
            Call LDF_Quit(1)
         End If
      End If

      l_X=MAB*nCD
      If (l_X.gt.l_X_) Then
         Call WarningMessage(2,SecNam//': Insufficient X dimension')
         Call LDF_Quit(1)
      End If
      Call Cho_dZero(X,l_X)

      If (LDF_nBasAux_Pair_wLD(CD).lt.1) Return

      l=MAB*max(MC,MD,MCD)
      Call GetMem('CVIFC2','Allo','Real',ip,l)
      ipC=1
      Call LDF_Compute2IndexIntegrals_12(C,AB,tau,l,Work(ip))
      Call dGeMM_('T','T',MAB,nCD,MC,
     &            1.0d0,Work(ip),max(MC,1),CCD(ipC),nCD,
     &            1.0d0,X,MAB)
      ipC=ipC+nCD*MC
      If (D.ne.C) Then
         Call LDF_Compute2IndexIntegrals_12(D,AB,tau,l,Work(ip))
         Call dGeMM_('T','T',MAB,nCD,MD,
     &               1.0d0,Work(ip),max(MD,1),CCD(ipC),nCD,
     &               1.0d0,X,MAB)
         ipC=ipC+nCD*MD
      End If
      If (MCD.gt.0) Then
         Call LDF_Compute2IndexIntegrals_22(AB,CD,tau,l,Work(ip))
         Call dGeMM_('N','T',MAB,nCD,MCD,
     &               1.0d0,Work(ip),MAB,CCD(ipC),nCD,
     &               1.0d0,X,MAB)
      End If
      Call GetMem('CVIFC2','Free','Real',ip,l)

      End
