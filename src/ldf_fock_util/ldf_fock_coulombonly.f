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
      Subroutine LDF_Fock_CoulombOnly(IntegralOption,
     &                                Timing,Mode,ThrPS,
     &                                Add,PackedD,PackedF,
     &                                nD,FactC,ip_D,ip_F)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Compute Coulomb contribution to Fock matrix using local
C              Density Fitting coefficients.
C
C     Args:    - IntegralOption (integer)
C                  111  => use two-electron integrals computed from
C                          LDF coefficients (for debugging)
C                  222  => use conventional two-electron integrals
C                          (for debugging)
C                  333  => use conventional or LDF two-electron
C                          integrals depending on positivity of the
C                          latter (for debugging)
C                  444  => use exact integral diagonal (blocks), all
C                          other blocks are LDF
C                  all other values => use production-level LDF code
C                  (INPUT)
C              - Timing (boolean) .True. if detailed timing (incl total)
C                should be printed from this routine (INPUT)
C              - Mode (integer) 1:robust, 2:nonrobust, 3:half-and-half
C                to be used (INPUT)
C              - ThrPS(2): prescreening thresholds (INPUT)
C                          ThrPS(1): integral prescreening threshold
C                          ThrPS(2): prescreening threshold on
C                                    contributions
C              - Add (boolean): .True. if Coulomb contribution should be
C                added to the Fock matrix (INPUT). Not available in
C                parallel execution!!
C              - PackedD (boolean): .True. if density matrices
C                are stored in lower triangular format; else quadratic
C                storage (INPUT)
C              - PackedF (boolean): .True. if Fock matrices
C                are stored in lower triangular format; else quadratic
C                storage (INPUT)
C              - nD: Number of densities/Fock matrices (INPUT)
C              - FactC(nD): scaling factor for each density (INPUT)
C              - ip_D(nD): pointers to nD density matrices (if PackedD:
C                lower triangular storage) (INPUT)
C              - ip_F(nD): pointers to nD Fock matrices (if PackedF:
C                lower triangular storage) (INPUT)
C
C              If (Add):  [NOT IMPLEMENTED IN PARALLEL]
C       (1)       F(uv) = F(uv) + FactC * sum_kl (uv|kl)*D(kl)
C              Else:
C       (2)       F(uv) = FactC * sum_kl (uv|kl)*D(kl)
C
C              where the integrals are given by the robust LDF
C              representation (Mode=1)
C
C              (uv|kl) = ([uv]|kl) + (uv|[kl]) - ([uv]|[kl])
C
C              or the non-robust representation (Mode=2)
C
C              (uv|kl) = ([uv]|[kl])
C
C              or the half-and-half representation (Mode=3)
C
C              (uv|kl) = 0.5*([uv]|kl) + 0.5*(uv|[kl])
C
C              The fitted products are given by
C
C              |[uv]) = sum_J C(uv,J) |J)
C
C              where the fitting functions |J) are centered on the atom
C              pair to which the product |uv) belongs.
C
C              In more detail (assuming robust fitting representation):
C
C              F(uv) = FactC * {sum_J C(uv,J)*W(J) + sum_K (uv|K)*V(K)}
C
C              W(J) = sum_kl (J|kl)*D(kl) - sum_K (J|K)*V(K)
C
C              V(K) = sum_kl C(kl,K)*D(kl)
C
C     NOTES:
C       - The density matrices are assumed to be symmetric.
C       - It is the complete sum over kl in Eqs(1,2)! This means that
C         the input density matrices should NOT be scaled by 2 on the
C         off-diagonal (scaling factors are handled internally).
C       - LDF information must be properly set up before calling this
C         routine (use subroutine LDF_X_Init).
C       - Two-center functions may or mat not be included in the fitting
C         basis for a given atom pair.
C       - The algorithm is integral driven: the integrals are computed
C         once.
C       - The fitting coefficients are read from disk twice (for Coulomb
C         [V] intermediates and for the final contractions to give the
C         Fock matrix). This is not a problem if they are fully
C         buffered, but could impose an I/O bottleneck. Since the number
C         of coefficients scales linearly with system size, this should
C         be less of a problem than the quadratic scaling of full
C         Cholesky vectors.
C       - The max number of processes for which the algorithm can
C         possibly scale well equals the number of LDF atom pairs.
C       - Integral prescreening info will be set up here if not done
C         prior to calling this routine.
C       - Add is not implemented in parallel !!!
C
C     ALGORITHM:
C=======================================================================
C
C
C     If (Mode=3): FactC := 0.5*FactC
C     Scale off-diagonal density blocks:
C                  D[u_A v_B] <-- 2*D[u_A v_B] iff A != B
C
C     Initialize V(J)=0 for all J
C     *Parallel loop over atom pairs AB (A>=B):
C        - V(J) += sum_[u_A v_B] C(u_A v_B,J)*D(u_A v_B)
C          where J belongs to AB [1C and 2C funcs].
C     *End parallel loop AB
C     Add V over nodes.
C
C     Initialize F(uv) for all significant uv.
C     Initialize W(J) for all J.
C
C     If (Mode=1 or Mode=3):
C        Compute 3-index contributions:
C        *Parallel loop over atom pairs AB (A>=B):
C           *Loop over atoms C
C              - Compute (u_A v_B|K) for K belonging to C [1C func]
C              - F(u_A v_B) += sum_K (u_A v_B|K)*V(K)
C              - W(K) += sum_[u_A v_B] (u_A v_B|K)*D(u_A v_B)
C           *End loop C
C           If (LDF2):
C              *loop over atom pairs CD
C                 - Compute (u_A v_B|K) for K belonging to CD [2C func]
C                 - F(u_A v_B) += sum_[K] (u_A v_B|K)*V(K)
C                 - W(K) += sum_[u_A v_B] (u_A v_B|K)*D(u_A v_B)
C              *End loop over CD
C           End If
C        *End parallel loop AB
C     End If
C
C     Compute 2-index contributions:
C     If (Mode=1 or Mode=2):
C        If (Mode=1):
C           - const=-1.0d0
C        Else:
C           - const=1.0d0
C        End If
C        *Parallel loop over atoms A:
C           *Loop over atoms B=1,A-1:
C              - Compute (J|K) for J in A and K in B [1C func]
C              - W(J) += const * sum_K (J|K)*V(K)
C              - W(K) += const * sum_J V(J)*(J|K)
C           *End loop B
C           - Compute (J|K) for J,K in A [1C func]
C           - W(J) += const* sum_K (J|K)*V(K)
C           If (LDF2):
C              *Loop over atom pairs CD (C>=D):
C                 - Compute (J|K) for J in A [1C func], K in CD [2Cfunc]
C                 - W(J) += const * sum_K (J | K)*V(K)
C                 - W(K) += const * sum_J V(J)*(J|K)
C              *End loop CD
C           End If
C        *End parallel loop A
C        If (LDF2):
C           *Parallel loop over atom pairs AB (A>=B):
C              *Loop over atom pairs CD=1,AB-1:
C                 - Compute (J|K) for J in AB and K in CD [2C func]
C                 - W(J) += const * sum_K (J|K)*V(K)
C                 - W(K) += const * sum_J V(J)*(J|K)
C              *End loop B
C              - Compute (J|K) for J,K in AB [2C func]
C              - W(J) += const* sum_K (J|K)*V(K)
C           *End loop AB
C        End If
C     End If
C     Add W over nodes
C
C     *Parallel loop over atom pairs AB (A>=B):
C        - Read coefficient C(u_A v_B,J) for J in AB [1C,2C func]
C        - F(u_A v_B) += sum_J C(u_a v_B,J)*W(J)
C     *End parallel loop AB
C     Add F over nodes
C
C     ALGORITHM DONE: Return F
C=======================================================================
C
C     Note that arrays V, W, and F are stored locally as O(N) arrays,
C     which should keep the communication bottleneck reasonable as the
C     system size grows.
C
C     Integral prescreening is based on the Cauchy-Schwarz inequality
C     and involves estimation of the max integral in a given block.
C     Integral prescreening is done at atom/atom pair block level as
C     well as on shell level.
C
C     Contribution prescreening is based on the submultiplicative
C     property of the Frobenius norm in combination with the Cauchy-
C     Schwarz inequality for positive (semi-) definite matrices.
C     Specifically, an upper bound to the norm of the matrix-vector
C     product
C
C     y = A*x
C
C     where A is a subblock of a positive (semi-) definite matrix,
C     is given by
C
C     ||y|| <= ||A||*||x||                  {submultiplicative property}
C            = sqrt[sum_ij A(i,j)**2]*sqrt[sum_j x(j)**2]
C           <= sqrt[sum_ij A(i,i)*A(j,j)]*sqrt[sum_j x(j)**2]
C                                                       {Cauchy-Schwarz}
C            = sqrt[sum_i A(i,i)]*sqrt[sum_j A(j,j)]*sqrt[sum_j x(j)**2]
C
C     This upper bound is used to avoid calculation of contributions
C     below a given threshold. Contribution prescreening is done at
C     atom/atom pair block level only, f.ex. for K on atom C
C
C     ||sum_K (AB|K)*V(K)|| <= sqrt[sum_[u_A v_B] (u_A v_B | u_A v_B)]
C                             *sqrt[sum_K (K|K)]
C                             *sqrt[sum_K V(K)**2]
C
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
      Integer IntegralOption
      Logical Timing
      Integer Mode
      Real*8  ThrPS(2)
      Logical Add
      Logical PackedD
      Logical PackedF
      Integer nD
      Real*8  FactC(nD)
      Integer ip_D(nD)
      Integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*20 SecNam
      Parameter (SecNam='LDF_Fock_CoulombOnly')

      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet

      Integer  LDF_nAtom
      External LDF_nAtom

#if defined (_DEBUGPRINT_)
      Logical  LDF_X_IsSet, LDF_TestBlockMatrix
      External LDF_X_IsSet, LDF_TestBlockMatrix
      Logical  DoTest
      Real*8   x, y
#endif

      Logical UseOldCode
      Logical IPI_set_here

      Real*8 tau(2)

      Integer i
      Integer nBas
      Integer iD, ip0, l
      Integer ip_DBlocks, l_DBlocks
      Integer ip_FBlocks, l_FBlocks
      Integer ip_VP, l_VP
      Integer ip_DNorm, l_DNorm
      Integer ip_CNorm, l_CNorm
      Integer ip_VNorm, l_VNorm
      Integer ip_FactC, l_FactC

      Real*8 tTotC1, tTotC2
      Real*8 tTotW1, tTotW2

      ! Start total timing.
      If (Timing) Then
         Write(6,'(/,84A1)') ('-',i=1,84)
         Call CWTime(tTotC1,tTotW1)
      End If

      UseOldCode=IntegralOption.eq.111 .or. IntegralOption.eq.222 .or.
     &           IntegralOption.eq.333
      If (UseOldCode) Then
         Call WarningMessage(0,
     &                   SecNam//': Using atom pair-driven (old) code!')
         Call xFlush(6)
         Call LDF_Fock_CoulombOnly0(IntegralOption,ThrPS(1),
     &                              Mode,Add,PackedD,PackedF,
     &                              nD,FactC,ip_D,ip_F)
         Go To 1 ! Return
      End If

      ! Return if nothing to do
      If (nD.lt.1) Go To 1

      ! Get number of basis functions (from localdf_bas.fh)
      nBas=nBas_Valence
      If (nBas.lt.1) Then
         Call WarningMessage(1,
     &                  SecNam//': nBas<1 -- Fock matrix NOT computed!')
         Write(6,'(A,I9)') 'nBas=',nBas
         Call xFlush(6)
         Go To 1  ! return
      End If

#if defined (_DEBUGPRINT_)
      If (.not.LDF_X_IsSet()) Then
         Call WarningMessage(2,SecNam//': LDF data not initialized!')
         Call LDF_Quit(1)
      End If
#endif

#if defined (_MOLCAS_MPP_)
      If (Add) Then
         If (nProcs.gt.1 .and. Is_Real_Par()) Then
            Write(6,'(A,A)') SecNam,
     &        ': >>Add<< feature not implemented in parallel execution!'
            Call LDF_NotImplemented()
         End If
      End If
#endif

      ! For half-and-half integral representation, scale FactC by 1/2.
      ! Save a copy of FactC (restored at the end)
      If (Mode.eq.3) Then
         l_FactC=nD
         Call GetMem('FactCBak','Allo','Real',ip_FactC,l_FactC)
         Call dCopy_(nD,FactC,1,Work(ip_FactC),1)
         Call dScal_(nD,0.5d0,FactC,1)
      Else
         l_FactC=0
         ip_FactC=0
      End If

      ! Set prescreening info (if not already done)
      If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
         Call LDF_SetIntegralPrescreeningInfo()
         IPI_set_here=.True.
      Else
         IPI_set_here=.False.
      End If
      tau(1)=max(ThrPS(1),0.0d0) ! integral prescreening threshold
      tau(2)=max(ThrPS(2),0.0d0) ! contribution prescreening threshold

      ! Initialize Fock matrices (if not Add)
      If (.not.Add) Then
         If (PackedF) Then
            l=nBas*(nBas+1)/2
         Else
            l=nBas**2
         End If
         Do iD=1,nD
            Call Cho_dZero(Work(ip_F(iD)),l)
         End Do
      End If

      ! Allocate and set blocked density matrices (atom pair blocks)
      l_DBlocks=nD
      Call GetMem('DBlk_P','Allo','Inte',ip_DBlocks,l_DBlocks)
      ip0=ip_DBlocks-1
#if defined (_DEBUGPRINT_)
      x=dble(NumberOfAtomPairs)
      y=dble(LDF_nAtom())*(dble(LDF_nAtom())+1.0d0)/2.0d0
      DoTest=int(x-y).eq.0
#endif
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('Den',iWork(ip0+iD))
         Call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,iWork(ip0+iD))
#if defined (_DEBUGPRINT_)
         If (DoTest) Then
            If (.not.LDF_TestBlockMatrix(iWork(ip0+iD),PackedD,
     &                                   Work(ip_D(iD)))) Then
               Call WarningMessage(2,
     &                            SecNam//': block matrix test failure')
               Write(6,'(A,I4,A,I9,3X,A,L1)')
     &         'Density matrix',iD,' at location',ip_D(iD),
     &         'Packed: ',PackedD
               Call LDF_Quit(1)
            End If
         End If
#endif
         Call LDF_ScaleOffdiagonalMatrixBlocks(iWork(ip0+iD),2.0d0)
      End Do

      ! Allocate and set blocked Fock matrices (atom pair blocks)
      l_FBlocks=nD
      Call GetMem('FBlk_P','Allo','Inte',ip_FBlocks,l_FBlocks)
      ip0=ip_FBlocks-1
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('Fck',iWork(ip0+iD))
         Call LDF_Full2Blocked(Work(ip_F(iD)),PackedF,iWork(ip0+iD))
      End Do

      ! Allocate and compute Frobenius norm of blocked density matrices
      l_DNorm=NumberOfAtomPairs*nD
      Call GetMem('DNorm','Allo','Real',ip_DNorm,l_DNorm)
      ip0=ip_DNorm
      Do iD=0,nD-1
         Call LDF_BlockMatrixNorm(iWork(ip_DBlocks+iD),ip0)
         ip0=ip0+NumberOfAtomPairs
      End Do

      ! Allocate Coulomb intermediates
      l_VP=nD
      Call GetMem('VP','Allo','Inte',ip_VP,l_VP)
      Do iD=0,nD-1
         Call LDF_AllocateAuxBasVector('CIn',iWork(ip_VP+iD))
      End Do

      ! Allocate array for norm of fitting coefficients
      l_CNorm=4*NumberOfAtomPairs
      Call GetMem('CNorm','Allo','Real',ip_CNorm,l_CNorm)

      ! Compute Coulomb intermediates
      ! V(J) = sum_uv C(uv,J)*D(uv)
      ! for each density matrix
      ! Compute Frobenius norm of fitting coefficients.
      Call LDF_ComputeCoulombIntermediates(Timing,nD,iWork(ip_DBlocks),
     &                                            iWork(ip_VP),ip_CNorm)

      ! Allocate and compute Frobenius norm of Coulomb intermediates
      l_VNorm=(LDF_nAtom()+NumberOfAtomPairs)*nD
      Call GetMem('VNorm','Allo','Real',ip_VNorm,l_VNorm)
      ip0=ip_VNorm
      Do iD=0,nD-1
         Call LDF_AuxBasVectorNorm(iWork(ip_VP+iD),ip0)
         ip0=ip0+LDF_nAtom()+NumberOfAtomPairs
      End Do

      ! Compute Coulomb contributions
      Call LDF_Fock_CoulombOnly_(IntegralOption.eq.444,Timing,Mode,tau,
     &                           nD,FactC,iWork(ip_DBlocks),
     &                           iWork(ip_VP),iWork(ip_FBlocks),
     &                           ip_CNorm,ip_DNorm,ip_VNorm)

      ! Get full storage (triangular or quadratic) Fock matrices from
      ! blocked ones.
      ip0=ip_FBlocks-1
      Do iD=1,nD
         Call LDF_Blocked2Full(iWork(ip0+iD),PackedF,Work(ip_F(iD)))
      End Do

      ! Deallocation
      Call GetMem('VNorm','Free','Real',ip_VNorm,l_VNorm)
      Call GetMem('CNorm','Free','Real',ip_CNorm,l_CNorm)
      Do iD=0,nD-1
         Call LDF_DeallocateAuxBasVector('CIn',iWork(ip_VP+iD))
      End Do
      Call GetMem('VP','Free','Inte',ip_VP,l_VP)
      Call GetMem('DNorm','Free','Real',ip_DNorm,l_DNorm)
      Do iD=0,nD-1
         Call LDF_DeallocateBlockMatrix('Fck',iWork(ip_FBlocks+iD))
      End Do
      Call GetMem('FBlk_P','Free','Inte',ip_FBlocks,l_FBlocks)
      Do iD=0,nD-1
         Call LDF_DeallocateBlockMatrix('Den',iWork(ip_DBlocks+iD))
      End Do
      Call GetMem('DBlk_P','Free','Inte',ip_DBlocks,l_DBlocks)

      ! Unset prescreening info (if set in this routine)
      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      ! Restore FactC
      If (l_FactC.gt.0) Then
         Call dCopy_(nD,Work(ip_FactC),1,FactC,1)
         Call GetMem('FactCBak','Free','Real',ip_FactC,l_FactC)
      End If

    1 Continue
      If (Timing) Then
         Call CWTime(tTotC2,tTotW2)
         Write(6,'(A,A,A,2(1X,F12.2),A)')
     &   'Total time spent in ',SecNam,':         ',
     &   tTotC2-tTotC1,tTotW2-tTotW1,' seconds'
         Write(6,'(84A1)') ('-',i=1,84)
         Call xFlush(6)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombOnly_(UseExactIntegralDiagonal,
     &                                 Timing,Mode,tau,
     &                                 nD,FactC,ip_DBlocks,ip_V,
     &                                 ip_FBlocks,
     &                                 ip_CNorm,ip_DNorm,ip_VNorm)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Compute Coulomb contributions to the Fock matrix using
C              Coulomb intermediates V. Integral-driven algorithm.
C
C     See LDF_Fock_CoulombOnly for more details.
C
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
      Logical UseExactIntegralDiagonal
      Logical Timing
      Integer Mode
      Real*8  tau(2)
      Integer nD
      Real*8  FactC(nD)
      Integer ip_DBlocks(nD)
      Integer ip_V(nD)
      Integer ip_FBlocks(nD)
      Integer ip_CNorm
      Integer ip_DNorm
      Integer ip_VNorm
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_integral_prescreening_info.fh"
#include "ldf_a2ap.fh"

      Character*21 SecNam
      Parameter (SecNam='LDF_Fock_CoulombOnly_')

      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer  LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom
      Integer  LDF_nBasAux_Pair_wLD
      External LDF_nAtom, LDF_nBas_Atom, LDF_nBasAux_Atom
      External LDF_nBasAux_Pair_wLD

      Integer ip_WBlkP, l_WBlkP
      Integer iD
      Integer nAtom
      Integer TaskListID
      Integer jAB, AB, A, B
      Integer CD, C
      Integer nuv, M, MA, MB, MAB, MCD
      Integer ip_Int, l_Int
      Integer ipD, ipV, ipF, ipW
      Integer ip_C, l_C, ipC
      Integer ip_tauW, l_tauW, ip_tauWA
      Integer ip, l
      Integer ip_VNrm, l_VNrm
      Integer ip_DNrm, l_DNrm

      Real*8  Const
      Real*8  tC1, tC2, tIC1, tIC2
      Real*8  tW1, tW2, tIW1, tIW2
      Real*8  tIC, tIW
      Real*8  tauW
      Real*8  GABCD, GCDAB

      Integer i, j
      Integer ip_W
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer A2AP
      Real*8  IAB
      Real*8  GA, GAB
      Real*8  VA, VAB
      Real*8  DAB
      Real*8  CAB_A, CAB_B, CAB_AB
      ip_W(i)=iWork(ip_WBlkP-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)
      IAB(i)=Work(ip_IDiag_Sm-1+i)
      GA(i)=Work(ip_GDiag_1C_Sm-1+i)
      GAB(i)=Work(ip_GDiag_2C_Sm-1+i)
      VA(i)=Work(ip_VNrm-1+i)
      VAB(i)=Work(ip_VNrm-1+nAtom+i)
      DAB(i)=Work(ip_DNrm-1+i)
      CAB_A(i)=Work(ip_CNorm+4*(i-1)+1)
      CAB_B(i)=Work(ip_CNorm+4*(i-1)+2)
      CAB_AB(i)=Work(ip_CNorm+4*(i-1)+3)

      ! Get number of atoms
      nAtom=LDF_nAtom()

      ! Set up prescreening info
      If (nD.eq.1) Then
         l_VNrm=0
         ip_VNrm=ip_VNorm
         l_DNrm=0
         ip_DNrm=ip_DNorm
      Else
         l_VNrm=nAtom+NumberOfAtomPairs
         Call GetMem('VNrm','Allo','Real',ip_VNrm,l_VNrm)
         Call Cho_dZero(Work(ip_VNrm),l_VNrm)
         Do iD=0,nD-1
            ip=ip_VNorm+(nAtom+NumberOfAtomPairs)*iD
            Do AB=0,nAtom+NumberOfAtomPairs-1
               Work(ip_VNrm+AB)=max(Work(ip_VNrm+AB),Work(ip+AB))
            End Do
         End Do
         l_DNrm=NumberOfAtomPairs
         Call GetMem('DNrm','Allo','Real',ip_DNrm,l_DNrm)
         Call Cho_dZero(Work(ip_DNrm),l_DNrm)
         Do iD=0,nD-1
            ip=ip_DNorm+NumberOfAtomPairs*iD
            Do AB=0,NumberOfAtomPairs-1
               Work(ip_DNrm+AB)=max(Work(ip_DNrm+AB),Work(ip+AB))
            End Do
         End Do
      End If
      l_tauW=nAtom+NumberOfAtomPairs
      Call GetMem('tauW','Allo','Real',ip_tauW,l_tauW)
      If (tau(2).gt.0.0d0) Then
         Call LDF_SetA2AP()
         Do A=1,nAtom
            l=A2AP(1,A)
            ip=A2AP(2,A)-1
            ip_tauWA=ip_tauW-1+A
            Work(ip_tauWA)=0.0d0
            Do jAB=1,l
               AB=iWork(ip+jAB)
               If (AP_Atoms(1,AB).eq.A) Then
                  Work(ip_tauWA)=max(Work(ip_tauWA),CAB_A(AB))
               Else If (AP_Atoms(2,AB).eq.A) Then
                  Work(ip_tauWA)=max(Work(ip_tauWA),CAB_B(AB))
               Else
                  Call WarningMessage(2,
     &                                 SecNam//': logical error [A2AP]')
                  Call LDF_Quit(1)
               End If
            End Do
            If (Work(ip_tauWA).gt.1.0d-16) Then
               tauW=tau(2)/Work(ip_tauWA)
            Else
               tauW=1.0d99
            End If
            Work(ip_tauWA)=tauW
         End Do
         Call LDF_UnsetA2AP()
         Do AB=1,NumberOfAtomPairs
            If (AP_2CFunctions(1,AB).gt.0) Then
               If (CAB_AB(AB).gt.1.0d-16) Then
                  Work(ip_tauW-1+nAtom+AB)=tau(2)/CAB_AB(AB)
               Else
                  Work(ip_tauW-1+nAtom+AB)=1.0d99
               End If
            Else
               Work(ip_tauW-1+nAtom+AB)=1.0d99
            End If
         End Do
      Else
         Call Cho_dZero(Work(ip_tauW),l_tauW)
      End If

      ! Allocate and initialize W intermediates
      l_WBlkP=nD
      Call GetMem('WBlkP','Allo','Inte',ip_WBlkP,l_WBlkP)
      Do iD=1,nD
         Call LDF_AllocateAuxBasVector('Win',iWork(ip_WBlkP-1+iD))
         Call LDF_ZeroAuxBasVector(ip_W(iD))
      End Do

C==========================
C     3-index contributions
C==========================

      If (Mode.eq.1 .or. Mode.eq.3) Then
         If (Timing) Then
            tIC=0.0d0
            tIW=0.0d0
            Call CWTime(tC1,tW1)
         End If
         Call Init_Tsk(TaskListID,NumberOfAtomPairs)
         Do While (Rsv_Tsk(TaskListID,AB))
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            nuv=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
            Do C=1,nAtom
               tauW=Work(ip_tauW-1+C)
               If (IAB(AB)*GA(C)*VA(C).ge.tau(2) .or.
     &             IAB(AB)*GA(C)*DAB(AB).ge.tauW) Then
                  ! Compute integrals (u_A v_B|J_C)
                  M=LDF_nBasAux_Atom(C)
                  l_Int=nuv*M
                  Call GetMem('Fck3Int1','Allo','Real',ip_Int,l_Int)
                  If (Timing) Call CWTime(tIC1,tIW1)
                  Call LDF_Compute3IndexIntegrals_1(AB,C,tau(1),l_Int,
     &                                              Work(ip_Int))
                  If (Timing) Then
                     Call CWTime(tIC2,tIW2)
                     tIC=tIC+(tIC2-tIC1)
                     tIW=tIW+(tIW2-tIW1)
                  End If
                  ! Compute Fock matrix contribution
                  ! F(u_A v_B)+=FactC*sum_J_C (u_A v_B|J_C)*V(J_C)
                  Do iD=1,nD
                     ipV=iWork(ip_V(iD)-1+C)
                     ipF=iWork(ip_FBlocks(iD)-1+AB)
                     Call dGeMV_('N',nuv,M,
     &                          FactC(iD),Work(ip_Int),nuv,
     &                          Work(ipV),1,1.0d0,Work(ipF),1)
                  End Do
                  ! Compute W contribution
                  ! W(J_C)+=sum_u_Av_B (u_A v_B|J_C)*D(u_A v_B)
                  Do iD=1,nD
                     ipD=iWork(ip_DBlocks(iD)-1+AB)
                     ipW=iWork(ip_W(iD)-1+C)
                     Call dGeMV_('T',nuv,M,
     &                          1.0d0,Work(ip_Int),nuv,
     &                          Work(ipD),1,1.0d0,Work(ipW),1)
                  End Do
                  Call GetMem('Fck3Int1','Free','Real',ip_Int,l_Int)
               End If
            End Do
            If (LDF2) Then
               ! Contributions from two-center aux functions
               Do CD=1,NumberOfAtomPairs
                  MCD=AP_2CFunctions(1,CD)
                  If (MCD.gt.0) Then
                     tauW=Work(ip_tauW-1+nAtom+CD)
                     If (IAB(AB)*GAB(CD)*VAB(CD).ge.tau(2) .or.
     &                   IAB(AB)*GAB(CD)*DAB(AB).ge.tauW) Then
                        ! Compute integrals (u_A v_B|J_CD)
                        l_Int=nuv*MCD
                        Call GetMem('Fck3Int2','Allo','Real',ip_Int,
     &                                                        l_Int)
                        If (Timing) Call CWTime(tIC1,tIW1)
                        Call LDF_Compute3IndexIntegrals_2(AB,CD,tau(1),
     &                                                    l_Int,
     &                                                    Work(ip_Int))
                        If (Timing) Then
                           Call CWTime(tIC2,tIW2)
                           tIC=tIC+(tIC2-tIC1)
                           tIW=tIW+(tIW2-tIW1)
                        End If
                        ! Compute Fock matrix contribution
                        ! F(u_A v_B)+=
                        ! FactC*sum_J_CD (u_A v_B|J_CD)*V(J_CD)
                        Do iD=1,nD
                           ipV=iWork(ip_V(iD)-1+nAtom+CD)
                           ipF=iWork(ip_FBlocks(iD)-1+AB)
                           Call dGeMV_('N',nuv,MCD,
     &                                FactC(iD),Work(ip_Int),nuv,
     &                                Work(ipV),1,1.0d0,Work(ipF),1)
                        End Do
                        ! Compute W contribution
                        ! W(J_CD)+=sum_u_Av_B (u_A v_B|J_CD)*D(u_A v_B)
                        Do iD=1,nD
                           ipD=iWork(ip_DBlocks(iD)-1+AB)
                           ipW=iWork(ip_W(iD)-1+nAtom+CD)
                           Call dGeMV_('T',nuv,MCD,
     &                                1.0d0,Work(ip_Int),nuv,
     &                                Work(ipD),1,1.0d0,Work(ipW),1)
                        End Do
                        Call GetMem('Fck3Int2','Free','Real',ip_Int,
     &                                                        l_Int)
                     End If
                  End If
               End Do
            End If
         End Do
         Call Free_Tsk(TaskListID)
         If (Timing) Then
            Call CWTime(tC2,tW2)
            Write(6,'(A,2(1X,F12.2),A)')
     &      'Time spent on 3-index contributions:              ',
     &      tC2-tC1,tW2-tW1,' seconds'
            Write(6,'(A,2(1X,F12.2),A)')
     &      '      - of which integrals required:              ',
     &      tIC,tIW,' seconds'
         End If
      End If

C==========================
C     2-index contributions
C==========================

      If (Mode.eq.1 .or. Mode.eq.2) Then
         If (Timing) Then
            tIC=0.0d0
            tIW=0.0d0
            Call CWTime(tC1,tW1)
         End If
         If (Mode.eq.1) Then
            Const=-1.0d0
         Else
            Const=1.0d0
         End If
         Call Init_Tsk(TaskListID,nAtom)
         Do While (Rsv_Tsk(TaskListID,A))
            MA=LDF_nBasAux_Atom(A)
            If (MA.gt.0) Then
               Do B=1,A-1
                  MB=LDF_nBasAux_Atom(B)
                  If (MB.gt.0) Then
                     If (GA(A)*GA(B)*VA(B).ge.Work(ip_tauW-1+A) .or.
     &                   GA(A)*GA(B)*VA(A).ge.Work(ip_tauW-1+B)) Then
                        ! Compute integrals (J_A|K_B)
                        l_Int=MA*MB
                        Call GetMem('Fck2Int11','Allo','Real',ip_Int,
     &                                                         l_Int)
                        If (Timing) Call CWTime(tIC1,tIW1)
                        Call LDF_Compute2IndexIntegrals_11(A,B,tau(1),
     &                                                     l_Int,
     &                                                     Work(ip_Int))
                        If (Timing) Then
                           Call CWTime(tIC2,tIW2)
                           tIC=tIC+(tIC2-tIC1)
                           tIW=tIW+(tIW2-tIW1)
                        End If
                        ! Compute W contribution
                        ! W(J_A)+=Const*sum_K_B (J_A|K_B)*V(K_B)
                        Do iD=1,nD
                           ipV=iWork(ip_V(iD)-1+B)
                           ipW=iWork(ip_W(iD)-1+A)
                           Call dGeMV_('N',MA,MB,
     &                                Const,Work(ip_Int),MA,
     &                                Work(ipV),1,1.0d0,Work(ipW),1)
                        End Do
                        ! Compute W contribution
                        ! W(K_B)+=Const*sum_J_A (J_A|K_B)*V(J_A)
                        Do iD=1,nD
                           ipV=iWork(ip_V(iD)-1+A)
                           ipW=iWork(ip_W(iD)-1+B)
                           Call dGeMV_('T',MA,MB,
     &                                Const,Work(ip_Int),MA,
     &                                Work(ipV),1,1.0d0,Work(ipW),1)
                        End Do
                        Call GetMem('Fck2Int11','Free','Real',ip_Int,
     &                                                         l_Int)
                     End If
                  End If
               End Do
               If (GA(A)*GA(A)*VA(A).ge.Work(ip_tauW-1+A)) Then
                  ! Compute integrals (J_A|K_A)
                  l_Int=MA**2
                  Call GetMem('Fck2Int11','Allo','Real',ip_Int,l_Int)
                  If (Timing) Call CWTime(tIC1,tIW1)
                  Call LDF_Compute2IndexIntegrals_11(A,A,tau(1),l_Int,
     &                                               Work(ip_Int))
                  If (Timing) Then
                     Call CWTime(tIC2,tIW2)
                     tIC=tIC+(tIC2-tIC1)
                     tIW=tIW+(tIW2-tIW1)
                  End If
                  ! Compute W contribution
                  ! W(J_A)+=Const*sum_K_A (J_A|K_A)*V(K_A)
                  Do iD=1,nD
                     ipV=iWork(ip_V(iD)-1+A)
                     ipW=iWork(ip_W(iD)-1+A)
                     Call dGeMV_('N',MA,MA,
     &                          Const,Work(ip_Int),MA,
     &                          Work(ipV),1,1.0d0,Work(ipW),1)
                  End Do
                  Call GetMem('Fck2Int11','Free','Real',ip_Int,l_Int)
               End If
               If (LDF2) Then
                  ! Two-center contributions
                  Do CD=1,NumberOfAtomPairs
                     MCD=AP_2CFunctions(1,CD)
                     If (MCD.gt.0) Then
                        If (GA(A)*GAB(CD)*VAB(CD).ge.Work(ip_tauW-1+A)
     &                      .or.
     &                      GA(A)*GAB(CD)*VA(A).ge.
     &                                         Work(ip_tauW-1+nAtom+CD))
     &                  Then
                           ! Compute integrals (J_A|K_CD)
                           l_Int=MA*MCD
                           Call GetMem('Fck2Int12','Allo','Real',ip_Int,
     &                                                            l_Int)
                           If (Timing) Call CWTime(tIC1,tIW1)
                           Call LDF_Compute2IndexIntegrals_12(A,CD,
     &                                                        tau(1),
     &                                                        l_Int,
     &                                                     Work(ip_Int))
                           If (Timing) Then
                              Call CWTime(tIC2,tIW2)
                              tIC=tIC+(tIC2-tIC1)
                              tIW=tIW+(tIW2-tIW1)
                           End If
                           ! Compute W contribution
                           ! W(J_A)+=Const*sum_K_CD (J_A|K_CD)*V(K_CD)
                           Do iD=1,nD
                              ipV=iWork(ip_V(iD)-1+nAtom+CD)
                              ipW=iWork(ip_W(iD)-1+A)
                              Call dGeMV_('N',MA,MCD,
     &                                   Const,Work(ip_Int),MA,
     &                                   Work(ipV),1,1.0d0,Work(ipW),1)
                           End Do
                           ! Compute W contribution
                           ! W(K_CD)+=Const*sum_J_A (J_A|K_CD)*V(J_A)
                           Do iD=1,nD
                              ipV=iWork(ip_V(iD)-1+A)
                              ipW=iWork(ip_W(iD)-1+nAtom+CD)
                              Call dGeMV_('T',MA,MCD,
     &                                   Const,Work(ip_Int),MA,
     &                                   Work(ipV),1,1.0d0,Work(ipW),1)
                           End Do
                           Call GetMem('Fck2Int12','Free','Real',ip_Int,
     &                                                            l_Int)
                        End If
                     End If
                  End Do
               End If
            End If
         End Do
         Call Free_Tsk(TaskListID)
         If (LDF2) Then
            Call Init_Tsk(TaskListID,NumberOfAtomPairs)
            Do While (Rsv_Tsk(TaskListID,AB))
               MAB=AP_2CFunctions(1,AB)
               If (MAB.gt.0) Then
                  Do CD=1,AB-1
                     MCD=AP_2CFunctions(1,CD)
                     If (MCD.gt.0) Then
                        GABCD=GAB(AB)*GAB(CD)*VAB(CD)
                        GCDAB=GAB(AB)*GAB(CD)*VAB(AB)
                        If (GABCD.ge.Work(ip_tauW-1+nAtom+AB) .or.
     &                      GCDAB.ge.Work(ip_tauW-1+nAtom+CD)) Then
                           ! Compute integrals (J_AB|K_CD)
                           l_Int=MAB*MCD
                           Call GetMem('Fck2Int22','Allo','Real',ip_Int,
     &                                                            l_Int)
                           If (Timing) Call CWTime(tIC1,tIW1)
                           Call LDF_Compute2IndexIntegrals_22(AB,CD,
     &                                                        tau(1),
     &                                                        l_Int,
     &                                                     Work(ip_Int))
                           If (Timing) Then
                              Call CWTime(tIC2,tIW2)
                              tIC=tIC+(tIC2-tIC1)
                              tIW=tIW+(tIW2-tIW1)
                           End If
                           ! Compute W contribution
                           ! W(J_AB)+=Const*sum_K_CD (J_AB|K_CD)*V(K_CD)
                           Do iD=1,nD
                              ipV=iWork(ip_V(iD)-1+nAtom+CD)
                              ipW=iWork(ip_W(iD)-1+nAtom+AB)
                              Call dGeMV_('N',MAB,MCD,
     &                                   Const,Work(ip_Int),MAB,
     &                                   Work(ipV),1,1.0d0,Work(ipW),1)
                           End Do
                           ! Compute W contribution
                           ! W(K_CD)+=Const*sum_J_AB (J_AB|K_CD)*V(J_AB)
                           Do iD=1,nD
                              ipV=iWork(ip_V(iD)-1+nAtom+AB)
                              ipW=iWork(ip_W(iD)-1+nAtom+CD)
                              Call dGeMV_('T',MAB,MCD,
     &                                   Const,Work(ip_Int),MAB,
     &                                   Work(ipV),1,1.0d0,Work(ipW),1)
                           End Do
                           Call GetMem('Fck2Int22','Free','Real',ip_Int,
     &                                                            l_Int)
                        End If
                     End If
                  End Do
                  ! Compute integrals (J_AB|K_AB)
                  If (GAB(AB)*GAB(AB)*VAB(AB).ge.
     &                                         Work(ip_tauW-1+nAtom+AB))
     &            Then
                     l_Int=MAB**2
                     Call GetMem('Fck2Int22','Allo','Real',ip_Int,l_Int)
                     If (Timing) Call CWTime(tIC1,tIW1)
                     Call LDF_Compute2IndexIntegrals_22(AB,AB,tau(1),
     &                                                  l_Int,
     &                                                  Work(ip_Int))
                     If (Timing) Then
                        Call CWTime(tIC2,tIW2)
                        tIC=tIC+(tIC2-tIC1)
                        tIW=tIW+(tIW2-tIW1)
                     End If
                     ! Compute W contribution
                     ! W(J_AB)+=Const*sum_K_AB (J_AB|K_AB)*V(K_AB)
                     Do iD=1,nD
                        ipV=iWork(ip_V(iD)-1+nAtom+AB)
                        ipW=iWork(ip_W(iD)-1+nAtom+AB)
                        Call dGeMV_('N',MAB,MAB,
     &                             Const,Work(ip_Int),MAB,
     &                             Work(ipV),1,1.0d0,Work(ipW),1)
                     End Do
                     Call GetMem('Fck2Int22','Free','Real',ip_Int,l_Int)
                  End If
               End If
            End Do
            Call Free_Tsk(TaskListID)
         End If
         If (Timing) Then
            Call CWTime(tC2,tW2)
            Write(6,'(A,2(1X,F12.2),A)')
     &      'Time spent on 2-index contributions:              ',
     &      tC2-tC1,tW2-tW1,' seconds'
            Write(6,'(A,2(1X,F12.2),A)')
     &      '      - of which integrals required:              ',
     &      tIC,tIW,' seconds'
         End If
      End If

      ! Deallocate prescreening info
      Call GetMem('tauW','Free','Real',ip_tauW,l_tauW)
      If (l_DNrm.gt.0) Then
         Call GetMem('DNrm','Free','Real',ip_DNrm,l_DNrm)
      End If
      If (l_VNrm.gt.0) Then
         Call GetMem('VNrm','Free','Real',ip_VNrm,l_VNrm)
      End If

#if defined (_MOLCAS_MPP_)
      ! Add W over nodes
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         If (Timing) Call CWTime(tC1,tW1)
         Do iD=1,nD
            Call LDF_P_AddAuxBasVector(ip_W(iD))
         End Do
         If (Timing) Then
            Call CWTime(tC2,tW2)
            Write(6,'(A,2(1X,F12.2),A)')
     &      'Parallel overhead for W intermediates:            ',
     &      tC2-tC1,tW2-tW1,' seconds'
         End If
      End If
#endif

C===================
C     Compute F+=C*W
C===================

      If (Timing) Call CWTime(tC1,tW1)
      Call Init_Tsk(TaskListID,NumberOfAtomPairs)
      Do While (Rsv_Tsk(TaskListID,AB))
         ! Read coefficients
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nuv=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         l_C=nuv*LDF_nBasAux_Pair_wLD(AB)
         Call GetMem('FckCoef','Allo','Real',ip_C,l_C)
         Call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l_C)
         ! Compute F(u_A v_B) += FactC*sum_J C(u_A v_B,J)*W(J)
         ipC=ip_C
         MA=LDF_nBasAux_Atom(A)
         If (MA.gt.0) Then
            Do iD=1,nD
               ipF=iWork(ip_FBlocks(iD)-1+AB)
               ipW=iWork(ip_W(iD)-1+A)
               Call dGeMV_('N',nuv,MA,
     &                    FactC(iD),Work(ipC),nuv,
     &                    Work(ipW),1,1.0d0,Work(ipF),1)
            End Do
            ipC=ipC+nuv*MA
         End If
         If (B.ne.A) Then
            MB=LDF_nBasAux_Atom(B)
            If (MB.gt.0) Then
               Do iD=1,nD
                  ipF=iWork(ip_FBlocks(iD)-1+AB)
                  ipW=iWork(ip_W(iD)-1+B)
                  Call dGeMV_('N',nuv,MB,
     &                       FactC(iD),Work(ipC),nuv,
     &                       Work(ipW),1,1.0d0,Work(ipF),1)
               End Do
               ipC=ipC+nuv*MB
            End If
         End If
         MAB=AP_2CFunctions(1,AB)
         If (MAB.gt.0) Then
            Do iD=1,nD
               ipF=iWork(ip_FBlocks(iD)-1+AB)
               ipW=iWork(ip_W(iD)-1+nAtom+AB)
               Call dGeMV_('N',nuv,MAB,
     &                    FactC(iD),Work(ipC),nuv,
     &                    Work(ipW),1,1.0d0,Work(ipF),1)
            End Do
         End If
         Call GetMem('FckCoef','Free','Real',ip_C,l_C)
      End Do
      Call Free_Tsk(TaskListID)
      If (Timing) Then
         Call CWTime(tC2,tW2)
         Write(6,'(A,2(1X,F12.2),A)')
     &   'Time spent computing C*W contribution:            ',
     &   tC2-tC1,tW2-tW1,' seconds'
      End If

      ! Deallocate W intermediates
      Do iD=0,nD-1
         Call LDF_DeallocateAuxBasVector('Win',iWork(ip_WBlkP+iD))
      End Do
      Call GetMem('WBlkP','Free','Inte',ip_WBlkP,l_WBlkP)

      ! Use exact integral diagonal blocks if requested.
      ! Computed by adding the correction
      ! dF(uv) = FactC * sum_kl { (uv|kl) - [uv|kl] }*D(kl)
      ! where |uv) and |kl) belong to the same atom pair
      ! and [uv|kl] is an LDF integral approximation.
      If (UseExactIntegralDiagonal) Then
         If (Timing) Call CWTime(tC1,tW1)
         Call LDF_Fock_CoulombOnly_XIDI(Mode,tau(1),nD,FactC,
     &                                  ip_DBlocks,ip_FBlocks)
         If (Timing) Then
            Call CWTime(tC2,tW2)
            Write(6,'(A,2(1X,F12.2),A)')
     &      'Time spent computing XIDI corrections:            ',
     &      tC2-tC1,tW2-tW1,' seconds'
         End If
      End If

#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         If (Timing) Call CWTime(tC1,tW1)
         Do iD=1,nD
            Call LDF_P_AddBlockMatrix(ip_FBlocks(iD))
         End Do
         If (Timing) Then
            Call CWTime(tC2,tW2)
            Write(6,'(A,2(1X,F12.2),A)')
     &      'Parallel overhead for F blocks:                   ',
     &      tC2-tC1,tW2-tW1,' seconds'
         End If
      End If
#endif

      If (Timing) Call xFlush(6)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC                                                          CCCCCCC
CCCCCCC                      OLD CODE                            CCCCCCC
CCCCCCC                                                          CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombOnly0(IntegralOption,tau,
     &                                 Mode,Add,PackedD,PackedF,
     &                                 nD,FactC,ip_D,ip_F)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Compute Coulomb contribution to Fock matrix using local
C              Density Fitting coefficients. Use only for debugging!
C              Poor performance!
C
C     Args:    - IntegralOption (integer)
C                  111  => use two-electron integrals computed from
C                          LDF coefficients (for debugging)
C                  222  => use conventional two-electron integrals
C                          (for debugging)
C                  333  => use conventional or LDF two-electron
C                          integrals depending on positivity of the
C                          latter (for debugging)
C                  all other values => use production-level LDF code
C                  (INPUT)
C              - tau (real*8) integral prescreening threshold (INPUT).
C                Only used with IntegralOption=111.
C              - Mode (integer) 1:robust,2:nonrobust,3:half-and-half
C              - Add (boolean): .True. if Coulomb contribution should be
C                added to the Fock matrix (INPUT). Not available in
C                parallel execution!!
C              - PackedD (boolean): .True. if density matrices
C                are stored in lower triangular format; else quadratic
C                storage (INPUT)
C              - PackedF (boolean): .True. if Fock matrices
C                are stored in lower triangular format; else quadratic
C                storage (INPUT)
C              - nD: Number of densities/Fock matrices (INPUT)
C              - FactC(nD): scaling factor for each density (INPUT)
C              - ip_D(nD): pointers to nD density matrices (if PackedD:
C                lower triangular storage) (INPUT)
C              - ip_F(nD): pointers to nD Fock matrices (if PackedF:
C                lower triangular storage) (INPUT)
C
C              If (Add):  [NOT IMPLEMENTED IN PARALLEL]
C       (1)       F(uv) = F(uv) + FactC * sum_kl (uv|kl)*D(kl)
C              Else:
C       (2)       F(uv) = FactC * sum_kl (uv|kl)*D(kl)
C
C              where the integrals are given by the robust LDF
C              representation
C
C              (uv|kl) = ([uv]|kl) + (uv|[kl]) - ([uv]|[kl])
C
C              or the non-robust representation
C
C              (uv|kl) = ([uv]|[kl])
C
C              or the half-and-half representation
C
C              (uv|kl) = 0.5*([uv]|kl) + 0.5*(uv|[kl])
C
C              and the fitted products are given by
C
C              |[uv]) = sum_J C(uv,J) |J)
C
C              where the fitting functions |J) are centered on the atom
C              pair to which the product |uv) belongs.
C
C              In more detail (assuming robust fitting representation):
C
C              F(uv) = FactC * {sum_J C(uv,J)*W(J) + sum_K (uv|K)*V(K)}
C
C              W(J) = sum_kl (J|kl)*D(kl) - sum_K (J|K)*V(K)
C
C              V(K) = sum_kl C(kl,K)*D(kl)
C
C     NOTES:
C       - The density matrices are assumed to be symmetric.
C       - It is the complete sum over kl in Eqs(1,2)! This means that
C         the input density matrices should NOT be scaled by 2 on the
C         off-diagonal (scaling factors are handled internally).
C       - LDF information must be properly set up before calling this
C         routine (use subroutine LDF_X_Init).
C       - Two-center functions may or mat not be included in the fitting
C         basis for a given atom pair.
C       - The algorithm is NOT integral driven. This means that 3-center
C         integrals (4-center integrals if 2-center fitting functions
C         are included) are computed several times, causing quite poor
C         performance.
C       - The max number of processes for which the algorithm can
C         possibly scale well equals the number of LDF atom pairs.
C       - No screening has been implemented.
C       - Add is not implemented in parallel !!!
C
C     Outline of the algorithm (robust representation):
C     =================================================
C
C     *If (Mode=3): FactC := 0.5*FactC
C
C     *Scale off-diagonal density blocks:
C      D[u_A v_B] <-- 2*D[u_A v_B] iff A != B
C
C     *Parallel loop over atom pairs AB (A>=B):
C        - V(J_AB) = sum_[u_A v_B] C(u_A v_B,J_AB)*D(u_A v_B)
C     *End loop AB
C     Add V over nodes.
C
C     *Parallel loop over atom pairs AB (A>=B):
C        *Loop over atom pairs CD (C>=D):
C           - F(u_A v_B) = F(u_A v_B)
C                        + sum_[K_CD] (u_A v_B | K_CD)*V(K_CD)
C           - W(J_AB) = W(J_AB)
C                     + sum_[k_C l_D] (J_AB | k_C l_D)*D(k_C l_D)
C           - W(J_AB) = W(J_AB)
C                     - sum_[K_CD] (J_AB | K_CD)*V(K_CD)
C        *End loop CD
C        - F(u_A v_B) = F(u_A v_B)
C                     + sum_[J_AB] C(u_A v_B,J_AB)*W(J_AB)
C     *End loope AB
C     Add F over nodes
C
C     Note that both arrays V and F are stored locally as O(N) arrays,
C     which should keep the communication bottleneck at a minimum as the
C     system size grows.
C
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
      Integer IntegralOption
      Real*8  tau
      Integer Mode
      Logical Add
      Logical PackedD
      Logical PackedF
      Integer nD
      Real*8  FactC(nD)
      Integer ip_D(nD)
      Integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*21 SecNam
      Parameter (SecNam='LDF_Fock_CoulombOnly0')

#if defined (_DEBUGPRINT_)
      Logical  LDF_X_IsSet, LDF_TestBlockMatrix
      External LDF_X_IsSet, LDF_TestBlockMatrix
      Integer  LDF_nAtom
      External LDF_nAtom
      Logical  DoTest
      Real*8   x, y
#endif

      Logical UsePartPermSym

      Integer nBas
      Integer iD, ip0, l
      Integer ip_DBlocks, l_DBlocks
      Integer ip_FBlocks, l_FBlocks
      Integer ip_VBlocks, l_VBlocks
      Integer ip_FactC, l_FactC

      ! Return if nothing to do
      If (nD.lt.1) Return

      ! Get number of basis functions (from localdf_bas.fh)
      nBas=nBas_Valence
      If (nBas.lt.1) Then
         Call WarningMessage(1,
     &                  SecNam//': nBas<1 -- Fock matrix NOT computed!')
         Write(6,'(A,I9)') 'nBas=',nBas
         Call xFlush(6)
         Return
      End If

#if defined (_DEBUGPRINT_)
      If (.not.LDF_X_IsSet()) Then
         Call WarningMessage(2,SecNam//': LDF data not initialized!')
         Call LDF_Quit(1)
      End If
#endif

#if defined (_MOLCAS_MPP_)
      If (Add) Then
         If (nProcs.gt.1 .and. Is_Real_Par()) Then
            Write(6,'(A,A)') SecNam,
     &        ': >>Add<< feature not implemented in parallel execution!'
            Call LDF_NotImplemented()
         End If
      End If
#endif

      ! For half-and-half integral representation, scale FactC by 1/2.
      ! Save a copy of FactC (restored at the end)
      If (Mode.eq.3) Then
         l_FactC=nD
         Call GetMem('FactCBak','Allo','Real',ip_FactC,l_FactC)
         Call dCopy_(nD,FactC,1,Work(ip_FactC),1)
         Call dScal_(nD,0.5d0,FactC,1)
      Else
         l_FactC=0
         ip_FactC=0
      End If

      ! Initialize Fock matrices (if not Add)
      If (.not.Add) Then
         If (PackedF) Then
            l=nBas*(nBas+1)/2
         Else
            l=nBas**2
         End If
         Do iD=1,nD
            Call Cho_dZero(Work(ip_F(iD)),l)
         End Do
      End If

      ! Allocate and set blocked density matrices (atom pair blocks)
      l_DBlocks=nD
      Call GetMem('DBlk_P','Allo','Inte',ip_DBlocks,l_DBlocks)
      ip0=ip_DBlocks-1
#if defined (_DEBUGPRINT_)
      x=dble(NumberOfAtomPairs)
      y=dble(LDF_nAtom())*(dble(LDF_nAtom())+1.0d0)/2.0d0
      DoTest=int(x-y).eq.0
#endif
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('Den',iWork(ip0+iD))
         Call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,iWork(ip0+iD))
#if defined (_DEBUGPRINT_)
         If (DoTest) Then
            If (.not.LDF_TestBlockMatrix(iWork(ip0+iD),PackedD,
     &                                   Work(ip_D(iD)))) Then
               Call WarningMessage(2,
     &                            SecNam//': block matrix test failure')
               Write(6,'(A,I4,A,I9,3X,A,L1)')
     &         'Density matrix',iD,' at location',ip_D(iD),
     &         'Packed: ',PackedD
               Call LDF_Quit(1)
            End If
         End If
#endif
         Call LDF_ScaleOffdiagonalMatrixBlocks(iWork(ip0+iD),2.0d0)
      End Do

      ! Allocate and set blocked Fock matrices (atom pair blocks)
      l_FBlocks=nD
      Call GetMem('FBlk_P','Allo','Inte',ip_FBlocks,l_FBlocks)
      ip0=ip_FBlocks-1
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('Fck',iWork(ip0+iD))
         Call LDF_Full2Blocked(Work(ip_F(iD)),PackedF,iWork(ip0+iD))
      End Do

      If (IntegralOption.eq.111) Then
         ! Compute Fock matrix using LDF integrals (debug)
         Call WarningMessage(0,
     &               SecNam//': Using integrals from LDF coefficients!')
         Call xFlush(6)
         UsePartPermSym=.True.
         If (Mode.eq.3) Then
            Call LDF_FVIFC(UsePartPermSym,Mode,max(tau,0.0d0),
     &                     nD,Work(ip_FactC),iWork(ip_DBlocks),
     &                     iWork(ip_FBlocks))
         Else
            Call LDF_FVIFC(UsePartPermSym,Mode,max(tau,0.0d0),
     &                     nD,FactC,iWork(ip_DBlocks),
     &                     iWork(ip_FBlocks))
         End If
      Else If (IntegralOption.eq.222) Then
         ! Compute Fock matrix using conventional integrals (debug)
         Call WarningMessage(0,
     &                        SecNam//': Using conventional integrals!')
         Call xFlush(6)
         UsePartPermSym=.True.
         Call LDF_FCI(UsePartPermSym,nD,FactC,iWork(ip_DBlocks),
     &                                        iWork(ip_FBlocks))
      Else If (IntegralOption.eq.333) Then
         ! Compute Fock matrix using conventional or LDF integrals
         ! depending on positivity of the latter (debug)
         Call WarningMessage(0,
     &                  SecNam//': Using PSD (LDF or conv.) integrals!')
         Call xFlush(6)
         UsePartPermSym=.True.
         If (Mode.eq.3) Then
            Call LDF_FTst(UsePartPermSym,Mode,max(tau,0.0d0),
     &                    nD,Work(ip_FactC),iWork(ip_DBlocks),
     &                    iWork(ip_FBlocks))
         Else
            Call LDF_FTst(UsePartPermSym,Mode,max(tau,0.0d0),
     &                    nD,FactC,iWork(ip_DBlocks),
     &                    iWork(ip_FBlocks))
         End If
      Else
         ! Allocate Coulomb intermediates
         l_VBlocks=nD
         Call GetMem('VBlk_P','Allo','Inte',ip_VBlocks,l_VBlocks)
         Do iD=0,nD-1
            Call LDF_AllocateBlockVector('CIn',iWork(ip_VBlocks+iD))
         End Do
         ! Compute Coulomb intermediates,
         ! V(J) = sum_uv C(uv,J)*D(uv)
         ! for each density matrix
         Call LDF_ComputeCoulombIntermediates0(nD,iWork(ip_DBlocks),
     &                                            iWork(ip_VBlocks))
         ! Compute Coulomb contributions
         Call LDF_Fock_CoulombOnly0_(Mode,nD,FactC,
     &                              iWork(ip_DBlocks),iWork(ip_VBlocks),
     &                              iWork(ip_FBlocks))
         ! Deallocate Coulomb intermediates
         Do iD=0,nD-1
            Call LDF_DeallocateBlockVector('CIn',iWork(ip_VBlocks+iD))
         End Do
         Call GetMem('VBlk_P','Free','Inte',ip_VBlocks,l_VBlocks)
      End If

      ! Get full storage (triangular or quadratic) Fock matrices from
      ! blocked ones.
      ip0=ip_FBlocks-1
      Do iD=1,nD
         Call LDF_Blocked2Full(iWork(ip0+iD),PackedF,Work(ip_F(iD)))
      End Do

      ! Restore FactC
      If (l_FactC.gt.0) Then
         Call dCopy_(nD,Work(ip_FactC),1,FactC,1)
         Call GetMem('FactCBak','Free','Real',ip_FactC,l_FactC)
      End If

      ! Deallocation
      Do iD=0,nD-1
         Call LDF_DeallocateBlockMatrix('Fck',iWork(ip_FBlocks+iD))
      End Do
      Call GetMem('FBlk_P','Free','Inte',ip_FBlocks,l_FBlocks)
      Do iD=0,nD-1
         Call LDF_DeallocateBlockMatrix('Den',iWork(ip_DBlocks+iD))
      End Do
      Call GetMem('DBlk_P','Free','Inte',ip_DBlocks,l_DBlocks)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombOnly0_(Mode,nD,FactC,
     &                                 ip_DBlocks,ip_VBlocks,ip_FBlocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Compute Coulomb contributions to the Fock matrix using
C              Coulomb intermediates V.
C
C     See LDF_Fock_CoulombOnly for an outline of the algorithm.
C
      Implicit None
      Integer Mode
      Integer nD
      Real*8  FactC(nD)
      Integer ip_DBlocks(nD)
      Integer ip_VBlocks(nD)
      Integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*22 SecNam
      Parameter (SecNam='LDF_Fock_CoulombOnly0_')

      Logical  Rsv_Tsk
      External Rsv_Tsk

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer ip_WBlkP, l_WBlkP
      Integer iD
      Integer TaskListID
      Integer AB, CD
      Integer ip, l
      Integer nuv, M
      Integer ipW, ipF


      Integer i, j
      Integer ip_WBlocks
      Integer AP_Atoms
      ip_WBlocks(i)=iWork(ip_WBlkP-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Allocate and initialize W intermediates
      l_WBlkP=nD
      Call GetMem('WBlk_P','Allo','Inte',ip_WBlkP,l_WBlkP)
      Do iD=1,nD
         Call LDF_AllocateBlockVector('Win',iWork(ip_WBlkP-1+iD))
         Call LDF_ZeroBlockVector(ip_WBlocks(iD))
      End Do

      If (Mode.eq.1 .or. Mode.eq.3) Then
         ! Parallel loop over atom pairs AB (A>=B)
         Call Init_Tsk(TaskListID,NumberOfAtomPairs)
         Do While (Rsv_Tsk(TaskListID,AB))
            ! Serial loop over atom pairs CD (C>=D)
            Do CD=1,NumberOfAtomPairs
               ! F(u_A v_B) = F(u_A v_B)
               !            + sum_[K_CD] (u_A v_B | K_CD)*V(K_CD)
               Call LDF_Fock_CoulombOnly0_1(nD,FactC,ip_VBlocks,
     &                                      ip_FBlocks,AB,CD)
               ! W(J_AB) = W(J_AB)
               !         + sum_[k_C l_D] (J_AB | k_C l_D)*D(k_C l_D)
               Call LDF_Fock_CoulombOnly0_2(nD,ip_DBlocks,
     &                                      iWork(ip_WBlkP),AB,CD)
               If (Mode.eq.1) Then
                  ! W(J_AB) = W(J_AB)
                  !         - sum_[K_CD] (J_AB | K_CD)*V(K_CD)
                  Call LDF_Fock_CoulombOnly0_3(-1.0d0,nD,ip_VBlocks,
     &                                         iWork(ip_WBlkP),AB,CD)
               End If
            End Do ! end serial loop over atom pairs CD
            ! F(u_A v_B) = F(u_A v_B)
            !            + sum_[J_AB] C(u_A v_B,J_AB)*W(J_AB)
            nuv=LDF_nBas_Atom(AP_Atoms(1,AB))
     &         *LDF_nBas_Atom(AP_Atoms(2,AB))
            M=LDF_nBasAux_Pair(AB)
            l=nuv*M
            Call GetMem('C_AB','Allo','Real',ip,l)
            Call LDF_CIO_ReadC(AB,Work(ip),l)
            Do iD=1,nD
               ipF=iWork(ip_FBlocks(iD)-1+AB)
               ipW=iWork(ip_WBlocks(iD)-1+AB)
               Call dGeMV_('N',nuv,M,
     &                    FactC(iD),Work(ip),nuv,
     &                    Work(ipW),1,1.0d0,Work(ipF),1)
            End Do
            Call GetMem('C_AB','Free','Real',ip,l)
         End Do ! end parallel loop over atom pairs AB
         Call Free_Tsk(TaskListID)
      Else If (Mode.eq.2) Then ! non-robust fitting
         ! Parallel loop over atom pairs AB (A>=B)
         Call Init_Tsk(TaskListID,NumberOfAtomPairs)
         Do While (Rsv_Tsk(TaskListID,AB))
            ! Serial loop over atom pairs CD (C>=D)
            Do CD=1,NumberOfAtomPairs
               ! W(J_AB) = W(J_AB)
               !         + sum_[K_CD] (J_AB | K_CD)*V(K_CD)
               Call LDF_Fock_CoulombOnly0_3(1.0d0,nD,ip_VBlocks,
     &                                      iWork(ip_WBlkP),AB,CD)
            End Do ! end serial loop over atom pairs CD
            ! F(u_A v_B) = F(u_A v_B)
            !            + sum_[J_AB] C(u_A v_B,J_AB)*W(J_AB)
            nuv=LDF_nBas_Atom(AP_Atoms(1,AB))
     &         *LDF_nBas_Atom(AP_Atoms(2,AB))
            M=LDF_nBasAux_Pair(AB)
            l=nuv*M
            Call GetMem('C_AB','Allo','Real',ip,l)
            Call LDF_CIO_ReadC(AB,Work(ip),l)
            Do iD=1,nD
               ipF=iWork(ip_FBlocks(iD)-1+AB)
               ipW=iWork(ip_WBlocks(iD)-1+AB)
               Call dGeMV_('N',nuv,M,
     &                    FactC(iD),Work(ip),nuv,
     &                    Work(ipW),1,1.0d0,Work(ipF),1)
            End Do
            Call GetMem('C_AB','Free','Real',ip,l)
         End Do ! end parallel loop over atom pairs AB
         Call Free_Tsk(TaskListID)
      Else
         Write(6,'(A,A,I6)') SecNam,': unknown Mode:',Mode
         Call LDF_NotImplemented()
      End If

#if defined (_MOLCAS_MPP_)
      ! Add F over nodes
      Do iD=1,nD
         Call LDF_P_AddBlockMatrix(ip_FBlocks(iD))
      End Do
#endif

      ! Deallocate W intermediates
      Do iD=0,nD-1
         Call LDF_DeallocateBlockVector('Win',iWork(ip_WBlkP+iD))
      End Do
      Call GetMem('WBlk_P','Free','Inte',ip_WBlkP,l_WBlkP)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombOnly0_1(nD,FactC,ip_VBlocks,ip_FBlocks,
     &                                   AB,CD)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Compute
C
C              F(u_A v_B) = F(u_A v_B)
C                         + FactC*sum_[K_CD] (u_A v_B | K_CD)*V(K_CD)
C
      Implicit None
      Integer nD
      Real*8  FactC(nD)
      Integer ip_VBlocks(nD)
      Integer ip_FBlocks(nD)
      Integer AB, CD
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer nuv, M
      Integer ip_Int, l_Int
      Integer iD
      Integer ipV, ipF

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Get row and column dimensions of integrals
      nuv=LDF_nBas_Atom(AP_Atoms(1,AB))
     &   *LDF_nBas_Atom(AP_Atoms(2,AB))
      M=LDF_nBasAux_Pair(CD)

      ! Return if nothing to do
      If (nuv.lt.1 .or. M.lt.1) Return

      ! Allocate integrals (u_A v_B | K_CD)
      l_Int=nuv*M
      Call GetMem('LDFFuvJ1','Allo','Real',ip_Int,l_Int)

      ! Compute integrals (u_A v_B | K_CD)
      Call LDF_ComputeIntegrals_uvJ_2P(AB,CD,l_Int,Work(ip_Int))

      ! Compute contributions
      Do iD=1,nD
         ipV=iWork(ip_VBlocks(iD)-1+CD)
         ipF=iWork(ip_FBlocks(iD)-1+AB)
         Call dGeMV_('N',nuv,M,
     &              FactC(iD),Work(ip_Int),nuv,
     &              Work(ipV),1,1.0d0,Work(ipF),1)
      End Do

      ! Deallocate integrals
      Call GetMem('LDFFuvJ1','Free','Real',ip_Int,l_Int)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombOnly0_2(nD,ip_DBlocks,ip_WBlocks,AB,CD)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Compute
C
C              W(J_AB) = W(J_AB)
C                      + sum_[k_C l_D] (k_C l_D|J_AB)*D(k_C l_D)
C
      Implicit None
      Integer nD
      Integer ip_DBlocks(nD)
      Integer ip_WBlocks(nD)
      Integer AB, CD
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer nkl, M
      Integer ip_Int, l_Int
      Integer iD
      Integer ipD, ipW

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Get row and column dimensions of integrals
      nkl=LDF_nBas_Atom(AP_Atoms(1,CD))
     &   *LDF_nBas_Atom(AP_Atoms(2,CD))
      M=LDF_nBasAux_Pair(AB)

      ! Return if nothing to do
      If (nkl.lt.1 .or. M.lt.1) Return

      ! Allocate integrals (k_C l_D | J_AB)
      l_Int=nkl*M
      Call GetMem('LDFFuvJ2','Allo','Real',ip_Int,l_Int)

      ! Compute integrals (k_C l_D | J_AB)
      Call LDF_ComputeIntegrals_uvJ_2P(CD,AB,l_Int,Work(ip_Int))

      ! Compute contributions
      Do iD=1,nD
         ipD=iWork(ip_DBlocks(iD)-1+CD)
         ipW=iWork(ip_WBlocks(iD)-1+AB)
         Call dGeMV_('T',nkl,M,
     &              1.0d0,Work(ip_Int),nkl,
     &              Work(ipD),1,1.0d0,Work(ipW),1)
      End Do

      ! Deallocate integrals
      Call GetMem('LDFFuvJ2','Free','Real',ip_Int,l_Int)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombOnly0_3(Const,nD,ip_VBlocks,ip_WBlocks,
     &                                  AB,CD)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Compute
C
C              W(J_AB) = W(J_AB)
C                      + Const * sum_[K_CD] (J_AB | K_CD)*V(K_CD)
C
      Implicit None
      Real*8  Const
      Integer nD
      Integer ip_VBlocks(nD)
      Integer ip_WBlocks(nD)
      Integer AB, CD
#include "WrkSpc.fh"

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Integer MAB, MCD
      Integer ip_Int, l_Int
      Integer iD
      Integer ipV, ipW

      ! Get row and column dimension of integrals
      MAB=LDF_nBasAux_Pair(AB)
      MCD=LDF_nBasAux_Pair(CD)

      ! Return if nothing to do
      If (MAB.lt.1 .or. MCD.lt.1) Return

      ! Allocate integrals (J_AB | K_CD)
      l_Int=MAB*MCD
      Call GetMem('LDFFJK','Allo','Real',ip_Int,l_Int)

      ! Compute integrals (J_AB | K_CD)
      Call LDF_ComputeIntegrals_JK_2P(AB,CD,l_Int,Work(ip_Int))

      ! Compute contributions
      Do iD=1,nD
         ipV=iWork(ip_VBlocks(iD)-1+CD)
         ipW=iWork(ip_WBlocks(iD)-1+AB)
         Call dGeMV_('N',MAB,MCD,
     &              Const,Work(ip_Int),MAB,
     &              Work(ipV),1,1.0d0,Work(ipW),1)
      End Do

      ! Deallocate integrals
      Call GetMem('LDFFJK','Free','Real',ip_Int,l_Int)

      End
