!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nD,FactC,ip_D,ip_F)
! Thomas Bondo Pedersen, October 2010.
!
! Purpose: Compute Coulomb contribution to Fock matrix using local
!          Density Fitting coefficients.
!
! Args:    - IntegralOption (integer)
!              111  => use two-electron integrals computed from
!                      LDF coefficients (for debugging)
!              222  => use conventional two-electron integrals
!                      (for debugging)
!              333  => use conventional or LDF two-electron
!                      integrals depending on positivity of the
!                      latter (for debugging)
!              444  => use exact integral diagonal (blocks), all
!                      other blocks are LDF
!              all other values => use production-level LDF code
!              (INPUT)
!          - Timing (boolean) .True. if detailed timing (incl total)
!            should be printed from this routine (INPUT)
!          - Mode (integer) 1:robust, 2:nonrobust, 3:half-and-half
!            to be used (INPUT)
!          - ThrPS(2): prescreening thresholds (INPUT)
!                      ThrPS(1): integral prescreening threshold
!                      ThrPS(2): prescreening threshold on
!                                contributions
!          - Add (boolean): .True. if Coulomb contribution should be
!            added to the Fock matrix (INPUT). Not available in
!            parallel execution!!
!          - PackedD (boolean): .True. if density matrices
!            are stored in lower triangular format; else quadratic
!            storage (INPUT)
!          - PackedF (boolean): .True. if Fock matrices
!            are stored in lower triangular format; else quadratic
!            storage (INPUT)
!          - nD: Number of densities/Fock matrices (INPUT)
!          - FactC(nD): scaling factor for each density (INPUT)
!          - ip_D(nD): pointers to nD density matrices (if PackedD:
!            lower triangular storage) (INPUT)
!          - ip_F(nD): pointers to nD Fock matrices (if PackedF:
!            lower triangular storage) (INPUT)
!
!          If (Add):  [NOT IMPLEMENTED IN PARALLEL]
!   (1)       F(uv) = F(uv) + FactC * sum_kl (uv|kl)*D(kl)
!          Else:
!   (2)       F(uv) = FactC * sum_kl (uv|kl)*D(kl)
!
!          where the integrals are given by the robust LDF
!          representation (Mode=1)
!
!          (uv|kl) = ([uv]|kl) + (uv|[kl]) - ([uv]|[kl])
!
!          or the non-robust representation (Mode=2)
!
!          (uv|kl) = ([uv]|[kl])
!
!          or the half-and-half representation (Mode=3)
!
!          (uv|kl) = 0.5*([uv]|kl) + 0.5*(uv|[kl])
!
!          The fitted products are given by
!
!          |[uv]) = sum_J C(uv,J) |J)
!
!          where the fitting functions |J) are centered on the atom
!          pair to which the product |uv) belongs.
!
!          In more detail (assuming robust fitting representation):
!
!          F(uv) = FactC * {sum_J C(uv,J)*W(J) + sum_K (uv|K)*V(K)}
!
!          W(J) = sum_kl (J|kl)*D(kl) - sum_K (J|K)*V(K)
!
!          V(K) = sum_kl C(kl,K)*D(kl)
!
! NOTES:
!   - The density matrices are assumed to be symmetric.
!   - It is the complete sum over kl in Eqs(1,2)! This means that
!     the input density matrices should NOT be scaled by 2 on the
!     off-diagonal (scaling factors are handled internally).
!   - LDF information must be properly set up before calling this
!     routine (use subroutine LDF_X_Init).
!   - Two-center functions may or mat not be included in the fitting
!     basis for a given atom pair.
!   - The algorithm is integral driven: the integrals are computed
!     once.
!   - The fitting coefficients are read from disk twice (for Coulomb
!     [V] intermediates and for the final contractions to give the
!     Fock matrix). This is not a problem if they are fully
!     buffered, but could impose an I/O bottleneck. Since the number
!     of coefficients scales linearly with system size, this should
!     be less of a problem than the quadratic scaling of full
!     Cholesky vectors.
!   - The max number of processes for which the algorithm can
!     possibly scale well equals the number of LDF atom pairs.
!   - Integral prescreening info will be set up here if not done
!     prior to calling this routine.
!   - Add is not implemented in parallel !!!
!
! ALGORITHM:
!=======================================================================
!
! If (Mode=3): FactC := 0.5*FactC
! Scale off-diagonal density blocks:
!              D[u_A v_B] <-- 2*D[u_A v_B] iff A != B
!
! Initialize V(J)=0 for all J
! *Parallel loop over atom pairs AB (A>=B):
!    - V(J) += sum_[u_A v_B] C(u_A v_B,J)*D(u_A v_B)
!      where J belongs to AB [1C and 2C funcs].
! *End parallel loop AB
! Add V over nodes.
!
! Initialize F(uv) for all significant uv.
! Initialize W(J) for all J.
!
! If (Mode=1 or Mode=3):
!    Compute 3-index contributions:
!    *Parallel loop over atom pairs AB (A>=B):
!       *Loop over atoms C
!          - Compute (u_A v_B|K) for K belonging to C [1C func]
!          - F(u_A v_B) += sum_K (u_A v_B|K)*V(K)
!          - W(K) += sum_[u_A v_B] (u_A v_B|K)*D(u_A v_B)
!       *End loop C
!       If (LDF2):
!          *loop over atom pairs CD
!             - Compute (u_A v_B|K) for K belonging to CD [2C func]
!             - F(u_A v_B) += sum_[K] (u_A v_B|K)*V(K)
!             - W(K) += sum_[u_A v_B] (u_A v_B|K)*D(u_A v_B)
!          *End loop over CD
!       End If
!    *End parallel loop AB
! End If
!
! Compute 2-index contributions:
! If (Mode=1 or Mode=2):
!    If (Mode=1):
!       - const=-1.0d0
!    Else:
!       - const=1.0d0
!    End If
!    *Parallel loop over atoms A:
!       *Loop over atoms B=1,A-1:
!          - Compute (J|K) for J in A and K in B [1C func]
!          - W(J) += const * sum_K (J|K)*V(K)
!          - W(K) += const * sum_J V(J)*(J|K)
!       *End loop B
!       - Compute (J|K) for J,K in A [1C func]
!       - W(J) += const* sum_K (J|K)*V(K)
!       If (LDF2):
!          *Loop over atom pairs CD (C>=D):
!             - Compute (J|K) for J in A [1C func], K in CD [2Cfunc]
!             - W(J) += const * sum_K (J | K)*V(K)
!             - W(K) += const * sum_J V(J)*(J|K)
!          *End loop CD
!       End If
!    *End parallel loop A
!    If (LDF2):
!       *Parallel loop over atom pairs AB (A>=B):
!          *Loop over atom pairs CD=1,AB-1:
!             - Compute (J|K) for J in AB and K in CD [2C func]
!             - W(J) += const * sum_K (J|K)*V(K)
!             - W(K) += const * sum_J V(J)*(J|K)
!          *End loop B
!          - Compute (J|K) for J,K in AB [2C func]
!          - W(J) += const* sum_K (J|K)*V(K)
!       *End loop AB
!    End If
! End If
! Add W over nodes
!
! *Parallel loop over atom pairs AB (A>=B):
!    - Read coefficient C(u_A v_B,J) for J in AB [1C,2C func]
!    - F(u_A v_B) += sum_J C(u_a v_B,J)*W(J)
! *End parallel loop AB
! Add F over nodes
!
! ALGORITHM DONE: Return F
!=======================================================================
!
! Note that arrays V, W, and F are stored locally as O(N) arrays,
! which should keep the communication bottleneck reasonable as the
! system size grows.
!
! Integral prescreening is based on the Cauchy-Schwarz inequality
! and involves estimation of the max integral in a given block.
! Integral prescreening is done at atom/atom pair block level as
! well as on shell level.
!
! Contribution prescreening is based on the submultiplicative
! property of the Frobenius norm in combination with the Cauchy-
! Schwarz inequality for positive (semi-) definite matrices.
! Specifically, an upper bound to the norm of the matrix-vector
! product
!
! y = A*x
!
! where A is a subblock of a positive (semi-) definite matrix,
! is given by
!
! ||y|| <= ||A||*||x||                  {submultiplicative property}
!        = sqrt[sum_ij A(i,j)**2]*sqrt[sum_j x(j)**2]
!       <= sqrt[sum_ij A(i,i)*A(j,j)]*sqrt[sum_j x(j)**2]  {Cauchy-Schwarz}
!        = sqrt[sum_i A(i,i)]*sqrt[sum_j A(j,j)]*sqrt[sum_j x(j)**2]
!
! This upper bound is used to avoid calculation of contributions
! below a given threshold. Contribution prescreening is done at
! atom/atom pair block level only, f.ex. for K on atom C
!
! ||sum_K (AB|K)*V(K)|| <= sqrt[sum_[u_A v_B] (u_A v_B | u_A v_B)]
!                         *sqrt[sum_K (K|K)]
!                         *sqrt[sum_K V(K)**2]

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif
implicit none
integer IntegralOption
logical Timing
integer Mode
real*8 ThrPS(2)
logical Add
logical PackedD
logical PackedF
integer nD
real*8 FactC(nD)
integer ip_D(nD)
integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

character*20 SecNam
parameter(SecNam='LDF_Fock_CoulombOnly')

logical LDF_IntegralPrescreeningInfoIsSet
external LDF_IntegralPrescreeningInfoIsSet

integer LDF_nAtom
external LDF_nAtom

#ifdef _DEBUGPRINT_
logical LDF_X_IsSet, LDF_TestBlockMatrix
external LDF_X_IsSet, LDF_TestBlockMatrix
logical DoTest
real*8 x, y
#endif

logical UseOldCode
logical IPI_set_here

real*8 tau(2)

integer i
integer nBas
integer iD, ip0, l
integer ip_DBlocks, l_DBlocks
integer ip_FBlocks, l_FBlocks
integer ip_VP, l_VP
integer ip_DNorm, l_DNorm
integer ip_CNorm, l_CNorm
integer ip_VNorm, l_VNorm
integer ip_FactC, l_FactC

real*8 tTotC1, tTotC2
real*8 tTotW1, tTotW2

! Start total timing.
if (Timing) then
  write(6,'(/,84A1)') ('-',i=1,84)
  call CWTime(tTotC1,tTotW1)
end if

UseOldCode = (IntegralOption == 111) .or. (IntegralOption == 222) .or. (IntegralOption == 333)
if (UseOldCode) then
  call WarningMessage(0,SecNam//': Using atom pair-driven (old) code!')
  call xFlush(6)
  call LDF_Fock_CoulombOnly0(IntegralOption,ThrPS(1),Mode,Add,PackedD,PackedF,nD,FactC,ip_D,ip_F)
  Go To 1 ! Return
end if

! Return if nothing to do
if (nD < 1) Go To 1

! Get number of basis functions (from localdf_bas.fh)
nBas = nBas_Valence
if (nBas < 1) then
  call WarningMessage(1,SecNam//': nBas<1 -- Fock matrix NOT computed!')
  write(6,'(A,I9)') 'nBas=',nBas
  call xFlush(6)
  Go To 1  ! return
end if

#ifdef _DEBUGPRINT_
if (.not. LDF_X_IsSet()) then
  call WarningMessage(2,SecNam//': LDF data not initialized!')
  call LDF_Quit(1)
end if
#endif

#ifdef _MOLCAS_MPP_
if (Add) then
  if ((nProcs > 1) .and. Is_Real_Par()) then
    write(6,'(A,A)') SecNam,': >>Add<< feature not implemented in parallel execution!'
    call LDF_NotImplemented()
  end if
end if
#endif

! For half-and-half integral representation, scale FactC by 1/2.
! Save a copy of FactC (restored at the end)
if (Mode == 3) then
  l_FactC = nD
  call GetMem('FactCBak','Allo','Real',ip_FactC,l_FactC)
  call dCopy_(nD,FactC,1,Work(ip_FactC),1)
  call dScal_(nD,0.5d0,FactC,1)
else
  l_FactC = 0
  ip_FactC = 0
end if

! Set prescreening info (if not already done)
if (.not. LDF_IntegralPrescreeningInfoIsSet()) then
  call LDF_SetIntegralPrescreeningInfo()
  IPI_set_here = .true.
else
  IPI_set_here = .false.
end if
tau(1) = max(ThrPS(1),0.0d0) ! integral prescreening threshold
tau(2) = max(ThrPS(2),0.0d0) ! contribution prescreening threshold

! Initialize Fock matrices (if not Add)
if (.not. Add) then
  if (PackedF) then
    l = nBas*(nBas+1)/2
  else
    l = nBas**2
  end if
  do iD=1,nD
    call Cho_dZero(Work(ip_F(iD)),l)
  end do
end if

! Allocate and set blocked density matrices (atom pair blocks)
l_DBlocks = nD
call GetMem('DBlk_P','Allo','Inte',ip_DBlocks,l_DBlocks)
ip0 = ip_DBlocks-1
#ifdef _DEBUGPRINT_
x = dble(NumberOfAtomPairs)
y = dble(LDF_nAtom())*(dble(LDF_nAtom())+1.0d0)/2.0d0
DoTest = int(x-y) == 0
#endif
do iD=1,nD
  call LDF_AllocateBlockMatrix('Den',iWork(ip0+iD))
  call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,iWork(ip0+iD))
#ifdef _DEBUGPRINT_
  if (DoTest) then
    if (.not. LDF_TestBlockMatrix(iWork(ip0+iD),PackedD,Work(ip_D(iD)))) then
      call WarningMessage(2,SecNam//': block matrix test failure')
      write(6,'(A,I4,A,I9,3X,A,L1)') 'Density matrix',iD,' at location',ip_D(iD),'Packed: ',PackedD
      call LDF_Quit(1)
    end if
  end if
#endif
  call LDF_ScaleOffdiagonalMatrixBlocks(iWork(ip0+iD),2.0d0)
end do

! Allocate and set blocked Fock matrices (atom pair blocks)
l_FBlocks = nD
call GetMem('FBlk_P','Allo','Inte',ip_FBlocks,l_FBlocks)
ip0 = ip_FBlocks-1
do iD=1,nD
  call LDF_AllocateBlockMatrix('Fck',iWork(ip0+iD))
  call LDF_Full2Blocked(Work(ip_F(iD)),PackedF,iWork(ip0+iD))
end do

! Allocate and compute Frobenius norm of blocked density matrices
l_DNorm = NumberOfAtomPairs*nD
call GetMem('DNorm','Allo','Real',ip_DNorm,l_DNorm)
ip0 = ip_DNorm
do iD=0,nD-1
  call LDF_BlockMatrixNorm(iWork(ip_DBlocks+iD),ip0)
  ip0 = ip0+NumberOfAtomPairs
end do

! Allocate Coulomb intermediates
l_VP = nD
call GetMem('VP','Allo','Inte',ip_VP,l_VP)
do iD=0,nD-1
  call LDF_AllocateAuxBasVector('CIn',iWork(ip_VP+iD))
end do

! Allocate array for norm of fitting coefficients
l_CNorm = 4*NumberOfAtomPairs
call GetMem('CNorm','Allo','Real',ip_CNorm,l_CNorm)

! Compute Coulomb intermediates
! V(J) = sum_uv C(uv,J)*D(uv)
! for each density matrix
! Compute Frobenius norm of fitting coefficients.
call LDF_ComputeCoulombIntermediates(Timing,nD,iWork(ip_DBlocks),iWork(ip_VP),ip_CNorm)

! Allocate and compute Frobenius norm of Coulomb intermediates
l_VNorm = (LDF_nAtom()+NumberOfAtomPairs)*nD
call GetMem('VNorm','Allo','Real',ip_VNorm,l_VNorm)
ip0 = ip_VNorm
do iD=0,nD-1
  call LDF_AuxBasVectorNorm(iWork(ip_VP+iD),ip0)
  ip0 = ip0+LDF_nAtom()+NumberOfAtomPairs
end do

! Compute Coulomb contributions
call LDF_Fock_CoulombOnly_(IntegralOption == 444,Timing,Mode,tau,nD,FactC,iWork(ip_DBlocks),iWork(ip_VP),iWork(ip_FBlocks), &
                           ip_CNorm,ip_DNorm,ip_VNorm)

! Get full storage (triangular or quadratic) Fock matrices from
! blocked ones.
ip0 = ip_FBlocks-1
do iD=1,nD
  call LDF_Blocked2Full(iWork(ip0+iD),PackedF,Work(ip_F(iD)))
end do

! Deallocation
call GetMem('VNorm','Free','Real',ip_VNorm,l_VNorm)
call GetMem('CNorm','Free','Real',ip_CNorm,l_CNorm)
do iD=0,nD-1
  call LDF_DeallocateAuxBasVector('CIn',iWork(ip_VP+iD))
end do
call GetMem('VP','Free','Inte',ip_VP,l_VP)
call GetMem('DNorm','Free','Real',ip_DNorm,l_DNorm)
do iD=0,nD-1
  call LDF_DeallocateBlockMatrix('Fck',iWork(ip_FBlocks+iD))
end do
call GetMem('FBlk_P','Free','Inte',ip_FBlocks,l_FBlocks)
do iD=0,nD-1
  call LDF_DeallocateBlockMatrix('Den',iWork(ip_DBlocks+iD))
end do
call GetMem('DBlk_P','Free','Inte',ip_DBlocks,l_DBlocks)

! Unset prescreening info (if set in this routine)
if (IPI_set_here) then
  call LDF_UnsetIntegralPrescreeningInfo()
end if

! Restore FactC
if (l_FactC > 0) then
  call dCopy_(nD,Work(ip_FactC),1,FactC,1)
  call GetMem('FactCBak','Free','Real',ip_FactC,l_FactC)
end if

1 continue
if (Timing) then
  call CWTime(tTotC2,tTotW2)
  write(6,'(A,A,A,2(1X,F12.2),A)') 'Total time spent in ',SecNam,':         ',tTotC2-tTotC1,tTotW2-tTotW1,' seconds'
  write(6,'(84A1)') ('-',i=1,84)
  call xFlush(6)
end if

end subroutine LDF_Fock_CoulombOnly
