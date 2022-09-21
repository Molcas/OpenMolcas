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

subroutine LDF_Fock_CoulombOnly0(IntegralOption,tau,Mode,Add,PackedD,PackedF,nD,FactC,ip_D,F)
! Thomas Bondo Pedersen, September 2010.
!
!          OLD CODE
!
! Purpose: Compute Coulomb contribution to Fock matrix using local
!          Density Fitting coefficients. Use only for debugging!
!          Poor performance!
!
! Args:    - IntegralOption (integer)
!              111  => use two-electron integrals computed from
!                      LDF coefficients (for debugging)
!              222  => use conventional two-electron integrals
!                      (for debugging)
!              333  => use conventional or LDF two-electron
!                      integrals depending on positivity of the
!                      latter (for debugging)
!              all other values => use production-level LDF code
!              (INPUT)
!          - tau (real*8) integral prescreening threshold (INPUT).
!            Only used with IntegralOption=111.
!          - Mode (integer) 1:robust,2:nonrobust,3:half-and-half
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
!          - F(*,nD): nD Fock matrices (if PackedF:
!            lower triangular storage) (INPUT/OUTPUT)
!
!          If (Add):  [NOT IMPLEMENTED IN PARALLEL]
!   (1)       F(uv) = F(uv) + FactC * sum_kl (uv|kl)*D(kl)
!          Else:
!   (2)       F(uv) = FactC * sum_kl (uv|kl)*D(kl)
!
!          where the integrals are given by the robust LDF
!          representation
!
!          (uv|kl) = ([uv]|kl) + (uv|[kl]) - ([uv]|[kl])
!
!          or the non-robust representation
!
!          (uv|kl) = ([uv]|[kl])
!
!          or the half-and-half representation
!
!          (uv|kl) = 0.5*([uv]|kl) + 0.5*(uv|[kl])
!
!          and the fitted products are given by
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
!   - The algorithm is NOT integral driven. This means that 3-center
!     integrals (4-center integrals if 2-center fitting functions
!     are included) are computed several times, causing quite poor
!     performance.
!   - The max number of processes for which the algorithm can
!     possibly scale well equals the number of LDF atom pairs.
!   - No screening has been implemented.
!   - Add is not implemented in parallel !!!
!
! Outline of the algorithm (robust representation):
! =================================================
!
! *If (Mode=3): FactC := 0.5*FactC
!
! *Scale off-diagonal density blocks:
!  D[u_A v_B] <-- 2*D[u_A v_B] iff A != B
!
! *Parallel loop over atom pairs AB (A>=B):
!    - V(J_AB) = sum_[u_A v_B] C(u_A v_B,J_AB)*D(u_A v_B)
! *End loop AB
! Add V over nodes.
!
! *Parallel loop over atom pairs AB (A>=B):
!    *Loop over atom pairs CD (C>=D):
!       - F(u_A v_B) = F(u_A v_B)
!                    + sum_[K_CD] (u_A v_B | K_CD)*V(K_CD)
!       - W(J_AB) = W(J_AB)
!                 + sum_[k_C l_D] (J_AB | k_C l_D)*D(k_C l_D)
!       - W(J_AB) = W(J_AB)
!                 - sum_[K_CD] (J_AB | K_CD)*V(K_CD)
!    *End loop CD
!    - F(u_A v_B) = F(u_A v_B)
!                 + sum_[J_AB] C(u_A v_B,J_AB)*W(J_AB)
! *End loope AB
! Add F over nodes
!
! Note that both arrays V and F are stored locally as O(N) arrays,
! which should keep the communication bottleneck at a minimum as the
! system size grows.

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
#ifdef _DEBUGPRINT_
use Constants, only: One
#endif
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IntegralOption, Mode, nD, ip_D(nD)
real(kind=wp), intent(in) :: tau
logical(kind=iwp), intent(in) :: Add, PackedD, PackedF
real(kind=wp), intent(inout) :: FactC(nD), F(*)
logical(kind=iwp) :: UsePartPermSym
integer(kind=iwp) :: nBas, iD, l, l_DBlocks, l_FBlocks, l_VBlocks, l_FactC
integer(kind=iwp), allocatable :: DBlocks(:), FBlocks(:), VBlocks(:)
real(kind=wp), allocatable :: FactCBak(:)
character(len=*), parameter :: SecNam = 'LDF_Fock_CoulombOnly0'
#ifdef _DEBUGPRINT_
logical(kind=iwp) :: DoTest
real(kind=wp) :: x, y
logical(kind=iwp), external :: LDF_X_IsSet, LDF_TestBlockMatrix
integer(kind=iwp), external :: LDF_nAtom
#endif
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

! Return if nothing to do
if (nD < 1) return

! Get number of basis functions (from localdf_bas.fh)
nBas = nBas_Valence
if (nBas < 1) then
  call WarningMessage(1,SecNam//': nBas<1 -- Fock matrix NOT computed!')
  write(u6,'(A,I9)') 'nBas=',nBas
  call xFlush(u6)
  return
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
    write(u6,'(A,A)') SecNam,': >>Add<< feature not implemented in parallel execution!'
    call LDF_NotImplemented()
  end if
end if
#endif

! For half-and-half integral representation, scale FactC by 1/2.
! Save a copy of FactC (restored at the end)
if (Mode == 3) then
  l_FactC = nD
  call mma_allocate(FactCBak,l_FactC,label='FactCBak')
  call dCopy_(nD,FactC,1,FactCBak,1)
  call dScal_(nD,Half,FactC,1)
else
  l_FactC = 0
end if

! Initialize Fock matrices (if not Add)
if (PackedF) then
  l = nBas*(nBas+1)/2
else
  l = nBas**2
end if
if (.not. Add) then
  F(1:nD*l) = Zero
end if

! Allocate and set blocked density matrices (atom pair blocks)
l_DBlocks = nD
call mma_allocate(DBlocks,l_DBlocks,label='DBlk_P')
#ifdef _DEBUGPRINT_
x = real(NumberOfAtomPairs,kind=wp)
y = real(LDF_nAtom(),kind=wp)*(real(LDF_nAtom(),kind=wp)+One)*Half
DoTest = int(x-y) == 0
#endif
do iD=1,nD
  call LDF_AllocateBlockMatrix('Den',DBlocks(iD))
  call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,DBlocks(iD))
# ifdef _DEBUGPRINT_
  if (DoTest) then
    if (.not. LDF_TestBlockMatrix(DBlocks(iD),PackedD,Work(ip_D(iD)))) then
      call WarningMessage(2,SecNam//': block matrix test failure')
      write(u6,'(A,I4,A,I9,3X,A,L1)') 'Density matrix',iD,' at location',ip_D(iD),'Packed: ',PackedD
      call LDF_Quit(1)
    end if
  end if
# endif
  call LDF_ScaleOffdiagonalMatrixBlocks(DBlocks(iD),Two)
end do

! Allocate and set blocked Fock matrices (atom pair blocks)
l_FBlocks = nD
call mma_allocate(FBlocks,l_FBlocks,label='FBlk_P')
do iD=1,nD
  call LDF_AllocateBlockMatrix('Fck',FBlocks(iD))
  call LDF_Full2Blocked(F((iD-1)*l+1),PackedF,FBlocks(iD))
end do

if (IntegralOption == 111) then
  ! Compute Fock matrix using LDF integrals (debug)
  call WarningMessage(0,SecNam//': Using integrals from LDF coefficients!')
  call xFlush(u6)
  UsePartPermSym = .true.
  if (Mode == 3) then
    call LDF_FVIFC(UsePartPermSym,Mode,max(tau,Zero),nD,FactCBak,DBlocks,FBlocks)
  else
    call LDF_FVIFC(UsePartPermSym,Mode,max(tau,Zero),nD,FactC,DBlocks,FBlocks)
  end if
else if (IntegralOption == 222) then
  ! Compute Fock matrix using conventional integrals (debug)
  call WarningMessage(0,SecNam//': Using conventional integrals!')
  call xFlush(u6)
  UsePartPermSym = .true.
  call LDF_FCI(UsePartPermSym,nD,FactC,DBlocks,FBlocks)
else if (IntegralOption == 333) then
  ! Compute Fock matrix using conventional or LDF integrals
  ! depending on positivity of the latter (debug)
  call WarningMessage(0,SecNam//': Using PSD (LDF or conv.) integrals!')
  call xFlush(u6)
  UsePartPermSym = .true.
  if (Mode == 3) then
    call LDF_FTst(UsePartPermSym,Mode,max(tau,Zero),nD,FactCBak,DBlocks,FBlocks)
  else
    call LDF_FTst(UsePartPermSym,Mode,max(tau,Zero),nD,FactC,DBlocks,FBlocks)
  end if
else
  ! Allocate Coulomb intermediates
  l_VBlocks = nD
  call mma_allocate(VBlocks,l_VBlocks,label='VBlk_P')
  do iD=1,nD
    call LDF_AllocateBlockVector('CIn',VBlocks(iD))
  end do
  ! Compute Coulomb intermediates,
  ! V(J) = sum_uv C(uv,J)*D(uv)
  ! for each density matrix
  call LDF_ComputeCoulombIntermediates0(nD,DBlocks,VBlocks)
  ! Compute Coulomb contributions
  call LDF_Fock_CoulombOnly0_(Mode,nD,FactC,DBlocks,VBlocks,FBlocks)
  ! Deallocate Coulomb intermediates
  do iD=1,nD
    call LDF_DeallocateBlockVector('CIn',VBlocks(iD))
  end do
  call mma_deallocate(VBlocks)
end if

! Get full storage (triangular or quadratic) Fock matrices from
! blocked ones.
do iD=1,nD
  call LDF_Blocked2Full(FBlocks(iD),PackedF,F((iD-1)*l+1))
end do

! Restore FactC
if (l_FactC > 0) then
  call dCopy_(nD,FactCBak,1,FactC,1)
  call mma_deallocate(FactCBak)
end if

! Deallocation
do iD=1,nD
  call LDF_DeallocateBlockMatrix('Fck',FBlocks(iD))
end do
call mma_deallocate(FBlocks)
do iD=1,nD
  call LDF_DeallocateBlockMatrix('Den',DBlocks(iD))
end do
call mma_deallocate(DBlocks)

end subroutine LDF_Fock_Coulombonly0
