!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ExRas(iCStart,nBaseQ,nBaseC,iQ_Atoms,nAtomsCC,Ax,Ay,Az,itristate,SmatRas,SmatPure,InCutOff,AOSum)

use qmstat_global, only: AvRed, BigT, c_orbene, Cordst, Cut_Ex1, Cut_Ex2, exrep2, iOrb, lExtr, lmax, MoAveRed, nCent, nPart, &
                         nRedMO, nState, outxyzRAS
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCStart, nBaseQ, nBaseC, iQ_Atoms, nAtomsCC, itristate
real(kind=wp), intent(out) :: Ax, Ay, Az, SmatRas(itristate), SmatPure(itristate)
logical(kind=iwp), intent(out) :: InCutOff
real(kind=wp), intent(inout) :: AOSum(nTri_Elem(nBaseQ))
integer(kind=iwp) :: i, ind, inwm, iS, js, kaunter, N, nDim1, nDimT, nGross, nHalf, nInsideCut
real(kind=wp) :: Addition, Cut_ExSq1, Cut_ExSq2, DH1, DH2, dist_sw, HighS, r2, r3, r3temp1, r3temp2
logical(kind=iwp) :: NearBy
real(kind=wp), allocatable :: ACC(:,:), ACCp(:,:), ACCt(:), ACCtp(:), AOAUX(:,:), AOAUXtri(:), AOint(:,:), AOG(:), AOMOOvl(:,:), &
                              AOMOOvlE(:,:), HalfE(:,:), HalfP(:), TEMP(:,:), V2(:,:)
logical(kind=iwp), allocatable :: Inside(:,:)
real(kind=wp), external :: Ddot_

!----------------------------------------------------------------------*
! Deduce how much the QM-molecule is translated from its position as   *
! defined in Seward.                                                   *
!----------------------------------------------------------------------*
Ax = Cordst(1,1)-outxyzRAS(1,1)
Ay = Cordst(2,1)-outxyzRAS(2,1)
Az = Cordst(3,1)-outxyzRAS(3,1)

! Make some initializations.

Cut_ExSq1 = Cut_Ex1**2
Cut_ExSq2 = Cut_Ex2**2
call mma_allocate(V2,nBaseC,iOrb(2),label='RotOrb')
call mma_allocate(AOint,nBaseQ,nBaseC,label='Sint')
nHalf = nBaseQ*iOrb(2)
nGross = nTri_Elem(nBaseQ)
call mma_allocate(HalfP,nBaseQ*iOrb(2),label='HalfPure')
!Jose****************************************************************
if (lExtr(8)) then
  call mma_allocate(AOMOOvl,nBaseQ,iOrb(2),label='qAOclMOOvl')
  call mma_allocate(AOMOOvlE,nBaseQ,iOrb(2),label='qAOclMOOvlE')
  call mma_allocate(AOAUX,nBaseQ,nBaseQ,label='AuxAOp')
  call mma_allocate(AOAUXtri,nGross,label='AuxAOpTri')
end if
!********************************************************************
InCutOff = .false.
SmatRas(:) = Zero
SmatPure(:) = Zero
nInsideCut = 0
if (MoAveRed) then
  nDim1 = nRedMO
  call mma_allocate(TEMP,nDim1,iOrb(2),label='TEMP')
else
  nDim1 = nBaseQ
end if
nDimT = nTri_Elem(nDim1)
call mma_allocate(HalfE,nDim1,iOrb(2),label='HalfOrbE')
call mma_allocate(AOG,nDimT,label='GammaAO')
call mma_allocate(ACC,nDim1,nDim1,label='Accumulate')
call mma_allocate(ACCp,nDim1,nDim1,label='AccumulateP')
call mma_allocate(ACCt,nDimT,label='AccumulateT')
call mma_allocate(ACCtp,nDimT,label='AccumulateTP')
ACC(:,:) = Zero
ACCp(:,:) = Zero

! Start loop over all solvent molecules.

call mma_allocate(Inside,iQ_Atoms,nAtomsCC,label='Inside')
do N=iCStart-1,nCent*(nPart-1),nCent

  ! Initialize

  dist_sw = huge(dist_sw)
  r3 = huge(r3)
  Inside(:,:) = .false.
  NearBy = .false.
  ! Loop over atoms.
  do inwm=1,iQ_Atoms
    r2 = (Cordst(1,N+1)-Cordst(1,inwm))**2+(Cordst(2,N+1)-Cordst(2,inwm))**2+(Cordst(3,N+1)-Cordst(3,inwm))**2
    dist_sw = min(dist_sw,r2)
    ! Distances for the inner cut-off. Also include the hydrogens.
    DH1 = (Cordst(1,N+2)-Cordst(1,inwm))**2+(Cordst(2,N+2)-Cordst(2,inwm))**2+(Cordst(3,N+2)-Cordst(3,inwm))**2
    DH2 = (Cordst(1,N+3)-Cordst(1,inwm))**2+(Cordst(2,N+3)-Cordst(2,inwm))**2+(Cordst(3,N+3)-Cordst(3,inwm))**2
    r3temp1 = min(DH1,DH2)
    r3temp2 = min(r3temp1,r2)
    r3 = min(r3,r3temp2)
    ! See if these atom-atom pairs inside cut-off.
    if (r2 < Cut_ExSq1) then
      Inside(inwm,1) = .true.
      NearBy = .true.
    end if
    if (DH1 < Cut_ExSq1) then
      Inside(inwm,2) = .true.
      NearBy = .true.
    end if
    if (DH2 < Cut_ExSq1) then
      Inside(inwm,3) = .true.
      NearBy = .true.
    end if
  end do

  ! Make some cut-off tests.

  if (.not. NearBy) cycle !If all distances larger than cut-off, jump to new solvent.
  ! Inner cut-off. Set flag to true then huge energy is added later. S*S matrix, however!
  if (r3 < Cut_ExSq2) InCutOff = .true.
  nInsideCut = nInsideCut+1

  ! Start integrating.
  call AOIntegrate(nBaseQ,nBaseC,Ax,Ay,Az,iQ_Atoms,nAtomsCC,AOint,V2,N,lmax,Inside)

  ! Transform overlaps from solvent AO to solvent MO.

  call Dgemm_('N','N',nBaseQ,iOrb(2),nBaseC,One,AOint,nBaseQ,V2,nBaseC,Zero,HalfP,nBaseQ)
  !Jose***************************************************************
  if (lExtr(8)) then
    call dcopy_(nHalf,HalfP,1,AOMOOvl,1)
    do i=1,iOrb(2)
      AOMOOvlE(:,i) = c_orbene(i)*AOMOOvl(:,i)
    end do
  end if
  !*******************************************************************

  ! If average natural orbital basis used, transform again. We also
  ! define Dim1 and Dim2 as dimensions of matrices depending on
  ! whether average natural orbitals have been used. This means that
  ! some vectors may be larger than necessary, but since no
  ! huge demand of memory is required, this should not cause problem.

  if (MoAveRed) then
    call Dgemm_('T','N',nRedMO,iOrb(2),nBaseQ,One,AvRed,nBaseQ,HalfP,nBaseQ,Zero,TEMP,nRedMO)
    call dcopy_(nDim1*iOrb(2),TEMP,1,HalfP,1)
  end if

  ! Hook on the orbital energy.

  do i=1,iOrb(2)
    ind = nDim1*(i-1)
    HalfE(:,i) = c_orbene(i)*HalfP(ind+1:ind+nDim1)
  end do

  ! Construct auxiliary matrix for the non-electrostatic operator
  ! in AO-basis for QM region. Also construct matrix of pure
  ! overlaps for the higher order terms.
  ! And accumulate.

  call Dgemm_('N','T',nDim1,nDim1,iOrb(2),One,HalfP,nDim1,HalfE,nDim1,One,ACC,nDim1)
  call Dgemm_('N','T',nDim1,nDim1,iOrb(2),One,HalfP,nDim1,HalfP,nDim1,One,ACCp,nDim1)

  !Jose*********************************
  if (lExtr(8)) then
    call Dgemm_('N','T',nBaseQ,nBaseQ,iOrb(2),exrep2,AOMOOvl,nBaseQ,AOMOOvlE,nBaseQ,Zero,AOAUX,nBaseQ)
    call SqToTri_Q(AOAUX,AOAUXtri,nBaseQ)
    AOSum(:) = AOSum+AOAUXTri
  end if
  !*************************************

  ! The end for this solvent molecule.

end do

! Now construct the matrix elements to the non-electrostatic
! operator in RASSI basis.

kaunter = 0
do iS=1,nState
  do jS=1,iS
    HighS = 0
    kaunter = kaunter+1
    ! Collect the relevant part of the transition density matrix.
    AOG(:) = BigT(1:nDimT,kaunter)
    ! Then transform according to theory.
    call SqToTri_Q(ACC,ACCt,nDim1)
    Addition = Ddot_(nDimT,AOG,1,ACCt,1)
    SmatRas(kaunter) = SmatRas(kaunter)+exrep2*Addition
    ! And include pure S*S for subsequent higher order overlap repulsion.
    call SqToTri_Q(ACCp,ACCtp,nDim1)
    HighS = Ddot_(nDimT,AOG,1,ACCtp,1)
    SmatPure(kaunter) = SmatPure(kaunter)+HighS
  end do
end do

! Deallocations.

call mma_deallocate(Inside)

call mma_deallocate(V2)
call mma_deallocate(AOint)
call mma_deallocate(HalfP)
call mma_deallocate(HalfE)
call mma_deallocate(AOG)
call mma_deallocate(ACC)
call mma_deallocate(ACCp)
call mma_deallocate(ACCt)
call mma_deallocate(ACCtp)
if (MoAveRed) call mma_deallocate(TEMP)
!Jose****************************************************************
if (lExtr(8)) then
  call mma_deallocate(AOMOOvl)
  call mma_deallocate(AOMOOvlE)
  call mma_deallocate(AOAUX)
  call mma_deallocate(AOAUXtri)
end if
!********************************************************************

return

end subroutine ExRas
