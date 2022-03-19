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

subroutine ExRas(iCStart,nBaseQ,nBaseC,nCnC_C,iQ_Atoms,nAtomsCC,Ax,Ay,Az,itristate,SmatRas,SmatPure,InCutOff,ipAOSum)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, r8

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iCStart, nBaseQ, nBaseC, nCnC_C(MxBasC), iQ_Atoms, nAtomsCC, itristate, ipAOSum
real(kind=wp) :: Ax, Ay, Az, SmatRas(itristate), SmatPure(itristate)
integer(kind=iwp) :: i, iAOMOOvl, iAOMOOvlE, iHalf, iHalfE, iHalfpar, ind, inwm, ipACC, ipACCp, ipACCt, ipACCtp, ipAOAUX, &
                     ipAOAUXtri, ipAOG, ipAOint, ipAOintpar, ipAux, ipAuxp, iS, iTEMP, iV2, js, k, kaunter, N, nAObaseSize, nDim1, &
                     nDim2, nDimP, nDimT, nGross, nHalf, nInsideCut, nV2size
real(kind=wp) :: Addition, CorTemp(3), Cut_ExSq1, Cut_ExSq2, DH1, DH2, dist_sw, HighS, r2, r3, r3temp1, r3temp2
logical(kind=iwp) :: InCutOff
logical(kind=iwp) :: NearBy
logical(kind=iwp), allocatable :: Inside(:,:)
real(kind=r8), external :: Ddot_

!----------------------------------------------------------------------*
! Deduce how much the QM-molecule is translated from its position as   *
! defined in Seward.                                                   *
!----------------------------------------------------------------------*
Ax = Cordst(1,1)-outxyzRAS(1,1)
Ay = Cordst(1,2)-outxyzRAS(1,2)
Az = Cordst(1,3)-outxyzRAS(1,3)

! Make some initializations.

Cut_ExSq1 = Cut_Ex1**2
Cut_ExSq2 = Cut_Ex2**2
nV2size = iOrb(2)*nBaseC
nAObaseSize = nBaseQ*nBaseC
call GetMem('RotOrb','Allo','Real',iV2,nV2size)
call GetMem('Sint','Allo','Real',ipAOint,nAObaseSize)
call GetMem('Sintpar','Allo','Real',ipAOintpar,nAObaseSize)
nHalf = nBaseQ*iOrb(2)
nGross = nTri_Elem(nBaseQ)
call GetMem('HalfTrans','Allo','Real',iHalfpar,nHalf)
call GetMem('HalfPure','Allo','Real',iHalf,nHalf)
call GetMem('HalfOrbE','Allo','Real',iHalfE,nHalf)
call GetMem('Auxiliary','Allo','Real',ipAux,nBaseQ**2)
call GetMem('AuxiliaryP','Allo','Real',ipAuxp,nBaseQ**2)
call GetMem('GammaAO','Allo','Real',ipAOG,nGross)
call GetMem('Accumulate','Allo','Real',ipACC,nBaseQ**2)
call GetMem('Accumulate','Allo','Real',ipACCp,nBaseQ**2)
call GetMem('AccumulateT','Allo','Real',ipACCt,nGross)
call GetMem('AccumulateTP','Allo','Real',ipACCtp,nGross)
call GetMem('TEMP','Allo','Real',iTEMP,nRedMO*iOrb(2))
call dcopy_(nBaseQ**2,[Zero],0,Work(ipACC),1)
call dcopy_(nBaseQ**2,[Zero],0,Work(ipACCp),1)
!Jose****************************************************************
if (lExtr(8)) then
  call GetMem('qAOclMOOvl','Allo','Real',iAOMOOvl,nHalf)
  call GetMem('qAOclMOOvlE','Allo','Real',iAOMOOvlE,nHalf)
  call GetMem('AuxAOp','Allo','Real',ipAOAUX,nBaseQ**2)
  call GetMem('AuxAOpTri','Allo','Real',ipAOAUXtri,nGross)
end if
!********************************************************************
InCutOff = .false.
do i=1,iTriState
  SmatRas(i) = 0
  SmatPure(i) = 0
end do
nInsideCut = 0
if (MoAveRed) then
  nDim1 = nRedMO
  nDim2 = iOrb(2)
else
  nDim1 = nBaseQ
  nDim2 = iOrb(2)
end if
nDimP = nDim1*nDim2
nDimT = nTri_Elem(nDim1)

! Start loop over all solvent molecules.

call mma_allocate(Inside,MxAt,3,label='Inside')
do N=iCStart-1,nCent*(nPart-1),nCent

  ! Initialize

  dist_sw = huge(dist_sw)
  r3 = huge(r3)
  do i=1,MxAt
    Inside(i,1) = .false.
    Inside(i,2) = .false.
    Inside(i,3) = .false.
  end do
  NearBy = .false.
  ! Loop over atoms.
  do inwm=1,iQ_Atoms
    do k=1,3
      CorTemp(k) = (Cordst(N+1,k)-Cordst(inwm,k))**2
    end do
    r2 = CorTemp(1)+CorTemp(2)+CorTemp(3)
    dist_sw = min(dist_sw,r2)
    DH1 = Zero !Distances for the inner cut-off. Also include the hydrogens.
    DH2 = Zero
    do k=1,3
      DH1 = DH1+(Cordst(N+2,k)-Cordst(inwm,k))**2
      DH2 = DH2+(Cordst(N+3,k)-Cordst(inwm,k))**2
    end do
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
  if (r3 < Cut_ExSq2) then !Inner cut-off. Set flag to true then huge energy is added later. S*S matrix, however!
    InCutOff = .true.
  end if
  nInsideCut = nInsideCut+1

  ! Start integrating.
  call AOIntegrate(iCStart,nBaseQ,nBaseC,Ax,Ay,Az,nCnC_C,iQ_Atoms,nAtomsCC,ipAOint,ipAOintpar,iV2,N,lmax,Inside)

  ! Transform overlaps from solvent AO to solvent MO.

  call Dgemm_('N','N',nBaseQ,iOrb(2),nBaseC,One,Work(ipAOint),nBaseQ,Work(iV2),nBaseC,Zero,Work(iHalf),nBaseQ)
  !Jose***************************************************************
  if (lExtr(8)) then
    call dcopy_(nHalf,Work(iHalf),1,Work(iAOMOOvl),1)
    call dcopy_(nHalf,[Zero],0,Work(iAOMOOvlE),1)
    do k=1,iOrb(2)
      ind = nBaseQ*(k-1)
      call DaxPy_(nBaseQ,c_orbene(k),Work(iAOMOOvl+ind),1,Work(iAOMOOvlE+ind),1)
    end do
  end if
  !*******************************************************************

  ! If average natural orbital basis used, transform again. We also
  ! define Dim1 and Dim2 as dimensions of matrices depending on
  ! whether average natural orbitals have been used. This means that
  ! some vectors may be larger than necessary, but since no
  ! huge demand of memory is required, this should not cause problem.

  if (MoAveRed) then
    call Dgemm_('T','N',nRedMO,iOrb(2),nBaseQ,One,Work(ipAvRed),nBaseQ,Work(iHalf),nBaseQ,Zero,Work(iTEMP),nRedMO)
    call dcopy_(nDim1*nDim2,Work(iTEMP),1,Work(iHalf),1)
  end if

  ! Hook on the orbital energy.

  call dcopy_(nDimP,[Zero],0,Work(iHalfE),1)
  do k=1,iOrb(2)
    ind = nDim1*(k-1)
    call DaxPy_(nDim1,c_orbene(k),Work(iHalf+ind),1,Work(iHalfE+ind),1)
  end do

  ! Construct auxiliary matrix for the non-electrostatic operator
  ! in AO-basis for QM region. Also construct matrix of pure
  ! overlaps for the higher order terms.

  call Dgemm_('N','T',nDim1,nDim1,iOrb(2),One,Work(iHalf),nDim1,Work(iHalfE),nDim1,Zero,Work(ipAUX),nDim1)
  call Dgemm_('N','T',nDim1,nDim1,iOrb(2),One,Work(iHalf),nDim1,Work(iHalf),nDim1,Zero,Work(ipAUXp),nDim1)

  ! Accumulate.

  call DaxPy_(nDim1**2,One,Work(ipAUX),1,Work(ipACC),1)
  call DaxPy_(nDim1**2,One,Work(ipAUXp),1,Work(ipACCp),1)

  !Jose*********************************
  if (lExtr(8)) then
    call Dgemm_('N','T',nBaseQ,nBaseQ,iOrb(2),exrep2,Work(iAOMOOvl),nBaseQ,Work(iAOMOOvlE),nBaseQ,Zero,Work(ipAOAUX),nBaseQ)
    call SqToTri_Q(Work(ipAOAUX),Work(ipAOAUXtri),nBaseQ)
    call DaxPy_(nGross,One,Work(ipAOAUXtri),1,Work(ipAOSum),1)
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
    call dCopy_(nDimT,Work(iBigT+nDimT*(kaunter-1)),1,Work(ipAOG),1)
    ! Then transform according to theory.
    call SqToTri_Q(Work(ipACC),Work(ipACCt),nDim1)
    Addition = Ddot_(nDimT,Work(ipAOG),1,Work(ipACCt),1)
    SmatRas(kaunter) = SmatRas(kaunter)+exrep2*Addition
    ! And include pure S*S for subsequent higher order overlap repulsion.
    call SqToTri_Q(Work(ipACCp),Work(ipACCtp),nDim1)
    HighS = Ddot_(nDimT,Work(ipAOG),1,Work(ipACCtp),1)
    SmatPure(kaunter) = SmatPure(kaunter)+HighS
  end do
end do

! Deallocations.

call mma_deallocate(Inside)

call GetMem('RotOrb','Free','Real',iV2,nV2size)
call GetMem('Sint','Free','Real',ipAOint,nAObaseSize)
call GetMem('Sintpar','Free','Real',ipAOintpar,nAObaseSize)
call GetMem('HalfTrans','Free','Real',iHalfpar,nHalf)
call GetMem('HalfPure','Free','Real',iHalf,nHalf)
call GetMem('HalfOrbE','Free','Real',iHalfE,nHalf)
call GetMem('Auxiliary','Free','Real',ipAux,nBaseQ**2)
call GetMem('AuxiliaryP','Free','Real',ipAuxp,nBaseQ**2)
call GetMem('GammaAO','Free','Real',ipAOG,nGross)
call GetMem('Accumulate','Free','Real',ipACC,nBaseQ**2)
call GetMem('Accumulate','Free','Real',ipACCp,nBaseQ**2)
call GetMem('AccumulateT','Free','Real',ipACCt,nGross)
call GetMem('AccumulateTP','Free','Real',ipACCtp,nGross)
call GetMem('TEMP','Free','Real',iTEMP,nRedMO*iOrb(2))
!Jose****************************************************************
if (lExtr(8)) then
  call GetMem('qAOclMOOvl','Free','Real',iAOMOOvl,nHalf)
  call GetMem('qAOclMOOvlE','Free','Real',iAOMOOvlE,nHalf)
  call GetMem('AuxAOp','Free','Real',ipAOAUX,nBaseQ**2)
  call GetMem('AuxAOpTri','Free','Real',ipAOAUXtri,nGross)
end if
!********************************************************************

return

end subroutine ExRas
