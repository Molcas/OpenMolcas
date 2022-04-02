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

subroutine ExScf(iCStart,nBaseQ,nBaseC,iQ_Atoms,nAtomsCC,Ax,Ay,Az,itri,Smat,SmatPure,InCutOff,AOSum)

use qmstat_global, only: c_orbene, Cordst, Cut_Ex1, Cut_Ex2, exrep2, iOrb, iPrint, lExtr, lmax, nCent, nPart, outxyz, V1
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCStart, nBaseQ, nBaseC, iQ_Atoms, nAtomsCC, itri
real(kind=wp), intent(out) :: Ax, Ay, Az, Smat(itri), SmatPure(itri)
logical(kind=iwp), intent(out) :: InCutOff
real(kind=wp), intent(inout) :: AOSum(*)
integer(kind=iwp) :: i, iAOAOTri, inwm, j, N, nInsideCut
real(kind=wp) :: Cut_ExSq1, Cut_ExSq2, DH1, DH2, dist_sw, r2, r3, r3temp1, r3temp2
logical(kind=iwp) :: NearBy
real(kind=wp), allocatable :: AOAUX(:,:), AOAUXtri(:), AOint(:,:), AOMOOvl(:,:), AOMOOvlE(:,:), AUX(:,:), AUXp(:,:), AUXtri(:), &
                              Inte(:,:), OvlMO(:,:), OvlMOE(:,:), V2(:,:)
logical(kind=iwp), allocatable :: Inside(:,:)

!----------------------------------------------------------------------*
! Deduce how much the QM-molecule is translated from its position as   *
! defined in Seward.                                                   *
!----------------------------------------------------------------------*
Ax = Cordst(1,1)-outxyz(1,1)
Ay = Cordst(2,1)-outxyz(2,1)
Az = Cordst(3,1)-outxyz(3,1)
!----------------------------------------------------------------------*
! Rotate solvent orbitals and make AO integration.                     *
!----------------------------------------------------------------------*
Cut_ExSq1 = Cut_Ex1**2
Cut_ExSq2 = Cut_Ex2**2
call mma_allocate(V2,nBaseC,iOrb(2),label='RotOrb')
call mma_allocate(AOint,nBaseQ,nBaseC,label='Sint')
call mma_allocate(OvlMO,iOrb(1),iOrb(2),label='OvlMO')
call mma_allocate(Inte,iOrb(1),nBaseC,label='Intermed')
call mma_allocate(OvlMOE,iOrb(1),iOrb(2),label='OvlMOene')
call mma_allocate(AUX,iOrb(1),iOrb(1),label='AUX')
call mma_allocate(AUXp,iOrb(1),iOrb(1),label='AUXp')
call mma_allocate(AUXtri,iTri,label='AUXtri')
!Jose****************************************************************
if (lExtr(8)) then
  iAOAOTri = nTri_Elem(nBaseQ)
  call mma_allocate(AOMOOvl,nBaseQ,iOrb(2),label='qAOclMOOvl')
  call mma_allocate(AOMOOvlE,nBaseQ,iOrb(2),label='qAOclMOOvlE')
  call mma_allocate(AOAUX,nBaseQ,nBaseQ,label='AuxAOp')
  call mma_allocate(AOAUXtri,iAOAOTri,label='AuxAOpTri')
end if
!********************************************************************
InCutOff = .false.
Smat(:) = Zero
SmatPure(:) = Zero
nInsideCut = 0
call mma_allocate(Inside,iQ_Atoms,nAtomsCC,label='Inside')
do N=iCStart-1,nCent*(nPart-1),nCent

  ! Initialize.

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
    ! Check if this atom-atom pair inside
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

  ! Now make the cut-off test.

  if (.not. NearBy) cycle
  ! Inner cut-off.
  if (r3 < Cut_ExSq2) InCutOff = .true.
  nInsideCut = nInsideCut+1

  ! Make the AO-AO overlap integration.

  call AOIntegrate(nBaseQ,nBaseC,Ax,Ay,Az,iQ_Atoms,nAtomsCC,AOint,V2,N,lmax,Inside)

  ! Transform to MO-MO overlap.

  call Dgemm_('T','N',iOrb(1),nBaseC,nBaseQ,One,V1,nBaseQ,AOint,nBaseQ,Zero,Inte,iOrb(1))
  call Dgemm_('N','N',iOrb(1),iOrb(2),nBaseC,One,Inte,iOrb(1),V2,nBaseC,Zero,OvlMO,iOrb(1))
  do i=1,iOrb(2)
    OvlMOE(:,i) = c_orbene(i)*OvlMO(:,i)
  end do

  !**Jose
  if (lExtr(8)) then
    call Dgemm_('N','N',nBaseQ,iOrb(2),nBaseC,One,AOint,nBaseQ,V2,nBaseC,Zero,AOMOOvl,nBaseQ)
    do i=1,iOrb(2)
      AOMOOvlE(:,i) = c_orbene(i)*AOMOOvl(:,i)
    end do
  end if
  !*******************

  ! If you are interested, print some bla bla bla.

  if (iPrint >= 29) then
    write(u6,*)
    write(u6,*) 'OVERLAP BETWEEN QM-SYSTEM AND SOLVENT MOLECULE ',N/nCent
    write(u6,*) 'QM-MO  SOLV-MO  OVERLAP'
    do i=1,iOrb(1)
      do j=1,iOrb(2)
        write(u6,8888) i,j,OvlMO(i,j)
      end do
    end do
  end if

  ! Construct the perturbation and accumulate pure overlap for
  ! subsequent construction of higher order term.

  call Dgemm_('N','T',iOrb(1),iOrb(1),iOrb(2),exrep2,OvlMO,iOrb(1),OvlMOE,iOrb(1),Zero,AUX,iOrb(1))
  call SqToTri_Q(AUX,AUXtri,iOrb(1))
  Smat(:) = Smat+AUXtri
  call Dgemm_('N','T',iOrb(1),iOrb(1),iOrb(2),One,OvlMO,iOrb(1),OvlMO,iOrb(1),Zero,AUXp,iOrb(1))
  call SqToTri_Q(AUXp,AUXtri,iOrb(1))
  SmatPure(:) = SmatPure+AUXtri

  !Jose*********************************
  if (lExtr(8)) then
    call Dgemm_('N','T',nBaseQ,nBaseQ,iOrb(2),exrep2,AOMOOvl,nBaseQ,AOMOOvlE,nBaseQ,Zero,AOAUX,nBaseQ)
    call SqToTri_Q(AOAUX,AOAUXtri,nBaseQ)
    AOSum(1:iAOAOTri) = AOSum(1:iAOAOTri)+AOAUXtri
  end if
  !*************************************

  ! This solvent molecule ends now!

end do

call mma_deallocate(Inside)

call mma_deallocate(V2)
call mma_deallocate(AOint)
call mma_deallocate(OvlMO)
call mma_deallocate(Inte)
call mma_deallocate(OvlMOE)
call mma_deallocate(AUX)
call mma_deallocate(AUXp)
call mma_deallocate(AUXtri)
!Jose****************************************************************
if (lExtr(8)) then
  call mma_deallocate(AOMOOvl)
  call mma_deallocate(AOMOOvlE)
  call mma_deallocate(AOAUX)
  call mma_deallocate(AOAUXtri)
end if
!********************************************************************

return

8888 format(I3,'    ',I3,'       ',F12.10)

end subroutine ExScf
