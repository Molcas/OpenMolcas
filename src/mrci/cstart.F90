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

subroutine CSTART(AREF,EREF,CI,ICI)

use mrci_global, only: ESHIFT, GFAC, IAD25S, ICPF, IDFREE, IDISKC, IDISKD, IDISKS, IREFX, IREST, IROOT, Lu_25, LUEIG, LUREST, &
                       MBUF, MXVEC, NCONF, NNEW, NREF, NRROOT, NSTOT, NVTOT, POTNUC
use guga_util_global, only: nCOP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: AREF(NREF,NREF), EREF(NREF)
real(kind=wp), intent(out) :: CI(NCONF)
integer(kind=iwp), intent(out) :: ICI(MBUF)
#include "Molcas.fh"
integer(kind=iwp) :: I, I1, I2, IAD25, ID, IR, IREF, ISTA, NN
real(kind=wp) :: GINV
integer(kind=iwp), allocatable :: ISTART(:)
real(kind=wp), allocatable :: Buf(:)

do I=1,MXVEC
  IDISKC(I) = -1
  IDISKS(I) = -1
end do
! FIRST, USE THE CI ARRAY TO STORE THE DIAGONAL ELEMENTS:
call mma_allocate(Buf,nCOP,label='Buf')
IAD25 = IAD25S
do I=1,NCONF,nCOP
  call dDAFILE(Lu_25,2,Buf,nCOP,IAD25)
  NN = min(nCOP,NCONF+1-I)
  CI(I:I+NN-1) = Buf(1:NN)
end do
call mma_deallocate(Buf)
! THESE ARE DIAGONAL ELEMENTS OF THE ELECTRONIC HAMILTONIAN.
! POTNUC SHOULD BE ADDED. IN ADDITION, WE USE AN ENERGY SHIFT.
! NOTE: DISPLACEMENT 1.0e-4 PROTECTS AGAINST DIVIDE ERRORS.
!PAM: Protect all diag elems. Needed in some weird cases.
! ENERGY SHIFT:
ESHIFT = EREF(1)
do I=1,NCONF
  CI(I) = CI(I)+POTNUC-ESHIFT+1.0e-4_wp
end do
call Add_Info('CI_DIAG2',CI(2),1,8)
! REPLACE REFERENCE ENERGIES:
do I=1,NREF
  IR = IREFX(I)
  CI(IR) = EREF(I)-ESHIFT-1.0e-4_wp
end do
if (ICPF == 1) then
  do IREF=1,NREF
    IR = IREFX(IREF)
    CI(IR) = GFAC*CI(IR)
  end do
  GINV = One/GFAC
  CI(:) = GINV*CI
end if
IDFREE = 0
IDISKD = 0
do ISTA=1,NCONF,MBUF
  NN = min(MBUF,(NCONF+1-ISTA))
  call dDAFILE(LUEIG,1,CI(ISTA),NN,IDFREE)
end do
! THEN, SET UP START CI VECTORS IN MCSF BASIS:
CI(:) = Zero
if (IREST == 0) then
  call mma_allocate(ISTART,MXROOT,label='ISTART')
  NNEW = IROOT(NRROOT)
  I1 = 1
  I2 = 1
  do I=1,NNEW
    if (I == IROOT(I1)) then
      ISTART(NNEW-NRROOT+I1) = I
      I1 = I1+1
    else
      ISTART(I2) = I
      I2 = I2+1
    end if
  end do
  if (NNEW > 1) then
    write(u6,*) ' THE FOLLOWING REFERENCE ROOTS ARE USED AS START VECTORS:'
    write(u6,'(12(A,I2))') ' ROOTS NR ',ISTART(1),(',',ISTART(I),I=2,NNEW-1),', AND ',ISTART(NNEW)
    if (NNEW > NRROOT) then
      write(u6,*) ' (THE FIRST EXTRA ROOT(S) WERE INCLUDED IN ORDER TO IMPROVE CONVERGENCE)'
    end if
  else
    write(u6,'(A,I2,A)') ' ROOT NR ',ISTART(1),' IS USED AS START VECTOR.'
  end if
  do I=1,NNEW
    ISTA = ISTART(I)
    IR = IREFX(ISTA)
    CI(IR) = One
    IDISKC(I) = IDFREE
    do ISTA=1,NCONF,MBUF
      NN = min(MBUF,(NCONF+1-ISTA))
      call PKVEC(NN,CI(ISTA),ICI)
      call iDAFILE(LUEIG,1,ICI,NN,IDFREE)
    end do
    CI(IR) = Zero
  end do
  call mma_deallocate(ISTART)
else
  ID = 0
  NNEW = NRROOT
  do I=1,NRROOT
    call dDAFILE(LUREST,2,CI,NCONF,ID)
    call CSFTRA('MCSF',CI,AREF)
    IDISKC(I) = IDFREE
    do ISTA=1,NCONF,MBUF
      NN = min(MBUF,(NCONF+1-ISTA))
      call PKVEC(NN,CI(ISTA),ICI)
      call iDAFILE(LUEIG,1,ICI,NN,IDFREE)
    end do
  end do
end if
NVTOT = NNEW
NSTOT = 0

return

end subroutine CSTART
