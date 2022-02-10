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

implicit real*8(A-H,O-Z)
dimension AREF(NREF,NREF), EREF(NREF), CI(NCONF), ICI(MBUF)
#include "SysDef.fh"
#include "mrci.fh"
dimension BUF(nCOP), ISTART(MXROOT)

do I=1,MXVEC
  IDISKC(I) = -1
  IDISKS(I) = -1
end do
! FIRST, USE THE CI ARRAY TO STORE THE DIAGONAL ELEMENTS:
IAD25 = IAD25S
do I=1,NCONF,nCOP
  call dDAFILE(Lu_25,2,BUF,nCOP,IAD25)
  NN = min(nCOP,NCONF+1-I)
  call DCOPY_(NN,BUF,1,CI(I),1)
end do
! THESE ARE DIAGONAL ELEMENTS OF THE ELECTRONIC HAMILTONIAN.
! POTNUC SHOULD BE ADDED. IN ADDITION, WE USE AN ENERGY SHIFT.
! NOTE: DISPLACEMENT 1.0d-4 PROTECTS AGAINST DIVIDE ERRORS.
!PAM: Protect all diag elems. Needed in some weird cases.
! ENERGY SHIFT:
ESHIFT = EREF(1)
do I=1,NCONF
  CI(I) = CI(I)+POTNUC-ESHIFT+1.0D-04
end do
call Add_Info('CI_DIAG2',CI(2),1,8)
! REPLACE REFERENCE ENERGIES:
do I=1,NREF
  IR = IREFX(I)
  CI(IR) = EREF(I)-ESHIFT-1.0D-04
end do
if (ICPF == 1) then
  do IREF=1,NREF
    IR = IREFX(IREF)
    CI(IR) = GFAC*CI(IR)
  end do
  GINV = 1.0d00/GFAC
  call DSCAL_(NCONF,GINV,CI,1)
end if
IDFREE = 0
IDISKD = 0
do ISTA=1,NCONF,MBUF
  NN = min(MBUF,(NCONF+1-ISTA))
  call dDAFILE(LUEIG,1,CI(ISTA),NN,IDFREE)
end do
! THEN, SET UP START CI VECTORS IN MCSF BASIS:
call DCOPY_(NCONF,[0.0d00],0,CI,1)
if (IREST == 0) then
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
    write(6,*) ' THE FOLLOWING REFERENCE ROOTS ARE USED AS START VECTORS:'
    call XFLUSH(6)
    write(6,'(12(A,I2))') ' ROOTS NR ',ISTART(1),(',',ISTART(I),I=2,NNEW-1),', AND ',ISTART(NNEW)
    call XFLUSH(6)
    if (NNEW > NRROOT) then
      write(6,*) ' (THE FIRST EXTRA ROOT(S) WERE INCLUDED IN ORDER TO IMPROVE CONVERGENCE)'
      call XFLUSH(6)
    end if
  else
    write(6,'(A,I2,A)') ' ROOT NR ',ISTART(1),' IS USED AS START VECTOR.'
    call XFLUSH(6)
  end if
  do I=1,NNEW
    ISTA = ISTART(I)
    IR = IREFX(ISTA)
    CI(IR) = 1.0d00
    IDISKC(I) = IDFREE
    do ISTA=1,NCONF,MBUF
      NN = min(MBUF,(NCONF+1-ISTA))
      call PKVEC(NN,CI(ISTA),ICI)
      call iDAFILE(LUEIG,1,ICI,NN,IDFREE)
    end do
    CI(IR) = 0.0d00
  end do
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
