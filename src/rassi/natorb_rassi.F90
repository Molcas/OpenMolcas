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

subroutine NATORB_RASSI(DMAT,TDMZZ,VNAT,OCC,EIGVEC)

use rassi_aux, only: iDisk_TDM
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Cntrl, only: nState, NrNATO
use cntrl, only: LuTDM
use Symmetry_Info, only: nSym => nIrrep
use rassi_data, only: NBSQ, NTDMZZ, NBST, NBTRI, NBMX, NBASF
use Constants, only: Zero, One, Half
use Definitions, only: wp, u6

implicit none
real*8 DMAT(NBSQ), TDMZZ(NTDMZZ), VNAT(NBSQ), OCC(NBST)
real*8 EIGVEC(NSTATE,NSTATE)
character(len=14) FNAME
character(len=8) KNUM, LABEL
integer, external :: ISFREEUNIT
real*8, external :: DDOT_
real*8 Dummy(1)
integer iDummy(7,8)
real*8, allocatable :: SZZ(:), VEC(:), VEC2(:), SCR(:), EIG(:)
integer NSZZ, NVEC, NVEC2, NSCR, NEIG, IRC, IOPT, ICMP, ISYLAB, LS, LV, LE, ISYM, NB, LS1, LV1, LE1, I, KEIG, J, IEMPTY, IDISK, &
        IGO, ID, INV, IOCC, ID1, ID2, ISCR, IJ, JI, I1, I2, ISTOCC, LUXXVEC, II
real*8 X, SumOcc

! ALLOCATE WORKSPACE AREAS.
NSZZ = NBTRI
NVEC = NBSQ
NVEC2 = NBMX**2
NSCR = NBMX**2
NEIG = NBST
call mma_allocate(SZZ,NSZZ)
call mma_allocate(VEC,NVEC)
call mma_allocate(VEC2,NVEC2)
call mma_allocate(SCR,NSCR)
call mma_allocate(EIG,NEIG)
! READ ORBITAL OVERLAP MATRIX.
IRC = -1
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICMP = 1
ISYLAB = 1
LABEL = 'MLTPL  0'
call RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '      *** ERROR IN SUBROUTINE  NATORB ***'
  write(u6,*) '      OVERLAP INTEGRALS ARE NOT AVAILABLE'
  write(u6,*)
  call ABEND()
end if
! DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
LS = 1
LV = 1
LE = 1
do ISYM=1,NSYM
  NB = NBASF(ISYM)
  call DCOPY_(NB**2,[Zero],0,VEC(LV),1)
  call DCOPY_(NB,[One],0,VEC(LV),NB+1)
  call JACOB(SZZ(LS),VEC(LV),NB,NB)
  ! SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
  LS1 = LS
  LV1 = LV
  LE1 = LE
  do I=1,NB
    EIG(LE1) = SZZ(LS1)
    X = One/sqrt(max(SZZ(LS1),1.0e-14_wp))
    call DSCAL_(NB,X,VEC(LV1),1)
    LS1 = LS1+I+1
    LV1 = LV1+NB
    LE1 = LE1+1
  end do
  LS = LS+(NB*(NB+1))/2
  LV = LV+NB**2
  LE = LE+NB
end do
call mma_deallocate(SZZ)

! VERY LONG LOOP OVER EIGENSTATES KEIG.
do KEIG=1,NRNATO

  call DCOPY_(NBSQ,[Zero],0,DMAT,1)
  ! DOUBLE LOOP OVER RASSCF WAVE FUNCTIONS, TRIANGULAR.
  do I=1,NSTATE
    do J=1,I
      ! PICK UP TRANSITION DENSITY MATRIX FOR THIS PAIR OF RASSCF STATES:
      ! WEIGHT WITH WHICH THIS TDM CONTRIBUTES IS EIGVEC(I,KEIG)*EIGVEC(J,KEIG).
      ! HOWEVER, WE ARE LOOPING TRIANGULARLY AND WILL RESTORE SYMMETRY BY
      ! ADDING TRANSPOSE AFTER DMAT HAS BEEN FINISHED, SO I=J IS SPECIAL CASE:
      X = EIGVEC(I,KEIG)*EIGVEC(J,KEIG)
      if (abs(X) > 1.0e-12_wp) then
        iEmpty = iDisk_TDM(I,J,2)
        if (iand(iEmpty,1) /= 0) then
          IDISK = iDisk_TDM(I,J,1)
          iOpt = 2
          iGo = 1
          call dens2file(TDMZZ,TDMZZ,TDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
          if (I == J) X = Half*X
          call DAXPY_(NTDMZZ,X,TDMZZ,1,DMAT,1)
        end if
      end if
    end do
  end do
  ! LOOP OVER SYMMETRY BLOCKS OF DMAT.
  ID = 1
  INV = 1
  IOCC = 0
  LV = 1
  LE = 1
  do ISYM=1,NSYM
    NB = NBASF(ISYM)
    ! TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
    ! BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
    ! SCALING WITH THE EIGENVECTORS OF THE OVERLAP MATRIX:
    call DGEMM_('N','N',NB,NB,NB,One,DMAT(ID),NB,VEC(LV),NB,Zero,SCR,NB)
    call DGEMM_('T','N',NB,NB,NB,One,VEC(LV),NB,SCR,NB,Zero,DMAT(ID),NB)
    ID1 = ID
    ID2 = ID
    do I=1,NB
      call DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
      call DSCAL_(NB,EIG(LE-1+I),DMAT(ID2),1)
      ID1 = ID1+1
      ID2 = ID2+NB
    end do
    ! SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
    ISCR = 1
    do I=1,NB
      do J=1,I
        IJ = I+NB*(J-1)
        JI = J+NB*(I-1)
        SCR(ISCR) = DMAT(ID-1+IJ)+DMAT(ID-1+JI)
        ISCR = ISCR+1
      end do
    end do
    ! DIAGONALIZE THE DENSITY MATRIX BLOCK:
    VEC2(:) = Zero
    call DCOPY_(NB,[One],0,VEC2,NB+1)
    call JACOB(SCR,VEC2,NB,NB)
    call JACORD(SCR,VEC2,NB,NB)
    ! JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
    II = 0
    do I=1,NB
      II = II+I
      OCC(IOCC+NB+1-I) = max(Zero,SCR(II))
    end do
    IOCC = IOCC+NB
    ! REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
    call DGEMM_('N','N',NB,NB,NB,One,VEC(LV),NB,VEC2,NB,Zero,SCR,NB)
    I1 = 1
    I2 = INV+NB**2
    do I=1,NB
      I2 = I2-NB
      call DCOPY_(NB,SCR(I1),1,VNAT(I2),1)
      I1 = I1+NB
    end do
    ID = ID+NB**2
    INV = INV+NB**2
    LV = LV+NB**2
    LE = LE+NB
  end do
  ! WRITE OUT THIS SET OF NATURAL ORBITALS. THE FILES WILL BE NAMED
  ! SIORB.1, SIORB.2, ...
  write(KNUM,'(I8)') KEIG
  KNUM = adjustl(KNUM)
  FNAME = 'SIORB.'//KNUM
  write(u6,'(A,I2)') ' NATURAL ORBITALS FOR EIGENSTATE NR ',KEIG
  write(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ID = ',FNAME
  write(u6,'(A)') ' OCCUPATION NUMBERS:'
  ISTOCC = 0
  do I=1,NSYM
    NB = NBASF(I)
    if (NB /= 0) then
      write(u6,'(A,I2)') ' SYMMETRY SPECIES:',I
      write(u6,'(1X,10F8.5)') (OCC(ISTOCC+J),J=1,NB)
    end if
    ISTOCC = ISTOCC+NB
  end do
  LuxxVec = 50
  LuxxVec = isfreeunit(LuxxVec)
  call WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,VNAT,OCC,Dummy,iDummy,'* NATURAL ORBITALS FROM RASSI EIGENSTATE NR '//trim(KNUM))
  SUMOCC = DDOT_(sum(NBASF),OCC,1,OCC,1)
  call ADD_INFO('NATORB',[SUMOCC],1,5)

  ! End of very long loop over eigenstates KEIG.
end do

write(u6,*) repeat('*',80)
call mma_deallocate(VEC)
call mma_deallocate(VEC2)
call mma_deallocate(SCR)
call mma_deallocate(EIG)

return

end subroutine NATORB_RASSI
