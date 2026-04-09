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

subroutine SONATORB_PLOT(DENS,FILEBASE,CHARTYPE,ASS,BSS)

use OneDat, only: sNoNuc, sNoOri
use Symmetry_Info, only: nIrrep
use rassi_data, only: NBASF, NBMX, NBSQ, NBST, NBTRI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: DENS(6,NBTRI)
character(len=*), intent(in) :: FILEBASE
character(len=8), intent(in) :: CHARTYPE
integer(kind=iwp), intent(in) :: ASS, BSS
integer(kind=iwp) :: I, I1, I2, ICMP, ID1, ID2, IDIR, iDummy(7,8), IEND, II, II2, IJ, INV, IOCC, IOPT, IRC, ISCR, ISTART, ISYLAB, &
                     ISYM, ITYPE, J, JI, LE, LE1, LS, LS1, LuXXVEC, LV, LV1, NB, NBMX2
real(kind=wp) :: Dummy(1)
character(len=25) :: FNAME
character(len=16) :: FNUM, KNUM, XNUM
character(len=8) :: LABEL
character :: CDIR
real(kind=wp), allocatable :: DMAT(:), EIG(:), OCC(:), SCR(:), SZZ(:), VEC(:), VEC2(:), VNAT(:)
integer(kind=iwp), external :: IsFreeUnit

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PLOTTING SECTION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Get the proper type of the property
ITYPE = 0
if (CHARTYPE == 'HERMSING') ITYPE = 1
if (CHARTYPE == 'ANTISING') ITYPE = 2
if (CHARTYPE == 'HERMTRIP') ITYPE = 3
if (CHARTYPE == 'ANTITRIP') ITYPE = 4
if (ITYPE == 0) then
  write(u6,*) 'RASSI/SONATORB internal error.'
  write(u6,*) 'Erroneous property type:',CHARTYPE
  call ABEND()
end if

NBMX2 = NBMX**2

! SZZ  - AO Overlap integral
! VEC  - AO Overlap eigenvectors
! EIG  - AO Overlap eigenvalues
! VEC2 - Eigenvectors of density matrix
! SCR  - Temporary for matrix multiplication
! NOTE: SCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
!       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
!       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
call mma_allocate(SZZ,NBTRI,Label='SZZ')
call mma_allocate(VEC,NBSQ,Label='VEC')
call mma_allocate(VEC2,NBMX2,Label='VEC2')
call mma_allocate(SCR,NBMX2,Label='SCR')
call mma_allocate(EIG,NBST,Label='EIG')
SZZ(:) = Zero
VEC(:) = Zero
VEC2(:) = Zero
SCR(:) = Zero
EIG(:) = Zero

call mma_allocate(VNAT,NBSQ,Label='VNAT')
VNAT(:) = Zero
call mma_allocate(OCC,NBST,Label='OCC')
OCC(:) = Zero

! READ ORBITAL OVERLAP MATRIX.
IRC = -1

! IOPT=6, origin and nuclear contrib not read
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICMP = 1
ISYLAB = 1
LABEL = 'MLTPL  0'
call RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '      *** ERROR IN SUBROUTINE  SONATORB ***'
  write(u6,*) '      OVERLAP INTEGRALS ARE NOT AVAILABLE'
  write(u6,*)
  call ABEND()
end if

! DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
LS = 1
LV = 1
LE = 1
do ISYM=1,nIrrep
  NB = NBASF(ISYM)
  call unitmat(VEC(LV:LV+NB**2-1),NB)
  call JACOB(SZZ(LS),VEC,NB,NB)
  ! SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
  LS1 = LS
  LV1 = LV
  LE1 = LE
  do I=1,NB
    EIG(LE1) = SZZ(LS1)
    VEC(LV1:LV1+NB-1) = VEC(LV1:LV1+NB-1)/sqrt(max(SZZ(LS1),1.0e-14_wp))
    LS1 = LS1+I+1
    LV1 = LV1+NB
    LE1 = LE1+1
  end do
  LS = LS+(NB*(NB+1))/2
  LV = LV+NB**2
  LE = LE+NB
end do

call mma_deallocate(SZZ)

call mma_allocate(DMAT,NBMX2,Label='DMAT')

if (ITYPE <= 2) then
  ISTART = 3
  IEND = 3
else
  ISTART = 1
  IEND = 3
end if

do IDIR=ISTART,IEND

  CDIR = '?'
  if (IDIR == 1) CDIR = 'X'
  if (IDIR == 2) CDIR = 'Y'
  if (IDIR == 3) CDIR = 'Z'

  INV = 1
  II2 = 0
  IOCC = 0
  LV = 1
  LE = 1
  do ISYM=1,nIrrep
    NB = NBASF(ISYM)
    if (NB == 0) cycle

    ! TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
    ! BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
    ! SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

    ! expand the triangular matrix for this symmetry to a square matrix
    DMAT(:) = Zero
    SCR(:) = Zero
    do J=1,NB
      do I=1,J
        II2 = II2+1
        IJ = NB*(J-1)+I
        JI = NB*(I-1)+J
        if (I /= J) then
          DMAT(IJ) = DENS(IDIR,II2)/Two
          DMAT(JI) = DENS(IDIR,II2)/Two
        else
          DMAT(IJ) = DENS(IDIR,II2)
        end if
      end do
    end do

    call DGEMM_('N','N',NB,NB,NB,One,DMAT,NB,VEC(LV),NB,Zero,SCR,NB)
    call DGEMM_('T','N',NB,NB,NB,One,VEC(LV),NB,SCR,NB,Zero,DMAT,NB)

    ID1 = 1
    ID2 = 1
    do I=1,NB
      call DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
      DMAT(ID2:ID2+NB-1) = EIG(LE-1+I)*DMAT(ID2:ID2+NB-1)
      ID1 = ID1+1
      ID2 = ID2+NB
    end do

    ! SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
    SCR(:) = Zero
    ISCR = 1
    do I=1,NB
      do J=1,I
        IJ = I+NB*(J-1)
        JI = J+NB*(I-1)
        ! simple averaging
        SCR(ISCR) = (DMAT(IJ)+DMAT(JI))/Two

        ! add a factor of two to convert spin -> sigma
        if (ITYPE >= 3) SCR(ISCR) = SCR(ISCR)*Two
        ISCR = ISCR+1
      end do
    end do

    ! DIAGONALIZE THE DENSITY MATRIX BLOCK:
    call unitmat(VEC2,NB)

    call JACOB(SCR,VEC2,NB,NB)
    call JACORD(SCR,VEC2,NB,NB)

    ! JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
    II = 0
    do I=1,NB
      II = II+I
      OCC(IOCC+NB+1-I) = SCR(II)
    end do
    IOCC = IOCC+NB

    ! REEXPRESS THE EIGENVALUES IN AO BASIS FUNCTIONS. REVERSE ORDER.
    call DGEMM_('N','N',NB,NB,NB,One,VEC(LV),NB,VEC2,NB,Zero,SCR,NB)
    I1 = 1
    I2 = INV+NB**2
    do I=1,NB
      I2 = I2-NB
      VNAT(I2:I2+NB-1) = SCR(I1:I1+NB-1)
      I1 = I1+NB
    end do
    INV = INV+NB**2
    LV = LV+NB**2
    LE = LE+NB

  end do

  ! WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
  if (ITYPE <= 2) then
    write(KNUM,'(I2.2,A,I2.2)') ASS,'.',BSS
  else
    write(KNUM,'(I2.2,A,I2.2,A,A)') ASS,'.',BSS,'.',CDIR
  end if
  write(FNUM,'(I8)') BSS
  FNUM = adjustl(FNUM)
  if (ASS /= BSS) then
    write(XNUM,'(I8,A)') ASS,'_'//trim(FNUM)
    FNUM = adjustl(XNUM)
  end if
  if (ITYPE > 2) FNUM = CDIR//trim(FNUM)

  FNAME = FILEBASE//'.'//trim(FNUM)
  if (ITYPE == 1) write(u6,'(A,A)') ' NATURAL ORBITALS FOR ',KNUM
  if (ITYPE == 2) write(u6,'(A,A)') ' ANTISING NATURAL ORBITALS FOR  ',KNUM
  if (ITYPE == 3) write(u6,'(A,A)') ' NATURAL SPIN ORBITALS FOR  ',KNUM
  if (ITYPE == 4) write(u6,'(A,A)') ' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

  write(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

  LuxxVec = 50
  LuxxVec = isfreeunit(LuxxVec)

  call WRVEC(FNAME,LUXXVEC,'CO',nIrrep,NBASF,NBASF,VNAT,OCC,Dummy,iDummy,'* DENSITY FOR PROPERTY TYPE '//CHARTYPE//KNUM)

  ! Test a few values
  !call ADD_INFO('SONATORB_PLOT',VNAT,1,4)

  ! ONLY FOR NATURAL ORBITALS
  if (ITYPE == 1) call ADD_INFO('SONATORB_NO_OCC',OCC,sum(NBASF),4)

end do

call mma_deallocate(DMAT)
call mma_deallocate(VEC)
call mma_deallocate(VEC2)
call mma_deallocate(SCR)
call mma_deallocate(EIG)
call mma_deallocate(VNAT)
call mma_deallocate(OCC)

end subroutine SONATORB_PLOT
