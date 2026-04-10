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

subroutine SONATORB_CPLOT(DENS,FILEBASE,CHARTYPE,ASS,BSS)

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri, sOpSiz
use rassi_aux, only: ipglob
use Symmetry_Info, only: nIrrep
use rassi_data, only: NBASF, NBMX, NBSQ, NBST, NBTRI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: DENS(6,NBTRI)
character(len=*) :: FILEBASE
character(len=8) :: CHARTYPE
integer(kind=iwp) :: ASS, BSS
integer(kind=iwp) :: I, I1, I1I, I2, ICMP, ID1, ID2, IDIR, IDUM(1), iDummy(7,8), IEND, II, II2, IJ, INV, INV2, IOCC, IOPT, IRC, &
                     ISCR, ISCRI, ISTART, ISYLAB, ISYM, ITYPE, J, JI, JOPT, LE, LE1, LS, LS1, LuXXVEC, LV, LV1, NB, NBMX2
real(kind=wp) :: Dummy(1), SUMI, SUMR
character(len=25) :: FNAME
character(len=16) :: FNUM, KNUM, XNUM
character(len=8) :: LABEL
character :: CDIR
real(kind=wp), allocatable :: DMAT(:), DMATI(:), EIG(:), OCC(:), SANG(:), SANGF(:), SANGTI(:), SANGTI2(:), SANGTR(:), SANGTR2(:), &
                              SCR(:), SCRI(:), SZZ(:), VEC(:), VEC2(:), VEC2I(:), VNAT(:), VNATI(:)
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
! NOTE: SCR COULD PROBABLY BE SOMETHING LIKE nTri_Elem(NBMX)
!       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
!       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
call mma_allocate(SZZ,NBTRI,Label='SZZ')
SZZ(:) = Zero
call mma_allocate(VEC,NBSQ,Label='VEC')
VEC(:) = Zero
call mma_allocate(VEC2,NBMX2,Label='VEC2')
VEC2(:) = Zero
call mma_allocate(VEC2I,NBMX2,Label='VEC2I')
VEC2I(:) = Zero
call mma_allocate(SCR,NBMX2,Label='SCR')
SCR(:) = Zero
call mma_allocate(SCRI,NBMX2,Label='SCRI')
SCRI(:) = Zero
call mma_allocate(EIG,NBST,Label='EIG')
EIG(:) = Zero

call mma_allocate(VNAT,NBSQ,Label='VNAT')
VNAT(:) = Zero
call mma_allocate(VNATI,NBSQ,Label='VNATI')
VNATI(:) = Zero
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
  call JACOB(SZZ(LS),VEC(LV),NB,NB)
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
  LS = LS+nTri_Elem(NB)
  LV = LV+NB**2
  LE = LE+NB
end do

call mma_deallocate(SZZ)

call mma_allocate(DMAT,NBMX2,Label='DMAT')
call mma_allocate(DMATI,NBMX2,Label='DMATI')

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

  !ccccccccccccccccccccccc
  !ccccccccccccccccccccccc
  !ccccccccccccccccccccccc
  !ccccccccccccccccccccccc
  ! read in ao matrix for angmom or mltpl
  call mma_allocate(SANG,NBTRI,Label='SANG')
  SANG(:) = Zero

  IRC = -1
  IOPT = ibset(ibset(0,sNoOri),sNoNuc)
  JOPT = ibset(0,sOpSiz)

  if ((ITYPE == 1) .or. (ITYPE == 3)) then
    ICMP = 1
    LABEL = 'MLTPL  0'
    call iRDONE(IRC,JOPT,LABEL,ICMP,IDUM,ISYLAB)
    call RDONE(IRC,IOPT,LABEL,ICMP,SANG,ISYLAB)

    if (IRC /= 0) then
      write(u6,*)
      write(u6,*) '      *** ERROR IN SUBROUTINE  SONATORB ***'
      write(u6,*) '      MLTPL0 INTEGRALS ARE NOT AVAILABLE'
      write(u6,*) '      IRC:',IRC
      write(u6,*)
      call ABEND()
    end if

  else if ((ITYPE == 2) .or. (ITYPE == 4)) then
    ICMP = 3
    LABEL = 'ANGMOM'
    call iRDONE(IRC,JOPT,LABEL,ICMP,IDUM,ISYLAB)
    call RDONE(IRC,IOPT,LABEL,ICMP,SANG,ISYLAB)

    if (IRC /= 0) then
      write(u6,*)
      write(u6,*) '      *** ERROR IN SUBROUTINE  SONATORB ***'
      write(u6,*) '      ANGMOM INTEGRALS ARE NOT AVAILABLE'
      write(u6,*) '      IRC:',IRC
      write(u6,*)
      call ABEND()
    end if

  end if

  !ccccccccccccccccccccccc
  !ccccccccccccccccccccccc
  !ccccccccccccccccccccccc
  !ccccccccccccccccccccccc
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
    DMATI(:) = Zero
    SCR(:) = Zero
    SCRI(:) = Zero

    do J=1,NB
      do I=1,J
        II2 = II2+1
        IJ = NB*(J-1)+I
        JI = NB*(I-1)+J
        if (I /= J) then
          DMAT(IJ) = DENS(IDIR,II2)/Two
          DMAT(JI) = DENS(IDIR,II2)/Two
          DMATI(IJ) = -DENS(IDIR+3,II2)/Two
          DMATI(JI) = DENS(IDIR+3,II2)/Two
        else
          DMAT(IJ) = DENS(IDIR,II2)
          DMATI(JI) = DENS(IDIR+3,II2)
        end if
      end do
    end do

    call DGEMM_('N','N',NB,NB,NB,One,DMAT,NB,VEC(LV),NB,Zero,SCR,NB)
    call DGEMM_('N','N',NB,NB,NB,One,DMATI,NB,VEC(LV),NB,Zero,SCRI,NB)

    call DGEMM_('T','N',NB,NB,NB,One,VEC(LV),NB,SCR,NB,Zero,DMAT,NB)
    call DGEMM_('T','N',NB,NB,NB,One,VEC(LV),NB,SCRI,NB,Zero,DMATI,NB)

    ID1 = 1
    ID2 = 1
    do I=1,NB
      call DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
      DMAT(ID2:ID2+NB-1) = EIG(LE-1+I)*DMAT(ID2:ID2+NB-1)
      call DSCAL_(NB,EIG(LE-1+I),DMATI(ID1),NB)
      DMATI(ID2:ID2+NB-1) = EIG(LE-1+I)*DMATI(ID2:ID2+NB-1)
      ID1 = ID1+1
      ID2 = ID2+NB
    end do

    ! SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
    SCR(:) = Zero
    SCRI(:) = Zero

    ISCR = 1
    ISCRI = 1
    do I=1,NB
      do J=1,I
        IJ = I+NB*(J-1)
        JI = J+NB*(I-1)
        ! simple averaging
        SCR(ISCR) = (DMAT(JI)+DMAT(IJ))/Two
        SCRI(ISCRI) = (DMATI(JI)-DMATI(IJ))/Two
        ! add a factor of two to convert spin -> sigma
        if (ITYPE >= 3) SCR(ISCR) = SCR(ISCR)*Two
        if (ITYPE >= 3) SCRI(ISCRI) = SCRI(ISCRI)*Two
        ISCR = ISCR+1
        ISCRI = ISCRI+1
      end do
    end do

    ! DIAGONALIZE THE DENSITY MATRIX BLOCK:
    VEC2(:) = Zero
    VEC2I(:) = Zero

    call CPLOT_DIAG(SCR,SCRI,NB,VEC2,VEC2I)

    ! LAPACK ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
    II = 0
    do I=1,NB
      II = II+I
      OCC(IOCC+NB+1-I) = SCR(II)
    end do
    IOCC = IOCC+NB

    ! REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
    call DGEMM_('N','N',NB,NB,NB,One,VEC(LV),NB,VEC2,NB,Zero,SCR,NB)
    call DGEMM_('N','N',NB,NB,NB,One,VEC(LV),NB,VEC2I,NB,Zero,SCRI,NB)

    I1 = 1
    I1I = 1
    I2 = INV+NB**2
    do I=1,NB
      I2 = I2-NB
      VNAT(I2:I2+NB-1) = SCR(I1:I1+NB-1)
      VNATI(I2:I2+NB-1) = SCRI(I1:I1+NB-1)
      I1 = I1+NB
      I1I = I1I+NB
    end do
    INV = INV+NB**2
    LV = LV+NB**2
    LE = LE+NB
  end do

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !CCCCCCC TESTING
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCC
  if (IPGLOB >= 4) then

    call mma_allocate(SANGF,NBMX**2,Label='SANGF')
    SANGF(:) = Zero
    call mma_allocate(SANGTR,NBMX**2,Label='SANGTR')
    call mma_allocate(SANGTI,NBMX**2,Label='SANGTI')
    SANGTR(:) = Zero
    SANGTI(:) = Zero
    call mma_allocate(SANGTR2,NBMX**2,Label='SANGTR2')
    call mma_allocate(SANGTI2,NBMX**2,Label='SANGTI2')
    SANGTR2(:) = Zero
    SANGTI2(:) = Zero

    INV = 0
    INV2 = 0
    II = 0
    SUMR = Zero
    SUMI = Zero

    do ISYM=1,nIrrep
      NB = NBASF(ISYM)
      if (NB /= 0) then

        ! Expand integrals for this symmetry to full storage
        SANGF(:) = Zero

        do J=1,NB
          do I=1,J
            IJ = NB*(J-1)+I-1
            JI = NB*(I-1)+J-1

            SANGF(1+JI) = SANG(1+II)

            if (I /= J) then
              if ((ITYPE == 2) .or. (ITYPE == 4)) then
                SANGF(1+IJ) = -SANG(1+II)
              else
                SANGF(1+IJ) = SANG(1+II)
              end if
            end if

            II = II+1

          end do
        end do

        if ((ITYPE == 1) .or. (ITYPE == 3)) then
          call DGEMM_('T','N',NB,NB,NB,One,SANGF,NB,VNAT(1+INV),NB,Zero,SANGTR,NB)
          call DGEMM_('T','N',NB,NB,NB,One,SANGF,NB,VNATI(1+INV),NB,Zero,SANGTI,NB)

          call DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,SANGTR,NB,Zero,SANGTR2,NB)
          call DGEMM_('T','N',NB,NB,NB,One,VNATI(1+INV),NB,SANGTI,NB,One,SANGTR2,NB)

          call DGEMM_('T','N',NB,NB,NB,-One,VNATI(1+INV),NB,SANGTR,NB,Zero,SANGTI2,NB)
          call DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,SANGTI,NB,One,SANGTI,NB)

        else if ((ITYPE == 2) .or. (ITYPE == 4)) then

          call DGEMM_('T','N',NB,NB,NB,One,SANGF,NB,VNAT(1+INV),NB,Zero,SANGTI,NB)
          call DGEMM_('T','N',NB,NB,NB,-One,SANGF,NB,VNATI(1+INV),NB,Zero,SANGTR,NB)

          call DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,SANGTR,NB,Zero,SANGTR2,NB)
          call DGEMM_('T','N',NB,NB,NB,One,VNATI(1+INV),NB,SANGTI,NB,One,SANGTR2,NB)

          call DGEMM_('T','N',NB,NB,NB,-One,VNATI(1+INV),NB,SANGTR,NB,Zero,SANGTI2,NB)
          call DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,SANGTI,NB,One,SANGTI2,NB)

        end if

        ! Sum over the trace
        do I=1,NB
          IJ = I+(I-1)*NB-1
          SUMR = SUMR+OCC(I+INV2)*SANGTR2(1+IJ)
          SUMI = SUMI+OCC(I+INV2)*SANGTI2(1+IJ)
        end do

      end if

      INV = INV+NB**2
      INV2 = INV2+NB

    end do

    write(u6,*) 'Ben P TEST for JA:'
    write(u6,*) 'REAL: ',SUMR
    write(u6,*) 'IMAG: ',SUMI

    call mma_deallocate(SANGF)
    call mma_deallocate(SANGTR)
    call mma_deallocate(SANGTI)
    call mma_deallocate(SANGTR2)
    call mma_deallocate(SANGTI2)
  end if ! IPGLOB >= 4

  call mma_deallocate(SANG)

  ! WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
  ! REAL PART
  if (ITYPE <= 2) then
    write(KNUM,'(I2.2,A,I2.2,A,A)') ASS,'.',BSS,'.','R'
  else
    write(KNUM,'(I2.2,A,I2.2,A,A,A,A)') ASS,'.',BSS,'.',CDIR,'.','R'
  end if
  write(FNUM,'(I8)') BSS
  FNUM = adjustl(FNUM)
  if (ASS /= BSS) then
    write(XNUM,'(I8,A)') ASS,'_'//trim(FNUM)
    FNUM = adjustl(XNUM)
  end if
  if (ITYPE > 2) FNUM = CDIR//trim(FNUM)

  FNAME = FILEBASE//'.'//trim(FNUM)//'.R'
  if (ITYPE == 1) write(u6,'(A,A)') ' NATURAL ORBITALS FOR ',KNUM
  if (ITYPE == 2) write(u6,'(A,A)') ' ANTISING NATURAL ORBITALS FOR  ',KNUM
  if (ITYPE == 3) write(u6,'(A,A)') ' NATURAL SPIN ORBITALS FOR  ',KNUM
  if (ITYPE == 4) write(u6,'(A,A)') ' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

  write(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

  LuxxVec = 50
  LuxxVec = isfreeunit(LuxxVec)

  call WRVEC(FNAME,LUXXVEC,'CO',nIrrep,NBASF,NBASF,VNAT,OCC,Dummy,iDummy,'* DENSITY FOR PROPERTY TYPE '//CHARTYPE//KNUM)

  ! IMAGINARY PART
  if (ITYPE <= 2) then
    write(KNUM,'(I2.2,A,I2.2,A,A)') ASS,'.',BSS,'.','I'
  else
    write(KNUM,'(I2.2,A,I2.2,A,A,A,A)') ASS,'.',BSS,'.',CDIR,'.','I'
  end if

  FNAME = FILEBASE//'.'//trim(FNUM)//'.I'
  if (ITYPE == 1) write(u6,'(A,A)') ' NATURAL ORBITALS FOR ',KNUM
  if (ITYPE == 2) write(u6,'(A,A)') ' ANTISING NATURAL ORBITALS FOR  ',KNUM
  if (ITYPE == 3) write(u6,'(A,A)') ' NATURAL SPIN ORBITALS FOR  ',KNUM
  if (ITYPE == 4) write(u6,'(A,A)') ' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

  write(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

  LuxxVec = 50
  LuxxVec = isfreeunit(LuxxVec)

  call WRVEC(FNAME,LUXXVEC,'CO',nIrrep,NBASF,NBASF,VNATI,OCC,Dummy,iDummy,'* DENSITY FOR PROPERTY TYPE '//CHARTYPE//KNUM)

  !Test a few values
  !call ADD_INFO('SONATORB_CPLOTR',VNAT,1,4)
  !call ADD_INFO('SONATORB_CPLOTI',VNATI,1,4)
  !call ADD_INFO('SONATORB_CPLOTO',OCC,1,4)

end do

call mma_deallocate(DMAT)
call mma_deallocate(DMATI)
call mma_deallocate(VEC)
call mma_deallocate(VEC2)
call mma_deallocate(VEC2I)
call mma_deallocate(SCR)
call mma_deallocate(SCRI)
call mma_deallocate(VNAT)
call mma_deallocate(VNATI)
call mma_deallocate(OCC)

end subroutine SONATORB_CPLOT
