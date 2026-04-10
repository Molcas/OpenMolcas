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

subroutine BINAT()

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: MUL, nIrrep
use OneDat, only: sNoNuc, sNoOri
use rassi_aux, only: iDisk_TDM
use rassi_data, only: NBASF, NBMX, NBSQ, NBST, NBTRI, NTDMZZ
use rassi_global_arrays, only: EIGVEC, JBNUM
use Cntrl, only: IBINA, IRREP, LuTDM, NBINA, NSTATE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, ICMP, IDISK, IDUMMY(2), IE, iEmpty, iGo, IJPAIR, IOFF_ISV(8), IOFF_SEV(8), IOFF_TDM(8), IOFF_VEC(8), IOPT, &
                     IRC, ISEL, ISV, ISYLAB, ISYM, ISYM1, ISYM2, ITD, ITD1, ITD2, IV, J, K, KEIG_BRA, KEIG_KET, L, LB, LE, LE1, &
                     LE2, LK, LS, LS1, LS2, LSYM12, LSYM_BRA, LSYM_KET, LUNIT, LV, LV1, LV2, NB, NB1, NB2, NBMIN
real(kind=wp) :: DUMMY(2), S_EV, SSEL, SUMSNG, SWAP, X
character(len=24) :: FNAME
character(len=21) :: TXT
character(len=16) :: KNUM
character(len=8) :: BNUM, LABEL
real(kind=wp), allocatable :: BRABNO(:), KETBNO(:), ONBAS(:), SAO(:), SCR(:), SEV(:), SNGV1(:), SNGV2(:), SVAL(:), TDMAO(:), &
                              TDMAT(:), UMAT(:), VTMAT(:)
integer(kind=iwp), external :: ISFREEUNIT
real(kind=wp), external :: DDOT_

! Tables of starting locations, created and used later
! IOFF_VEC, IOFF_SEV, IOFF_TDM, IOFF_ISV

! Nr of basis functions, total
NBSQ = sum(NBASF(1:nIrrep)**2)
!============================================================
! START BY CREATING A SET OF ORTHONORMAL VECTORS:
call mma_allocate(ONBAS,NBSQ,Label='ONBAS')
! EIGENVALUES OF OVERLAP MATRIX:
call mma_allocate(SEV,NBST,Label='SEV')
! READ ORBITAL OVERLAP MATRIX.
call mma_allocate(SAO,NBTRI,Label='SAO')
IRC = -1
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICMP = 1
ISYLAB = 1
LABEL = 'MLTPL  0'
call RDONE(IRC,IOPT,LABEL,ICMP,SAO,ISYLAB)
if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '      *** ERROR IN SUBROUTINE  BINAT ***'
  write(u6,*) '      OVERLAP INTEGRALS ARE NOT AVAILABLE'
  write(u6,*)
  call ABEND()
end if

! LOOP OVER SYMMETRY BLOCKS
! DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
LS = 1
LV = 1
LE = 1
do ISYM=1,nIrrep
  NB = NBASF(ISYM)
  call unitmat(ONBAS(LV:LV+NB**2-1),NB)
  call JACOB(SAO(LS),ONBAS(LV),NB,NB)
  ! SORT IN ORDER OF DECREASING EIGENVALUES.
  LS1 = LS
  do I=1,NB-1
    ISEL = I
    SSEL = SAO(LS1)
    LS2 = LS1
    do J=I+1,NB
      LS2 = LS2+J
      if (SAO(LS2) > SSEL) then
        ISEL = J
        SSEL = SAO(LS2)
      end if
    end do
    if (ISEL > I) then
      LS2 = LS-1+nTri_Elem(ISEL)
      SWAP = SAO(LS2)
      SAO(LS2) = SAO(LS1)
      SAO(LS1) = SWAP
      do K=1,NB
        SWAP = ONBAS(LV-1+K+NB*(ISEL-1))
        ONBAS(LV-1+K+NB*(ISEL-1)) = ONBAS(LV-1+K+NB*(I-1))
        ONBAS(LV-1+K+NB*(I-1)) = SWAP
      end do
    end if
    LS1 = LS1+I+1
  end do
  ! SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
  LS1 = LS
  LV1 = LV
  LE1 = LE
  do I=1,NB
    S_EV = SAO(LS1)
    SEV(LE1) = S_EV
    if (S_EV > 1.0e-14_wp) then
      ONBAS(LV1:LV1+NB-1) = ONBAS(LV1:LV1+NB-1)/sqrt(S_EV)
    else
      ONBAS(LV1:LV1+NB-1) = Zero
    end if
    LS1 = LS1+I+1
    LV1 = LV1+NB
    LE1 = LE1+1
  end do
  LS = LS+nTri_Elem(NB)
  LV = LV+NB**2
  LE = LE+NB
end do
call mma_deallocate(SAO)
! Starting at ONBAS there is now symmetry blocks of CMO arrays
! describing orthonormal vectors. In case the AO overlap matrix is
! (almost) singular, one or more vectors at the end of each symmetry
! block will be null vectors.
!============================================================

! Left and right singular vectors, used temporarily
! in calls to SVD routine. Also temporary, singular values.
call mma_allocate(UMAT,NBMX**2,Label='UMAT')
call mma_allocate(VTMAT,NBMX**2,Label='VTMAT')
call mma_allocate(SVAL,NBMX,Label='SVAL')
! Final bra and ket singular value array:
call mma_allocate(SNGV1,NBST,Label='SNGV1')
! An extra copy for singular values in a different order
! used until proper GV support for binatural orbitals.
call mma_allocate(SNGV2,NBST,Label='SNGV2')
! The transition density matrix, symmetry-blocked
! Symmetry blocks may combine different symmetries, but the size
! is certainly less than or equal to NBSQ.
call mma_allocate(TDMAT,NBSQ,Label='TDMAT')
! Same, read buffer
call mma_allocate(TDMAO,NBSQ,Label='TDMAO')
! Temporary intermediate in matrix multiplies
! (Also used as temporary when transposing some TDMAO matrices)
call mma_allocate(SCR,NBSQ,Label='SCR')
! The BRA and KET binatural orbitals:
call mma_allocate(BRABNO,NBSQ,Label='BRABNO')
call mma_allocate(KETBNO,NBSQ,Label='KETBNO')

! A long loop over eigenstate pairs:
do IJPAIR=1,NBINA
  ! Requested state pairs for computation: (OBSOLETE)
  KEIG_BRA = IBINA(1,IJPAIR)
  KEIG_KET = IBINA(2,IJPAIR)
  ! Get symmetries, via jobiph number for the states:
  LSYM_BRA = IRREP(JBNUM(KEIG_BRA))
  LSYM_KET = IRREP(JBNUM(KEIG_KET))
  ! Combined symmetry:
  LSYM12 = MUL(LSYM_BRA,LSYM_KET)
  ! For relating left and right symmetry blocks, offset tables are
  ! needed for the singular values and for the TDM.
  ITD = 0
  ISV = 0
  do ISYM1=1,nIrrep
    IOFF_TDM(ISYM1) = ITD
    IOFF_ISV(ISYM1) = ISV
    ISYM2 = MUL(ISYM1,LSYM12)
    ITD = ITD+NBASF(ISYM1)*NBASF(ISYM2)
    ISV = ISV+NBASF(ISYM1)
  end do
  TDMAT(:) = Zero
  ! DOUBLE LOOP OVER RASSCF WAVE FUNCTIONS
  do I=1,NSTATE
    if (IRREP(JBNUM(I)) /= LSYM_BRA) cycle
    do J=1,NSTATE
      if (IRREP(JBNUM(J)) /= LSYM_KET) cycle
      ! PICK UP TRANSITION DENSITY MATRIX FOR THIS PAIR OF RASSCF STATES:
      ! WEIGHT WITH WHICH THEY CONTRIBUTE IS EIGVEC(I,KEIG_BRA)*EIGVEC(J,KEIG_KET).
      X = EIGVEC(KEIG_BRA,i)*EIGVEC(KEIG_KET,j)
      IDISK = iDisk_TDM(J,I,1)
      IEMPTY = iDisk_TDM(J,I,2)
      iOpt = 2
      iGo = 1
      if (btest(iEMPTY,0)) then
        if (I > J) then
          call dens2file(TDMAO,TDMAO,TDMAO,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
        else
          ! Pick up conjugate TDM array, and transpose it into TDMAO.
          call dens2file(SCR,SCR,SCR,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
          ! Loop over the receiving side:
          do ISYM1=1,nIrrep
            ISYM2 = MUL(ISYM1,LSYM12)
            NB1 = NBASF(ISYM1)
            NB2 = NBASF(ISYM2)
            do K=1,NB1
              do L=1,NB2
                TDMAO(IOFF_TDM(ISYM1)+K+NB1*(L-1)) = SCR(IOFF_TDM(ISYM2)+L+NB2*(K-1))
              end do
            end do
          end do
        end if
        TDMAT(:) = TDMAT(:)+X*TDMAO(:)
      end if
    end do
  end do
  ! TDMAT() now contains the transition density matrix in AO basis for
  ! the eigenstates.
  ! ------------------------------------------------------------------

  ! LOOP OVER SYMMETRY BLOCKS OF TDMAT.
  ! On the ket side, ISYM2 is not looping sequentially so we need
  ! tables of offsets:
  IV = 0
  IE = 0
  do ISYM=1,nIrrep
    IOFF_VEC(ISYM) = IV
    IOFF_SEV(ISYM) = IE
    IV = IV+NBASF(ISYM)**2
    IE = IE+NBASF(ISYM)
  end do
  SNGV1(:) = Zero
  SNGV2(:) = Zero
  ITD = 0
  do ISYM1=1,nIrrep
    ISYM2 = MUL(ISYM1,LSYM12)
    NB1 = NBASF(ISYM1)
    NB2 = NBASF(ISYM2)
    LV1 = 1+IOFF_VEC(ISYM1)
    LV2 = 1+IOFF_VEC(ISYM2)
    LE1 = 1+IOFF_SEV(ISYM1)
    LE2 = 1+IOFF_SEV(ISYM2)
    ! TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
    ! BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
    ! SCALING WITH THE EIGENVECTORS OF THE OVERLAP MATRIX:
    call DGEMM_('N','N',NB1,NB2,NB2,One,TDMAT(1+ITD),NB1,ONBAS(LV2),NB2,Zero,SCR,NB1)
    call DGEMM_('T','N',NB1,NB2,NB1,One,ONBAS(LV1),NB1,SCR,NB1,Zero,TDMAT(1+ITD),NB1)
    ITD1 = ITD
    do I=1,NB1
      call DSCAL_(NB2,SEV(LE1-1+I),TDMAT(1+ITD1),NB1)
      ITD1 = ITD1+1
    end do
    ITD2 = ITD
    do I=1,NB2
      TDMAT(1+ITD2:NB1+ITD2) = SEV(LE2-1+I)*TDMAT(1+ITD2:NB1+ITD2)
      ITD2 = ITD2+NB1
    end do

    ! SVD DECOMPOSITION OF THIS MATRIX BLOCK:
    call FULL_SVD(NB1,NB2,TDMAT(1+ITD),UMAT,VTMAT,SVAL)
    ! On return, UMAT has dimension (NB1,NB1),
    ! On return, VTMAT has dimension (NB2,NB2), transpose storage
    NBMIN = min(NB1,NB2)
    ! REEXPRESS THE SINGULAR VECTORS USING AO BASIS:
    LB = 1+IOFF_VEC(ISYM1)
    LK = 1+IOFF_VEC(ISYM2)
    call DGEMM_('N','N',NB1,NB1,NB1,One,ONBAS(LV1),NB1,UMAT,NB1,Zero,BRABNO(LB),NB1)
    call DGEMM_('N','T',NB2,NB2,NB2,One,ONBAS(LV2),NB2,VTMAT,NB2,Zero,KETBNO(LK),NB2)

    ! Move the singular values into their proper places:
    SNGV1(IOFF_ISV(ISYM1)+1:IOFF_ISV(ISYM1)+NBMIN) = SVAL(1:NBMIN)
    SNGV2(IOFF_ISV(ISYM2)+1:IOFF_ISV(ISYM1)+NBMIN) = SVAL(1:NBMIN)
    ITD = ITD+NB1*NB2
  end do

  ! WRITE OUT THIS SET OF BI-NATURAL ORBITALS. THE FILES WILL BE NAMED
  ! BIORB.x_y, where x,y are KEIG_BRA and KEIG_KET.
  ! The BRA and KET orbitals will be written as alpha and beta, respectively,
  ! and the singular values will be written as "occupation numbers".

  write(u6,*) ' Binatural singular values for the transition from'
  write(u6,*) ' ket eigenstate KEIG_KET to bra eigenstate KEIG_BRA'
  write(u6,'(1x,I2,a,i2)') KEIG_BRA,' <-- ',KEIG_KET
  do I=1,nIrrep
    NB = NBASF(I)
    if (NB /= 0) then
      write(u6,'(A,I2)') ' SYMMETRY SPECIES:',I
      LS = 1+IOFF_SEV(I)
      write(u6,'(1X,10F8.5)') (SNGV1(LS-1+J),J=1,NB)
    end if
  end do

  write(BNUM,'(I8)') KEIG_BRA
  BNUM = adjustl(BNUM)
  write(KNUM,'(I8)') KEIG_KET
  KNUM = adjustl(KNUM)
  TXT = trim(BNUM)//' <-- '//trim(KNUM)
  KNUM = trim(BNUM)//'_'//trim(KNUM)

  FNAME = 'BIORB.'//KNUM
  write(u6,'(A,A)') ' Orbitals are written onto file id = ',FNAME
  LUNIT = ISFREEUNIT(50)
  call WRVEC_(FNAME,LUNIT,'CO',1,nIrrep,NBASF,NBASF,BRABNO,KETBNO,SNGV1,SNGV2,DUMMY,DUMMY,IDUMMY, &
              '* Binatural orbitals from transition '//trim(TXT),0)
  close(LUNIT)
  SUMSNG = DDOT_(sum(NBASF),SNGV1,1,SNGV2,1)
  call ADD_INFO('BINAT',[SUMSNG],1,5)

  ! End of very long loop over eigenstate pairs.
end do

write(u6,*) repeat('*',80)
call mma_deallocate(ONBAS)
call mma_deallocate(SEV)
call mma_deallocate(SCR)
call mma_deallocate(UMAT)
call mma_deallocate(VTMAT)
call mma_deallocate(SVAL)
call mma_deallocate(SNGV1)
call mma_deallocate(SNGV2)
call mma_deallocate(TDMAT)
call mma_deallocate(TDMAO)
call mma_deallocate(BRABNO)
call mma_deallocate(KETBNO)

end subroutine BINAT
