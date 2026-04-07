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

subroutine mkSODIAG(UMATR,UMATI,NSS)

use rassi_aux, only: ipglob
use rassi_data, only: NBST
use Cntrl, only: IFCURD, SODIAG, SODIAGNSTATE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, cZero, cOne, cm_s, hPlanck, gElectron, mBohr
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NSS
real(kind=wp) :: UMATR(NSS,NSS), UMATI(NSS,NSS)

integer(kind=iwp) :: I, I2, IC, IDIR, IJ, INFO, IOPT, ISTATE, J, J2, JSTATE, L, LCWORK, LUMAXES, N
real(kind=wp) :: AXI, AXR, AYI, AYR, AZI, AZR, GE, GTENS(3), IDENTMAT(3,3), MAXES(3,3), MAXES2(3,3), MU_BOHR, SXI, SXR, SYI, SYR, &
                 SZI, SZR
character(len=11) :: FILEBASE, FILEBASEL
complex(kind=wp), allocatable :: BPTST(:,:), DEIGVEC(:,:), H_ZEE(:,:), PROP(:,:,:), PROP2(:,:,:), ZOUT(:,:), ZWORK(:)
real(kind=wp), allocatable :: DEIGVAL(:), DMATTMP(:), EIGVECI(:,:), EIGVECR(:,:), LMATI(:,:,:,:), LMATR(:,:,:,:), &
                              MUMAT2I(:,:,:,:), MUMAT2R(:,:,:,:), RWORK(:), SMATI(:,:,:,:), SMATR(:,:,:,:), UWI(:), UWR(:)
integer(kind=iwp), external :: IsFreeUnit

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Matrices

! PROP, PROP2 Storage for expectation values (all directions)
! UMATR,I     Transformation matrix (real, imaginary parts)
! NSS         Number of spin states
! DEIGVAL     Storage for eigenvalues (REAL!)
! DEIGVEC     Storage for eigenvectors
! BPTST       Storage for some testing

ge = -gElectron
MU_BOHR = mBohr*(cm_s*hPlanck) ! in cm-1/T

write(u6,*)
write(u6,*)
write(u6,*) '***********************************************'
write(u6,*) '* STARTING SODIAG *****************************'
write(u6,*) '***********************************************'
write(u6,*)
write(u6,*)

N = SODIAGNSTATE
write(u6,*) 'Number of states to be diagonalized: ',N
write(u6,*) 'STATES:'
do I=1,N
  write(u6,*) SODIAG(I)
end do

call mma_allocate(LMATR,SODIAGNSTATE,SODIAGNSTATE,3,3,Label='LMATR')
call mma_allocate(LMATI,SODIAGNSTATE,SODIAGNSTATE,3,3,Label='LMATI')
call mma_allocate(SMATR,SODIAGNSTATE,SODIAGNSTATE,3,3,Label='SMATR')
call mma_allocate(SMATI,SODIAGNSTATE,SODIAGNSTATE,3,3,Label='SMATI')
call mma_allocate(MUMAT2R,SODIAGNSTATE,SODIAGNSTATE,3,3,Label='MUMAT2R')
call mma_allocate(MUMAT2I,SODIAGNSTATE,SODIAGNSTATE,3,3,Label='MUMAT2I')

call DCOPY_(9*N**2,[Zero],0,LMATR,1)
call DCOPY_(9*N**2,[Zero],0,LMATI,1)
call DCOPY_(9*N**2,[Zero],0,SMATR,1)
call DCOPY_(9*N**2,[Zero],0,SMATI,1)
call DCOPY_(9*N**2,[Zero],0,MUMAT2R,1)
call DCOPY_(9*N**2,[Zero],0,MUMAT2I,1)

call mma_allocate(DMATTMP,3*(NBST*(NBST+1)),Label='DMATTMP')

!> identity mat
call unitmat(IDENTMAT,3)

call mma_allocate(DEIGVAL,SODIAGNSTATE,Label='DEIGVAL')
call mma_allocate(RWORK,3*SODIAGNSTATE-2,Label='RWORK')
call mma_allocate(BPTST,SODIAGNSTATE,SODIAGNSTATE,Label='BPTST')
call mma_allocate(DEIGVEC,SODIAGNSTATE,SODIAGNSTATE,Label='DEIGVEC')
call mma_allocate(H_ZEE,SODIAGNSTATE,SODIAGNSTATE,Label='H_ZEE')
call mma_allocate(PROP,3,SODIAGNSTATE,SODIAGNSTATE,Label='PROP')
call mma_allocate(PROP2,3,SODIAGNSTATE,SODIAGNSTATE,Label='PROP2')
call mma_allocate(ZOUT,SODIAGNSTATE,SODIAGNSTATE,Label='ZOUT')
call mma_allocate(ZWORK,2*SODIAGNSTATE-1,Label='ZWORK')

! First, we calculate the expectation values of
!  (L+ge*S)x (L+ge*S)y (L+ge*S)z
! These get stored in PROP()

! Only work with one triangle - this is a hermitian matrix
do J=1,N
  do I=1,J

    ISTATE = SODIAG(I)
    JSTATE = SODIAG(J)
    write(u6,*) 'States: ',ISTATE,JSTATE

    iOpt = 0
    call SONATORBM('ANTISING',UMATR,UMATI,ISTATE,JSTATE,NSS,iOpt,IDENTMAT,DMATTMP)

    IC = -1
    iOpt = 1
    call SONATORBM_INT(DMATTMP,'ANGMOM',IC,'ANTISING',ISTATE,JSTATE,iOpt,IDENTMAT,AXR,AYR,AZR,AXI,AYI,AZI)

    iOpt = 1
    call SONATORBM('HERMTRIP',UMATR,UMATI,ISTATE,JSTATE,NSS,iOpt,IDENTMAT,DMATTMP)

    IC = 1
    iOpt = 0
    call SONATORBM_INT(DMATTMP,'MLTPL  0',IC,'HERMTRIP',ISTATE,JSTATE,iOpt,IDENTMAT,SXR,SYR,SZR,SXI,SYI,SZI)

    ! The first index of PROP is the direction
    PROP(1,I,J) = -cmplx(AXR+ge*SXR,AXI+ge*SXI,kind=wp)
    PROP(2,I,J) = -cmplx(AYR+ge*SYR,AYI+ge*SYI,kind=wp)
    PROP(3,I,J) = -cmplx(AZR+ge*SZR,AZI+ge*SZI,kind=wp)

    if (I /= J) then
      PROP(1,J,I) = conjg(PROP(1,I,J))
      PROP(2,J,I) = conjg(PROP(2,I,J))
      PROP(3,J,I) = conjg(PROP(3,I,J))
    end if

  end do
end do

!call ADD_INFO('SODIAG_PROP',PROP,3*SODIAGNSTATE**2,4)

! Calculate the atens as in single_aniso
call ATENS_RASSI(PROP,N,GTENS,MAXES,IPGLOB)

do l=1,3
  do i=1,N
    do j=1,N
      PROP2(l,i,j) = sum(PROP(:,i,j)*maxes(:,l))
    end do
  end do
end do

call atens_RASSI(PROP2,N,GTENS,MAXES2,2)

! Diagonalize along each direction
! LOOP OVER THE DIRECTIONS
do IDIR=1,3

  ! apply the magnetic field along the main iDir axis
  DEIGVAL(:) = Zero

  DEIGVEC(:,:) = cZero

  H_ZEE(:,:) = -MU_BOHR*PROP2(IDIR,:,:)

  if (IPGLOB >= 4) then
    write(u6,*) 'BP: H_ZEE:'
    write(u6,*) H_ZEE
    write(u6,*) 'BP: PROP:'
    write(u6,*) 'L+2S in main axis direction (after rotation)',IDIR
    do I=1,N
      do J=1,N

        write(u6,*) 'I,J',SODIAG(I),SODIAG(J),PROP2(IDIR,I,J)

      end do
    end do

    write(u6,*)
    write(u6,*)
    write(u6,*)
    write(u6,*) 'Direction: ',IDIR
    write(u6,*) 'H_ZEE Before Diagonalization'

    do I=1,N
      do J=1,N
        ISTATE = SODIAG(I)
        JSTATE = SODIAG(J)
        write(u6,*) ISTATE,JSTATE,H_ZEE(I,J)
      end do
    end do
  end if

  ! DIAGONALIZE
  lcwork = (2*n-1); info = 0
  call zheev_('V','U',n,h_zee,n,deigval,zwork,lcwork,rwork,info)

  !> put eigenvectors in deigvec
  call zcopy(n**2,h_zee,1,deigvec,1)

  if (INFO /= 0) then
    write(u6,*) 'DIAGONALIZATION FAILED! ERROR: ',INFO
    call ABEND()
  end if

  if (IPGLOB >= 4) then
    write(u6,*) 'EIGENVALUES OF L+ge*S in direction: ',IDIR
    write(u6,*) DEIGVAL
    write(u6,*) 'EIGENVECTORS OF L+ge*S:'
    do I=1,N
      do J=1,N
        write(u6,*) I,J,DEIGVEC(I,J)
      end do
    end do

    !CCCCCCCCCCCCCCC
    ! Test the eigenvectors
    !CCCCCCCCCCCCCCC
    write(u6,*) 'TESTING EIGENVECTORS'

    call ZGEMM('C','N',N,N,N,cOne,DEIGVEC,N,DEIGVEC,N,cZero,BPTST,N)

    write(u6,*) 'V*V: Should be real unit matrix'
    do I=1,N
      do J=1,N
        ISTATE = SODIAG(I)
        JSTATE = SODIAG(J)
        write(u6,*) ISTATE,JSTATE,BPTST(I,J)
      end do
    end do

    !CCCCCCCCCCCCCC
    ! END testing the eigenvectors
    !CCCCCCCCCCCCCCC
  end if ! IPGLOB >= 4

  ! CALL SPIN_PHASE FROM THE OTHER CODE
  !      SUBROUTINE SPIN_PHASE(IPGLOB,DIPSO2,GMAIN,DIM,ZIN,ZOUT)
  ZOUT(:,:) = cZero

  call SPIN_PHASE_RASSI(2,PROP2,GTENS,N,DEIGVEC,ZOUT)

  ! EXPAND EIGENVECTORS TO SEPARATE R,I MATRICES AND
  ! AS A PART OF AN IDENTITY MATRIX
  call mma_allocate(EIGVECR,NSS,NSS,Label='EIGVECR')
  call mma_allocate(EIGVECI,NSS,NSS,Label='EIGVECI')

  call unitmat(EigVecR,NSS)
  EigVecI(:,:) = Zero

  do I=1,N
    do J=1,N
      I2 = SODIAG(I)
      J2 = SODIAG(J)

      IJ = NSS*(J2-1)+I2
      write(u6,*) I2,J2,IJ
      EIGVECR(I2,J2) = real(ZOUT(I,J))
      EIGVECI(I2,J2) = aimag(ZOUT(I,J))
    end do
  end do

  !write(u6,*) 'EIGENVECTORS SPLIT INTO REAL/IMAG PARTS'
  !do I=1,NSS
  !  do J=1,NSS
  !    IJ = NSS*(J-1)+I
  !    write(u6,*) I,J,EIGVECR(I,J)
  !    write(u6,*) I,J,EIGVECI(I,J)
  !  end do
  !end do

  ! Multiply by SO eigenvectors to get new UW matrix

  call mma_allocate(UWR,NSS**2,Label='UWR')
  call mma_allocate(UWI,NSS**2,Label='UWI')

  call DGEMM_('N','N',NSS,NSS,NSS,One,UMATR,NSS,EIGVECR,NSS,Zero,UWR,NSS)
  call DGEMM_('N','N',NSS,NSS,NSS,One,UMATI,NSS,EIGVECR,NSS,Zero,UWI,NSS)
  call DGEMM_('N','N',NSS,NSS,NSS,One,UMATR,NSS,EIGVECI,NSS,One,UWI,NSS)
  call DGEMM_('N','N',NSS,NSS,NSS,-One,UMATI,NSS,EIGVECI,NSS,One,UWR,NSS)

  ! REDO USING SONATORB_MIX

  do I=1,SODIAGNSTATE
    do J=1,SODIAGNSTATE

      ISTATE = SODIAG(I)
      JSTATE = SODIAG(J)

      IJ = J*(J-1)/2+I

      write(u6,*) 'State: ',ISTATE,JSTATE

      ! file name for the spin density orb file
      if (IDIR == 1) FILEBASE = 'SODISDENS.X'
      if (IDIR == 2) FILEBASE = 'SODISDENS.Y'
      if (IDIR == 3) FILEBASE = 'SODISDENS.Z'
      if (IDIR == 1) FILEBASEL = 'SODILDENS.X'
      if (IDIR == 2) FILEBASEL = 'SODILDENS.Y'
      if (IDIR == 3) FILEBASEL = 'SODILDENS.Z'

      ! For L, mix the AO integrals, leave the density alone
      !    -> Call SONATORB then SONATORBM_INT
      ! For S, leave AO integrals alone, mix density matrices
      !    -> Call SONATORBM, SONATORBM_INT

      ! store antising density in LDMATTMP
      iOpt = 0
      call SONATORBM('ANTISING',UWR,UWI,ISTATE,JSTATE,NSS,iOpt,IDENTMAT,DMATTMP)

      ! Expectation values of L -> LMAT{R,I}
      IC = -1
      iOpt = 1
      call SONATORBM_INT(DMATTMP,'ANGMOM',IC,'ANTISING',ISTATE,JSTATE,iOpt,MAXES,LMATR(I,J,IDIR,1),LMATR(I,J,IDIR,2), &
                         LMATR(I,J,IDIR,3),LMATI(I,J,IDIR,1),LMATI(I,J,IDIR,2),LMATI(I,J,IDIR,3))

      ! Plot for generation of current density
      if (IFCURD) call SONATORB_CPLOT(DMATTMP,FILEBASEL,'ANTISING',ISTATE,JSTATE)

      ! store hermtrip density in LDMATTMP
      iOpt = 1
      call SONATORBM('HERMTRIP',UWR,UWI,ISTATE,JSTATE,NSS,iOpt,MAXES,DMATTMP)

      ! Expectation values of S -> SMAT{R,I}
      IC = 1
      iOpt = 0
      call SONATORBM_INT(DMATTMP,'MLTPL  0',IC,'HERMTRIP',ISTATE,JSTATE,iOpt,IDENTMAT,SMATR(I,J,IDIR,1),SMATR(I,J,IDIR,2), &
                         SMATR(I,J,IDIR,3),SMATI(I,J,IDIR,1),SMATI(I,J,IDIR,2),SMATI(I,J,IDIR,3))

      ! plot the rotated density
      call SONATORB_PLOT(DMATTMP,FILEBASE,'HERMTRIP',ISTATE,JSTATE)

    end do
  end do

  ! write the magnetic axes to a file
  write(u6,*) 'Writing magnetic axes to file SODIAG.MAXES'
  LUMAXES = 54
  LUMAXES = IsFreeUnit(LUMAXES)
  call Molcas_Open(LUMAXES,'SODIAG.MAXES')
  write(LUMAXES,*) MAXES
  close(LUMAXES)

  call mma_deallocate(UWI)
  call mma_deallocate(UWR)
  call mma_deallocate(EIGVECI)
  call mma_deallocate(EIGVECR)

end do ! end loop over directions

! Final output
call DCOPY_(9*SODIAGNSTATE**2,LMATR,1,MUMAT2R,1)
call DCOPY_(9*SODIAGNSTATE**2,LMATI,1,MUMAT2I,1)
call DSCAL_(9*SODIAGNSTATE**2,-One,MUMAT2R,1)
call DSCAL_(9*SODIAGNSTATE**2,-One,MUMAT2I,1)
call DAXPY_(9*SODIAGNSTATE**2,-ge,SMATR,1,MUMAT2R,1)
call DAXPY_(9*SODIAGNSTATE**2,-ge,SMATI,1,MUMAT2I,1)

! we are only interested in the 1,1 state
write(u6,*) '-----------------------------------------------------'
write(u6,*) 'Final output after diagonalizing along all directions'
write(u6,'(A9,A15,A15,A15)') '','X','Y','Z'
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Mux)',MUMAT2R(1,1,1,1),MUMAT2R(1,1,2,1),MUMAT2R(1,1,3,1)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Mux)',MUMAT2I(1,1,1,1),MUMAT2I(1,1,2,1),MUMAT2I(1,1,3,1)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Muy)',MUMAT2R(1,1,1,2),MUMAT2R(1,1,2,2),MUMAT2R(1,1,3,2)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Muy)',MUMAT2I(1,1,1,2),MUMAT2I(1,1,2,2),MUMAT2I(1,1,3,2)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Muz)',MUMAT2R(1,1,1,3),MUMAT2R(1,1,2,3),MUMAT2R(1,1,3,3)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Muy)',MUMAT2I(1,1,1,3),MUMAT2I(1,1,2,3),MUMAT2I(1,1,3,3)
write(u6,*)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Lx)',LMATR(1,1,1,1),LMATR(1,1,2,1),LMATR(1,1,3,1)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Lx)',LMATI(1,1,1,1),LMATI(1,1,2,1),LMATI(1,1,3,1)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Ly)',LMATR(1,1,1,2),LMATR(1,1,2,2),LMATR(1,1,3,2)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Ly)',LMATI(1,1,1,2),LMATI(1,1,2,2),LMATI(1,1,3,2)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Lz)',LMATR(1,1,1,3),LMATR(1,1,2,3),LMATR(1,1,3,3)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Ly)',LMATI(1,1,1,3),LMATI(1,1,2,3),LMATI(1,1,3,3)
write(u6,*)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Sx)',SMATR(1,1,1,1),SMATR(1,1,2,1),SMATR(1,1,3,1)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Sx)',SMATI(1,1,1,1),SMATI(1,1,2,1),SMATI(1,1,3,1)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Sy)',SMATR(1,1,1,2),SMATR(1,1,2,2),SMATR(1,1,3,2)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Sy)',SMATI(1,1,1,2),SMATI(1,1,2,2),SMATI(1,1,3,2)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Re(Sz)',SMATR(1,1,1,3),SMATR(1,1,2,3),SMATR(1,1,3,3)
write(u6,'(A9,F20.12,F20.12,F20.12)') 'Im(Sy)',SMATI(1,1,1,3),SMATI(1,1,2,3),SMATI(1,1,3,3)
write(u6,*) '-----------------------------------------------------'

!call ADD_INFO('SODIAG_MUMAT2R',MUMAT2R,9*N*N,4)
!call ADD_INFO('SODIAG_MUMAT2I',MUMAT2R,9*N*N,4)
!call ADD_INFO('SODIAG_LMATR',LMATR,9*N*N,4)
!call ADD_INFO('SODIAG_LMATI',LMATI,9*N*N,4)
!call ADD_INFO('SODIAG_SMATR',SMATR,9*N*N,4)
!call ADD_INFO('SODIAG_SMATI',SMATI,9*N*N,4)

write(u6,*)
write(u6,*) '***********************************************'
write(u6,*) '* ENDING SODIAG *******************************'
write(u6,*) '***********************************************'

call mma_deallocate(DMATTMP)
call mma_deallocate(DEIGVAL)
call mma_deallocate(RWORK)
call mma_deallocate(BPTST)
call mma_deallocate(DEIGVEC)
call mma_deallocate(H_ZEE)
call mma_deallocate(PROP)
call mma_deallocate(PROP2)
call mma_deallocate(ZOUT)
call mma_deallocate(ZWORK)
call mma_deallocate(LMATR)
call mma_deallocate(LMATI)
call mma_deallocate(SMATR)
call mma_deallocate(SMATI)
call mma_deallocate(MUMAT2R)
call mma_deallocate(MUMAT2I)

end subroutine mkSODIAG
