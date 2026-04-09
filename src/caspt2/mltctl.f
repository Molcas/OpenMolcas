************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MLTCTL(HEFF,EIGVEC,U0)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, Half, One
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: TERSE, USUAL, VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSTATE, IfChol, IFDW, IFRMS, IFXMS, JMS,
     &                         MSTATE,  ENERGY
      IMPLICIT None
      real(kind=wp), intent(inout):: HEFF(NSTATE,NSTATE)
      real(kind=wp), intent(out):: EIGVEC(NSTATE,NSTATE)
      real(kind=wp), Intent(in):: U0(Nstate,Nstate)

      integer(kind=iwp) LAXITY, I, IEND, II0, IJ, ISTA, J, NHTRI, NUMAT
      integer(kind=iwp), external:: Cho_X_GetTol
      CHARACTER(LEN=8) INLAB
      character(len=3) variant
      real(kind=wp), allocatable :: Utmp(:,:), UMAT(:,:), HTRI(:)
      real(kind=wp) DSHIFT

      IF(IPRGLB.GE.TERSE) THEN
        CALL CollapseOutput(1,'Multi-State CASPT2 section:')
        WRITE(u6,'(20A4)')('****',I=1,20)
        WRITE(u6,*)' MULTI-STATE CASPT2 SECTION'
        IF(IPRGLB.GE.USUAL) THEN
          WRITE(u6,'(20A4)')('----',I=1,20)
          WRITE(u6,*)
        END IF
      END IF

C Write out the effective Hamiltonian, for use in e.g. RASSI:
      INLAB='HEFF'
      CALL put_darray(INLAB,HEFF,NSTATE**2)

C Analyze the effective Hamiltonian:
      DSHIFT=Zero
      IF(HEFF(1,1).LE.-100.0D0) THEN
        DSHIFT=-DBLE(INT(-HEFF(1,1)))
      END IF
      DO I=1,NSTATE
        HEFF(I,I)=HEFF(I,I)-DSHIFT
      END DO

      IF(IPRGLB.GE.TERSE .and. DSHIFT.NE.Zero) THEN
        WRITE(u6,*)
     &  ' Output diagonal energies have been shifted. Add ',DSHIFT
       END IF

      IF((IPRGLB.GE.VERBOSE).OR.JMS) THEN
        WRITE(u6,*)' Effective Hamiltonian matrix (Asymmetric):'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          WRITE(u6,*)
          WRITE(u6,'(1x,5I16)')(MSTATE(I),I=ISTA,IEND)
          DO J=1,NSTATE
            WRITE(u6,'(1x,I3,3X,5F16.8)')
     &            MSTATE(J),(HEFF(J,I),I=ISTA,IEND)
          END DO
        END DO
      END IF

C Diagonalize:
C Use a symmetrized matrix, in triangular storage:
      NUMAT=NSTATE**2
      NHTRI=(NUMAT+NSTATE)/2
      CALL mma_allocate(UMAT,NSTATE,NSTATE,LABEL='UMAT')
      CALL mma_allocate(HTRI,NHTRI,LABEL='HTRI')
      IJ=0
      DO I=1,NSTATE
        DO J=1,I
          IJ=IJ+1
          HTRI(IJ)=Half*(HEFF(I,J)+HEFF(J,I))
        END DO
      END DO
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(u6,*)
        WRITE(u6,*)' Effective Hamiltonian matrix (Symmetric):'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          WRITE(u6,*)
          WRITE(u6,'(1x,5I16)')(MSTATE(I),I=ISTA,IEND)
          DO I=ISTA,NSTATE
            II0=(I*(I-1))/2
            WRITE(u6,'(1x,I3,3X,5F16.8)')
     &            MSTATE(I),(HTRI(II0+J),J=ISTA,MIN(I,IEND))
          END DO
        END DO
      END IF
      UMAT(:,:)=Zero
      CALL DCOPY_(NSTATE,[One],0,UMAT,NSTATE+1)
      CALL NIDiag(HTRI,UMAT,NSTATE,NSTATE)
      CALL JACORD(HTRI,UMAT,NSTATE,NSTATE)
      DO I=1,NSTATE
        ENERGY(I)=DSHIFT+HTRI((I*(I+1))/2)
        EIGVEC(:,I)=UMAT(:,I)
      END DO
      CALL mma_deallocate(UMAT)
      CALL mma_deallocate(HTRI)

      IF(IPRGLB.GE.TERSE) THEN
        If (IFRMS) Then
          variant = 'RMS'
        Else if (IFXMS.and.IFDW) then
          variant = 'XDW'
        Else if (IFXMS) then
          variant = 'XMS'
        Else if (IFDW) then
          variant = 'DW '
        Else
          variant = 'MS '
        End If
          WRITE(u6,*)
          WRITE(u6,'(6X,A,A)')' Total ',trim(variant)//
     &      '-CASPT2 energies:'
          DO I=1,NSTATE
            Call PrintResult(u6,'(6x,A,I3,5X,A,F16.8)',trim(variant)//
     &      '-CASPT2 Root',I,'Total energy:',ENERGY(I),1)
          END DO
      END IF

      IF(IPRGLB.GE.USUAL) THEN
        WRITE(u6,*)
        WRITE(u6,'(6X,A)')' Eigenvectors:'
        DO ISTA=1,NSTATE,5
          IEND=MIN(ISTA+4,NSTATE)
          DO J=1,NSTATE
            WRITE(u6,'(6x,5F16.8)')(EIGVEC(J,I),I=ISTA,IEND)
          END DO
          WRITE(u6,*)
        END DO
        if (IFXMS.or.IFRMS) then
* Transform eigenvectors into the original input basis
          call mma_allocate(Utmp,Nstate,Nstate,Label='Utmp')
          call dgemm_('N','N',Nstate,Nstate,Nstate,
     &                One,U0,Nstate,eigvec,Nstate,
     &                Zero,Utmp,Nstate)
          WRITE(u6,'(6X,A)')' In terms of the input states:'
          DO ISTA=1,NSTATE,5
            IEND=MIN(ISTA+4,NSTATE)
            DO J=1,NSTATE
              WRITE(u6,'(6x,5F16.8)')(Utmp(J,I),I=ISTA,IEND)
            END DO
            WRITE(6,*)
          END DO
          call mma_deallocate(Utmp)
        end if
        CALL CollapseOutput(0,'Multi-State CASPT2 section:')
        WRITE(u6,*)
      END IF

* Restore original effective Hamiltonian
      DO I=1,NSTATE
        HEFF(I,I)=HEFF(I,I)+DSHIFT
      END DO

* In automatic verification calculations, the precision is lower
* in case of Cholesky calculation.
      LAXITY=8
      IF(IfChol) LAXITY=Cho_X_GetTol(LAXITY)
      Call Add_Info('E_MSPT2',ENERGY,nState,LAXITY)

      END SUBROUTINE MLTCTL
