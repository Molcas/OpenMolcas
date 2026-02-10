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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE TRDNS2D(IVEC,JVEC,DPT2,NDPT2,SCAL)

      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, only: imag_shift, sigma_p_epsilon
      use caspt2_global, only: do_grad
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_global, only: LUSBT, LISTS
      use EQSOLV, only: iDBMat
      use stdalloc, only: mma_allocate, mma_deallocate
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: nSym, nInDep, nISup, nASup

      IMPLICIT NONE

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      integer(kind=iwp), intent(in) :: IVEC,JVEC,NDPT2
      real(kind=wp), intent(inout) :: DPT2(NDPT2)
      real(kind=wp), intent(in) :: SCAL

      real(kind=wp), ALLOCATABLE:: BD(:), ID(:)
#ifdef _MOLCAS_MPP_
      real(kind=wp), ALLOCATABLE:: VEC1(:), VEC2(:)
#endif
      integer(kind=iwp) ICASE, ISYM, jD, lg_v1, lg_v2, NAS, NIN, NIS,
     &                  nVec

C Add to the diagonal blocks of transition density matrix,
C    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>,
C i.e. inactive/inactive, active/active, and virt/virt
C submatrices.IVEC, JVEC stands for the 1st-order perturbed
C CASPT2 wave functions in vectors nr IVEC, JVEC on LUSOLV.

C Inact/Inact and Virt/Virt blocks:
      DO ICASE=1,13
C       if (icase/=12 .and. icase.ne.13) cycle ! H
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN==0) CYCLE
          NIS=NISUP(ISYM,ICASE)
          NVEC=NIN*NIS
          IF(NVEC==0) CYCLE
          !! lg_V1: T+lambda
          !! lg_V2: T
          !! IVEC = iVecX
          !! JVEC = iVecR
          CALL RHS_ALLO(NIN,NIS,lg_V1)
          CALL RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
          IF(IVEC==JVEC) THEN
            lg_V2=lg_V1
          ELSE
            CALL RHS_ALLO(NIN,NIS,lg_V2)
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
            If (do_grad) Then
              CALL RHS_SCAL(NIN,NIS,lg_V1,SCAL)
              if (sigma_p_epsilon /= Zero) then
                !! derivative of the numerator
                nAS = nASUP(iSym,iCase)
                Call mma_allocate(BD,nAS,Label='BD')
                Call mma_allocate(ID,nIS,Label='ID')
                jD = iDBMat(iSym,iCase)
                Call dDaFile(LUSBT,2,BD,nAS,jD)
                Call dDaFile(LUSBT,2,ID,nIS,jD)
                Call CASPT2_ResD(3,nIN,nIS,lg_V2,lg_V1,BD,ID)
                Call mma_deallocate(BD)
                Call mma_deallocate(ID)
              end if
              CALL RHS_DAXPY(NIN,NIS,One,lg_V2,lg_V1)
              CALL RHS_READ_SR(lg_V2,ICASE,ISYM,IVEC)
            End If
          END IF

CSVC: DIADNS can currently not handle pieces of RHS, so pass the
C full array in case we are running in parallel
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            IF (KING()) THEN
              ! copy global array to local buffer
              CALL mma_allocate(VEC1,NVEC,Label='VEC1')
              CALL GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
              IF(IVEC==JVEC) THEN
                CALL DIADNS(ISYM,ICASE,VEC1,VEC1,DPT2,LISTS)
              ELSE
                CALL mma_allocate(VEC2,NVEC,Label='VEC2')
                CALL GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
                CALL DIADNS(ISYM,ICASE,VEC1,VEC2,DPT2,LISTS)
                CALL mma_deallocate(VEC2)
              END IF
              CALL mma_deallocate(VEC1)
            END IF
            CALL GASYNC()
          ELSE
#endif
            CALL DIADNS(ISYM,ICASE,GA_Arrays(lg_V1)%A,
     &                             GA_Arrays(lg_V2)%A,DPT2,LISTS)
#ifdef _MOLCAS_MPP_
          END IF
#endif
          If (do_grad .and. (imag_shift /=Zero
     *                  .or. sigma_p_epsilon /= Zero)) Then
            !! for sigma-p CASPT2, derivative of the denominator
            nAS = nASUP(iSym,iCase)
            Call mma_allocate(BD,nAS,Label='BD')
            Call mma_allocate(ID,nIS,Label='ID')
            jD = iDBMat(iSym,iCase)
            Call dDaFile(LUSBT,2,BD,nAS,jD)
            Call dDaFile(LUSBT,2,ID,nIS,jD)
C
            CALL RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
            Call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,BD,ID)
C
            Call DScal_(NDPT2,-One,DPT2,1)
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              IF (KING()) THEN
                ! copy global array to local buffer
                CALL mma_allocate(VEC1,NVEC,Label='VEC1')
                CALL GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
                IF(IVEC==JVEC) THEN
                  CALL DIADNS(ISYM,ICASE,VEC1,VEC1,DPT2,LISTS)
                ELSE
                  CALL mma_allocate(VEC2,NVEC,Label='VEC2')
                  CALL GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
                  CALL DIADNS(ISYM,ICASE,VEC1,VEC2,DPT2,LISTS)
                  CALL mma_deallocate(VEC2)
                END IF
                CALL mma_deallocate(VEC1)

              END IF
              CALL GASYNC()
            ELSE
#endif
              CALL DIADNS(ISYM,ICASE,GA_Arrays(lg_V1)%A,
     &                               GA_Arrays(lg_V2)%A,DPT2,LISTS)
#ifdef _MOLCAS_MPP_
            END IF
#endif
            Call DScal_(NDPT2,-One,DPT2,1)
            Call mma_deallocate(BD)
            Call mma_deallocate(ID)
          End IF

          CALL RHS_FREE(lg_V1)
          IF(IVEC/=JVEC) CALL RHS_FREE(lg_V2)
        End Do
      End Do
#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par().and.do_grad) call gadsum(DPT2,NDPT2)
#endif

      END SUBROUTINE TRDNS2D
