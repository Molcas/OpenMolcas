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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE FMAT_CASPT2(FIFA,nFIFA,FIMO,NFIMO,DREF,NDREF,
     &                       HONE,nHONE)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, Half, One, Two
      use caspt2_global, only: LUINTM
      use caspt2_module, only: NSYM, NORB, NISH, NOSH, NAES, NoMx, NoTri
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT None
      integer(kind=iwp), intent(in):: nFIFA, NFIMO, NDREF, nHONE
      real(kind=wp), intent(inout):: FIFA(nFIFA), FIMO(NFIMO)
      real(kind=wp), intent(in)::DREF(NDREF), HONE(nHONE)

      integer(kind=iwp) IAD2M(3,36*36)
      integer(kind=iwp) NDIM2M, IDISK, IFSTA, ISYR, NBR, NB3, NBNB,
     &                  ISYS, IS3RS, ISYP, NIP, NOP, NAESP, ISYQ, IS3PQ,
     &                  ISADDR, IDISK1, IT, ITABS, ITU, IU, IUABS,
     &                  nInts, nBUF
      real(kind=wp) DTU
      real(kind=wp), Allocatable:: BUF(:)

C COMPUTE THE INACTIVE FOCK MATRIX IN MO BASIS.
C COMPUTE THE ACTIVE FOCK MATRIX IN MO BASIS.
C DEFINITIONS ARE:
C FIMO(P,Q)= H(P,Q)+ SUM( (2*(PQ,KK)-(PK,QK)) OVER INACTIVE K)
C WHERE H(P,Q) IS THE CORE FOCK MATRIX OBTAINED FROM THE
C TRANSFORMATION PART, I.E., THE CI 1-ELECTRON HAMILTONIAN.
C SIMILARLY,
C FAMO(P,Q)= (1/2)*SUM( (2*(PQ,TU)-(PU,QT))*DREF(T,U) OVER ACTIVE T,U)
C THIS ROUTINE USES DIRECTLY THE TRANSFORMED INTEGRALS TO SET UP
C FIMO AND FAMO.
C CODED 92-12-04 BY MALMQVIST FOR CASPT2, MOLCAS-3 VERSION.

c Inactive and active Fock matrices:
      FIMO(:)=Zero
      FIFA(:)=Zero ! ued temporarily as FAMO

c notri=Size of an array with symmetry-blocked triangular
c submatrices, using non-frozen, non-deleted MO indices.
c NBUF=Max size of a LUINTM buffer.
      NBUF=MAX(NOMX**2,notri)
      CALL mma_allocate(BUF,NBUF,Label='BUF')

      NDIM2M=(NSYM*(NSYM+1))/2
      IDISK=0
      CALL IDAFILE(LUINTM,2,IAD2M,3*36*36,IDISK)

      Call Do_Loops(1) !     Do Coulomb contributions
      Call Do_Loops(2) !     Do exchange contributions

      CALL mma_deallocate(BUF)

* FIMO comes from contractions over inactive orbitals, while FAMO from
* contractions over active orbitals and therefore are summed up together
* here into FIFA.

      FIMO(1:notri) = FIMO(:) + HONE(:)
      FIFA(1:notri) = FIFA(:) + FIMO(:)

      Contains

      Subroutine Do_Loops(icase)
      integer(kind=iwp), intent(in):: iCase

      IFSTA=1
      DO ISYR=1,NSYM
        NBR=NORB(ISYR)
        IF(NBR.EQ.0) CYCLE
        NB3=(NBR**2+NBR)/2
        NBNB=NBR**2
        ISYS=ISYR
        IS3RS=(ISYR**2+ISYR)/2
        DO ISYP=1,NSYM
          NIP=NISH(ISYP)
          NOP=NOSH(ISYP)
          NAESP=NAES(ISYP)
          IF(NOP.EQ.0) CYCLE
          ISYQ=ISYP
          IS3PQ=(ISYP*(ISYP+1))/2
          ISADDR=IS3RS+NDIM2M*(IS3PQ-1)
          IDISK1=IAD2M(iCase,ISADDR)
          IF (IDISK1.EQ.0) THEN
             If (iCase==1) Then
               WRITE(u6,*)' FMAT: COULOMB INTEGRAL BUFFER MISSING!'
             Else
               WRITE(u6,*)' FMAT: EXCH-1  INTEGRAL BUFFER MISSING!'
             End If
             WRITE(u6,'(1X,A,4I4)')'SYMMETRY BLOCK:',ISYP,ISYQ,ISYR,ISYS
             CALL ABEND()
          END IF
          If (iCase==1) Then
             nInts=NB3
          ELSE
             nInts=NBNB
          END IF

          DO IT=1,NOP
            DO IU=1,IT

C Process here contributions to FIMO or FAMO.

              IF(IT.LE.NIP .AND. IT.EQ.IU) THEN
                CALL DDAFILE(LUINTM,2,BUF,nInts,IDISK1)
                If (iCase==1) Then
C ADD 2* BUFFER INTO CORRECT POSITION OF  FIMO.
                   CALL DAXPY_(NB3,Two,BUF,1,FIMO(IFSTA),1)
                ELSE
                   CALL TRIANG(NBR,BUF)
C ADD -1* BUFFER INTO CORRECT POSITION OF  FIMO.
                   CALL DAXPY_(NB3,-One,BUF,1,FIMO(IFSTA),1)
                End If
              ELSE IF(IT.GT.NIP .AND. IT.LE.NOP .AND.
     &                IU.GT.NIP .AND. IU.LE.NOP) THEN
                ITABS=NAESP+(IT-NIP)
                IUABS=NAESP+(IU-NIP)
                ITU=(ITABS*(ITABS-1))/2 + IUABS
                DTU=DREF(ITU)
                IF(IT.EQ.IU) DTU=half*DTU
                CALL DDAFILE(LUINTM,2,BUF,nInts,IDISK1)
                If (iCase==1) Then
                   CALL DAXPY_(NB3,Two*DTU,BUF,1,FIFA(IFSTA),1)
                Else
                   CALL TRIANG(NBR,BUF)
                   CALL DAXPY_(NB3,-DTU,BUF,1,FIFA(IFSTA),1)
                End If
              ELSE
                CALL DDAFILE(LUINTM,0,BUF,nInts,IDISK1)
              END IF

            END DO
          END DO

        END DO
        IFSTA=IFSTA+NB3
      END DO
      End Subroutine Do_Loops

      END SUBROUTINE FMAT_CASPT2
