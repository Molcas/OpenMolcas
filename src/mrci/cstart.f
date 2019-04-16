************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CSTART(AREF,EREF,CI,ICI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AREF(NREF,NREF),EREF(NREF),CI(NCONF),ICI(MBUF)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION BUF(nCOP),ISTART(MXROOT)
*
      CALL QENTER('CSTART')
      DO  5 I=1,MXVEC
        IDISKC(I)=-1
        IDISKS(I)=-1
5     CONTINUE
C FIRST, USE THE CI ARRAY TO STORE THE DIAGONAL ELEMENTS:
      IAD25=IAD25S
      DO 10 I=1,NCONF,nCOP
        CALL dDAFILE(Lu_25,2,BUF,nCOP,IAD25)
        NN=MIN(nCOP,NCONF+1-I)
        CALL DCOPY_(NN,BUF,1,CI(I),1)
10    CONTINUE
C THESE ARE DIAGONAL ELEMENTS OF THE ELECTRONIC HAMILTONIAN.
C POTNUC SHOULD BE ADDED. IN ADDITION, WE USE AN ENERGY SHIFT.
C NOTE: DISPLACEMENT 1.0d-4 PROTECTS AGAINST DIVIDE ERRORS.
CPAM: Protect all diag elems. Needed in some weird cases.
C ENERGY SHIFT:
      ESHIFT=EREF(1)
      DO 25 I=1,NCONF
        CI(I)=CI(I)+POTNUC-ESHIFT+1.0D-04
25    CONTINUE
      Call Add_Info('CI_DIAG2',CI(2),1,8)
C REPLACE REFERENCE ENERGIES:
      DO 20 I=1,NREF
        IR=IREFX(I)
        CI(IR)=EREF(I)-ESHIFT-1.0D-04
20    CONTINUE
      IF(ICPF.EQ.1) THEN
        DO 30 IREF=1,NREF
          IR=IREFX(IREF)
          CI(IR)=GFAC*CI(IR)
30      CONTINUE
        GINV=1.0D00/GFAC
        CALL DSCAL_(NCONF,GINV,CI,1)
      END IF
      IDFREE=0
      IDISKD=0
      DO 31 ISTA=1,NCONF,MBUF
        NN=MIN(MBUF,(NCONF+1-ISTA))
        CALL dDAFILE(LUEIG,1,CI(ISTA),NN,IDFREE)
31    CONTINUE
C THEN, SET UP START CI VECTORS IN MCSF BASIS:
      CALL DCOPY_(NCONF,[0.0D00],0,CI,1)
      IF(IREST.EQ.0) THEN
        NNEW=IROOT(NRROOT)
        I1=1
        I2=1
        DO 35 I=1,NNEW
          IF(I.EQ.IROOT(I1)) THEN
            ISTART(NNEW-NRROOT+I1)=I
            I1=I1+1
          ELSE
            ISTART(I2)=I
            I2=I2+1
          END IF
35      CONTINUE
        IF(NNEW.GT.1) THEN
          WRITE(6,*)
     *    ' THE FOLLOWING REFERENCE ROOTS ARE USED AS START VECTORS:'
      CALL XFLUSH(6)
          WRITE(6,'(12(A,I2))') ' ROOTS NR ',ISTART(1),
     *                           (',',ISTART(I),I=2,NNEW-1),
     *                          ', AND ',ISTART(NNEW)
      CALL XFLUSH(6)
          IF(NNEW.GT.NRROOT) THEN
            WRITE(6,*)' (THE FIRST EXTRA ROOT(S) WERE INCLUDED IN'//
     *         'ORDER TO IMPROVE CONVERGENCE)'
      CALL XFLUSH(6)
          END IF
        ELSE
          WRITE(6,'(A,I2,A)') ' ROOT NR ',ISTART(1),
     *                        ' IS USED AS START VECTOR.'
      CALL XFLUSH(6)
        END IF
        DO 40 I=1,NNEW
          ISTA=ISTART(I)
          IR=IREFX(ISTA)
          CI(IR)=1.0D00
          IDISKC(I)=IDFREE
          DO 41 ISTA=1,NCONF,MBUF
            NN=MIN(MBUF,(NCONF+1-ISTA))
            CALL PKVEC(NN,CI(ISTA),ICI)
            CALL iDAFILE(LUEIG,1,ICI,NN,IDFREE)
41        CONTINUE
          CI(IR)=0.0D00
40      CONTINUE
      ELSE
        ID=0
        NNEW=NRROOT
        DO 50 I=1,NRROOT
          CALL dDAFILE(LUREST,2,CI,NCONF,ID)
          CALL CSFTRA('MCSF',CI,AREF)
          IDISKC(I)=IDFREE
          DO 51 ISTA=1,NCONF,MBUF
            NN=MIN(MBUF,(NCONF+1-ISTA))
            CALL PKVEC(NN,CI(ISTA),ICI)
            CALL iDAFILE(LUEIG,1,ICI,NN,IDFREE)
51        CONTINUE
50      CONTINUE
      END IF
      NVTOT=NNEW
      NSTOT=0
      CALL QEXIT('CSTART')
      RETURN
      END
