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
* Copyright (C) 2004, Teodoro Laino                                    *
************************************************************************
C     ******************************************************************C
      SUBROUTINE COMP_NAC(ISTATE, JSTATE, SCR, nSCR,ISY12, IOFF,LCI1)
C***********************************************************************C
C COMP_NAC : This routine was created to compute NonAdiabatic Couplings C
C                                                                       C
C                      d                                                C
C    NAC =   <Psi_1 | ---- Psi_2>                                       C
C                      dR                                               C
C                                                                       C
C    where Psi_1 and Psi_2 are two MCSCF WaveFunctions                  C
C                                                                       C
C Author:     Teodoro Laino.  Scuola Normale Superiore di Pisa (Italy)  C
C Date:       27.02.2004  Lund (Sweden)                                 C
C                                                                       C
C The original idea to implement this NAC computation was due to        C
C Per-Ake Malmqvist on January 2000                                     C
C                                                                       C
C On Input:                                                             C
C             SCR:  Transition density matrice in AO basis.             C
C                                                                       C
C***********************************************************************C
      Use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep
      IMPLICIT REAL*8 (A-H,O-Z)

#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='COMP_NAC')
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "disp.fh"
#include "diff.fh"
#include "symmul.fh"
* Arguments
      REAL*8 SCR(nSCR)
      DIMENSION IOFF(*)
      Integer IndGrd(0:7)
      Logical, External :: TF
*
* Main Loop on all  geometrical displacements to perform
* calculation of NonAdiabatic Couplings
*
*
      NTEMP = 3 * nAlCnt
      CALL GETMEM('NACS','ALLO','REAL',NNACS,NTEMP)
      CALL DCOPY_(NTEMP, [0.D0], 0, WORK(NNACS), 1)
      idcnt = 0
      Nprop = 0
      Job   = Jstate
      Do ICnttp=1,nCnttp
         Do ICnt=1,dbsc(iCnttp)%nCntr
            IdCnt=IdCnt+1
            Do IdCar=1,3
               Nprop=Nprop+1
*
* Compute idisp to retrieve data from disk...
*
* Differentiated symmetry-unique center IDCNT
* Derivative wrt component IDCAR=1,2,3 (d/dx,d/dy,d/dz)
* INDDSP(IDCNT,IIRREP) is the number of displacements in
* earlier center/irrep. Thus it is an offset.
               Call iCopy(nIrrep,[0],0,IndGrd,1)
               lOper=0
               Do iIrrep=0,nIrrep-1
                  nDisp = IndDsp(iDcnt,iIrrep)
* First set NDISP=ordering number of this displacement.
* Then loop over directions d/dx,d/dy,d/dz
                  Do iCar=1,3
                     iComp = 2**(iCar-1)
                     If ( TF(IdCnt,iIrrep,iComp)) Then
                        ndisp=ndisp+1
* NDISP is now the ordering number of this displacement.
                        If (iDCar.eq.icar) Then
                           lOper=lOper+2**iIrrep
                           IndGrd(iIrrep) = nDisp
                        End If
                     End If
                  End Do
               End Do
               If (lOper.ne.0) Then
* For the displacement represented by this symmetry-unique
* center IDCNT and this component IDCAR, the differentiation
* operator has components with irreps that have been marked
* with '1' in LOPER, regarded as a flag array.
* INDGRD(IIRREP) will be zero, except for those irreps, and
* will then contain the ordering number of the displacement.
                  Do iIrrep = 0, nIrrep-1
                     iSmLbl = 2**iIrrep
                     If ((iAnd(2**iIrrep,lOper).ne.0).and.
     &                    (MUL(iIrrep+1,isy12).eq.1))  Then
                        idisp = indgrd(iirrep)
                        CALL COMP_NAC_IDISP(Job,
     &                                      idisp,
     &                                      iIrrep+1,
     &                                      isy12,
     &                                      SCR,
     &                                      Work(LCI1),
     &                                      Prop,
     &                                      Ioff)
                        Work(NNACS-1+idisp) = Prop
                     End If
                  End Do
               End If
            End Do
         End Do
      End Do
* Print Some Stuff related to the computation of NACs
      IF (IPGLOB.ge.TERSE) THEN
      WRITE(6,'(/,"NONADIABATIC COUPLINGS BETWEEN STATE",'//
     &          'I5,"AND STATE",I5," .",/)')ISTATE, JSTATE
      Do I = 1, nAlCnt*3
         WRITE(6,'(I5,F15.9)')I,Work(NNACS-1+I)
      End Do
      END IF
*
      CALL GETMEM('NACS','FREE','REAL',NNACS,NTEMP)
*
      RETURN
C     ******************************************************************C
      END
C     ******************************************************************C


C     ******************************************************************C
      SUBROUTINE COMP_NAC_IDISP(JOB, IDISP, ISYMP, ISY12, TDMZZ,
     &                          CISTATJ, PROPVAL, IOFF )
C***********************************************************************C
C COMP_NAC : This routine was created to compute NonAdiabatic Couplings C
C                                                                       C
C                      d                                                C
C    NAC =   <Psi_1 | ---- Psi_2>                                       C
C                      dR                                               C
C                                                                       C
C    where Psi_1 and Psi_2 are two MCSCF WaveFunctions                  C
C                                                                       C
C Author:     Teodoro Laino.  Scuola Normale Superiore di Pisa (Italy)  C
C Date:       27.02.2004  Lund (Sweden)                                 C
C                                                                       C
C The original idea to implement this NAC computation was due to        C
C Per-Ake Malmqvist on January 2000                                     C
C                                                                       C
C Files needed by this routine:                                         C
C                                                                       C
C          File             Label        Description                    C
C                                                                       C
C     -   MCKINT                                                        C
C                       1) OVRGRDA - "ANTISYMMETRIC DERIVATIVE"         C
C                       2) KAPPA   - "ORBITAL ROTATION MATRIX "         C
C                       3) CI      - "CI ARRAY DERIVATIVES    "         C
C                                                                       C
C***********************************************************************C
      IMPLICIT REAL*8 (A-H,O-Z)

#include "Molcas.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "Files.fh"
* Arguments
      DIMENSION TDMZZ(*), CISTATJ(*), IOFF(*)
* Local Variables
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='COMP_NAC')
      CHARACTER*8 LABEL, STYPE
* Subroutine Statements


*
* Reading information from MCKINT file...
*
* Let's Compute the dimension of the 1-electron array
*
      ISUM  = 0
      DO IS=1,NSYM
       JS=MUL(IS,ISYMP)
       NBI=NBASF(IS)
       NBJ=NBASF(JS)
       NBIJ=NBI*NBJ
       ISUM=ISUM+NBIJ
      END DO
      NTEMP=ISUM
      NSCR=(NBST*(NBST+1))/2
*
* OVRGRDA...
*
      LABEL ='OVRGRDA '
      STYPE ='ANTI    '
      CALL GETMEM('OVRGRDA','ALLO','REAL',N_OVRGRDA,NTEMP)
      CALL RDMGRD(JOB,IDISP,LABEL,STYPE,ISYMP,NTEMP,WORK(N_OVRGRDA))
*
* KAPPA...
*
      LABEL ='KAPPA   '
      CALL GETMEM('KAPPA','ALLO','REAL',N_KAPPA,NTEMP)
      CALL RDMGRD(JOB,IDISP,LABEL,STYPE,ISYMP,NTEMP,WORK(N_KAPPA))
*
* CI VECTORS...
*
      LABEL ='CI      '
      CALL GETMEM('DCIVEC','ALLO','REAL',N_DCIVEC,NCONF1)
      CALL RDMCCI(JOB,IDISP,LABEL,ISYMP,NCONF1,WORK(N_DCIVEC))
*
* Compute NACs w.r.t. the IDISP displacement derivative
*
      CALL GETMEM('XMATRIX','ALLO','REAL',N_XMAT,NTEMP)
* Summ up OVRGRDA AND KAPPA contribution in XMATRIX
*
      DO I = 1, NTEMP
         WORK(N_XMAT-1+I) = 0.5D0 *  WORK(N_OVRGRDA-1+I) +
     &        WORK(N_KAPPA-1+I)
      END DO
*
      IINT  = 1
      PSUM  = 0.0D0
      ITYPE = 2
      DO  ISY1 = 1, NSYM
         NB1=NBASF(ISY1)
         IF (NB1.NE.0) THEN
            DO ISY2= 1, ISY1
               I12=MUL(ISY1,ISY2)
C     IF(IAND(2**(I12-1),ISCHK).NE.0) THEN
               NB2=NBASF(ISY2)
               IF(NB2.NE.0) THEN
                  NB12=NB1*NB2
                  IF(ISY1.EQ.ISY2) NB12=(NB12+NB1)/2
                  IF(I12.EQ.ISY12) THEN
                     IPOS=IOFF(ISY1)+1+NSCR*(ITYPE-1)
                     PSUM=PSUM+
     &               DDOT_(NB12,WORK(N_XMAT-1+IINT),1,TDMZZ(IPOS),1)
* write(6,*)'ISY1, ISY2, PSUM',ISY1, ISY2, PSUM
                  END IF
                  IINT=IINT+NB12
               END IF
C     END IF
            END DO
         END IF
      END DO
* Compute CI contribution to NACs
      IF(IPGLOB.ge.DEBUG) THEN
      write(6,*)
      write(6,*)'PSUM, CIcon',PSUM,
     &  DDOT_(NCONF1, CISTATJ, 1, WORK(N_DCIVEC), 1)
      END IF
*      do ii=1,nconf1
*         write(6,*)CISTATJ(ii),WORK(N_DCIVEC+ii-1)
*      end do
      PSUM = PSUM + DDOT_(NCONF1, CISTATJ, 1, WORK(N_DCIVEC), 1)
      PROPVAL = PSUM
*
* Deallocate Vectors
*
      CALL GETMEM('XMATRIX','FREE','REAL',N_XMAT,    NTEMP)
      CALL GETMEM('DCIVEC', 'FREE','REAL',N_DCIVEC,  NCONF1)
      CALL GETMEM('KAPPA',  'FREE','REAL',N_KAPPA,   NTEMP)
      CALL GETMEM('OVRGRDA','FREE','REAL',N_OVRGRDA, NTEMP)
*
* Now you can leave...
*
      RETURN
      END
