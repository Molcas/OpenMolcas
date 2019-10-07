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
      SUBROUTINE SODIAG(UMATR, UMATI, NSS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "Morsel.fh"
#include "Struct.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "rasdef.fh"
#include "jobin.fh"
#include "symmul.fh"
#include "constants.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SODIAG')

C subroutine arguments
      REAL*8 UMATR(NSS,NSS),UMATI(NSS,NSS)


      COMPLEX*16 PROP(3,SODIAGNSTATE,SODIAGNSTATE)
      COMPLEX*16 PROP2(3,SODIAGNSTATE,SODIAGNSTATE)
      REAL*8 DEIGVAL(SODIAGNSTATE)
      COMPLEX*16 DEIGVEC(SODIAGNSTATE,SODIAGNSTATE)
      COMPLEX*16 BPTST(SODIAGNSTATE,SODIAGNSTATE)

      REAL*8 GTENS(3),MAXES(3,3),MAXES2(3,3)
      COMPLEX*16 H_ZEE(SODIAGNSTATE,SODIAGNSTATE)
      COMPLEX*16 ZOUT(SODIAGNSTATE,SODIAGNSTATE)

      REAL*8 RWORK(3*SODIAGNSTATE-2)
      COMPLEX*16 ZWORK(2*SODIAGNSTATE-1)

      REAL*8 IDENTMAT(3,3)
      REAL*8 LMATR(SODIAGNSTATE,SODIAGNSTATE,3,3)
      REAL*8 LMATI(SODIAGNSTATE,SODIAGNSTATE,3,3)
      REAL*8 SMATR(SODIAGNSTATE,SODIAGNSTATE,3,3)
      REAL*8 SMATI(SODIAGNSTATE,SODIAGNSTATE,3,3)
      REAL*8 MUMAT2R(SODIAGNSTATE,SODIAGNSTATE,3,3)
      REAL*8 MUMAT2I(SODIAGNSTATE,SODIAGNSTATE,3,3)

      REAL*8 MU_BOHR

C For creating the filename of the ORB file
      CHARACTER*11 FILEBASE
      CHARACTER*11 FILEBASEL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Matrices
C
C PROP, PROP2 Storage for expectation values (all directions)
C UMATR,I     Transformation matrix (real, imaginary parts)
C NSS         Number of spin states
C DEIGVAL     Storage for eigenvalues (REAL!)
C DEIGVEC     Storage for eigenvectors
C BPTST       Storage for some testing

      ge=-(CONST_ELECTRON_G_FACTOR_)
      MU_BOHR=CONST_BOHR_MAGNETON_IN_SI_*
     &        (CONV_AU_TO_CM1_/CONV_AU_TO_KJ_/1.0D3) ! in cm-1/T

      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) "***********************************************"
      WRITE(6,*) "* STARTING SODIAG *****************************"
      WRITE(6,*) "***********************************************"
      WRITE(6,*)
      WRITE(6,*)

      N=SODIAGNSTATE
      WRITE(6,*) "Number of states to be diagonalized: ",N
      WRITE(6,*) "STATES:"
      DO I = 1,N
        WRITE(6,*) IWORK(LSODIAG-1+I)
      END DO

      CALL DCOPY_(9*N**2,[0.0d0],0,LMATR,1)
      CALL DCOPY_(9*N**2,[0.0d0],0,LMATI,1)
      CALL DCOPY_(9*N**2,[0.0d0],0,SMATR,1)
      CALL DCOPY_(9*N**2,[0.0d0],0,SMATI,1)
      CALL DCOPY_(9*N**2,[0.0d0],0,MUMAT2R,1)
      CALL DCOPY_(9*N**2,[0.0d0],0,MUMAT2I,1)

      CALL GETMEM('DMATTMPA','ALLO','REAL',LDMATTMP,3*(NBST*(NBST+1)))

      !> identity mat
      IDENTMAT(:,:)=0.0D0
      FORALL (I=1:3) IDENTMAT(I,I)=1.0D0

C First, we calculate the expectation values of
C  (L+ge*S)x (L+ge*S)y (L+ge*S)z
C These get stored in PROP()

C Only work with one triangle - this is a hermitian matrix
      DO J=1,N
      DO I=1,J

        ISTATE=IWORK(LSODIAG-1+I)
        JSTATE=IWORK(LSODIAG-1+J)
        WRITE(6,*) "States: ",ISTATE,JSTATE

        CALL SONATORB('ANTISING',UMATR,UMATI,
     &       ISTATE,JSTATE,NSS,WORK(LDMATTMP))

        IC=-1
        iOpt=1
        CALL SONATORBM_INT(WORK(LDMATTMP),'ANGMOM  ',IC,'ANTISING',
     &                    ISTATE,JSTATE,NSS,iOpt,IDENTMAT,
     &                    AXR,AYR,AZR,AXI,AYI,AZI)


        CALL SONATORBM('HERMTRIP',UMATR,UMATI,
     &       ISTATE,JSTATE,NSS,IDENTMAT,
     &       WORK(LDMATTMP))

        IC=1
        iOpt=0
        CALL SONATORBM_INT(WORK(LDMATTMP),'MLTPL  0',IC,'HERMTRIP',
     &                    ISTATE,JSTATE,NSS,iOpt,IDENTMAT,
     &                    SXR,SYR,SZR,SXI,SYI,SZI)


c The first index of PROP is the direction
        PROP(1,I,J)=-1.0d0*DCMPLX(AXR+ge*SXR,AXI+ge*SXI)
        PROP(2,I,J)=-1.0d0*DCMPLX(AYR+ge*SYR,AYI+ge*SYI)
        PROP(3,I,J)=-1.0d0*DCMPLX(AZR+ge*SZR,AZI+ge*SZI)

        IF(I.NE.J) THEN
          PROP(1,J,I)=DCONJG(PROP(1,I,J))
          PROP(2,J,I)=DCONJG(PROP(2,I,J))
          PROP(3,J,I)=DCONJG(PROP(3,I,J))
        END IF


      END DO
      END DO

C      CALL ADD_INFO("SODIAG_PROP",PROP,3*SODIAGNSTATE**2,4)


c Calculate the atens as in single_aniso
      CALL ATENS_RASSI(PROP,N,GTENS,MAXES,IPGLOB)


      do l=1,3
         do i=1,N
            do j=1,N
      PROP2(l,i,j)=(0.d0,0.d0)
               do k=1,3
      PROP2(l,i,j) = PROP2(l,i,j) + PROP(k,i,j) * maxes(k,l)
               enddo
            enddo
         enddo
      enddo


      call atens_RASSI(PROP2, N, GTENS, MAXES2, 2)

c Diagonalize along each direction
C LOOP OVER THE DIRECTIONS
      DO IDIR=1,3



c  apply the magnetic field along the main iDir axis
      do I=1,N
        DEIGVAL(I)=0.d0

        do J=1,N
          H_ZEE(I,J)=(0.d0,0.d0)
          DEIGVEC(I,J)=(0.d0,0.d0)
        enddo
      enddo

      do I=1,N
        do J=1,N
          H_ZEE(I,J)=H_ZEE(I,J) - MU_BOHR * PROP2(IDIR,I,J)
        enddo
      enddo

      IF(IPGLOB.GE.DEBUG) THEN
        WRITE(6,*) "BP: H_ZEE: "
        WRITE(6,*) H_ZEE
        WRITE(6,*) "BP: PROP: "
        WRITE(6,*) "L+2S in main axis direction (after rotation)",IDIR
        DO I=1,N
        DO J=1,N

          WRITE(6,*) "I,J",IWORK(LSODIAG-1+I),IWORK(LSODIAG-1+J),
     &                PROP2(IDIR, I,J)

        END DO
        END DO

        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*) "Direction: ",IDIR
        WRITE(6,*) "H_ZEE Before Diagonalization"

        DO I=1,N
        DO J=1,N
          ISTATE=IWORK(LSODIAG-1+I)
          JSTATE=IWORK(LSODIAG-1+J)
          WRITE(6,*) ISTATE,JSTATE,H_ZEE(I,J)
        END DO
        END DO
      END IF

c DIAGONALIZE
      lcwork = (2*n-1); info = 0
      call zheev('V','U',n,h_zee,n,deigval,zwork,lcwork,
     &           rwork,info)

      !> put eigenvectors in deigvec
      call zcopy(n**2,h_zee,1,deigvec,1)

      IF(INFO.NE.0) THEN
        WRITE(6,*) "DIAGONALIZATION FAILED! ERROR: ",INFO
        CALL ABEND()
      END IF


      IF(IPGLOB.GE.DEBUG) THEN
        WRITE(6,*) "EIGENVALUES OF L+ge*S in direction: ",IDIR
        WRITE(6,*) DEIGVAL
        WRITE(6,*) "EIGENVECTORS OF L+ge*S: "
        DO I=1,N
        DO J=1,N
          WRITE(6,*) I,J,DEIGVEC(I,J)
        END DO
        END DO


CCCCCCCCCCCCCCCC
C Test the eigenvectors
CCCCCCCCCCCCCCCC
      WRITE(6,*) "TESTING EIGENVECTORS"


        CALL ZGEMM('C','N',N,N,N,
     &           (1.0d0,0.0d0),DEIGVEC,N,DEIGVEC,N,(0.0d0,0.0d0),
     &            BPTST,N)


        WRITE(6,*) "V*V: Should be real unit matrix"
        DO I=1,N
        DO J=1,N
          ISTATE=IWORK(LSODIAG-1+I)
          JSTATE=IWORK(LSODIAG-1+J)
          WRITE(6,*) ISTATE,JSTATE,BPTST(I,J)
        END DO
        END DO

CCCCCCCCCCCCCCC
C END testing the eigenvectors
CCCCCCCCCCCCCCCC
      END IF ! IPGLOB >= DEBUG


C CALL SPIN_PHASE FROM THE OTHER CODE (at the bottom of this file)
c      SUBROUTINE SPIN_PHASE(IPGLOB,DIPSO2,GMAIN,DIM,ZIN,ZOUT)
      do i=1,N
        do j=1,N
        ZOUT(i,j) = (0.d0,0.d0)
        enddo
      enddo

      call SPIN_PHASE_RASSI(2,PROP2,GTENS,N,DEIGVEC,ZOUT)


c EXPAND EIGENVECTORS TO SEPARATE R,I MATRICES AND
c AS A PART OF AN IDENTITY MATRIX
      CALL GETMEM('SODEIGR','ALLO','REAL',LEIGVECR,NSS**2)
      CALL GETMEM('SODEIGI','ALLO','REAL',LEIGVECI,NSS**2)

      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LEIGVECR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LEIGVECI),1)

      DO I=1,NSS
      DO J=1,NSS
        IJ=NSS*(J-1)+I
        WORK(LEIGVECI-1+IJ)=0.0d0
        IF(I.EQ.J) THEN
          WORK(LEIGVECR-1+IJ)=1.0d0
        ELSE
          WORK(LEIGVECR-1+IJ)=0.0d0
        END IF
      END DO
      END DO

      DO I=1,N
      DO J=1,N
        I2=IWORK(LSODIAG-1+I)
        J2=IWORK(LSODIAG-1+J)

        IJ=NSS*(J2-1)+I2
        WRITE(6,*) I2,J2,IJ
        WORK(LEIGVECR-1+IJ)=REAL(ZOUT(I,J))
        WORK(LEIGVECI-1+IJ)=AIMAG(ZOUT(I,J))
      END DO
      END DO

c      WRITE(6,*) "EIGENVECTORS SPLIT INTO REAL/IMAG PARTS"
c      DO I=1,NSS
c      DO J=1,NSS
c        IJ=NSS*(J-1)+I
c        WRITE(6,*) I,J,WORK(LEIGVECR-1+IJ)
c        WRITE(6,*) I,J,WORK(LEIGVECI-1+IJ)
c      END DO
c      END DO


c Multiply by SO eigenvectors to get new UW matrix

      CALL GETMEM('SONUWR','ALLO','REAL',LUWR,NSS**2)
      CALL GETMEM('SONUWI','ALLO','REAL',LUWI,NSS**2)

      CALL DGEMM_('N','N',NSS,NSS,NSS,1.0d0,
     &           UMATR,NSS,WORK(LEIGVECR),NSS,0.0d0,WORK(LUWR),NSS)
      CALL DGEMM_('N','N',NSS,NSS,NSS,1.0d0,
     &           UMATI,NSS,WORK(LEIGVECR),NSS,0.0d0,WORK(LUWI),NSS)
      CALL DGEMM_('N','N',NSS,NSS,NSS,1.0d0,
     &           UMATR,NSS,WORK(LEIGVECI),NSS,1.0d0,WORK(LUWI),NSS)
      CALL DGEMM_('N','N',NSS,NSS,NSS,-1.0d0,
     &           UMATI,NSS,WORK(LEIGVECI),NSS,1.0d0,WORK(LUWR),NSS)


c REDO USING SONATORB_MIX

      DO I=1,SODIAGNSTATE
      DO J=1,SODIAGNSTATE

        ISTATE=IWORK(LSODIAG-1+I)
        JSTATE=IWORK(LSODIAG-1+J)

        IJ=J*(J-1)/2+I

        WRITE(6,*) "State: ",ISTATE,JSTATE

c file name for the spin density orb file
        IF(IDIR.EQ.1) FILEBASE='SODISDENS.X'
        IF(IDIR.EQ.2) FILEBASE='SODISDENS.Y'
        IF(IDIR.EQ.3) FILEBASE='SODISDENS.Z'
        IF(IDIR.EQ.1) FILEBASEL='SODILDENS.X'
        IF(IDIR.EQ.2) FILEBASEL='SODILDENS.Y'
        IF(IDIR.EQ.3) FILEBASEL='SODILDENS.Z'


C For L, mix the AO integrals, leave the density alone
C    -> Call SONATORB then SONATORBM_INT
C For S, leave AO integrals alone, mix density matrices
c    -> Call SONATORBM, SONATORBM_INT


c store antising density in LDMATTMP
        CALL SONATORB('ANTISING',WORK(LUWR),WORK(LUWI),
     &       ISTATE,JSTATE,NSS,WORK(LDMATTMP))


c Expectation values of L -> LMAT{R,I}
        IC=-1
        iOpt=1
        CALL SONATORBM_INT(WORK(LDMATTMP),'ANGMOM  ',IC,'ANTISING',
     &                    ISTATE,JSTATE,NSS,iOpt,MAXES,
     &                    LMATR(I,J,IDIR,1),LMATR(I,J,IDIR,2),
     &                    LMATR(I,J,IDIR,3),
     &                    LMATI(I,J,IDIR,1),LMATI(I,J,IDIR,2),
     &                    LMATI(I,J,IDIR,3))

c Plot for generation of current density
        IF(IFCURD) THEN
          CALL SONATORB_CPLOT(WORK(LDMATTMP),FILEBASEL,'ANTISING',
     &                        ISTATE,JSTATE)
        END IF


c store hermtrip density in LDMATTMP
        CALL SONATORBM('HERMTRIP',WORK(LUWR),WORK(LUWI),
     &       ISTATE,JSTATE,NSS,MAXES,WORK(LDMATTMP))

c Expectation values of S -> SMAT{R,I}
        IC=1
        iOpt=0
        CALL SONATORBM_INT(WORK(LDMATTMP),'MLTPL  0',IC,'HERMTRIP',
     &                     ISTATE,JSTATE,NSS,iOpt,IDENMAT,
     &                     SMATR(I,J,IDIR,1),SMATR(I,J,IDIR,2),
     &                     SMATR(I,J,IDIR,3),
     &                     SMATI(I,J,IDIR,1),SMATI(I,J,IDIR,2),
     &                     SMATI(I,J,IDIR,3))

c plot the rotated density
        CALL SONATORB_PLOT(WORK(LDMATTMP),FILEBASE,'HERMTRIP',
     &                     ISTATE,JSTATE)

      END DO
      END DO

c write the magnetic axes to a file
      write(6,*) "Writing magnetic axes to file SODIAG.MAXES"
      LUMAXES=54
      LUMAXES=IsFreeUnit(LUMAXES)
      CALL Molcas_Open(LUMAXES,'SODIAG.MAXES')
      WRITE(LUMAXES,*) MAXES
      CLOSE(LUMAXES)

      CALL GETMEM('SONUWR','FREE','REAL',LUWR,NSS**2)
      CALL GETMEM('SONUWI','FREE','REAL',LUWI,NSS**2)
      CALL GETMEM('SODEIGR','FREE','REAL',LEIGVECR,2*N**2)
      CALL GETMEM('SODEIGI','FREE','REAL',LEIGVECI,2*N**2)

      END DO ! end loop over directions


c Final output
      CALL DCOPY_(9*SODIAGNSTATE**2,LMATR,1,MUMAT2R,1)
      CALL DCOPY_(9*SODIAGNSTATE**2,LMATI,1,MUMAT2I,1)
      CALL DSCAL_(9*SODIAGNSTATE**2,-1.0d0,MUMAT2R,1)
      CALL DSCAL_(9*SODIAGNSTATE**2,-1.0d0,MUMAT2I,1)
      CALL DAXPY_(9*SODIAGNSTATE**2,-1.0d0*ge,SMATR,1,MUMAT2R,1)
      CALL DAXPY_(9*SODIAGNSTATE**2,-1.0d0*ge,SMATI,1,MUMAT2I,1)

C we are only interested in the 1,1 state
      WRITE(6,*) '-----------------------------------------------------'
      WRITE(6,*) "Final output after diagonalizing along all directions"
      WRITE(6,'(A9,A15,A15,A15)') '','X','Y','Z'
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Mux)',
     &             MUMAT2R(1,1,1,1),MUMAT2R(1,1,2,1),MUMAT2R(1,1,3,1)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Mux)',
     &             MUMAT2I(1,1,1,1),MUMAT2I(1,1,2,1),MUMAT2I(1,1,3,1)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Muy)',
     &             MUMAT2R(1,1,1,2),MUMAT2R(1,1,2,2),MUMAT2R(1,1,3,2)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Muy)',
     &             MUMAT2I(1,1,1,2),MUMAT2I(1,1,2,2),MUMAT2I(1,1,3,2)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Muz)',
     &             MUMAT2R(1,1,1,3),MUMAT2R(1,1,2,3),MUMAT2R(1,1,3,3)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Muy)',
     &             MUMAT2I(1,1,1,3),MUMAT2I(1,1,2,3),MUMAT2I(1,1,3,3)
      WRITE(6,*)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Lx)',
     &             LMATR(1,1,1,1),LMATR(1,1,2,1),LMATR(1,1,3,1)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Lx)',
     &             LMATI(1,1,1,1),LMATI(1,1,2,1),LMATI(1,1,3,1)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Ly)',
     &             LMATR(1,1,1,2),LMATR(1,1,2,2),LMATR(1,1,3,2)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Ly)',
     &             LMATI(1,1,1,2),LMATI(1,1,2,2),LMATI(1,1,3,2)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Lz)',
     &             LMATR(1,1,1,3),LMATR(1,1,2,3),LMATR(1,1,3,3)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Ly)',
     &             LMATI(1,1,1,3),LMATI(1,1,2,3),LMATI(1,1,3,3)
      WRITE(6,*)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Sx)',
     &             SMATR(1,1,1,1),SMATR(1,1,2,1),SMATR(1,1,3,1)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Sx)',
     &             SMATI(1,1,1,1),SMATI(1,1,2,1),SMATI(1,1,3,1)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Sy)',
     &             SMATR(1,1,1,2),SMATR(1,1,2,2),SMATR(1,1,3,2)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Sy)',
     &             SMATI(1,1,1,2),SMATI(1,1,2,2),SMATI(1,1,3,2)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Re(Sz)',
     &             SMATR(1,1,1,3),SMATR(1,1,2,3),SMATR(1,1,3,3)
      WRITE(6,'(A9,F20.12,F20.12,F20.12)')'Im(Sy)',
     &             SMATI(1,1,1,3),SMATI(1,1,2,3),SMATI(1,1,3,3)
      WRITE(6,*) '-----------------------------------------------------'

C      CALL ADD_INFO("SODIAG_MUMAT2R",MUMAT2R,9*N*N,4)
C      CALL ADD_INFO("SODIAG_MUMAT2I",MUMAT2R,9*N*N,4)
C      CALL ADD_INFO("SODIAG_LMATR",LMATR,9*N*N,4)
C      CALL ADD_INFO("SODIAG_LMATI",LMATI,9*N*N,4)
C      CALL ADD_INFO("SODIAG_SMATR",SMATR,9*N*N,4)
C      CALL ADD_INFO("SODIAG_SMATI",SMATI,9*N*N,4)

      WRITE(6,*)
      WRITE(6,*) "***********************************************"
      WRITE(6,*) "* ENDING SODIAG *******************************"
      WRITE(6,*) "***********************************************"

      CALL GETMEM('DMATTMPA','FREE','REAL',LDMATTMP,3*(NBST*(NBST+1)))

      RETURN
      END





      SUBROUTINE SPIN_PHASE_RASSI(IPGLOB,DIPSO2,GMAIN,DIM,ZIN,ZOUT)
C
C     The RASSI program gives a random phase to the spin-orbit functions.
C
C     This routine performs a simple check with the obtained spin functions,
C     in order to determine the phase of the spin functions.
C     IF the phase is not the same, then the spin functions will be multiplied
C     with the correspondind coefficient that sets the same phase to all spin
C     eigenfunctions
C

      INTEGER  l,i,j,i1,i2,NPAR,ms1,ms2,DIM, IPGLOB

      REAL*8  GMAIN(3), ALFA(DIM)

      COMPLEX*16  PHS(3,DIM,DIM), ZIN(DIM,DIM), DIPSO2(3,DIM,DIM),
     & Spin2(3,DIM,DIM), PHSA(DIM,DIM), PHSA2(DIM,DIM),
     & ZOUT(DIM,DIM)

#include "spin.fh"

C Determine the Parity:
      NPAR=MOD(DIM,2)

C  Change the basis of the magnetic moment matrices from the RASSI functions to the
C  effective spin eigenfunctions-- eigenfunctions of the Mu_Z
C
      do L=1,3
       do I=1,DIM
        do J=1,DIM
      PHS(L,I,J)  =(0.d0,0.d0)
      SPIN2(L,I,J)=(0.d0,0.d0)
      PHSA(I,J)   =(0.d0,0.d0)
      PHSA2(I,J)  =(0.d0,0.d0)
        enddo
       enddo
      enddo

      do l=1,3
        do i=1,DIM
          do j=1,DIM
            do i1=1,DIM
              do i2=1,DIM
      PHS(l,i,j)=PHS(l,i,j)+DIPSO2(l,i1,i2)*DCONJG(ZIN(i1,i))*
     & ZIN(i2,j)
              enddo
            enddo
          enddo
        enddo
      enddo


CC Rewrite the Spin m.e. in a new basis:

      i=0
      do ms1=(DIM-NPAR)/2,-(DIM-NPAR)/2,-1
      if((ms1.eq.0).AND.(NPAR.EQ.0)) goto 18
      i=i+1
      j=0
      do ms2=(DIM-NPAR)/2,-(DIM-NPAR)/2,-1
      if((ms2.eq.0).AND.(NPAR.EQ.0)) goto 17
      j=j+1
      Spin2(1,i,j)=Spin(1,DIM,ms1,ms2)
      Spin2(2,i,j)=Spin(2,DIM,ms1,ms2)
      Spin2(3,i,j)=Spin(3,DIM,ms1,ms2)
17    continue
      enddo
18    continue
      enddo

      do i=1,DIM
         do j=1,DIM
      if(Spin2(1,i,j).EQ.(0.0d0,0.0d0)) goto 20
      PHSA2(i,j)=PHS(1,i,j)/(GMAIN(1)*Spin2(1,i,j))
      PHSA(i,j)=(0.0d0,-1.0d0)*LOG(PHSA2(i,j))
20    continue
         enddo
      enddo

      do I=1,DIM
      ALFA(I)=0.0d0
      enddo

      do I=2,DIM
      ALFA(I)=ALFA(I-1)+DBLE(PHSA(I-1,I))
      enddo

      do i=1,DIM
        do J=1,DIM
      ZOUT(J,I)=(0.0d0,0.0d0)
      ZOUT(J,I)=EXP((0.0d0,-1.0d0)*ALFA(i))*ZIN(j,i)
        enddo
      enddo

      if(IPGLOB.GT.2) then
      write(6,'(/)')
      write(6,'( 5x,a)')  'MAGNETIC MOMENT MATRIX ELEMENTS IN THE '//
     &'BASIS OF SPIN EIGENFUNCTIONS -- PHS(ic,i,j)'
      write(6,*)
      do l=1,3
      write(6,'(5x,a,i2)') 'PROJECTION=', l
      write(6,*)
        do i=1,DIM
        write(6,'(16(2F12.8,2x))') (PHS(l,i,j), j=1,DIM)
        enddo
      enddo
      write(6,'(5x,a)') 'Spin2(1,i,j)'
      do i=1,DIM
      write(6,'(3x,16(2F12.6,2x))') (Spin2(1,i,j),j=1,DIM)
      enddo
      write(6,'(5x,a)') 'Spin2(2,i,j)'
      do i=1,DIM
      write(6,'(3x,16(2F12.6,2x))') (Spin2(2,i,j),j=1,DIM)
      enddo
      write(6,'(5x,a)') 'PHSA2(i,j)'
      do i=1,DIM
      write(6,'(3x,16(2F12.6,2x))') (PHSA2(i,j),j=1,DIM)
      enddo
      write(6,'(5x,a)') 'PHSA(i,j)'
      do i=1,DIM
      write(6,'(3x,16(2F12.6,2x))') (PHSA(i,j),j=1,DIM)
      enddo
      write(6,'(5x,a)') 'ALFA'
      write(6,'(3x,16(2F12.6,2x))') (ALFA(j),j=1,DIM)
      write(6,'(5x,a)') 'ZOUT'
      do j=1,DIM
      write(6,'(3x,16(2F12.6,2x))') (ZOUT(i,j),i=1,DIM)
      enddo
      endif

      RETURN
      END




      SUBROUTINE ATENS_RASSI(moment, dim, gtens, maxes, IPGLOB)

      IMPLICIT NONE
      INTEGER dim,ic1,ic2,i,j,k,l,IPGLOB,info

      REAL*8  A_TENS_TERM(3,3),W(3),MAIN(3),Z(3,3),maxes(3,3),gtens(3)
      REAL*8 Det_gtens, diff12,diff23,ZR(3,3)
      COMPLEX*16 moment(3,dim,dim),AC_TENS(3,3),A_TEMP(3,3,dim,dim)
c------------------------------------------------
c  dim    -- size of the matrices
c  moment -- matrix of size (3,dim,dim) of the moment (magnetic, spin or angular)
c  gtens  -- array of size (3) keeping the main values of the A tensor ( dsqrt(main_values) )
c  maxes  -- array of size (3,3) keeping the main axes of the A tensor writen in
c            the right coordinate system (Determinant = +1)
c  IPGLOB -- the print level of the subroutine
c----------------------------------------------
C
C   initialization:
C

      do I=1,3
       do J=1,3
      AC_TENS(I,J)=(0.d0,0.d0)
      A_TENS_TERM(I,J)=0.d0
        do K=1,dim
         do L=1,dim
      A_TEMP(I,J,K,L)=(0.d0,0.d0)
         enddo
        enddo
       enddo
      enddo

      do ic1=1,3
         do ic2=1,3
            do i=1,dim
               do j=1,dim
                  do k=1,dim
      A_temp(ic1,ic2,i,j) = A_temp(ic1,ic2,i,j)+moment(ic1,i,k)*
     & moment(ic2,k,j)
                  enddo
               enddo
            enddo
         enddo
      enddo

      IF(IPGLOB.GE.3) THEN
      write(6,'(/)')
      write(6,'(5X,A)') 'BPMOMENT(ic1,ic2):'
      write(6,*)
      do ic1=1,3
      do i=1,dim
      do j=1,dim
      write(6,*) moment(ic1,i,j)
      end do
      end do
      end do

      write(6,'(/)')
      write(6,'(5X,A)') 'A_TEMP(ic1,ic2):'
      write(6,*)
      do ic1=1,3
      do ic2=1,3
      do i=1,dim
      do j=1,dim
      write(6,*) A_temp(ic1,ic2,i,j)
      end do
      end do
      end do
      end do
      END IF

      do ic1=1,3
         do ic2=1,3
            do k=1,dim
      Ac_tens(ic1,ic2)=Ac_tens(ic1,ic2)+A_temp(ic1,ic2,k,k)
            enddo
         enddo
      enddo
      do ic1=1,3
       do ic2=1,3
      A_TENS_TERM(ic1,ic2)=DBLE(Ac_tens(ic1,ic2)+Ac_tens(ic2,ic1))/2.d0
         enddo
      enddo
      do ic1=1,3
       do ic2=1,3
      A_TENS_TERM(ic1,ic2)=12.d0*A_TENS_TERM(ic1,ic2)/DBLE(dim**3-dim)
         enddo
      enddo

c      if(IPGLOB.GT.2) then
      write(6,'(/)')
      write(6,'(5X,A)') 'A_TENS_TERM(ic1,ic2):'
      write(6,*)
      do ic1=1,3
c      write(6,'(5X,3(2F14.7,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
      write(6,'(5X,3(2F21.14,3x))') (A_TENS_TERM(ic1,ic2),ic2=1,3)
      enddo
c      endif
C
C   Diagonalization of A_tens - g tensors
C
      do I=1,3
      main(I)=0.0D0
      w(I)=0.0D0
       do J=1,3
      z(I,J)=0.0D0
       enddo
      enddo
      info=0

      call DIAG_R2_RASSI(A_TENS_TERM,3,info,w,z)
      if(INFO.NE.0) goto 199
      if((w(1).LT.0.D0).AND.(w(2).LT.0.D0).AND.(w(3).LT.0.D0)) then
      write(6,'(2x,A)') 'ALL EIGENVALUES OF THE A-TENSOR ARE NEGATIVE'
      write(6,'(2X,A)') 'THIS IS A VERY UNUSUAL SITUATION. PLEASE'//
     & 'CHECK MANUALLY '
      write(6,'(2x,A)') 'THE FOLLOWING PART OF THE PSEUDOSPIN SECTION'
      write(6,'(2x,A)') 'MUST BE DISREGARDED. THE RESULTS ARE NOT' //
     & 'TRUSTABLE.'
      goto 199
      endif
c

      IF(IPGLOB.GE.3) THEN
      write(6,*)
      write(6,'(4x,A)') 'A_TENS_TERM TENSOR:'
      write(6,'(65a)') ('-',i=1,56),'|'
      write(6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|',
     & 'MAIN MAGNETIC AXES','|', 'x , y , z  -- initial Cartesian axes'
      write(6,'(57a,3x,a)') ('-',i=1,19),'|',('-',i=1,36),'|',
     & 'Xm, Ym, Zm -- main magnetic axes'
      write(6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z',
     & '|'
      write(6,'(65a)') ('-',i=1,19),'|',('-',i=1,4),'|',('-',i=1,31),
     & '|'
      write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',w(1),' | Xm |',
     & (z(j,1),j=1,3),'|'
      write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',w(2),' | Ym |',
     & (z(j,2),j=1,3),'|'
      write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',w(3),' | Zm |',
     & (z(j,3),j=1,3),'|'
      write(6,'(65a)') ('-',i=1,56),'|'
      END IF

      do I=1,3
      if(W(I).LT.0.d0) then
      W(I)=0.1d-14
      endif
      MAIN(i)=sqrt(W(i))
      enddo

      if(IPGLOB.GT.2) write(6,'(5x,a,3F9.5)') 'EIGenValues after DSPEV:'
     & , (W(I),I=1,3)

C  Check the sign of the coordinate system. if CS is Left-handed,
C  then change it to RIGHT-handed
      Det_gtens=0.d0
      do I=1,3
       do J=1,3
      ZR(I,J)=0.0D0
      ZR(I,J)=Z(I,J)
       enddo
      enddo
      Det_gtens=ZR(1,1)*(ZR(2,2)*ZR(3,3)-ZR(2,3)*ZR(3,2))
     &         -ZR(1,2)*(ZR(2,1)*ZR(3,3)-ZR(2,3)*ZR(3,1))
     &         +ZR(1,3)*(ZR(2,1)*ZR(3,2)-ZR(2,2)*ZR(3,1))
      if(Det_gtens.LT.0.0d0) then
      do i=1,3
      Z(i,1)=-Z(i,1)
      enddo
      if(IPGLOB.GT.2) write(6,'(a)') 'The original coordinate system '//
     & 'was LEFT-handed. It has been changed to the RIGHT-handed'
      endif
      diff12=0.d0
      diff23=0.d0
      diff12=MAIN(2)-MAIN(1)
      diff23=MAIN(3)-MAIN(2)
      if(IPGLOB.GT.2) then
      write(6,'(5x,a,3F19.15)') 'diff12 = ', diff12
      write(6,'(5x,a,3F19.15)') 'diff23 = ', diff23
      endif

      do i=1,3
      gtens(i)=0.0d0
         do j=1,3
         maxes(i,j)=0.0d0
         enddo
      enddo
c set the main Z axis:
      if(diff12.gt.diff23) then
      gtens(3)=MAIN(1)
      gtens(2)=MAIN(2)
      gtens(1)=MAIN(3)
         if(Z(3,1).ge.0.0d0) then
           do i=1,3
           maxes(i,3)=Z(i,1)
           maxes(i,1)=Z(i,3)
           enddo
         elseif(Z(3,1).lt.0.0d0) then
           do i=1,3
           maxes(i,3)=-1.d0*Z(i,1)
           maxes(i,1)=Z(i,3)
           enddo
         endif
           maxes(1,2)=maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
           maxes(2,2)=maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
           maxes(3,2)=maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)

      elseif(diff23.gt.diff12) then
      gtens(3)=MAIN(3)
      gtens(2)=MAIN(2)
      gtens(1)=MAIN(1)
         if(Z(3,3).ge.0.0d0) then
           do i=1,3
           maxes(i,3)=Z(i,3)
           maxes(i,1)=Z(i,1)
           enddo
         elseif(Z(3,3).lt.0.0d0) then
           do i=1,3
           maxes(i,3)=-1.d0*Z(i,3)
           maxes(i,1)=Z(i,1)
           enddo
         endif
           maxes(1,2)=maxes(2,3)*maxes(3,1)-maxes(2,1)*maxes(3,3)
           maxes(2,2)=maxes(1,1)*maxes(3,3)-maxes(1,3)*maxes(3,1)
           maxes(3,2)=maxes(1,3)*maxes(2,1)-maxes(1,1)*maxes(2,3)
      endif

      if(IPGLOB.GT.2) then
      write(6,*)
      write(6,'(20X,A)') 'A-TENSOR:'
      write(6,*)
      write(6,'(10X,A,10X,3(F11.5,2X))') '|  xx    xy    xz  |',
     & (A_TENS_TERM(1,ic2),ic2=1,3)
      write(6,'(10X,A,10X,3(F11.5,2X))') '|  yx    yy    yz  |',
     & (A_TENS_TERM(2,ic2),ic2=1,3)
      write(6,'(10X,A,10X,3(F11.5,2X))') '|  zx    zy    zz  |',
     & (A_TENS_TERM(3,ic2),ic2=1,3)
      endif

      if(IPGLOB.GT.2) then
      write(6,*)
      write(6,'(4x,A)') 'g TENSOR:'
      write(6,'(65a)') ('-',i=1,56),'|'
      write(6,'(4x,A,4x,A,13x,A,5x,a,3x,a)') 'MAIN VALUES','|',
     & 'MAIN MAGNETIC AXES','|', 'x , y , z  -- initial Cartesian axes'
      write(6,'(57a,3x,a)') ('-',i=1,19),'|',('-',i=1,36),'|',
     & 'Xm, Ym, Zm -- main magnetic axes'
      write(6,'(19x,a,4x,a,5x,a,9x,a,9x,a,5x,a)') '|','|','x','y','z',
     & '|'
      write(6,'(65a)') ('-',i=1,19),'|',('-',i=1,4),'|',('-',i=1,31),
     & '|'
      write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gX = ',gtens(1),' | Xm |',
     & (maxes(j,1),j=1,3),'|'
      write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),' | Ym |',
     & (maxes(j,2),j=1,3),'|'
      write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),' | Zm |',
     & (maxes(j,3),j=1,3),'|'
      write(6,'(65a)') ('-',i=1,56),'|'
C      Call Add_Info('GTENS_MAIN',gtens,3,5)
      endif

 199  continue

      RETURN
      END



      Subroutine DIAG_R2_RASSI(MATRIX,NBTOT,INFO,W1,Z1)
C
C   THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A REAL SQUARE
C   MATRIX WITH THE DIMENSION NBTOT. THE EIGENVALUES OF THE DIAGONALIZATION
C   ARE DIRECTED INTO W1 AND THE REAL EIGENVECTORS ARE WRITTEN TO Z1.
C

      IMPLICIT NONE
      INTEGER   INFO,I,J,NBTOT
      REAL*8 AP(Nbtot*(Nbtot+1)/2), WORK(3*Nbtot),W1(NBTOT),
     & W(Nbtot),Z(Nbtot,Nbtot),Z1(NBTOT,NBTOT),MATRIX(NBTOT,NBTOT)

C initializations
      INFO=0
      Do j=1,Nbtot
        Do i=1,j
       AP(i+(j-1)*j/2)=0.0D0
        enddo
      enddo
      do i=1,3*Nbtot
      work(i)=0.0d0
      enddo
      Do I=1,Nbtot
      W1(I)=0.0D0
      W(i)=0.d0
       Do J=1,Nbtot
      Z1(J,I)=0.0D0
      Z(j,i)=0.0D0
       enddo
      enddo
      Do j=1,Nbtot
        Do i=1,j
       AP(i+(j-1)*j/2)=MATRIX(i,j)
        enddo
       enddo

      call dspev_('V','U',Nbtot,AP,W,Z,Nbtot,WORK,INFO)
      Do I=1,Nbtot
      W1(I)=W(I)
       Do J=1,Nbtot
      Z1(J,I)=Z(J,I)
       enddo
      enddo
      Return
      End
