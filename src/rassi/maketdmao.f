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
* Copyright (C) 2021, Rulin Feng                                       *
************************************************************************
*       ***************************************************
*                Get TDM in AO basis from SO states
*       ****************************************************
*        This routine is modified from sonatorbm.f to give the
*        full transition density matrix (TDM) in atomic orbital
*        basis, which is not a nbtri sized matrix as in sonatorbm
*        but a full nbst**2 sized matrix.
*

      SUBROUTINE MAKETDMAO(CHARTYPE,
     &                     USOR,USOI,ASS,BSS,NSS,
     &                     iOpt,ROTMAT,DENSOUT)
      use rassi_aux, only : idisk_TDM
      use rassi_global_arrays, only: JBNUM
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='MAKETDMAO')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION USOR(NSS,NSS),USOI(NSS,NSS)
      DIMENSION DENSOUT(6,nbst**2)
      Dimension IZMR(3),IZMI(3)
      Dimension IZMR2(3),IZMI2(3)
      DIMENSION IOFF(8)
      CHARACTER*8 CHARTYPE!,LABEL
      Dimension ROTMAT(3,3)
      INTEGER ASS,BSS
c      INTEGER ASF,BSF
c      INTEGER, DIMENSION(1)::SIZ
c      INTEGER jm,count_i
      nbsts=nbst**2
c      count_i = 0

c VV: dummy initialization
      CGY=-1
      CGX=-1
      CG0=-1
c Initialize
      Do I=1,6
        CALL DCOPY_(nbsts,[0.0D00],0,DENSOUT(I,:),1)
      Enddo
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

C The following creates an array that is used to
C map a specific spin state to the corresponding
C spin-free state and to its spin
C (see prprop.f and others)

      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)

      ISS=0
      DO ISF=1,NSTATE
        JOB=JBNUM(ISF)
        MPLET=MLTPLT(JOB)

        DO MSPROJ=-MPLET+1,MPLET-1,2
          ISS=ISS+1
          IWORK(LMAPST-1+ISS)=ISF
          IWORK(LMAPSP-1+ISS)=MPLET
          IWORK(LMAPMS-1+ISS)=MSPROJ
        END DO
      END DO


c Allocate some arrays
c LSDMXR, etc      DM/TDM for this iteration
c LSDMXR2, etc     Accumulated DM/TDM
c LTMPR,I          Temporary array for U*AU multiplication
c LTDMZZ           DM/TDM as read from file
c LSCR             Scratch for expansion of LTDMZZ
      CALL GETMEM('TSDMXR','ALLO','REAL',LSDMXR,nbsts)
      CALL GETMEM('TSDMYR','ALLO','REAL',LSDMYR,nbsts)
      CALL GETMEM('TSDMZR','ALLO','REAL',LSDMZR,nbsts)
      CALL GETMEM('TSDMXI','ALLO','REAL',LSDMXI,nbsts)
      CALL GETMEM('TSDMYI','ALLO','REAL',LSDMYI,nbsts)
      CALL GETMEM('TSDMZI','ALLO','REAL',LSDMZI,nbsts)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMXR),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMYR),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMZR),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMXI),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMYI),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMZI),1)
      IZMR(1)=LSDMXR
      IZMR(2)=LSDMYR
      IZMR(3)=LSDMZR
      IZMI(1)=LSDMXI
      IZMI(2)=LSDMYI
      IZMI(3)=LSDMZI

      CALL GETMEM('TSDMXR2','ALLO','REAL',LSDMXR2,nbsts)
      CALL GETMEM('TSDMYR2','ALLO','REAL',LSDMYR2,nbsts)
      CALL GETMEM('TSDMZR2','ALLO','REAL',LSDMZR2,nbsts)
      CALL GETMEM('TSDMXI2','ALLO','REAL',LSDMXI2,nbsts)
      CALL GETMEM('TSDMYI2','ALLO','REAL',LSDMYI2,nbsts)
      CALL GETMEM('TSDMZI2','ALLO','REAL',LSDMZI2,nbsts)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMXR2),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMYR2),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMZR2),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMXI2),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMYI2),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMZI2),1)
      IZMR2(1)=LSDMXR2
      IZMR2(2)=LSDMYR2
      IZMR2(3)=LSDMZR2
      IZMI2(1)=LSDMXI2
      IZMI2(2)=LSDMYI2
      IZMI2(3)=LSDMZI2

      CALL GETMEM('TSDMTMPR','ALLO','REAL',LTMPR,nbsts)
      CALL GETMEM('TSDMTMPI','ALLO','REAL',LTMPI,nbsts)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LTMPR),1)
      CALL DCOPY_(nbsts,[0.0D00],0,WORK(LTMPI),1)

      CALL GETMEM('TDMSCR','ALLO','REAL',LSCR,nbsts)
c zeroed inside the loop
c      CALL DCOPY_(nbsts,0.0D00,0,WORK(LSCR),1)

      CALL GETMEM('TDMZZ','ALLO','REAL',LTDMZZ,NTDMZZ)
      CALL DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZ),1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MAIN LOOP OVER KSF/LSF
C WRITTEN AS IN PRPROP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CORRESPONDING SPIN-FREE STATES OF THE
C REQUESTED SPIN STATES
c      ASF=IWORK(LMAPST-1+ASS)
c      BSF=IWORK(LMAPST-1+BSS)

      DO KSS=1,NSS
       KSF=IWORK(LMAPST-1+KSS)
       MPLETK=IWORK(LMAPSP-1+KSS)
       MSPROJK=IWORK(LMAPMS-1+KSS)

       DO LSS=1,NSS
        LSF=IWORK(LMAPST-1+LSS)
        MPLETL=IWORK(LMAPSP-1+LSS)
        MSPROJL=IWORK(LMAPMS-1+LSS)

        JOB1=JBNUM(KSF)
        JOB2=JBNUM(LSF)
        LSYM1=IRREP(JOB1)
        LSYM2=IRREP(JOB2)
        ISY12=MUL(LSYM1,LSYM2)

C SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS
        Call mk_IOFF(IOFF,nSYM,NBASF,ISY12)

c These are going to be zero, so head them off at the pass
        IF(ITYPE.LE.2
     &     .AND.(MPLETK.NE.MPLETL.OR.MSPROJK.NE.MSPROJL)) GOTO 2200


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C Transition density matrices, TDMZZ, in AO basis.
C WDMZZ similar, but WE-reduced 'triplet' densities.
C TDMZZ will store either, depending on the type
        CALL DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZ),1)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C IDTDM: TOC array for transition 1-matrices
c TDMZZ is stored on disk from i = 1, NSTATE j=1, i
c so swap if needed
        iEmpty=iDisk_TDM(KSF,LSF,2)
        IDISK=iDisk_TDM(KSF,LSF,1)
        iOpt=2
        IF (ITYPE.GE.3) Then
           iGo=4
           CALL dens2file(Work(LTDMZZ),Work(LTDMZZ),Work(LTDMZZ),
     &                    nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
C NOTE-the TD matrix as read in has an incorrect sign
           CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)
        Else
           iGo=1
           CALL dens2file(Work(LTDMZZ),Work(LTDMZZ),Work(LTDMZZ),
     &                    nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
        END IF


c Anti-hermitian properties need a little fixing
        IF((ITYPE.EQ.2.OR.ITYPE.EQ.4).AND.(KSF.LE.LSF))
     &          CALL DSCAL_(NTDMZZ,-1.0d0,WORK(LTDMZZ),1)


C CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
C AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSCR),1)
cccccccccccccc
c nbsq and NTDMZZ are the same, but they are the
c sizes of the symmetry-adapted matrices
c what we need is a NBST**2-sized matrix without symmetry
c Thus we expand the final TDM in C1 symmetry
c Leave the zero matrix elements as they are
ccccccccccccccc

cc Code seems to work for both diagonal symmetry blocks and NON-DIAGONAL symmetry blocks
cc This is strange, but WORK(TDMZZ) stores the matrix elements transposed in its
cc symmetry blocks, e.g., see mk_twdm.f
cc I don't know what's the reason, assuming it must be correct,
cc the following scheme is to change the orders of lower-triangular
cc elements of WORK(TDMZZ), putting them in WORK(TDMZZSCR) for later expansion into
cc C1 matrix.
c        Call GETMEM('TDMZZSCR','ALLO','REAL',LTDMZZSCR,NTDMZZ)
c        Call DCOPY_(NTDMZZ,[0.0D00],0,WORK(LTDMZZSCR),1)
c        ITD=0
c        IOF=0
c        Do ISY1=1,NSYM
c          NB1=NBASF(ISY1)
c          Do ISY2=1,NSYM
c            NB2=NBASF(ISY2)
c            ISY12_ma=MUL(ISY1,ISY2)
c            IF(ISY12_ma.EQ.ISY12) then
c              IF(ISY1.GT.ISY2) THEN
c                Do J=1,NB2
c                  Do I=1,NB1
c                    ITD=IOF+I+NB1*(J-1)
c                    IJ=IOF+J+NB2*(I-1)
c                    TDM=WORK(LTDMZZ-1+IJ)
c                    WORK(LTDMZZSCR-1+ITD)=TDM
c                  Enddo
c                Enddo
c              ELSE
c                Do J=1,NB2
c                  Do I=1,NB1
c                    ITD=IOF+J+NB2*(I-1)
c                    IJ=IOF+I+NB1*(J-1)
c                    TDM=WORK(LTDMZZ-1+IJ)
c                    WORK(LTDMZZSCR-1+ITD)=TDM
c                  Enddo
c                Enddo
c              Endif
c              IOF=IOF+NB1*NB2
c            Endif
c          Enddo
c        Enddo
cc Expand into C1
        ITD=0
        NB1_i=0
        NB1_f=0
        Do ISY1=1,NSYM
          NB2_i=0
          NB2_f=0
          NB1=NBASF(ISY1)
          NB1_f=NB1_i+NB1
          Do ISY2=1,NSYM
            ISY12_ma=MUL(ISY1,ISY2)
            NB2=NBASF(ISY2)
            NB2_f=NB2_i+NB2
            If(ISY12_ma.EQ.ISY12) then
              Do J=NB2_i+1,NB2_f
                Do I=NB1_i+1,NB1_f
                  ITD=ITD+1
c                  TDM=WORK(LTDMZZSCR-1+ITD)
                  TDM=WORK(LTDMZZ-1+ITD)
                  IJ=I+NBST*(J-1)
                  WORK(LSCR-1+IJ)=WORK(LSCR-1+IJ)+TDM
                Enddo
              Enddo
            Endif
            NB2_i=NB2_i+NB2
          Enddo
          NB1_i=NB1_i+NB1
        Enddo
c        Call GETMEM('TDMZZSCR','FREE','REAL',LTDMZZSCR,NTDMZZ)
***********************************
c        write(6,*) 'In maketdmao'
c        LABEL(1:8)='MLTPL  1'
c        IRC=-1
c        ICMP=3
c        ISYLAB=1
c        Call GETMEM('MSq','ALLO','REAL',LDIPs,nbsts)
c        Call GETMEM('BUFF1','ALLO','REAL',LBUFF1,nbsts)
c        IOPT=ibset(0,sOpSiz)
c        Call IRDONE(IRC,IOPT,LABEL,ICMP,SIZ,ISYLAB)
c        write(6,*) 'SIZ', SIZ
c        write(6,*) 'ISYLAB', ISYLAB
c        Call GETMEM('MLTPL  1','ALLO','REAL',LDIP,SIZ)
c        IOPT=ibset(ibset(0,sNoOri),sNoNuc)
c        Call RDONE(IRC,IOPT,LABEL,ICMP,WORK(LDIP),ISYLAB)
c        Do I=0,NTDMZZ-1
c          write(6,*) 'TDMZZ', WORK(LTDMZZ+I)
c        Enddo
c          write(6,*) '*********************'
c        Do I=0,SIZ(1)-1
c          write(6,*) 'DIP', WORK(LDIP+i)
c        Enddo
c        Call DESYM_SONTO(WORK(LDIP),SIZ,WORK(LDIPs),ISYLAB)
c        Call print_matrix('LSCR ',nbst,nbsts,WORK(LSCR))
c        Call print_matrix('LDIPs ',nbst,nbsts,WORK(LDIPs))
c        call DGEMM_('N','T',nbst,nbst,nbst,1.0D0,WORK(LDIPs),nbst,
c     &              WORK(LSCR),nbst,0.0D0,WORK(LBUFF1),nbst)
c        Call print_matrix('LBUFF ',nbst,nbsts,WORK(LBUFF1))
cc Trace the BUFF1 matrix
c        Transition_Dipole=0.0D0
c        do i=0,nbst-1
c            Transition_Dipole=Transition_Dipole+WORK(LBUFF1+i*nbst+i)
c        enddo
c        write(6,*) 'Transition_Dipole in maketdmao', Transition_Dipole
c        write(6,*)
*************************************

c ie, see how AMFI is processed in soeig.f
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMXR),1)
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMYR),1)
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMZR),1)
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMXI),1)
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMYI),1)
        CALL DCOPY_(nbsts,[0.0D00],0,WORK(LSDMZI),1)

        IF(ITYPE.GE.3) THEN
          S1=0.5D0*DBLE(MPLETK-1)
          SM1=0.5D0*DBLE(MSPROJK)
          S2=0.5D0*DBLE(MPLETL-1)
          SM2=0.5D0*DBLE(MSPROJL)
          FACT=1.0D0/SQRT(DBLE(MPLETK))
          IF(MPLETK.EQ.MPLETL-2) FACT=-FACT

          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX=SQRT(0.5D0)*(CGM-CGP)
          CGY=SQRT(0.5D0)*(CGM+CGP)
        END IF

        IF((ITYPE.EQ.1.OR.ITYPE.EQ.2)
     &          .AND.MPLETK.EQ.MPLETL
     &          .AND.MSPROJK.EQ.MSPROJL) THEN
          CALL DAXPY_(nbsts,1.0d0,WORK(LSCR),1,WORK(LSDMZR),1)
        ELSE IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
          If (iOpt.eq.1) Then
          CALL DAXPY_(nbsts,CGX*ROTMAT(1,1),WORK(LSCR),1,WORK(LSDMXR),1)
          CALL DAXPY_(nbsts,CGY*ROTMAT(2,1),WORK(LSCR),1,WORK(LSDMXI),1)
          CALL DAXPY_(nbsts,CG0*ROTMAT(3,1),WORK(LSCR),1,WORK(LSDMXR),1)

          CALL DAXPY_(nbsts,CGX*ROTMAT(1,2),WORK(LSCR),1,WORK(LSDMYR),1)
          CALL DAXPY_(nbsts,CGY*ROTMAT(2,2),WORK(LSCR),1,WORK(LSDMYI),1)
          CALL DAXPY_(nbsts,CG0*ROTMAT(3,2),WORK(LSCR),1,WORK(LSDMYR),1)

          CALL DAXPY_(nbsts,CGX*ROTMAT(1,3),WORK(LSCR),1,WORK(LSDMZR),1)
          CALL DAXPY_(nbsts,CGY*ROTMAT(2,3),WORK(LSCR),1,WORK(LSDMZI),1)
          CALL DAXPY_(nbsts,CG0*ROTMAT(3,3),WORK(LSCR),1,WORK(LSDMZR),1)
          Else
          CALL DAXPY_(nbsts,CGX,WORK(LSCR),1,WORK(LSDMXR),1)
          CALL DAXPY_(nbsts,CGY,WORK(LSCR),1,WORK(LSDMYI),1)
          CALL DAXPY_(nbsts,CG0,WORK(LSCR),1,WORK(LSDMZR),1)
          End If
        END IF

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SPINORBIT
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Sign of the left-hand imaginary part is handled
c when doing DAXPY
        URR=USOR(LSS,BSS)
        UIR=USOI(LSS,BSS)
        URL=USOR(KSS,ASS)
        UIL=USOI(KSS,ASS)

        DO IDIR=1,3
          CALL DCOPY_(nbsts,[0.0D00],0,WORK(LTMPR),1)
          CALL DCOPY_(nbsts,[0.0D00],0,WORK(LTMPI),1)

C right side
          CALL DAXPY_(nbsts,       URR,WORK(IZMR(IDIR)),1,WORK(LTMPR),1)
          CALL DAXPY_(nbsts,-1.0d0*UIR,WORK(IZMI(IDIR)),1,WORK(LTMPR),1)
          CALL DAXPY_(nbsts,       UIR,WORK(IZMR(IDIR)),1,WORK(LTMPI),1)
          CALL DAXPY_(nbsts,       URR,WORK(IZMI(IDIR)),1,WORK(LTMPI),1)

C left side
         CALL DAXPY_(nbsts,       URL,WORK(LTMPR),1,WORK(IZMR2(IDIR)),1)
         CALL DAXPY_(nbsts,       UIL,WORK(LTMPI),1,WORK(IZMR2(IDIR)),1)
         CALL DAXPY_(nbsts,       URL,WORK(LTMPI),1,WORK(IZMI2(IDIR)),1)
         CALL DAXPY_(nbsts,-1.0d0*UIL,WORK(LTMPR),1,WORK(IZMI2(IDIR)),1)
        END DO
cccccccccccccccccccccc
c END SPINORBIT STUFF
cccccccccccccccccccccc

C END MAIN LOOP OVER STATES (KSS,LSS)
 2200   CONTINUE

       END DO
      END DO


C Store this density to DENSOUT
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
        DO I=1,nbsts
          DENSOUT(1,I)=WORK(LSDMXR2-1+I)
          DENSOUT(2,I)=WORK(LSDMYR2-1+I)
          DENSOUT(3,I)=WORK(LSDMZR2-1+I)
          DENSOUT(4,I)=WORK(LSDMXI2-1+I)
          DENSOUT(5,I)=WORK(LSDMYI2-1+I)
          DENSOUT(6,I)=WORK(LSDMZI2-1+I)
        END DO
      ELSE
        DO I=1,nbsts
          DENSOUT(1,I)=WORK(LSDMZR2-1+I)
          DENSOUT(2,I)=WORK(LSDMZR2-1+I)
          DENSOUT(3,I)=WORK(LSDMZR2-1+I)
          DENSOUT(4,I)=WORK(LSDMZI2-1+I)
          DENSOUT(5,I)=WORK(LSDMZI2-1+I)
          DENSOUT(6,I)=WORK(LSDMZI2-1+I)
        END DO
      END IF

c Free memory
      CALL GETMEM('TDMSCR','FREE','REAL',LSCR,nbsts)
      CALL GETMEM('TDMZZ','FREE','REAL',LTDMZZ,NTDMZZ)

      CALL GETMEM('TSDMXR','FREE','REAL',LSDMXR,nbsts)
      CALL GETMEM('TSDMYR','FREE','REAL',LSDMYR,nbsts)
      CALL GETMEM('TSDMZR','FREE','REAL',LSDMZR,nbsts)
      CALL GETMEM('TSDMXI','FREE','REAL',LSDMXI,nbsts)
      CALL GETMEM('TSDMYI','FREE','REAL',LSDMYI,nbsts)
      CALL GETMEM('TSDMZI','FREE','REAL',LSDMZI,nbsts)

      CALL GETMEM('TSDMTMPR','FREE','REAL',LTMPR,nbsts)
      CALL GETMEM('TSDMTMPI','FREE','REAL',LTMPI,nbsts)

      CALL GETMEM('TSDMXR2','FREE','REAL',LSDMXR2,nbsts)
      CALL GETMEM('TSDMYR2','FREE','REAL',LSDMYR2,nbsts)
      CALL GETMEM('TSDMZR2','FREE','REAL',LSDMZR2,nbsts)
      CALL GETMEM('TSDMXI2','FREE','REAL',LSDMXI2,nbsts)
      CALL GETMEM('TSDMYI2','FREE','REAL',LSDMYI2,nbsts)
      CALL GETMEM('TSDMZI2','FREE','REAL',LSDMZI2,nbsts)

      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

      RETURN
      END
