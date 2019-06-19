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
      Subroutine ONCMO(CMO1,CMO2)
      use general_data, only :
     &    nSym, nBAs, nDel, nActEl, nDelt, nSSH, nDel, nOrb
      use rasscf_data, only :
     &    nSec, nTOT3, nTOT4, Tot_Nuc_Charge, nFr, nIn, nOrbT
      implicit none
      real*8, intent(in) :: CMO1(*)
      real*8, intent(out) :: CMO2(*)
#include "warnings.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
      integer :: DDOT_
      integer :: iPRLEV, nBM, NOM, NSBUF, iSYM, nB, nO, LSBUF, LSMAT,
     &    LSCTMP, LOVL, i_Rc, i_Opt, i_Component, i_SymLbl,
     &    xMol_Charge, nDSAVe, NNEGSS, ip_SBUF, ip_CMO, NNEW,
     &    IOLD, IPOLD, IPNEW, NREMOV, ND, NS, NDNEW, NSNEW
      real*8 :: XNRM2, XSCL
      Parameter (ROUTINE='ONCMO   ')
      Call qEnter('ONCMO')

C Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      NBM=0
      NOM=0
      NSBUF=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       NO=NB-NDEL(ISYM)
       NBM=MAX(NBM,NB)
       NOM=MAX(NOM,NO)
       NSBUF=NSBUF+(NB*(NB+1))/2
      END DO
      IF(NOM.eq.0) GOTO 900
      CALL GETMEM('SBUF','ALLO','REAL',LSBUF,NSBUF+4)
      CALL GETMEM('SMAT','ALLO','REAL',LSMAT,NBM**2)
      CALL GETMEM('SCTMP','ALLO','REAL',LSCTMP,NBM)
      CALL GETMEM('OVL','ALLO','REAL',LOVL,NBM)
* Read overlap matrix SMAT:
      i_Rc=0
      i_Opt=2
      i_Component=1
      i_SymLbl=1
      Call RdOne(i_Rc,i_Opt,'Mltpl  0',i_Component,WORK(LSBUF),i_SymLbl)
      If ( i_Rc.ne.0 ) Then
        Write(LF,*)' RASSCF is trying to orthonormalize orbitals but'
        Write(LF,*)' could not read overlaps from ONEINT. Something'
        Write(LF,*)' is wrong with the file, or possibly with the'
        Write(LF,*)' program. Please check.'
        Call Quit(_RC_IO_ERROR_READ_)
      End If
*
*---- Print out nuclear charge
*
      Tot_Nuc_Charge=Work(LSBUF+NSBUF+3)
      xMol_Charge=Tot_Nuc_Charge-DBLE(2*(NFR+NIN)+NACTEL)
      Call put_dscalar('Total Charge    ',xMol_Charge)
      IF(IPRLEV.GE.USUAL) THEN
        Write(LF,*)
        Write(LF,'(6x,A,f8.2)') 'Total molecular charge',xMol_Charge
      End If
*
* Orthonormalize symmetry blocks:
*
      NDSAVE=NDELT
      NNEGSS=0
      ip_SBUF=LSBUF
      ip_CMO=1
      Do  iSym=1,nSym
       NB=nBas(iSym)
       NO=NB-nDel(iSym)
       If ( NB.gt.0 ) then
         IF (NO.GT.0) THEN
           Call SQUARE(WORK(ip_SBUF),WORK(LSMAT),1,NB,NB)
* NNEW=Nr of already orthonormal new CMO''s
           NNEW=0
           DO IOLD=1,NO
            IPOLD=ip_CMO+NB*(IOLD-1)
            IPNEW=ip_CMO+NB*NNEW
           IF(IPNEW.LT.IPOLD)CALL DCOPY_(NB,CMO1(IPOLD),1,CMO2(IPNEW),1)
  10        CONTINUE
            CALL DGEMM_('N','N',NB,1,NB,1.0D0,WORK(LSMAT),NB,
     &                 CMO2(IPNEW),NB,0.0D0,WORK(LSCTMP),NB)
            IF(NNEW.GT.0) THEN
              CALL DGEMM_('T','N',NNEW,1,NB,1.0D0,CMO2(ip_CMO),NB,
     &                  WORK(LSCTMP),NB,0.0D0,WORK(LOVL),NNEW)
              CALL DGEMM_('N','N',NB,1,NNEW,-1.0D0,CMO2(ip_CMO),NB,
     &                  WORK(LOVL),NNEW,1.0D0,CMO2(IPNEW),NB)
            END IF
            XNRM2=DDOT_(NB,WORK(LSCTMP),1,CMO2(IPNEW),1)
            XSCL=1.0D0/SQRT(XNRM2)
            IF (XNRM2.GT.1.0D-10) THEN
             CALL DSCAL_(NB,XSCL,CMO2(IPNEW),1)
             IF(XNRM2.LT.0.2D0) GOTO 10
             NNEW=NNEW+1
            END IF
           END DO
           NREMOV=NO-NNEW
           IF (NREMOV.GT.0) THEN
             ND=NDEL(ISYM)
             NS=NSSH(ISYM)
             NDNEW=NB-NNEW
             NSNEW=NS+ND-NDNEW
             IF(NSNEW.GE.0) THEN
              IF(IPRLEV.GE.TERSE) THEN
               Call WarningMessage(1,'ONCMO Warning')
               Write(LF,*)' * Exact or very near linear dependence '
               Write(LF,*)' * forces RASSCF to delete additional '//
     &                     'orbitals.'
               Write(LF,*)' *                  Symmetry block:',ISYM
               Write(LF,*)' * Earlier number of deleted orbs =',ND
               Write(LF,*)' *     New number of deleted orbs =',NDNEW
              END IF
             ELSE
              Write(LF,*)' **** ONCMO Error **************************'
              Write(LF,*)' * Exact or very near linear dependence '
              Write(LF,*)' * forces RASSCF to stop execution.'
              Write(LF,*)' *                  Symmetry block:',ISYM
              Write(LF,*)' * Effective nr of orthonormal orbs =',NNEW
              Write(LF,*)' *   Earlier number of deleted orbs =',ND
              Write(LF,*)' * Earlier number of secondary orbs =',NS
              Write(LF,*)' *       New number of deleted orbs =',NDNEW
              Write(LF,*)' *     New number of secondary orbs =',NSNEW
              NNEGSS=NNEGSS+1
             END IF
             NDEL(ISYM)=NDNEW
             NSSH(ISYM)=NSNEW
             NORB(ISYM)=NORB(ISYM)-NREMOV
             NDELT=NDELT+NREMOV
             NSEC =NSEC -NREMOV
             NORBT=NORBT-NREMOV
           END IF
         END IF
         ip_SBUF=ip_SBUF+(NB*NB+NB)/2
         ip_CMO=ip_CMO+NB*NB
       End If
      End Do
      IF(NNEGSS.GT.0) CALL QUIT(_RC_GENERAL_ERROR_)

      IF (NDSAVE.NE.NDELT) THEN
       NTOT3=0
       NTOT4=0
       DO ISYM=1,NSYM
          NTOT3=NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
          NTOT4=NTOT4+NORB(ISYM)**2
       END DO
      END IF

      CALL GETMEM('SBUF','FREE','REAL',LSBUF,NSBUF)
      CALL GETMEM('SMAT','FREE','REAL',LSMAT,NBM)
      CALL GETMEM('SCTMP','FREE','REAL',LSCTMP,NBM)
      CALL GETMEM('OVL','FREE','REAL',LOVL,NBM)
*
 900  CONTINUE
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Leaving ',ROUTINE
      END IF
      Call qExit('ONCMO')
*
      Return
      End
