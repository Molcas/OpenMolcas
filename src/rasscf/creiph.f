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
      SUBROUTINE CREIPH
C
C     RASSCF program: version IBM 3090: input section
C
C     Purpose: to set up an result file (JOBIPH = 15) where the most
C              important results of the calculation are collected as
C              an interphase to other programs. JOBIPH is a direct
C              acces file where the first record is a table of content.
C              The address to this record is zero. The table of content
C              has presently the length 15 and is stored in
C              COMMON/INTAUX/ as IADR15(15).
C     Subsections of JOBIPH:
C     IADR15(1):  List of input data written here with WR_RASSCF_Info
C     IADR15(2):  Molecular orbitals first written here and then updated
C                 at the end of SXCTL in each MCSCF itration.
C                 To be used as starting orbitals for restarts|
C                 Modified finally by NEWORB to pseudo-nat orbitals
C                 Average occupation numbers appended.
C     IADR15(3):  density matrices for active orbitals
C                 (4 matrices per root)
C                 D  : one-body density matrix
C                 DS : spin density matrix
C                 P  : symmetric two-body density matrix
C                 PA : antsymmetric two-body density matrix
C     IADR15(4):  The CI-vectors corresponding to the LROOTS lowest
C                 roots written in CICTL after each iteration.
C     IADR15(5):  The Fock matrix for the occupied orbitals.
C                 written in FOCKOC at the end of the calculation.
C     IADR15(6):  Energies given in array ENER(10,50) in COMMON/RELAUX/.
C     IADR15(7):  Convergence parameters in array CONV(6,50) in RELAUX.
C     IADR15(8):  Not used
C     IADR15(9):  Canonical MO's constructed and written in FCKPT2 at
C                 the end of the calculation.
C     IADR15(10): Inactive Fock matrix for CASPT2 in MO basis
C                 Note: Frozen and deleted orbitals not included,
C                 Modified to include the CASPT2 matrix FP 900608 (BOR)
C     IADR15(11): The diagonal of the CASPT Fock matrix.
C     IADR15(12): Natural orbitals for the final wave functions.
C                 Note: These orbitals do not correspond to the same
C                 CI wave function and should only be used to compute
C                 properties.
C                 One set of orbitals is stored for each root.
C                 Each set of NO's is followed by a record
C                 containing the corresponding occupation numbers.
C     IADR15(13): One-body spin density matrices (one for each root)
C     IADR15(14): Spin orbitals for the final wave functions.
C                 Note: These orbitals do not correspond to the same
C                 CI wave function and should only be used to compute
C                 spin properties.
C                 One set of orbitals for each root in an average
C                 calculation. Stored sequentially (891219).
C                 Each set of NO's is followed by a record
C                 containing the corresponding occupation numbers.
C     IADR15(15): First unwritten disk address, if positive;
C                 If = -1, it means a new JOBIPH structure is used,
C                 same as before, except:
C     IADR15(16): First unwritten disk address (new layout),
C     IADR15(17): Effective Hamiltonian, written here by MS-CASPT2,
C                 LROOTS**2 elements, Fortran standard layout.
C     IADR15(18): Table that translates between levels and absolute
C                 orbital indices, setup by setsxci.
C     IADR15(19)--IADR15(30): Presently unused.
C     ********** IBM 3090 MOLCAS Release 90 02 22 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "general.fh"
      Common /IDSXCI/ IDXCI(mxAct),IDXSX(mxAct)

      Call qEnter('CREIPH')

      DO I=1,15
       IADR15(I)=0
      END DO
C Dummy write table of contents.
C New layout scheme; length is 30 integers.
      IAD15=0
      CALL IDAFILE(JOBIPH,1,IADR15,30,IAD15)
      IADR15(1)=IAD15
C
C     DUMMY WRITE THE REMAINING RECORDS TO OBTAIN ADDRESSES
C
      CALL WR_RASSCF_Info(JOBIPH,1,iAD15,NACTEL,ISPIN,NSYM,LSYM,
     &            NFRO,NISH,NASH,NDEL,NBAS,MxSym,
     &            NAME,LENIN8*mxOrb,NCONF,HEADER,144,
     &            TITLE,4*18*mxTit,POTNUC,LROOTS,NROOTS,
     &            IROOT,MxRoot,NRS1,NRS2,NRS3,
     &            NHOLE1,NELEC3,IPT2,WEIGHT)
      IADR15(2)=IAD15
C
      CALL DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
      CALL DDAFILE(JOBIPH,0,DUM,NTOT,IAD15)
      IADR15(3)=IAD15
      Do i = 1,lRoots
        Call DDafIle(JOBIPH,0,Dum,NACPAR,IAD15)
        Call DDafIle(JOBIPH,0,Dum,NACPAR,IAD15)
        Call DDafIle(JOBIPH,0,Dum,NACPR2,IAD15)
        Call DDafIle(JOBIPH,0,Dum,NACPR2,IAD15)
      End Do
      IADR15(4)=IAD15
      Do i = 1,lRoots
        Call DDafIle(JOBIPH,0,Dum,nConf,IAD15)
      End Do
      IADR15(5)=IAD15
      NFOCK=0
      DO ISYM=1,NSYM
       NOO=NISH(ISYM)+NASH(ISYM)
       NFOCK=NFOCK+NOO**2
      END DO
      CALL DDAFILE(JOBIPH,0,DUM,NFOCK,IAD15)
      IADR15(6)=IAD15
      CALL DDAFILE(JOBIPH,0,DUM,mxRoot*mxIter,IAD15)
      IADR15(7)=IAD15
      CALL DDAFILE(JOBIPH,0,DUM,6*mxIter,IAD15)
      IADR15(8)=IAD15
      IADR15(9)=IAD15
      CALL DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
      IADR15(10)=IAD15
      CALL DDAFILE(JOBIPH,0,DUM,NTOT3,IAD15)
      CALL DDAFILE(JOBIPH,0,DUM,NTOT3,IAD15)
      IADR15(11)=IAD15
      CALL DDAFILE(JOBIPH,0,DUM,NORBT,IAD15)
      IADR15(12)=IAD15
      Do i = 1,lRoots
        CALL DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
        CALL DDAFILE(JOBIPH,0,DUM,NTOT,IAD15)
      End Do
      IADR15(13)=IAD15
      IADR15(14)=IAD15
      Do i=1,lroots
        CALL DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
        CALL DDAFILE(JOBIPH,0,DUM,NTOT,IAD15)
      End Do
C New layout scheme:
      IADR15(15)=-1
      CALL GETMEM('HEFF','ALLO','REAL',LHEFF,LROOTS**2)
      DO J=1,LROOTS
       DO I=1,LROOTS
        WORK(LHEFF-1+I+LROOTS*(J-1))=0.0D0
       END DO
       WORK(LHEFF-1+J+LROOTS*(J-1))=1.0D12
      END DO
      IADR15(17)=IAD15
      CALL DDAFILE(JOBIPH,1,WORK(LHEFF),LROOTS**2,IAD15)
      CALL GETMEM('HEFF','FREE','REAL',LHEFF,LROOTS**2)
      IADR15(18)=IAD15
CSVC: translates levels to orbital index (L2ACT in caspt2)
      CALL IDAFILE(JOBIPH,1,IDXSX,mxAct,IAD15)
CSVC: translates orbital index to levels (LEVEL in caspt2)
      CALL IDAFILE(JOBIPH,1,IDXCI,mxAct,IAD15)
      DO I=19,30
       IADR15(I)=0
      END DO
C First unused disk address at IADR(16):
      IADR15(16)=IAD15
C
C     Write table of content: IADR15
C
      IAD15=0
      CALL IDAFILE(JOBIPH,1,IADR15,30,IAD15)
C
      Call qExit('CREIPH')
      RETURN
      END
