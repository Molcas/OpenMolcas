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
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE ORBCTL(CMO)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "chocaspt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      REAL*8 CMO(NCMO)
      INTEGER ISYM
* PAM Feb 2015 NTORB, LTORB are in Include/caspt2.fh!
*      INTEGER NTORB, LTORB
      INTEGER I1,I2,LORBE
      INTEGER IDISK
      REAL*8  OCC_DUM(1)

C Calculate transformation matrix to PT2 orbitals, defined as those
C that have standard Fock matrix FIFA diagonal within inactive,
C active, and secondary subblocks.
      CALL QENTER('ORBCTL')

c Determine PT2 orbitals, and transform CI coeffs.
      IF(IPRGLB.GE.DEBUG) THEN
       WRITE(6,*)' ORBCTL calling MKRPTORB...'
      END IF
* The CMO coefficient array is changed by ortonormal transformations of
* each of the inactive,Ras1,Ras2,Ras3,and secondary (i.e. virtual) orbitals
* in each symmetry. The transformation matrices are stored as a sequence of
* square matrices in TORB. The transformation is such that each of the
* diagonal blocks of the Fock matrix FIFA is diagonalized.
* MKRPTORB will at the same time transform each of the CI arrays on file
* such that, with new CMO vectors, they still represent the original
* wave function.

* The CI arrays are on file with unit number LUCIEX. There is NSTATE
* CI arrays, stored sequentially. The original set starts at disk address
* IDCIEX, the transformed ones are written after IDTCEX.
      CALL MKRPTORB(WORK(LFIFA),WORK(LTORB),CMO)
      IF(IPRGLB.GE.DEBUG) THEN
       WRITE(6,*)' ORBCTL back from MKRPTORB.'
      END IF

* Use the transformation matrices to change the HONE, FIMO, and FIFA arrays:
      CALL TRANSFOCK(WORK(LTORB),WORK(LHONE),1)
      CALL TRANSFOCK(WORK(LTORB),WORK(LFIMO),1)

* When doing XMS, FAMO refers only to the last state, therefore it's wrong!
* However, we never use it anywhere else...
      ! CALL TRANSFOCK(WORK(LTORB),WORK(LFAMO),1)
*****

      CALL TRANSFOCK(WORK(LTORB),WORK(LFIFA),1)

* When doing XMS, DREF refers to the last state considered and it is not the
* state average density, therefore it's wrong to transform it!
* However, it is never used again in this part, and next time it is used, it
* is actually recomputed for the right place.
      ! CALL TRANSDREF(WORK(LTORB),WORK(LDREF))
*****

* DREF is not really used for anything important in MKEPS, this is why we don't
* care that we pass the wrong one in...
      CALL MKEPS(WORK(LFIFA),WORK(LDREF))

      IF(IPRGLB.GE.DEBUG) THEN
       WRITE(6,*)' ORBCTL back from TRANSFOCK.'
      END IF

C Save new MO coeffs, and the transformation matrices:
      IDISK=IAD1M(2)
      CALL DDAFILE(LUONEM,1,WORK(LCMO),NCMO,IDISK)
      IAD1M(4)=IEOF1M
      IDISK=IAD1M(4)
      CALL DDAFILE(LUONEM,1,WORK(LTORB),NTORB,IDISK)
      IEOF1M=IDISK

c Print new orbitals. First, form array of orbital energies.
      CALL GETMEM('ORBE','ALLO','REAL',LORBE,NBAST)
      I1=1
      I2=1
      DO ISYM=1,NSYM
        IF(NFRO(ISYM).GT.0) THEN
          CALL DCOPY_(NFRO(ISYM),[0.0D0],0,WORK(LORBE-1+I2),1)
          I2=I2+NFRO(ISYM)
        END IF
        IF(NORB(ISYM).GT.0) THEN
          CALL DCOPY_(NORB(ISYM),EPS(I1),1,WORK(LORBE-1+I2),1)
          I1=I1+NORB(ISYM)
          I2=I2+NORB(ISYM)
        END IF
        IF(NDEL(ISYM).GT.0) THEN
          CALL DCOPY_(NDEL(ISYM),[0.0D0],0,WORK(LORBE-1+I2),1)
          I2=I2+NDEL(ISYM)
        END IF
      END DO
c Then call utility routine PRIMO.
      IF ( IPRGLB.GE.VERBOSE ) THEN
        WRITE(6,*) ' The internal wave function representation has'//
     &             ' been changed to use quasi-canonical orbitals:'
        WRITE(6,*) ' those which diagonalize the Fock matrix within'//
     &             ' inactive-inactive,'
        WRITE(6,*) ' active-active and virtual-virtual submatrices.'
        IF(.NOT. PRORB) THEN
          WRITE(6,*)' On user''s request, the quasi-canonical orbitals'
          WRITE(6,*)' will not be printed.'
          GOTO 99
        END IF
      END IF

      IF ( IPRGLB.GE.VERBOSE) THEN
C Print orbitals. Different options:
        IF ( OUTFMT.EQ.'LONG    ' ) THEN
          THRENE=2.0d0**31
          THROCC=-2.0d0**31
        ELSE IF ( OUTFMT.EQ.'DEFAULT ' ) THEN
          THRENE=5.0d+00
          THROCC=5.0d-04
        END IF
        CALL PRIMO(' Quasi-canonical orbitals',.FALSE.,.TRUE.,
     &              THROCC,THRENE,NSYM,NBAS,NBAS,NAME,
     &              WORK(LORBE),OCC_DUM,CMO,-1)
      END IF

  99  CONTINUE
      CALL GETMEM('ORBE','FREE','REAL',LORBE,NBAST)

      CALL QEXIT('ORBCTL')

      RETURN
      END
