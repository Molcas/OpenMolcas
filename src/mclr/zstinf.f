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
      SUBROUTINE ZSTINF_MCLR(IPRNT)
      use Str_Info
*
* Set up common block /STINF/ from information in /STINP/
*
*=========
* Input
*=========
* Information in /STINP/ and /ORBINP/
*
*======================
* Output ( in /STINF/ )
*======================
* ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
*                    or annihilating (ISTAC(ITYP,1)) an electron
*                    from a string of type  ITYP . A zero indicates
*                    that this mapping is not included
*                    Only strings having the same ISTTP index are
*                    mapped
* NOCTYP(ITYP) : Number of occupation classes for given type
*
*
* NSTFTP(ITYP) : Number of strings of this type
*
* INUMAP(ITYP) : Mapping of string type to next more general type
* INDMAP(ITYP) : Mapping of string type to next more restricted type
*
*   / \           Zero order space                         !
*    !            Double excitations from reference space  !  Down
* Up !            single excitation from reference space   !
*    !            reference space                         \ /
*
#include "detdim.fh"
#include "orbinp_mclr.fh"
*
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
* ******************************************************************
* Mappings between strings with the same type ISTTP index , +/- 1 el
* ******************************************************************
      ISTAC(:,:)=0
      DO 90 ITYP = 1, NSTTYP-1
        IF(NELEC(ITYP+1).EQ.NELEC(ITYP)-1) THEN
          ISTAC(ITYP,1) = ITYP+1
          ISTAC(ITYP+1,2) = ITYP
        END IF
90    CONTINUE
*
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Type - type mapping array ISTAC '
        WRITE(6,*) ' =============================== '
        CALL IWRTMA(ISTAC,NSTTYP,2,MXPSTT,2)
      END IF
* **************************************************
*. Number of occupation classes and strings per type
* **************************************************
      DO 200 ITYP = 1,NSTTYP
        NOCTYP(ITYP) = (MXRS1(ITYP)-MNRS1(ITYP)+1)
     *               * (MXRS3(ITYP)-MNRS3(ITYP)+1)
200   CONTINUE
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Number of occupation classes per type '
        WRITE(6,*) ' ===================================== '
        CALL IWRTMA(NOCTYP,1,NSTTYP,1,NSTTYP)
      END IF
*
      DO 300 ITYP = 1, NSTTYP
        NSTFTP(ITYP) = NUMST3(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),
     *                        NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP))
300   CONTINUE
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Number of strings per  type '
        WRITE(6,*) ' =========================== '
        CALL IWRTMA(NSTFTP,1,NSTTYP,1,NSTTYP)
      END IF
* *****************************************************************
*. Mappings between strings containing the same number of electrons
* *****************************************************************
      INUMAP(:)=0
      INDMAP(:)=0
*. Mapping to and from zero order space
*     Note: some lines are commented out here since IARTP and IBRTP
*           have never been defined. (R. Lindh 2006)
C     INUMAP(IARTP(3,5)) = IAZTP
C     IF(IARTP(3,4).NE.0)  INUMAP(IARTP(3,4)) = IAZTP+1
C     IF(IARTP(3,3).NE.0)  INUMAP(IARTP(3,3)) = IAZTP+2
*
C     INUMAP(IBRTP(3,5)) = IBZTP
C     IF(IBRTP(3,4).NE.0)  INUMAP(IBRTP(3,4)) = IBZTP+1
C     IF(IBRTP(3,3).NE.0)  INUMAP(IBRTP(3,3)) = IBZTP+2
*
C     NAEL = NELEC(IAZTP)
C     INDMAP(IAZTP) = IARTP(3,5)
C     IF(NAEL.GE.1) INDMAP(IAZTP+1) = IARTP(3,4)
C     IF(NAEL.GE.2) INDMAP(IAZTP+2) = IARTP(3,3)
*
C     NBEL = NELEC(IBZTP)
C     INDMAP(IBZTP) = IBRTP(3,5)
C     IF(NBEL.GE.1) INDMAP(IBZTP+1) = IBRTP(3,4)
C     IF(NBEL.GE.2) INDMAP(IBZTP+2) = IBRTP(3,3)
*.Number of electrons compared to reference
      DO 450 IDEL = -4,2
      DO 430 IEX = 1,2
* Up mappings
C       IF(IARTP(IEX,IDEL+5).NE.0) THEN
C         INUMAP(IARTP(IEX,IDEL+5)) = IARTP(IEX+1,IDEL+5)
C       END IF
C       IF(IBRTP(IEX,IDEL+5).NE.0) THEN
C         INUMAP(IBRTP(IEX,IDEL+5)) = IBRTP(IEX+1,IDEL+5)
C       END IF
* Down mappings
C       IF(IARTP(IEX+1,IDEL+5).NE.0) THEN
C         INDMAP(IARTP(IEX+1,IDEL+5)) = IARTP(IEX,IDEL+5)
C       END IF
C       IF(IBRTP(IEX+1,IDEL+5).NE.0) THEN
C         INDMAP(IBRTP(IEX+1,IDEL+5)) = IBRTP(IEX,IDEL+5)
C       END IF
430   CONTINUE
450   CONTINUE
*
      IF(NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Up mappings of string types '
        CALL IWRTMA(INUMAP,1,NSTTYP,1,NSTTYP)
        WRITE(6,*) ' Down mappings of string types '
        CALL IWRTMA(INDMAP,1,NSTTYP,1,NSTTYP)
      END IF
*
*
      RETURN
      END
