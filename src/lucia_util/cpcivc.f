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
      SUBROUTINE cpcivc(CIVec,nCIVEC,ifile, mxrec, isym, iway,lrec)
*
* Copies the CI-vector between Molcas Rasscf and Lucia enviroment
* IWAY = 1: from Molcas to Lucia (from core to disk unit ifile).
* IWAY = 2: from Lucia to Molcas (from disk unit ifile to core).
*
      implicit real*8 (a-h,o-z)
      integer nCIVEC, ifile, mxrec, isym, iway
      integer lrec(mxrec)
      real*8 CIVec(nCIVec)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "io_util.fh"
*
*   ========================
*      Find nrec and lrec
*   ========================
*
      CALL blkfo_min(isym, nrec, lrec)
      IDISK(IFILE)=0
*
*   =============================
*      Write CI-vector to disc
*   =============================
*
      IF (iway.eq.1) then
#ifdef _DEBUGPRINT_
         ioff = 1
         write(6,*) 'CI-vector put to disk:'
         DO IREC = 1,NREC
            IF(LREC(IREC) .GE. 0) THEN
               call wrtmat(CIVec(ioff),1,lrec(irec),
     &            1,lrec(irec))
            ioff = ioff + lrec(irec)
            ENDIF
         ENDDO
#endif
         CALL todscn(CIVec, nrec, lrec,-1, ifile)
         CALL itods([-1],1,-1,ifile)
*
*   ==============================
*      Read CI-vector from disc
*   ==============================
*
      ELSE
         CALL frmdscn(CIVec, nrec, -1, ifile)
      ENDIF
*
      END
*
      SUBROUTINE cpsivc(ifile, mxrec, vec,lrec)
      use rasscf_lucia
*
* Copies the Sigma-vector between Molcas Rasscf and Lucia enviroment
*
      implicit real*8 (a-h,o-z)
      dimension lrec(mxrec)
      dimension vec(mxrec)
#include "rasdim.fh"
#include "general.fh"
#include "io_util.fh"
*
*   ========================
*      Find nrec and lrec
*   ========================
*
      CALL blkfo_min(stsym, nrec, lrec)
      IDISK(IFILE)=0
*
*   =================================
*      Read Sigma-vector from disc
*   =================================
*
      CALL frmdscn(vec, nrec, -1, ifile)
*
      END

      SUBROUTINE CP_ONE_INT(W1,NDIM)
      use GLBBAS, only: INT1, INT1O
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer nDIM
      Real*8 W1(NDIM)
#include "mxpdim.fh"
#include "orbinp.fh"

      INT1(:)=0.0D0
      INT1(1:NDIM)=W1(1:NDIM)
      INT1O(:)=0.0D0
      INT1O(:)=INT1(:)

      End
