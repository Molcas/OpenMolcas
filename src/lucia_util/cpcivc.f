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
      SUBROUTINE cpcivc(ifile, mxrec, isym, iway,lrec)
*
* Copies the CI-vector between Molcas Rasscf and Lucia enviroment
* IWAY = 1: from Molcas to Lucia (from core to disk unit ifile).
* IWAY = 2: from Lucia to Molcas (from disk unit ifile to core).
*
      implicit real*8 (a-h,o-z)
      dimension lrec(mxrec)
#include "rasscf_lucia.fh"
#include "rasdim.fh"
#include "WrkSpc.fh"
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
      ioff = 0
      IF (iway.eq.1) then
         NTEST = 0
         IF (NTEST .GE. 100) THEN
            write(6,*) 'Jesper: CI-vector put to disk:'
            DO IREC = 1,NREC
               IF(LREC(IREC) .GE. 0) THEN
                  call wrtmat(work(c_pointer+ioff),1,lrec(irec),
     &               1,lrec(irec))
               ioff = ioff + lrec(irec)
               ENDIF
            ENDDO
         ENDIF

         CALL todscn(work(c_pointer), nrec, lrec,
     &        -1, ifile)
         CALL itods(-1,1,-1,ifile)
*
*   ==============================
*      Read CI-vector from disc
*   ==============================
*
      ELSE
         CALL frmdscn(work(c_pointer), nrec, -1, ifile)
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE cpsivc(ifile, mxrec, vec,lrec)
*
* Copies the Sigma-vector between Molcas Rasscf and Lucia enviroment
*
      implicit real*8 (a-h,o-z)
      dimension lrec(mxrec)
      dimension vec(mxrec)
#include "rasdim.fh"
#include "rasscf_lucia.fh"
#include "general.fh"
#include "io_util.fh"
*
*   ========================
*      Find nrec and lrec
*   ========================
*
      CALL blkfo_min(lsym, nrec, lrec)
      IDISK(IFILE)=0
*
*   =================================
*      Read Sigma-vector from disc
*   =================================
*
      CALL frmdscn(vec, nrec, -1, ifile)
*
      RETURN
      END

      SUBROUTINE CP_ONE_INT(W1,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W1(NDIM)
#include "mxpdim.fh"
#include "orbinp.fh"
#include "WrkSpc.fh"
#include "rasscf_lucia.fh"

      CALL DCOPY_(NTOOB**2,0.0D0,0,Work(kint1_pointer),1)
      CALL DCOPY_(NDIM,W1,1,Work(kint1_pointer),1)
      CALL DCOPY_(NTOOB**2,Work(kint1_pointer),1,Work(kint1o_pointer),1)

      Return
      End
