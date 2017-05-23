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
      SUBROUTINE HMAT(C,HC,HH,HD,NDIM,NDIMH,NTRIAL)
C
C RASSCF program: version IBM-3090: SX section
C
C Purpose: Calculation of the Davidson Hamiltonian.
C The sigma vector HC is computed in SIGVEC
C
C ********** IBM-3090 Release 88 09 08 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),HC(*),HH(*)
      DIMENSION HD(NDIM)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "wadr.fh"
#include "WrkSpc.fh"
C
      CALL QENTER('HMAT')
C
C COMPUTE SIGMA VECTORS
C
      IST=1+NDIMH*NDIM
      CALL GETMEM('SXXX','ALLO','REAL',LXX,NLX)
      CALL GETMEM('SXC3','ALLO','REAL',LC,NSXS)
      CALL SIGVEC(C(IST),HC(IST),HD,WORK(LBM),WORK(LSXN),
     *            WORK(LG),WORK(LH),WORK(LDIA),
     *            WORK(LF1),WORK(LF2),WORK(LXX),
     *            WORK(LC),NTRIAL)
      CALL GETMEM('XXXX','FREE','REAL',LXX,NLX)
      CALL GETMEM('XXXX','FREE','REAL',LC,NSXS)
C
C Compute a new block of the HH matrix
C
      IJ=(NDIMH+NDIMH**2)/2
      L1=NDIMH+1
      NDIMH=NDIMH+NTRIAL
      DO I=L1,NDIMH
       JST=1
       DO J=1,I
        IJ=IJ+1
        HH(IJ)=DDOT_(NDIM,HC(IST),1,C(JST),1)
        JST=JST+NDIM
       END DO
       IST=IST+NDIM
      END DO
C
      CALL QEXIT('HMAT')
      RETURN
      END
