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
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE T_ROW_TO_H(T,H,K,TKK)
*
* Set H integrals
*
*    Column K : H(P,K) = T(P,K)/T(K,K), P.NE.K
*    Other Columns     = 0
* - and return T_{kk} in TKK
*
*
* Jeppe Olsen, Jan 98
* For rotation of CI vectors
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "mxpdim.fh"
#include "glbbas.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "lucinp.fh"
*. Input ( in blocked form)
      DIMENSION T(*)
*. Output ( also in blocked form)
      DIMENSION H(*)
*
      KSM = ISMFSO(K)
      KOFF = IBSO(KSM)
      KREL = K - KOFF + 1
      NK = NTOOBS(KSM)

*
      ZERO = 0.0D0
      CALL SETVEC(H,ZERO,NTOOB**2)
*
      IOFF = IFRMR(IWORK(KPGINT1A(1)),1,KSM)
      CALL COPVEC(T(IOFF+(KREL-1)*NK),H(IOFF+(KREL-1)*NK),NK)
      TKK = H(IOFF-1+(KREL-1)*NK+KREL)
      IF(TKK .NE. 0.0D0) THEN
        FAC = 1.0D0/TKK
        CALL SCALVE(H(IOFF+(KREL-1)*NK),FAC,NK)
C       H(IOFF-1+(K-1)*NK+K) = H(IOFF-1+(K-1)*NK+K) -1.0D0
        H(IOFF-1+(KREL-1)*NK+KREL) = 0.0D0
      ELSE
C       TKK = 1.0D0
        TKK = 0.0D0
      END IF
*
*
      RETURN
      END
