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
      SUBROUTINE GTJK_MCLR(RJ,RK)
      use Arrays, only: Int2
C
C     PURPOSE: GET ALL INTEGRALS COULOMB AND EXCHANGE INTEGRALS
C              WITH THE CHARGE DISTRIBUTION JK
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "detdim.fh"
#include "Pointers.fh"
#include "orbinp_mclr.fh"
#include "glbbas_mclr.fh"
#include "WrkSpc.fh"
      DIMENSION RJ(NACOB,NACOB),RK(NACOB,nACOB)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
C
C     FORM THE COULOMB (RJ) AND EXCHANGE (RK) INTEGRAL MATRICES FROM
C     THE TWO-ELECTRON INTEGRAL LIST
C
      DO NT=1,NACOB
         DO  NU=1,NT
            NTUK=itri(itri(nt,nt),itri(nu,nu))
            RJ(NT,NU)=INT2(NTUK)
            RJ(NU,NT)=INT2(NTUK)

            NTUJ=itri(itri(nt,nu),itri(nu,nt))
            RK(NT,NU)=INT2(NTUJ)
            RK(NU,NT)=INT2(NTUJ)
         End Do
      End Do
C
      RETURN
      END
