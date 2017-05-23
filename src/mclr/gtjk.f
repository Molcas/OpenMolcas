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
C
C
      SUBROUTINE GTJK_MCLR(RJ,RK)
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

      DO  NT=1,NACOB
      DO  NU=1,NT
       NTUK=itri(itri(nt,nt),itri(nu,nu))-1
       RJ(NT,NU)=Work(K2INT+NTUK)
       RJ(NU,NT)=Work(K2INT+NTUK)
C
       NTUJ=itri(itri(nt,nu),itri(nu,nt))-1
       RK(NT,NU)=Work(K2INT+NTUJ)
       RK(NU,NT)=Work(K2INT+NTUJ)
      End Do
      End Do
C
C     EXIT
C
      RETURN
      END
