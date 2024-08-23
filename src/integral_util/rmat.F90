!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!   K.P.
!   Dieses File enthaelt Inputdaten fuer die Modifikation der
!   Ein-Elektronen-Integrale fuer die R-Matrix.
!                                                                      *
!***********************************************************************
!                                                                      *
!     The following section should be activated if a NAG library is
!     present. At that point I suggest that it gets it's own include
!     file. -RL-
!
!     Parameter(INTparm=200)
!     REAL*8  wrk1(4*INTparm)
!     Integer Iwrk1(INTparm)
!                                                                      *
!***********************************************************************
!                                                                      *
!     RmatR    : radius of the R-matrix sphere (bohr)
!     Epsabs   : absolute precision of numerical radial integration
!     Epsrel   : relative precision of numerical radial integration
!     qCoul    : effective charge of the target molecule
!     dipol(3) : effective dipole moment of the target molecule
!     dipol1   : abs sum of dipol(i)
!     epsq     : minimal value of qCoul and/or dipol1 to be considered
!     bParm    : Bloch term parameter
!
Module RMat
      Real*8         RmatR,Epsabs,Epsrel,qCoul,Epsq,bParm,dipol(3),     &
     &               Dipol1
!
!     keyr: option to the dqag quadpack routine
!
      Integer        keyr
!
!     Quadpack: logical flag to indicate use of QUADPACK routines
!     NagInt  : logical flag to indicate use of NAG routines
!     testint : logical flag to indicate use of both QUADPACK and NAG
!               routines
!     RMat_On : Logical flag to signal that R-matrix type integrals
!               are to be computed.
!
      Logical        Quadpack,nagint,testint,RMat_On
      Logical        RMat_Type_Integrals

!   K.P.
!   Dieses File enthaelt Feld fuer die theta/phi  Integration der
!   Ein-Elektronen-Integrale fuer die R-Matrix.

Integer lgamma,n_gam,m_gam
Integer, Parameter :: lgammax=15
Real*8 gammath(-2:2*lgammax+3,-2:2*lgammax+4),                    &
       gammaph(-2:2*lgammax+3,-2:2*lgammax+4)

!   K.P.
!   Dieses File enthaelt Inputdaten fuer die Modifikation der
!   Ein-Elektronen-Integrale fuer die R-Matrix.
Real*8 expsum
Integer l,lcost,lsint,lcosf,lsinf

End Module RMat
