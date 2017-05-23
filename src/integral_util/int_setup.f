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
      Subroutine Int_Setup(iSD,nSkal,iS,jS,kS,lS,
     &                     Coor,Shijij,
     &                     iAngV,iCmpV,iShelV,iShllV,iAOV,iStabs)
      Implicit Real*8 (a-h,o-z)
*
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
*
#include "nsd.fh"
#include "FMM.fh"
      Integer iSD(0:nSD,nSkal)
*
      Real*8  Coor(3,4)
      Integer iAngV(4),iCmpV(4),iShelV(4),iShllV(4),iAOV(4),iStabs(4),
     &        jQuad(4)
      Logical Shijij
*
      iCnttp=iSD(13,iS)
      kCnttp=iSD(13,kS)
*
      If (AuxCnttp(iCnttp)) Then
         call dcopy_(3,Work(iSD(8,jS)),1,Coor(1,1),1)
      Else
         call dcopy_(3,Work(iSD(8,iS)),1,Coor(1,1),1)
      End If
      call dcopy_(3,Work(iSD(8,jS)),1,Coor(1,2),1)
*
      If (AuxCnttp(kCnttp)) Then
         call dcopy_(3,Work(iSD(8,lS)),1,Coor(1,3),1)
      Else
         call dcopy_(3,Work(iSD(8,kS)),1,Coor(1,3),1)
      End If
      call dcopy_(3,Work(iSD(8,lS)),1,Coor(1,4),1)
*
      Shijij=(iSD(0,iS).eq.iSD(0,kS).and.iSD(10,iS).eq.iSD(10,kS))
     &       .and.
     &       (iSD(0,jS).eq.iSD(0,lS).and.iSD(10,jS).eq.iSD(10,lS))
*
      jQuad(1)=iS
      jQuad(2)=jS
      jQuad(3)=kS
      jQuad(4)=lS
      Do iQuad = 1, 4
         iSkal=jQuad(iQuad)
         iAngV(iQuad)  = iSD( 1,iSkal)
         iCmpV(iQuad)  = iSD( 2,iSkal)
         iAOV(iQuad)   = iSD( 7,iSkal)
         iStabs(iQuad) = iSD(10,iSkal)
         iShelV(iQuad) = iSD(11,iSkal)
         iShllV(iQuad) = iSD( 0,iSkal)
      End Do
CMAW start
*
*  For the FMM coulomb integrals <AB(r1)|1/r12|CD(r2)>
*  Here we flag the integral routines that we only want to compute
*  the short-range non-multipole component of integrals over this
*  shell quartet if midpoint(A,B) is sufficiently far from
*  midpoint(C,D) for numerical stability.
*  Note that midpoint(X,Y) corresponds to the multipole expansion
*  centre of an XY AO-pair, regardless of exponents.
*
      FMM_shortrange = .False.
      If (DoFMM) Then
         D = 0.0d0
         DO i = 1, 3
            P = (Coor(i,1) + Coor(i,2))/2.0d0    ! AB shell-pair
            Q = (Coor(i,3) + Coor(i,4))/2.0d0    ! CD shell-pair
            D = D + (P-Q)*(P-Q)
         End Do
         IF (D .gt. RPQMIN*RPQMIN) FMM_shortrange = .True.
      End If
CMAW end
*
      Return
      End
