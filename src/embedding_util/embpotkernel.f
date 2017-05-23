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
* Copyright (C) Thomas Dresselhaus                                     *
************************************************************************
      SubRoutine EmbPotKernel(expA,nPrimA,expB,nPrimB,Zeta,ZInv,rKappa,
     &                 P,Final,nZeta,nIC,nComp,lA,lB,coordA,coordB,nRys,
     &                 Array,nArr,CCoor,nOrdOp,lOper,iChO,
     &                 iStabM,nStabM,
     &                 PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine to calculate integrals over an embedding      *
*         potential on a grid.                                         *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : GetMem                                                  *
*                                                                      *
*     Author: Thomas Dresselhaus                                       *
*                                                                      *
************************************************************************

      Implicit None

****** Includes
#include "WrkSpc.fh"

      ! Holds the data which is read in in this subroutine
#include "embpotdata.fh"

******* Variables ******************************************************

      ! Number of primitives
      Integer nPrimA, nPrimB

      ! Exponents of primitives
      Real*8 expA(nPrimA), expB(nPrimB)

      ! Angular momentum
      Integer lA, lB

      ! Output: <A|sum_(gridpoints){potValue*weight}|B>
      Real*8 Final(nPrimA, nPrimB,(lA+1)*(lA+2)/2,(lB+1)*(lB+2)/2)

      ! Orbital coordinates
      Real*8 coordA(3), coordB(3)


      ! Other input variables (unused, unknown)
      Integer nZeta,nIC,nComp,nRys,nArr,nOrdOp,nStabM,nGrid,iAddPot
      Real*8 Zeta(nZeta), ZInv(nZeta),
     &       rKappa(nZeta), P(nZeta,3), CCoor(3,nComp),
     &       Array(nZeta*nArr)
      Integer iStabM(0:nStabM-1),lOper(nComp),iChO(nComp)
      Real*8 PtChrg

****** Local Variables

      ! Index
      Integer i, j, a, b, m
      Integer ix, iy, iz

      Integer nShellA, nShellB

      ! Positions in the work array
      Integer posRadA, posRadB, posSphA, posSphB

      ! distance of the current grid point to A/B
      Real*8 dRA(3), dRB(3)

      Real*8 prefactor

      ! Function return value
      Real*8 gaussRad

****** Initialization ***************************************************

      nShellA = (lA+1)*(lA+2)/2
      nShellB = (lB+1)*(lB+2)/2

      ! Intermediate arrays, radial and angular part of gaussians
      Call GetMem('radA','ALLO','REAL',posRadA,nPrimA)
      Call GetMem('radB','ALLO','REAL',posradB,nPrimB)
      Call GetMem('sphA','ALLO','REAL',posSphA,nShellA)
      Call GetMem('sphB','ALLO','REAL',posSphB,nShellB)


****** Calculation ******************************************************

      ! Loop over primitives and angular momenta, init result var
      do a=1, nPrimA
         do b=1, nPrimB
            do i=1, nShellA
               do j=1, nShellB
                  Final(a,b,i,j) = 0.0d0
               end do
            end do
         end do
      end do

      ! Now loop over grid points first
      do m=1, nEmbGridPoints
       ! Calculate distances from gaussian centers to current grid point
       do i=1, 3
        dRA(i) = Work(posEmbGridCoord+(m-1)*3+i-1) - coordA(i)
        dRB(i) = Work(posEmbGridCoord+(m-1)*3+i-1) - coordB(i)
       end do

       prefactor = Work(posEmbPotVal+m-1) *
     &             Work(posEmbWeight+m-1)

       ! Precalculate radial factors
       do a=1, nPrimA
         Work(posRadA+a-1) = gaussRad(expA(a), dRA)
       end do

       do b=1, nPrimB
         Work(posRadB+b-1) = gaussRad(expB(b), dRB)
       end do

       ! Precalculate spherial factors
        ! shell A
       i = 0
       do ix = lA, 0, -1
          do iy = lA - ix , 0, -1
             iz = lA - ix - iy
             i = i + 1
             Work(posSphA+i-1) = dRA(1)**ix * dRA(2)**iy * dRA(3)**iz
          end do
       end do
        ! shell B
       i = 0
       do ix = lB, 0, -1
          do iy = lB - ix , 0, -1
             iz = lB - ix - iy
             i = i + 1
             Work(posSphB+i-1) = dRB(1)**ix * dRB(2)**iy * dRB(3)**iz
          end do
       end do

       ! Put them together
       do a=1, nPrimA
         do b=1, nPrimB
            do i=1, nShellA
               do j=1, nShellB
                  Final(a,b,i,j) = Final(a,b,i,j) +
     &                             (prefactor *
     &                              Work(posRadA+a-1) *
     &                              Work(posRadB+b-1) *
     &                              Work(posSphA+i-1) *
     &                              Work(posSphB+j-1))
               end do
            end do
         end do
       end do

       if (embDebug) then
         if (mod(m,587).eq.0) then
          write(6,*) "-------------------------------------------------"
          write(6,*) "m=",m
          write(6,*) "A=", coordA(1), coordA(2), coordA(3)
          write(6,*) "r=",Work(posEmbGridCoord+(m-1)*3),
     &               Work(posEmbGridCoord+(m-1)*3+1),
     &               Work(posEmbGridCoord+(m-1)*3+2)
          write(6,*) "dRA=", dRA(1), dRA(2), dRA(3)
          write(6,*) "expA(1)=", expA(1), "; gaussRad(...)=",
     &                gaussRad(expA(1), dRA)
          write(6,*) "-------------------------------------------------"
          write(6,*) "m=", m, "prefactor=", prefactor
          write(6,*) "radA(1)=",Work(posRadA),"; sphA(1)=",Work(posSphA)
     &        ,", radB(2)=",Work(posRadB+1),"; sphB(2)=",Work(posSphB+1)
          write(6,*) "-------------------------------------------------"
         end if
       end if

      end do ! m, grid point

      if (embDebug) then
         write(6,*) "-------------------------------------------------"
       do a=1, nPrimA
         do b=1, nPrimB
            do i=1, nShellA
               do j=1, nShellB
                  write(6,*) Final(a,b,i,j)
               end do
            end do
         end do
       end do
!         write(*,*) "Final(1,1,1,1)=", Final(1,1,1,1),
!     &              "; Final(1,2,1,2)=", Final(1,2,1,2)
         write(6,*) "-------------------------------------------------"
         write(6,*) "-------------------------------------------------"
      end if



****** Done. Tidy up. ****************************************************

      Call GetMem('radA','FREE','REAL',posRadA,nPrimA)
      Call GetMem('radB','FREE','REAL',posradB,nPrimB)
      Call GetMem('sphA','FREE','REAL',posSphA,nShellA)
      Call GetMem('sphB','FREE','REAL',posSphB,nShellB)

      Return
! Avoid unused argument warnings
      If (.False.) Then
        Call unused_real_array(Zeta)
        Call unused_real_array(ZInv)
        Call unused_real_array(rKappa)
        Call unused_real_array(P)
        Call unused_integer(nIC)
        Call unused_integer(nRys)
        Call unused_real_array(Array)
        Call unused_real_array(CCoor)
        Call unused_integer(nOrdOp)
        Call unused_integer_array(lOper)
        Call unused_integer_array(iCho)
        Call unused_integer_array(iStabM)
        Call unused_real(PtChrg)
        Call unused_integer(nGrid)
        Call unused_integer(iAddPot)
      End If
      End


**************************************************************************
*     Returns the radial part of the value of a GTO with given exponent, *
*     centered at the origin.                                            *
**************************************************************************
      Real*8 Function gaussRad(alpha, r)

      Implicit None

      Real*8 alpha
      Real*8 r(3)

      Real*8 rSquare
      Integer i

      rSquare=0.0d0
      do i=1, 3
        rSquare = rSquare + r(i)**2
      end do

      ! Radial part
      gaussRad = exp(-alpha*rSquare)

      Return
      End
