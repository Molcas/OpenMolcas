!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Thomas Dresselhaus                                     *
!***********************************************************************

subroutine EmbPotKernel(expA,nPrimA,expB,nPrimB,Zeta,ZInv,rKappa,P,finalval,nZeta,nIC,nComp,lA,lB,coordA,coordB,nRys,Array,nArr, &
                        CCoor,nOrdOp,lOper,iChO,iStabM,nStabM,PtChrg,nGrid,iAddPot)
!***********************************************************************
!                                                                      *
! Object: kernel routine to calculate integrals over an embedding      *
!         potential on a grid.                                         *
!                                                                      *
! Called from: OneEl                                                   *
!                                                                      *
! Calling    : GetMem                                                  *
!                                                                      *
!     Author: Thomas Dresselhaus                                       *
!                                                                      *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
! Number of primitives
integer(kind=iwp), intent(in) :: nPrimA, nPrimB
! Exponents of primitives
real(kind=wp), intent(in) :: expA(nPrimA), expB(nPrimB)
! Angular momentum
integer(kind=iwp), intent(in) :: lA, lB
! Orbital coordinates
real(kind=wp), intent(in) :: coordA(3), coordB(3)
! Output: <A|sum_(gridpoints){potValue*weight}|B>
real(kind=wp), intent(out) :: finalval(nPrimA,nPrimB,(lA+1)*(lA+2)/2,(lB+1)*(lB+2)/2)
! Other input variables (unused, unknown)
integer(kind=iwp), intent(in) :: nZeta, nIC, nComp, nRys, nArr, nOrdOp, nStabM, nGrid, iAddPot, iStabM(0:nStabM-1), lOper(nComp), &
                                 iChO(nComp)
real(kind=wp), intent(in) :: Zeta(nZeta), ZInv(nZeta), rKappa(nZeta), P(nZeta,3), CCoor(3,nComp), Array(nZeta*nArr), PtChrg

!***** Includes
#include "WrkSpc.fh"

! Holds the data which is read in in this subroutine
#include "embpotdata.fh"

!***** Local Variables

! Index
integer(kind=iwp) :: i, j, a, b, m, ix, iy, iz, nShellA, nShellB
! Positions in the work array
integer(kind=iwp) :: posRadA, posRadB, posSphA, posSphB
! distance of the current grid point to A/B
real(kind=wp) :: dRA(3), dRB(3)
real(kind=wp) :: prefactor
! Function return value
real(kind=wp), external :: gaussRad

!***** Initialization ***************************************************

nShellA = (lA+1)*(lA+2)/2
nShellB = (lB+1)*(lB+2)/2

! Intermediate arrays, radial and angular part of gaussians
call GetMem('radA','ALLO','REAL',posRadA,nPrimA)
call GetMem('radB','ALLO','REAL',posradB,nPrimB)
call GetMem('sphA','ALLO','REAL',posSphA,nShellA)
call GetMem('sphB','ALLO','REAL',posSphB,nShellB)

!***** Calculation ******************************************************

! Init result var
finalval(:,:,:,:) = Zero

! Now loop over grid points first
do m=1,nEmbGridPoints
  ! Calculate distances from gaussian centers to current grid point
  do i=1,3
    dRA(i) = Work(posEmbGridCoord+(m-1)*3+i-1)-coordA(i)
    dRB(i) = Work(posEmbGridCoord+(m-1)*3+i-1)-coordB(i)
  end do

  prefactor = Work(posEmbPotVal+m-1)*Work(posEmbWeight+m-1)

  ! Precalculate radial factors
  do a=1,nPrimA
    Work(posRadA+a-1) = gaussRad(expA(a),dRA)
  end do

  do b=1,nPrimB
    Work(posRadB+b-1) = gaussRad(expB(b),dRB)
  end do

  ! Precalculate spherial factors
  ! shell A
  i = 0
  do ix=lA,0,-1
    do iy=lA-ix,0,-1
      iz = lA-ix-iy
      i = i+1
      Work(posSphA+i-1) = dRA(1)**ix*dRA(2)**iy*dRA(3)**iz
    end do
  end do
  ! shell B
  i = 0
  do ix=lB,0,-1
    do iy=lB-ix,0,-1
      iz = lB-ix-iy
      i = i+1
      Work(posSphB+i-1) = dRB(1)**ix*dRB(2)**iy*dRB(3)**iz
    end do
  end do

  ! Put them together
  do a=1,nPrimA
    do b=1,nPrimB
      do i=1,nShellA
        do j=1,nShellB
          finalval(a,b,i,j) = finalval(a,b,i,j)+(prefactor*Work(posRadA+a-1)*Work(posRadB+b-1)*Work(posSphA+i-1)*Work(posSphB+j-1))
        end do
      end do
    end do
  end do

  if (embDebug) then
    if (mod(m,587) == 0) then
      write(u6,*) '-------------------------------------------------'
      write(u6,*) 'm=',m
      write(u6,*) 'A=',coordA(1),coordA(2),coordA(3)
      write(u6,*) 'r=',Work(posEmbGridCoord+(m-1)*3),Work(posEmbGridCoord+(m-1)*3+1),Work(posEmbGridCoord+(m-1)*3+2)
      write(u6,*) 'dRA=',dRA(1),dRA(2),dRA(3)
      write(u6,*) 'expA(1)=',expA(1),'; gaussRad(...)=',gaussRad(expA(1),dRA)
      write(u6,*) '-------------------------------------------------'
      write(u6,*) 'm=',m,'prefactor=',prefactor
      write(u6,*) 'radA(1)=',Work(posRadA),'; sphA(1)=',Work(posSphA),', radB(2)=',Work(posRadB+1),'; sphB(2)=',Work(posSphB+1)
      write(u6,*) '-------------------------------------------------'
    end if
  end if

end do ! m, grid point

if (embDebug) then
  write(u6,*) '-------------------------------------------------'
  do a=1,nPrimA
    do b=1,nPrimB
      do i=1,nShellA
        do j=1,nShellB
          write(u6,*) finalval(a,b,i,j)
        end do
      end do
    end do
  end do
  !write(u6,*) 'Final(1,1,1,1)=',Final(1,1,1,1),'; Final(1,2,1,2)=',Final(1,2,1,2)
  write(u6,*) '-------------------------------------------------'
  write(u6,*) '-------------------------------------------------'
end if

!***** Done. Tidy up. ****************************************************

call GetMem('radA','FREE','REAL',posRadA,nPrimA)
call GetMem('radB','FREE','REAL',posradB,nPrimB)
call GetMem('sphA','FREE','REAL',posSphA,nShellA)
call GetMem('sphB','FREE','REAL',posSphB,nShellB)

return
! Avoid unused argument warnings
if (.false.) then
  call unused_real_array(Zeta)
  call unused_real_array(ZInv)
  call unused_real_array(rKappa)
  call unused_real_array(P)
  call unused_integer(nIC)
  call unused_integer(nRys)
  call unused_real_array(Array)
  call unused_real_array(CCoor)
  call unused_integer(nOrdOp)
  call unused_integer_array(lOper)
  call unused_integer_array(iCho)
  call unused_integer_array(iStabM)
  call unused_real(PtChrg)
  call unused_integer(nGrid)
  call unused_integer(iAddPot)
end if

end subroutine EmbPotKernel

!*********************************************************************
! Returns the radial part of the value of a GTO with given exponent, *
! centered at the origin.                                            *
!*********************************************************************
real(kind=wp) function gaussRad(alpha,r)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: alpha, r(3)

real(kind=wp) :: rSquare
integer(kind=iwp) :: i

rSquare = Zero
do i=1,3
  rSquare = rSquare+r(i)**2
end do

! Radial part
gaussRad = exp(-alpha*rSquare)

return

end function gaussRad
