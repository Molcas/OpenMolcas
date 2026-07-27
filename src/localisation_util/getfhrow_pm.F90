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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine GetFHrow_PM(nAtoms,nOrb2Loc,PA,fsdim,FHrow_k,k,CMO,nBasis)
! Purpose: compute the full FHrow_k of the Pipek-Mezey functional w.r.t. elements of the kappa matrix

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Eight
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc, fsdim, k, nBasis
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms), CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(out) :: FHrow_k(fsdim)
integer(kind=iwp) :: a, ab, b, c, cd, d, iAtom
real(kind=wp) :: d_ac, d_ad, d_bc, d_bd, diffnorm, Q_aa, Q_ab, Q_ac, Q_ad, Q_bb, Q_bc, Q_bd, Q_cc, Q_dd
real(kind=wp), allocatable :: NumHess(:,:)
logical(kind=iwp), parameter :: DebugHess = .false.
real(kind=wp), external :: DDot_

Q_aa = Zero
Q_bb = Zero
Q_cc = Zero
Q_dd = Zero

Q_ab = Zero
Q_ac = Zero
Q_ad = Zero
Q_bd = Zero
Q_bc = Zero

d_ac = One
d_ad = One
d_bc = One
d_bd = One

diffnorm = One
FHrow_k(:) = Zero
ab = 0
cd = 0

do a=1,nOrb2Loc-1
  do b=a+1,nOrb2Loc

    cd = 0
    ab = ab+1 ! compound index rows
    if (ab == k) then
      ! get off diagonal elements
      do c=1,nOrb2Loc-1
        do d=c+1,nOrb2Loc

          ! reset Kronecker deltas
          d_ac = One
          d_ad = One
          d_bc = One
          d_bd = One

          ! evaluate Kronecker deltas
          if (a /= c) d_ac = Zero
          if (a /= d) d_ad = Zero
          if (b /= c) d_bc = Zero
          if (b /= d) d_bd = Zero

          cd = cd+1 ! compound index columns

          do iAtom=1,nAtoms
            Q_aa = PA(a,a,iAtom)
            Q_bb = PA(b,b,iAtom)
            Q_cc = PA(c,c,iAtom)
            Q_dd = PA(d,d,iAtom)

            Q_ab = PA(a,b,iAtom)
            Q_ac = PA(a,c,iAtom)
            Q_ad = PA(a,d,iAtom)
            Q_bc = PA(b,c,iAtom)
            Q_bd = PA(b,d,iAtom)

            FHrow_k(cd) = FHrow_k(cd)+Eight*Q_ab*(d_ac*Q_ad-d_ad*Q_ac-d_bc*Q_bd+d_bd*Q_bc)+Two*d_ac*Q_bd*(Two*Q_aa-Q_bb-Q_dd)- &
                          Two*d_ad*Q_bc*(Two*Q_aa-Q_bb-Q_cc)-Two*d_bc*Q_ad*(Two*Q_bb-Q_aa-Q_dd)+Two*d_bd*Q_ac*(Two*Q_bb-Q_aa-Q_cc)
          end do

          !write(u6,'(A,4(I2),4X,A,2(I4))') 'a,b,c,d=',a,b,c,d,'ab,cd=',ab,cd
        end do
      end do
    else
      cycle
    end if
  end do
end do

if (DebugHess) then

  write(u6,*)
  write(u6,*) 'In GetFHrow_PM'
  write(u6,*) '--------------'

  ! check numerically
  call mma_allocate(NumHess,fsdim,fsdim,Label='NumHess')

  call GetNumHess_PM(CMO,nOrb2Loc,nBasis,fsdim,NumHess,.false.)
  call RecPrt('Analytical','',FHrow_k,fsdim,1)
  call RecPrt('Numerical','',NumHess(k,:),fsdim,1)
  call RecPrt('Analytical-Numerical','',FHrow_k(:)-NumHess(k,:),fsdim,1)
  diffnorm = sqrt(DDot_(fsdim*2,NumHess(k,:)-FHrow_k,1,NumHess(k,:)-FHrow_k,1))/real(fsdim,kind=wp)**2
  if (diffnorm > 1.0e-4_wp) then
    write(u6,*) 'ERROR: Numerical and analytical FHrow_k deviate too much, diffnorm=',diffnorm
    call Abend()
  else
    write(u6,*) 'Numerical and analytical FHrow_k agree: diffnorm= ',diffnorm
  end if
  call mma_deallocate(NumHess)
  write(u6,*)
end if

end subroutine GetFHrow_PM
