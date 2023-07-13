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

subroutine Dissoc(xyz,nCntr,mCntr,rMss,Dist,B,lWrite,Label,dB,ldB)
!***********************************************************************
!                                                                      *
!     Object: To evaluate the B matrix elements of an internal         *
!             coordinate corresponding to a dissociation of two        *
!             parts of a molecule.                                     *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCntr, mCntr
real(kind=wp), intent(in) :: xyz(3,nCntr+mCntr), rMss(nCntr+mCntr)
real(kind=wp), intent(out) :: Dist, B(3,nCntr+mCntr)
logical(kind=iwp), intent(in) :: lWrite, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dB(3,nCntr+mCntr,3,nCntr+mCntr)
integer(kind=iwp) :: i, iCntr, ix, j, jCntr, jx
real(kind=wp) :: dRdri, dRdrj, Fact, Facti, Factj, R(3,2), RM(2), Sgn, Signi, Signj

RM(:) = Zero
R(:,:) = Zero

#ifdef _DEBUGPRINT_
write(u6,*) ' nCntr,mCntr=',nCntr,mCntr
call RecPrt(' Masses',' ',rMss,nCntr+mCntr,1)
call RecPrt(' xyz',' ',xyz,3,nCntr+mCntr)
#endif
do iCntr=1,nCntr+mCntr
  i = 1
  if (iCntr > nCntr) i = 2
  ! Sum up mass of fragment
  RM(i) = RM(i)+rMss(iCntr)
  ! Compute center of mass of the fragment
  R(:,i) = R(:,i)+rMss(iCntr)*xyz(:,iCntr)
end do
#ifdef _DEBUGPRINT_
call RecPrt('RM',' ',RM,1,2)
#endif

! Evaluate center of mass of the two parts and the distance between

R(:,1) = R(:,1)/RM(1)
R(:,2) = R(:,2)/RM(2)
Dist = Zero
do ix=1,3
  Dist = Dist+(R(ix,1)-R(ix,2))**2
end do
#ifdef _DEBUGPRINT_
call RecPrt(' Center of mass of fragments',' ',R,3,2)
#endif

Dist = sqrt(Dist)

if (lWrite) write(u6,'(1X,A,A,2(F10.6,A))') Label,' : Dissociation distance=',Dist,'/bohr',Dist*Angstrom,'/Angstrom'

! Compute the B-matrix

do iCntr=1,nCntr+mCntr
  if (iCntr <= nCntr) then
    Sgn = One
    i = 1
  else
    Sgn = -One
    i = 2
  end if
  do ix=1,3
    if (xyz(ix,iCntr) /= Zero) then
      Fact = Sgn*rMss(iCntr)/RM(i)
    else
      Fact = Zero
    end if
    B(ix,iCntr) = Fact*(R(ix,1)-R(ix,2))/Dist
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('B',' ',B,3,nCntr+mCntr)
#endif

! Compute the Cartesian derivative of the B-Matrix

if (ldB) then
  dB(:,:,:,:) = Zero
  do iCntr=1,nCntr+mCntr
    if (iCntr <= nCntr) then
      Signi = One
      i = 1
    else
      Signi = -One
      i = 2
    end if
    Facti = Signi*rMss(iCntr)/RM(i)

    do jCntr=1,nCntr+mCntr
      if (jCntr <= nCntr) then
        Signj = One
        j = 1
      else
        Signj = -One
        j = 2
      end if
      Factj = Signj*rMss(jCntr)/RM(j)

      do ix=1,3
        if (xyz(ix,iCntr) /= Zero) then
          dRdri = Facti
        else
          dRdri = Zero
        end if
        do jx=1,3
          if (xyz(jx,jCntr) /= Zero) then
            dRdrj = Factj
          else
            dRdrj = Zero
          end if

          if (ix == jx) then
            dB(ix,iCntr,jx,jCntr) = (dRdri*dRdrj-B(ix,iCntr)*B(jx,jCntr))/Dist
          else
            dB(ix,iCntr,jx,jCntr) = (-B(ix,iCntr)*B(jx,jCntr))/Dist
          end if
        end do
      end do

    end do
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('dB',' ',dB,3*(nCntr+mCntr),3*(nCntr+mCntr))
# endif
end if

return

end subroutine Dissoc
