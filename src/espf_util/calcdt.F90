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

subroutine CalcDT(nMult,nGrdPt,natom,nAtQM,IsMM,Coord,Grid,T,TT,DT,DTT,DTTTT,DTTT)
! The dependency of T with respect to grid points is taken
! into account if the GEPOL grid is used (iGrdTyp == 2)
!
! 1) dT    = dT/dx + Sum_k dT/dq_k * dq_k/dx
! 2) dTT   = (dTt.T+Tt.dT)
! 3) dTTTT = (dTt.T+Tt.dT)((TtT)**-1)
! 4) dTT   = d((TtT)**-1) = -((TtT)**-1)(dTt.T+Tt.dT)((TtT)**-1)
! 5) dTTT = d((TtT)**-1)*Tt + ((TtT)**-1)*dTt

use Constants, only: Zero, One, Three
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMult, nGrdPt, natom, nAtQM, IsMM(natom)
real(kind=wp), intent(in) :: Coord(3,natom), Grid(3,nGrdPt), T(nMult,nGrdPt), TT(nMult,nMult)
real(kind=wp), intent(out) :: DT(nMult,nGrdPt,3,nAtQM), DTT(nMult,nMult,3,nAtQM), DTTTT(nMult,nMult,3,nAtQM), &
                              DTTT(nMult,nGrdPt,3,nAtQM)
integer(kind=iwp) :: iAt, iMlt, iPL, iPnt, iQM, jMlt, jPnt, jQM, kMlt, nOrd
real(kind=wp) :: DG(3,3), dIJ, R, R2, R3, R5, X, Y, Z
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()

nOrd = nMult/nAtQM

! Step 1

iQM = 0
do iPnt=1,nGrdPt
  iQM = 0
  do iAt=1,natom
    if (IsMM(iAt) == 1) cycle
    iQM = iQM+1
    X = Grid(1,iPnt)-Coord(1,iAt)
    Y = Grid(2,iPnt)-Coord(2,iAt)
    Z = Grid(3,iPnt)-Coord(3,iAt)
    R = sqrt(X*X+Y*Y+Z*Z)
    R2 = R*R
    R3 = R2*R
    R5 = R2*R3
    do jQM=1,nAtQM
      DG(:,:) = Zero
      ! This is commented out since the derivative of the GEPOL grid points look ugly
      !if (iGrdTyp == 2) DG(:,:) = DGrid(iPnt,jQM,:,:)
      dIJ = Zero
      if (iQM == jQM) dIJ = One
      if (nOrd == 1) then
        DT(iQM,iPnt,1,jQM) = X/R3*(dIJ-DG(1,1))+Y/R3*(-DG(2,1))+Z/R3*(-DG(3,1))
        DT(iQM,iPnt,2,jQM) = X/R3*(-DG(1,2))+Y/R3*(dIJ-DG(2,2))+Z/R3*(-DG(3,2))
        DT(iQM,iPnt,3,jQM) = X/R3*(-DG(1,3))+Y/R3*(-DG(2,3))+Z/R3*(dIJ-DG(3,3))
      else
        DT(4*iQM-3,iPnt,1,jQM) = X/R3*(dIJ-DG(1,1))+Y/R3*(-DG(2,1))+Z/R3*(-DG(3,1))
        DT(4*iQM-3,iPnt,2,jQM) = X/R3*(-DG(1,2))+Y/R3*(dIJ-DG(2,2))+Z/R3*(-DG(3,2))
        DT(4*iQM-3,iPnt,3,jQM) = X/R3*(-DG(1,3))+Y/R3*(-DG(2,3))+Z/R3*(dIJ-DG(3,3))
        DT(4*iQM-2,iPnt,1,jQM) = (three*X*X-R2)/R5*(dIJ-DG(1,1))+(three*X*Y)/R5*(-DG(2,1))+(three*X*Z)/R5*(-DG(3,1))
        DT(4*iQM-2,iPnt,2,jQM) = (three*X*X-R2)/R5*(-DG(1,2))+(three*X*Y)/R5*(dIJ-DG(2,2))+(three*X*Z)/R5*(-DG(3,2))
        DT(4*iQM-2,iPnt,3,jQM) = (three*X*X-R2)/R5*(-DG(1,3))+(three*X*Y)/R5*(-DG(2,3))+(three*X*Z)/R5*(dIJ-DG(3,3))
        DT(4*iQM-1,iPnt,1,jQM) = (three*Y*X)/R5*(dIJ-DG(1,1))+(three*Y*Y-R2)/R5*(-DG(2,1))+(three*Y*Z)/R5*(-DG(3,1))
        DT(4*iQM-1,iPnt,2,jQM) = (three*Y*X)/R5*(-DG(1,2))+(three*Y*Y-R2)/R5*(dIJ-DG(2,2))+(three*Y*Z)/R5*(-DG(3,2))
        DT(4*iQM-1,iPnt,3,jQM) = (three*Y*X)/R5*(-DG(1,3))+(three*Y*Y-R2)/R5*(-DG(2,3))+(three*Y*Z)/R5*(dIJ-DG(3,3))
        DT(4*iQM,iPnt,1,jQM) = (three*Z*X)/R5*(dIJ-DG(1,1))+(three*Z*Y)/R5*(-DG(2,1))+(three*Z*Z-R2)/R5*(-DG(3,1))
        DT(4*iQM,iPnt,2,jQM) = (three*Z*X)/R5*(-DG(1,2))+(three*Z*Y)/R5*(dIJ-DG(2,2))+(three*Z*Z-R2)/R5*(-DG(3,2))
        DT(4*iQM,iPnt,3,jQM) = (three*Z*X)/R5*(-DG(1,3))+(three*Z*Y)/R5*(-DG(2,3))+(three*Z*Z-R2)/R5*(dIJ-DG(3,3))
      end if
    end do
  end do
end do
if (iQM /= nAtQM) then
  write(u6,'(A)') ' Error in espf/calcdt: iQM != nAtQM !!!'
  call Abend()
end if
if (iPL >= 4) then
  do iMlt=1,nMult
    do jPnt=1,nGrdPt
      write(u6,'(A,i4,i4)') ' DT for iMlt and jPnt: ',iMlt,jPnt
      do iQM=1,nAtQM
        write(u6,'(i4,3F20.7)') iQM,DT(iMlt,jPnt,:,iQM)
      end do
    end do
  end do
end if

! Step 2

DTT(:,:,:,:) = Zero
do iMlt=1,nMult
  do jMlt=1,nMult
    do iQM=1,nAtQM
      do iPnt=1,nGrdPt
        DTT(iMlt,jMlt,:,iQM) = DTT(iMlt,jMlt,:,iQM)+T(iMlt,iPnt)*DT(jMlt,iPnt,:,iQM)+DT(iMlt,iPnt,:,iQM)*T(jMlt,iPnt)
      end do
    end do
    if (iPL >= 4) then
      write(u6,'(A,i4,i4)') 'DTT1 for iMlt and jMlt: ',iMlt,jMlt
      do iQM=1,nAtQM
        write(u6,'(i4,3F20.7)') iQM,DTT(iMlt,jMlt,:,iQM)
      end do
    end if
  end do
end do

! Step 3

DTTTT(:,:,:,:) = Zero
do iMlt=1,nMult
  do jMlt=1,nMult
    do iQM=1,nAtQM
      do kMlt=1,nMult
        DTTTT(iMlt,jMlt,:,iQM) = DTTTT(iMlt,jMlt,:,iQM)+DTT(iMlt,kMlt,:,iQM)*TT(kMlt,jMlt)
      end do
    end do
    if (iPL >= 4) then
      write(u6,'(A,i4,i4)') 'DTTTT for iMlt and jMlt: ',iMlt,jMlt
      do iQM=1,nAtQM
        write(u6,'(i4,3F20.7)') iQM,DTTTT(iMlt,jMlt,:,iQM)
      end do
    end if
  end do
end do

! Step 4

DTT(:,:,:,:) = Zero
do iMlt=1,nMult
  do jMlt=1,nMult
    do iQM=1,nAtQM
      do kMlt=1,nMult
        DTT(iMlt,jMlt,:,iQM) = DTT(iMlt,jMlt,:,iQM)-TT(iMlt,kMlt)*DTTTT(kMlt,jMlt,:,iQM)
      end do
    end do
    if (iPL >= 4) then
      write(u6,'(A,i4,i4)') 'DTT for iMlt and jMlt: ',iMlt,jMlt
      do iQM=1,nAtQM
        write(u6,'(i4,3F20.7)') iQM,DTT(iMlt,jMlt,:,iQM)
      end do
    end if
  end do
end do

! Step 5

DTTT(:,:,:,:) = Zero
do iMlt=1,nMult
  do iPnt=1,nGrdPt
    do iQM=1,nAtQM
      do jMlt=1,nMult
        DTTT(iMlt,iPnt,:,iQM) = DTTT(iMlt,iPnt,:,iQM)+TT(iMlt,jMlt)*DT(jMlt,iPnt,:,iQM)+DTT(iMlt,jMlt,:,iQM)*T(jMlt,iPnt)
      end do
    end do
    if (iPL >= 4) then
      write(u6,'(A,i4,i4)') 'DTTT for iMlt and iPnt: ',iMlt,iPnt
      do iQM=1,nAtQM
        write(u6,'(i4,3F20.7)') iQM,DTTT(iMlt,iPnt,:,iQM)
      end do
    end if
  end do
end do

return

end subroutine CalcDT
