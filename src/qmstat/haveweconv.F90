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

subroutine HaveWeConv(iCNum,iCStart,iQ_Atoms,Indma,DT,FFp,xyzMyI,Egun,Energy,NVarv,JaNej,Haveri)

use qmstat_global, only: Cordst, Enelim, iPrint, itMax, nCent, nPart, nPol, Pol, PolLim
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCNum, iCStart, iQ_Atoms, Indma, NVarv
real(kind=wp), intent(inout) :: DT(3,nPol*nPart), Egun
real(kind=wp), intent(in) :: FFp(nPol*nPart,3), Energy
real(kind=wp), intent(out) :: xyzMyI(3)
logical(kind=iwp), intent(out) :: JaNej, Haveri
integer(kind=iwp) :: i, imin, j, k, kmin, l
real(kind=wp) :: Diff, Diffab, dist, distmin, Dtil, Egtest

!----------------------------------------------------------------------*
! With the new and the old induced dipoles, check if we have converged.*
! We also have energy check.                                           *
!----------------------------------------------------------------------*
JaNej = .true.
Haveri = .false.
Diffab = Zero
xyzMyi(:) = Zero
do i=1+(nPol*iCnum),IndMa
  k = i-((i-1)/nPol)*nPol
  do l=1,3
    Dtil = FFp(i,l)*Pol(k)
    Diff = abs(DT(l,i)-Dtil)
    if (Diff > Diffab) Diffab = Diff
    xyzMyi(l) = xyzMyi(l)+Dtil
    ! This is the quantity that has
    ! changed during the iteration and that through FFp
    ! includes the effect of the polarization of the
    ! QM-molecule. It enters the iteration above, unless
    ! we have converged.
    DT(l,i) = Dtil
  end do
end do
Egtest = Egun-Energy
Egun = Energy
if (nVarv >= itMax) then !itMax is from input or default.
  write(u6,*)
  write(u6,*) '  No convergence for the induced dipoles.'
  write(u6,*) '  Difference remaining after ',nVarv,' iterations: ',Diffab
  Haveri = .true.
  iPrint = 10
  do j=icstart,nPart*nCent,nCent
    distmin = 1.0e4_wp
    kmin = 0
    imin = 0
    do i=1,iq_atoms
      do k=0,nCent-1
        dist = sqrt((Cordst(1,i)-Cordst(1,j+k))**2+(Cordst(2,i)-Cordst(2,j+k))**2+(Cordst(3,i)-Cordst(3,j+k))**2)
        if (dist < distmin) then
          distmin = dist
          imin = i
          kmin = k
        end if
      end do
    end do
    write(u6,*) 'solv.',j,'iq_atom',imin,'center',kmin+1,'dist',distmin
  end do
  write(u6,*)
else
  if (abs(egtest) > Enelim) JaNej = .false.
  if (Diffab > PolLim) JaNej = .false.
end if

return

end subroutine HaveWeConv
