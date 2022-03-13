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

subroutine HaveWeConv(iCNum,iCStart,iQ_Atoms,Indma,iDT,FFp,xyzMyI,Egun,Energy,NVarv,JaNej,Haveri)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
dimension iDT(3)
dimension FFp(npol*npart,3), xyzMyI(3)
logical JaNej, Haveri

!----------------------------------------------------------------------*
! With the new and the old induced dipoles, check if we have converged.*
! We also have energy check.                                           *
!----------------------------------------------------------------------*
JaNej = .true.
Haveri = .false.
Diffab = 0
xyzMyi(1) = 0
xyzMyi(2) = 0
xyzMyi(3) = 0
do i=1+(nPol*iCnum),IndMa
  k = i-((i-1)/nPol)*nPol
  do l=1,3
    Dtil = FFp(i,l)*Pol(k)
    Diff = abs(Work(iDT(l)+i-1)-Dtil)
    if (Diff > Diffab) Diffab = Diff
    xyzMyi(l) = xyzMyi(l)+Dtil
    ! This is the quantity that has
    ! changed during the iteration and that through FFp
    ! includes the effect of the polarization of the
    ! QM-molecule. It enters the iteration above, unless
    ! we have converged.
    Work(iDT(l)+i-1) = Dtil
  end do
end do
Egtest = Egun-Energy
Egun = Energy
if (nVarv >= itMax) then !itMax is from input or default.
  write(6,*)
  write(6,*) '  No convergence for the induced dipoles.'
  write(6,*) '  Difference remaining after ',nVarv,' iterations: ',Diffab
  Haveri = .true.
  iPrint = 10
  do j=icstart,npart*ncent,ncent
    distmin = 1000.0d0
    kmin = 0
    imin = 0
    do i=1,iq_atoms
      do k=0,ncent-1
        dist = sqrt((cordst(i,1)-cordst(j+k,1))**2+(cordst(i,2)-cordst(j+k,2))**2+(cordst(i,3)-cordst(j+k,3))**2)
        if (dist < distmin) then
          distmin = dist
          imin = i
          kmin = k
        end if
      end do
    end do
    write(6,*) 'solv.',j,'iq_atom',imin,'center',kmin+1,'dist',distmin
  end do
  write(6,*)
else
  if (abs(egtest) > Enelim) JaNej = .false.
  if (Diffab > PolLim) JaNej = .false.
end if

return

end subroutine HaveWeConv
