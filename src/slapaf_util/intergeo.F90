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

subroutine intergeo(FileName,Enrg,Crd,Grd,nAtm,nIter)
!---------------------------------*
! Add geometry optimization info  *
!   to the Molden inputfile       *
!---------------------------------*

use Symmetry_Info, only: nIrrep
use Phase_Info
use Slapaf_Info, only: Cx, nStab

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
#include "angstr.fh"
#include "periodic_table.fh"
real*8 Crd(3,nAtm,nIter), Enrg(nIter), Grd(3,nAtm,nIter)
integer, allocatable :: icoset2(:,:,:), jStab2(:,:), nStab2(:)
character*(*) FileName
real*8, allocatable :: Cx_p(:,:), Charge(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Pick information for centers and pseudo centers

call Get_iScalar('Unique atoms',msAtom)
call Get_iScalar('Pseudo atoms',msAtom_p)

call mma_Allocate(Charge,msAtom+msAtom_p,Label='Charge')
call Get_dArray('Nuclear charge',Charge,msAtom)

if (msAtom_p > 0) then
  call Get_dArray('Pseudo charge',Charge(msAtom+1),msAtom_p)
  call mma_allocate(Cx_p,3,msAtom_p,Label='Cx_p')
  call Get_dArray('Pseudo Coordinates',Cx_p,3*msAtom_p)
else
  call mma_allocate(Cx_p,3,1,Label='Cx_p')
  Cx_p(:,:) = 0.0d0
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (msAtom > nAtm) then
  call WarningMessage(2,'Error in InterGEO')
  write(6,*) 'msAtom > nAtm'
  write(6,*) 'msAtom=',msAtom
  write(6,*) 'nAtm=',nAtm
  call Abend()
end if

Lu_Molden = 19
call molcas_open(Lu_Molden,FileName)
write(Lu_Molden,*) '[Molden Format]'
write(Lu_Molden,*) '[N_GEO]'
write(Lu_Molden,*) nIter
write(Lu_Molden,*) '[GEOCONV]'
write(Lu_Molden,*) 'energy'

iEner = 0
do iIter=1,nIter
  write(Lu_Molden,'(E24.17)') Enrg(iIter)
  iEner = iEner+1
end do

write(Lu_Molden,*) 'max-force'
iGx = 0
do iIter=1,nIter
  grmax = 0.0d0
  do ndc=1,msAtom
    grx = abs(Grd(1,ndc,iIter))
    gry = abs(Grd(2,ndc,iIter))
    grz = abs(Grd(3,ndc,iIter))
    if (grx > grmax) grmax = grx
    if (gry > grmax) grmax = gry
    if (grz > grmax) grmax = grz
    iGx = iGx+3
  end do
  write(Lu_Molden,'(F12.7)') grmax
end do

write(Lu_Molden,*) 'rms-force'
iGx = 0
do iIter=1,nIter
  grtot = 0.0d0
  ngrad = 0
  do ndc=1,msAtom
    grx = Grd(1,ndc,iIter)
    gry = Grd(2,ndc,iIter)
    grz = Grd(3,ndc,iIter)
    do i=0,nIrrep/nStab(ndc)-1
      grtot = grtot+grx*grx+gry*gry+grz*grz
      ngrad = ngrad+1
    end do
    iGx = iGx+3
  end do
  write(Lu_Molden,'(F12.7)') sqrt(grtot)/dble(ngrad)
end do

! This part disabled because gv refuses to open the file
#ifdef write_molden_steps
write(Lu_Molden,*) 'max-step'
do iIter=1,nIter-1
  stepmax = 0.0d0
  do ndc=1,msAtom
    dx = Crd(1,ndc,iIter+1)-Crd(1,ndc,iIter)
    dy = Crd(2,ndc,iIter+1)-Crd(2,ndc,iIter)
    dz = Crd(3,ndc,iIter+1)-Crd(3,ndc,iIter)
    if (dx > stepmax) stepmax = dx
    if (dy > stepmax) stepmax = dy
    if (dz > stepmax) stepmax = dz
  end do
  write(Lu_Molden,'(F12.7)') stepmax
end do

write(Lu_Molden,*) 'rms-step'
do iIter=1,nIter-1
  step = 0.0d0
  do ndc=1,msAtom
    dx = Crd(1,ndc,iIter+1)-Crd(1,ndc,iIter)
    dy = Crd(2,ndc,iIter+1)-Crd(2,ndc,iIter)
    dz = Crd(3,ndc,iIter+1)-Crd(3,ndc,iIter)
    do i=0,nIrrep/nStab(ndc)-1
      step = step+dx*dx+dy*dy+dz*dz
    end do
  end do
  write(Lu_Molden,'(F12.7)') sqrt(step)/dble(ngrad)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up the desymmetrization of the coordinates

ixyz = 0
ixyz_p = 0
MaxDCR = 0
call mma_allocate(icoset2,[0,7],[0,7],[1,msAtom+msAtom_p],label='icoset2')
call mma_allocate(jStab2,[0,7],[1,msAtom+msAtom_p],label='jStab2')
call mma_allocate(nStab2,[1,msAtom+msAtom_p],label='nStab2')
do ndc=1,msAtom+msAtom_p
  if (ndc <= msAtom) then
    ixyz = ixyz+1
    iChxyz = iChAtm(Cx(1,ixyz,1))
  else
    ixyz_p = ixyz_p+1
    iChxyz = iChAtm(Cx_p(1,ixyz_p))
  end if
  call Stblz(iChxyz,nStab2(ndc),jStab2(0,ndc),MaxDCR,iCoSet2(0,0,ndc))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
nAt = 0
write(Lu_Molden,*) '[GEOMETRIES] (XYZ)'
do ndc=1,msAtom+msAtom_p
  do i=0,nIrrep/nStab2(ndc)-1
    nAt = nAt+1
  end do
end do

do iIter=1,nIter
  write(Lu_Molden,'(I4)') nAt
  write(Lu_Molden,*) Enrg(iIter)
  do ndc=1,msAtom+msAtom_p
    if (ndc <= msAtom) then
      x = Crd(1,ndc,iIter)
      y = Crd(2,ndc,iIter)
      z = Crd(3,ndc,iIter)
    else
      x = Cx_p(1,ndc-msAtom)
      y = Cx_p(2,ndc-msAtom)
      z = Cx_p(3,ndc-msAtom)
    end if
    do i=0,nIrrep/nStab2(ndc)-1
      iFacx = iPhase(1,icoset2(i,0,ndc))
      iFacy = iPhase(2,icoset2(i,0,ndc))
      iFacz = iPhase(3,icoset2(i,0,ndc))
      x1 = angstr*x*dble(iFacx)
      y1 = angstr*y*dble(iFacy)
      z1 = angstr*z*dble(iFacz)
      LbAtom = int(charge(ndc))
      write(Lu_Molden,102) pTab(LbAtom),x1,y1,z1
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
write(Lu_Molden,*) '[FORCES]'

do iIter=1,nIter
  write(Lu_Molden,'(A,1X,I4)') 'point',iIter
  write(Lu_Molden,'(I4)') nAt
  do ndc=1,msAtom+msAtom_p
    if (ndc <= msAtom) then
      x = Grd(1,ndc,iIter)
      y = Grd(2,ndc,iIter)
      z = Grd(3,ndc,iIter)
    else
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
    end if
    do i=0,nIrrep/nStab2(ndc)-1
      iFacx = iPhase(1,icoset2(i,0,ndc))
      iFacy = iPhase(2,icoset2(i,0,ndc))
      iFacz = iPhase(3,icoset2(i,0,ndc))
      x1 = x/angstr*dble(iFacx)
      y1 = y/angstr*dble(iFacy)
      z1 = z/angstr*dble(iFacz)
      write(Lu_Molden,103) x1,y1,z1
    end do
  end do
end do
call mma_deallocate(nStab2)
call mma_deallocate(jStab2)
call mma_deallocate(icoset2)

close(Lu_Molden)
!                                                                      *
!***********************************************************************
!                                                                      *
if (allocated(Cx_p)) call mma_deallocate(Cx_p)
call mma_deallocate(Charge)
!                                                                      *
!***********************************************************************
!                                                                      *
return

102 format(A2,3(3x,F12.7))
103 format(3(3x,F12.7))

end subroutine intergeo
