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

subroutine AOIntegrate(iCStart,nBaseQ,nBaseC,Ax,Ay,Az,nCnC_C,iQ_Atoms,nAtomsCC,ipAOint,ipAOintpar,iV2,N,lmax,Inside)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "integral.fh"
#include "WrkSpc.fh"
#include "lenin.fh"
dimension V2(MxBasC,MxOrb_C)
dimension nCnC_C(MxBasC)
dimension Sint(MxBas,MxBasC), SintPar(MxBas,MxBasC), Rot(3,3)
dimension Inside(MxAt,3)
character Snack*30, BsLbl*(LENIN8*MxBasC)
logical PrEne, PrOcc, Inside
dimension Dummy(1)

!----------------------------------------------------------------------*
! Call Transrot. There we compute the rotation matrix for the classical*
! water under consideration. Used later.                               *
!----------------------------------------------------------------------*
call TransRot(Cordst,N+1,Rot,Dx,Dy,Dz,Ax,Ay,Az)
if (iPrint >= 17) then
  write(6,*)
  write(6,*) 'ROTATION MATRIX, Molecule ',N/nCent
  write(6,*) (Rot(1,k),k=1,3)
  write(6,*) (Rot(2,k),k=1,3)
  write(6,*) (Rot(3,k),k=1,3)
end if
!----------------------------------------------------------------------*
! Call OrbRot2. Given the rotation matrix (Rot) and the original MO-   *
! coefficients, we transform them to new MO-coefficients. V2 is on     *
! input the original MO-coefficients (stored in V3), and on output the *
! rotated.                                                             *
!----------------------------------------------------------------------*
do iOrS=1,iOrb(2) !Collect original MO-coeff.
  do iBaS=1,nBaseC
    V2(iBaS,iOrS) = V3(iBaS,iOrS)
  end do
end do
call OrbRot2(Rot,V2,iQn,iOrb(2),lMax,nCnC_C)
kaunt = 0
do iMO=1,iOrb(2) !Store the rotated in vector for later convenience.
  do iBa=1,nBaseC
    Work(iV2+kaunt) = V2(iBa,iMO)
    kaunt = kaunt+1
  end do
end do
if (iPrint >= 25) then !Optional print-out.
  PrOcc = .false.
  PrEne = .false.
  write(snack,'(A,I3)') 'Rotated orbitals for water ',N/ncent
  call GetMem('PrCMO','Allo','Real',ipPPP,nBaseC*iOrb(2))
  kauntadetta = 0
  do i=1,iOrb(2)
    do j=1,nBaseC
      Work(ipPPP+kauntadetta) = V2(j,i)
      kauntadetta = kauntadetta+1
    end do
  end do
  call NameRun('WRUNFIL')
  call Get_cArray('Unique Basis Names',BsLbl,LENIN8*nBaseC)
  call Primo(Snack,PrOcc,PrEne,Dummy(1),Dummy(1),1,[nBaseC],iOrb(2),BsLbl,Dummy,Dummy,Work(ipPPP),3)
  call GetMem('PrCMO','Free','Real',ipPPP,nBaseC*iOrb(2))
end if
do m=1,lMax !New basis function origo definied.
  x = 0
  y = 0
  z = 0
  do j=1,3
    x = x+Rot(1,j)*SavOri(j,m)
    y = y+Rot(2,j)*SavOri(j,m)
    z = z+Rot(3,j)*SavOri(j,m)
  end do
  CasOri(1,m) = x+Dx
  CasOri(2,m) = y+Dy
  CasOri(3,m) = z+Dz
end do
!----------------------------------------------------------------------*
! Compute overlap between the contracted basis functions on the water  *
! molecule presently studied and the QM-molecule.                      *
!----------------------------------------------------------------------*
do i=1,nBaseQ
  do j=1,nBaseC
    Sint(i,j) = 0
    SintPar(i,j) = 0
  end do
end do
call ContractOvl(Sint,SintPar,nBaseQ,nBaseC,N,nCent,iEl,iQ_Atoms,nAtomsCC,iPrint,Inside)
! To be able to use the fast matrix multiplication routine DGEMM_,
! we have to put the Sint (and Sintpar) matrices in vector form.
! In the future we might 'cut out the middle-man' and already
! above put the overlap matrix in vector shape.
kaunt = 0
do iC=1,nBaseC
  do iQ=1,nBaseQ
    Work(ipAOint+kaunt) = Sint(iQ,iC)
    kaunt = kaunt+1
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(iCStart)
  call Unused_integer(ipAOintpar)
end if

end subroutine AOIntegrate
