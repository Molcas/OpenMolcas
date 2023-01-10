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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine Drvetc(nGrad)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradients of one-electron matrices.*
!                                                                      *
!             Written by Anders Bernhardsson for electric field        *
!             Gradients                                                *
!             October  97                                              *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use MckDat, only: sLength
use Basis_Info, only: dbsc, nBas, nCnttp
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
integer(kind=iwp) :: iCar, iCnt, iCnttp, idCar, idcnt, idum, iIrrep, iopt, irc, isym, loper, nDens
real(kind=wp) :: Ccoor(3)
character(len=8) :: Lbl
real(kind=wp), allocatable :: D0(:), EG(:), Temp(:)
integer(kind=iwp), external :: irrfnc
external :: ElGrd, ElGrddot, ElMem

Ccoor(:) = Zero

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nTri_Elem(nBas(iIrrep))
end do

call mma_Allocate(D0,nDens,Label='D0')
call Get_D1ao_Var(D0,nDens)

call mma_allocate(EG,3*nGrad,Label='EG')
call Dot1El2(ElGrddot,ElMem,EG,3*nGrad,.true.,CCoor,D0,1)
call mma_deallocate(D0)

EG(:) = -EG(:)

call mma_allocate(Temp,3*nGrad,Label='Temp')
Temp(:) = Zero

call Drvel1(Temp)

EG(:) = EG(:)+Temp(:)

Lbl = 'NUCELGR'
idum = 1
iopt = ibset(0,sLength)
irc = 3*nGrad
call dWrMCk(irc,iopt,LBL,idum,Temp,idum)
if (irc /= 0) call SysAbendMsg('drvect','error during write in dwrmck',' ')

idum = 1
iopt = ibset(0,sLength)
irc = 3*nGrad
Lbl = 'DOTELGR'
call dWrMCk(irc,iopt,LBL,idum,EG,idum)
if (irc /= 0) call SysAbendMsg('drvect','error during write in dwrmck',' ')
call mma_deallocate(EG)
call mma_deallocate(Temp)

! needed in RASSI
loper = 0
do iCar=1,3

  isym = irrfnc(2**(icar-1)) ! nropr(ichbas(1+iCar))
  write(Lbl,'(a,i2)') 'ELEC ',iCar
  idcnt = 0
  do iCnttp=1,nCnttp
    do iCnt=1,dbsc(iCnttp)%nCntr
      idcnt = idcnt+1
      do idCar=1,3
        call Cnt1El2(ELGRD,ELMEM,Lbl,idcnt,idcar,loper,One,.true.,lbl,0,isym,icar,1)
      end do
    end do
  end do

end do

return

end subroutine Drvetc
