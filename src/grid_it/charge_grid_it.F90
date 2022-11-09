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
! Copyright (C) 2011, Francesco Aquilante                              *
!***********************************************************************

subroutine Charge_GRID_IT(nSym,nBas,CMO,nCMO,OCCN,iDoIt,long_prt)

!***********************************************************************
!
!  Author : F. Aquilante
!
!
!   Purpose: Compute Mulliken charges for each MO separately.
!            The analysis is performed ONLY for the occupied MOs
!            specified in GRID_IT input (this info is stored in iDoIt)
!
!   Note:  this functionality was requested by some Turbomole users
!          recently converted to MolCas. Its scientific value is not
!          too high in my opinion, therefore this subroutine is simply
!          a hack of the existing CHARGE_ and thus infinitely far from
!          efficient coding.
!
!                                            Toulouse, 28 Nov 2011
!
!***********************************************************************

use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nCMO, iDoIt(*)
real(kind=wp), intent(in) :: CMO(nCMO), OCCN(*)
logical(kind=iwp), intent(in) :: long_prt
#include "Molcas.fh"
integer(kind=iwp) :: iCase, iComp, iOpt, iOrb, iRc, iSyLbl, iSym, jOcc, MxTyp, nNUC, nTot1
character(len=8) :: Label
real(kind=wp), allocatable :: QQ(:), S(:), Xocc(:)
character(len=LenIn8), allocatable :: UBName(:)

MxTyp = 0
nTot1 = 0
do iSym=1,nSym
  MxTyp = MxTyp+nBas(iSym)
  nTot1 = nTot1+nBas(iSym)*(nBas(iSym)+1)/2
end do
call mma_allocate(UBName,MxTyp,label='UBName')
call Get_cArray('Unique Basis Names',UBName,LenIn8*MxTyp)
call mma_allocate(Xocc,MxTyp,label='XOCC')
call Get_iScalar('Unique atoms',nNUC)
call mma_allocate(QQ,MxTyp*nNuc,label='QQ')
call mma_allocate(S,nTot1,label='Ovrlp')
iRc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(iRc,iOpt,Label,iComp,S,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'charge_grid_it: iRc from Call RdOne not 0'
  write(u6,*) 'Label = ',Label
  write(u6,*) 'iRc = ',iRc
  call Abend()
end if
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,'(A)') '         **************************'
call CollapseOutput(1,'       Charges per occupied MO ')
write(u6,'(A)') '         **************************'
write(u6,*)
write(u6,*)
write(u6,*)

Xocc(:) = Zero

iCase = 2
jOcc = 1
do iSym=1,nSym
  do iOrb=1,nBas(iSym)

    if ((IdoIt(jOcc) == 1) .and. (OCCN(jOcc) > Zero)) then

      write(u6,'(A,I4,A,I1,A,F6.4)') '          MO:',iOrb,'      Symm.: ',iSym,'      Occ. No.: ',OCCN(jOcc)

      Xocc(jOcc) = OCCN(jOcc)

      QQ(:) = Zero
      call One_CHARGE(NSYM,NBAS,UBName,CMO,Xocc,S,iCase,long_prt,MxTyp,QQ,nNuc)
      Xocc(jOcc) = Zero
    end if

    jOcc = jOcc+1
  end do
end do

call mma_deallocate(Xocc)
call mma_deallocate(QQ)
call mma_deallocate(S)
call mma_deallocate(UBName)

return

end subroutine Charge_GRID_IT
