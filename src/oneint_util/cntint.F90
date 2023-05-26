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
! Copyright (C) 1991,2008, Roland Lindh                                *
!***********************************************************************

subroutine CntInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: to compute contact integrals.                                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, February '91                 *
!             Modified from D1Int January 2008.                        *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: ia, ib, iIC, ipArr, ipAxyz, ipBxyz, iPrint, iRout, na, nb, nip
character(len=80) :: Label

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(nOrdOp)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 150
iPrint = nPrint(iRout)

rFinal(:,:,:,:) = Zero

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1)
ipArr = nip
na = (la+1)*(la+2)/2
nb = (lb+1)*(lb+2)/2
nip = nip+nZeta*na*nb
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'CntInt: nip-1 > nArr*nZeta')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nArr,nZeta=',nArr,nZeta
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In CntInt: A',' ',A,1,3)
  call RecPrt(' In CntInt: RB',' ',RB,1,3)
  call RecPrt(' In CntInt: CoorO',' ',CoorO,1,3)
  call RecPrt(' In CntInt: P',' ',P,nZeta,3)
  write(u6,*) ' In CntInt: la,lb=',la,lb
end if

! Compute the contact terms.

call Contact(Zeta,P,nZeta,A,Array(ipAxyz),la,RB,Array(ipBxyz),lb,CoorO,lOper,iCho,nIC,Array(ipArr),rFinal,iStabM,nStabM,nComp, &
             rKappa)

if (iPrint >= 99) then
  do iIC=1,nIC
    do ia=1,nTri_Elem1(la)
      do ib=1,nTri_Elem1(lb)
        write(Label,'(A,I2,A,I2,A)') 'Contact term(',ia,',',ib,')'
        call RecPrt(Label,' ',rFinal(:,ia,ib,iIC),1,nZeta)
      end do
    end do
  end do
end if

return

end subroutine CntInt
