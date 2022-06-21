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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine D1Int( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute the 1-electron Darwin contact term.               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, February '91                 *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
character*80 Label
! Statement function
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 150
iPrint = nPrint(iRout)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1)
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'D1Int: nip-1 > nArr*nZeta')
  write(6,*) 'nip=',nip
  write(6,*) 'nArr,nZeta=',nArr,nZeta
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In D1Int: A',' ',A,1,3)
  call RecPrt(' In D1Int: RB',' ',RB,1,3)
  call RecPrt(' In D1Int: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In D1Int: P',' ',P,nZeta,3)
  write(6,*) ' In D1Int: la,lb=',la,lb
end if

! Compute the contact terms.

call Darwin(Zeta,P,nZeta,A,Array(ipAxyz),la,RB,Array(ipBxyz),lb,final,iStabM,nStabM,nComp,rKappa)

if (iPrint >= 99) then
  do ia=1,nElem(la)
    do ib=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') 'Darwin contact(',ia,',',ib,')'
      call RecPrt(Label,' ',final(1,1,ia,ib),nZeta,nComp)
    end do
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_real_array(ZInv)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
  call Unused_integer_array(iChO)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine D1Int
