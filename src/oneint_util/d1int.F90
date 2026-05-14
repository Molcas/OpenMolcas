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

use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: ipAxyz, ipBxyz, nip
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ia, ib
character(len=80) :: Label
#endif

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(ZInv)
unused_var(CoorO)
unused_var(nOrdOp)
unused_var(lOper)
unused_var(iChO)
unused_var(PtChrg)
unused_var(iAddPot)

nip = 1
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+1)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+1)
if (nip-1 > nArr*nZeta) then
  call WarningMessage(2,'D1Int: nip-1 > nArr*nZeta')
  write(u6,*) 'nip=',nip
  write(u6,*) 'nArr,nZeta=',nArr,nZeta
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In D1Int: A',' ',A,1,3)
call RecPrt(' In D1Int: RB',' ',RB,1,3)
call RecPrt(' In D1Int: P',' ',P,nZeta,3)
write(u6,*) ' In D1Int: la,lb=',la,lb
#endif

! Compute the contact terms.

call Darwin(Zeta,P,nZeta,A,Array(ipAxyz),la,RB,Array(ipBxyz),lb,rFinal,iStabM,nStabM,nComp,rKappa)

#ifdef _DEBUGPRINT_
do ia=1,nTri_Elem1(la)
  do ib=1,nTri_Elem1(lb)
    write(Label,'(A,I2,A,I2,A)') 'Darwin contact(',ia,',',ib,')'
    call RecPrt(Label,' ',rFinal(:,ia,ib,:),nZeta,nComp)
  end do
end do
#endif

end subroutine D1Int
