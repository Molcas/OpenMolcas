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
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

subroutine Newq(q,nInter,nIter,dq,H,g,error,B,RHS,Scrt1,nScrt1,Beta,nFix,iP,Energy,Step_Trunc,Thr_RS)
!***********************************************************************
!                                                                      *
! Object: Driver for optimization procedures.                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!***********************************************************************

use Slapaf_Parameters, only: iOptC, UpMeth, Line_Search, E_Delta

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
!logical Print
real*8 q(nInter,nIter+1), dq(nInter,nIter), g(nInter,nIter), error(nInter,nIter+1), B((nIter+1)*(nIter+1)), &
       RHS(nIter+1), Scrt1(nScrt1), Energy(nIter), H(nInter,nInter)
integer iP(nIter)
character*1 Step_Trunc
real*8, allocatable :: t_q(:), t_g(:), t_dq(:)

Lu = 6
iRout = 113
iPrint = nPrint(iRout)
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(Lu,*) ' Newq: nIter,Beta=',nIter,Beta
call RecPrt(' Newq (Enter): q',' ',q,nInter,nIter+1)
call RecPrt(' Newq (Enter): dq',' ',dq,nInter,nIter)
call RecPrt(' Newq (Enter): g',' ',g,nInter,nIter)
iSave = nPrint(21)
nPrint(21) = 99
call DiagMtrx(H,nInter,ix)
nPrint(21) = iSave
#endif

E_Delta = Zero

! Perform first a linear search for the last two points
! to find minimum along the direction of the last step.

! The new point will temporarily replace the last point!

if (iPrint >= 6) write(Lu,*)
if (Line_Search) then
  if (nIter >= 2) then
    call mma_allocate(t_q,nInter,Label='t_q')
    call mma_allocate(t_g,nInter,Label='t_g')
    call mma_allocate(t_dq,nInter,Label='t_dq')

    call dcopy_(nInter,dq(1,nIter-1),1,t_dq,1)
    call dcopy_(nInter,q(1,nIter),1,t_q,1)
    call dcopy_(nInter,g(1,nIter),1,t_g,1)

    call LnSrch(Energy,q,dq,g,nInter,nIter,E_Delta)
  else
    if (iPrint >= 6) write(Lu,*) '-- First iteration no line search'
  end if
else
  if (iPrint >= 6) write(Lu,*) '-- Line search is disabled'
end if
if (iPrint >= 6) write(Lu,*)

! Invoke the quadratic optimization procedure

! iOptC=0001   : quasi Newton-Raphson
! iOptC=0010   : C1-DIIS
! iOptC=0100   : C2-DIIS
! iOptC=1000   : Rational Function

call dcopy_(nInter,[Zero],0,Scrt1,1)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iOptC == 0) then

  ! No update of the geometry

  call FZero(dq,nInter)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (iand(iOptC,1) == 1) then

  ! quasi Newton-Raphson

  UpMeth = ' qNR  '
  call QNR(nInter,nIter,dq,H,g)
  Beta_New = sqrt(DDot_(nInter,dq(1,nIter),1,dq(1,nIter),1))
  if (Beta_New > Beta) then
    call DScal_(nInter,Beta/Beta_New,dq(1,nIter),1)
    Step_Trunc = '*'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (iand(iOptC,2) == 2) then

  ! C1-DIIS

  UpMeth = 'c1DIIS'
  MinWdw = 5
  call C1DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,nFix,iP,MinWdw)
  Beta_New = sqrt(DDot_(nInter,dq(1,nIter),1,dq(1,nIter),1))
  if (Beta_New > Beta) then
    call DScal_(nInter,Beta/Beta_New,dq(1,nIter),1)
    Step_Trunc = '*'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (iand(iOptC,4) == 4) then

  ! C2-DIIS

  UpMeth = 'c2DIIS'
  call C2DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,Scrt1,nScrt1,nFix,iP)
  Beta_New = sqrt(DDot_(nInter,dq(1,nIter),1,dq(1,nIter),1))
  if (Beta_New > Beta) then
    call DScal_(nInter,Beta/Beta_New,dq(1,nIter),1)
    Step_Trunc = '*'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (iand(iOptC,8) == 8) then

  ! Rational Function Optimization (RFO)

  UpMeth = '  RF  '

  if (iand(iOptC,128) /= 128) then

    ! TS optimization

    if (iand(iOptC,512) /= 512) then

      ! Restricted Step Partitioned RFO

      call RS_P_RFO(H,g(1,nIter),nInter,dq(1,nIter),UpMeth,E_Delta,Beta,Step_Trunc)
    else

      ! Restricted Step Image RFO

      call RS_I_RFO(H,g(1,nIter),nInter,dq(1,nIter),UpMeth,E_Delta,Beta,Step_Trunc,Thr_RS)
    end if

  else

    ! Restricted Step RFO

    call RS_RFO(H,g(1,nIter),nInter,dq(1,nIter),UpMeth,E_Delta,Beta,Step_Trunc,Thr_RS)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  call WarningMessage(2,'Error in NewQ')
  write(6,*) ' Newq: Illegal setting of iOptC!'
  write(6,*) '  iOptC=',iOptC
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of a line search restore some data and add the replacements.

if (Line_Search .and. (nIter >= 2)) then
  if (iPrint >= 99) then
    call RecPrt(' Newq: q ',' ',q,nInter,nIter+1)
    call RecPrt(' Newq: dq',' ',dq,nInter,nIter)
    call RecPrt(' Newq: g ',' ',g,nInter,nIter)
  end if
  call dcopy_(nInter,q(1,nIter),1,q(1,nIter+1),1)
  call DaXpY_(nInter,One,dq(1,nIter),1,q(1,nIter+1),1)
  call dcopy_(nInter,t_q,1,q(1,nIter),1)
  call dcopy_(nInter,q(1,nIter+1),1,dq(1,nIter),1)
  call DaXpY_(nInter,-One,q(1,nIter),1,dq(1,nIter),1)

  call dcopy_(nInter,t_dq,1,dq(1,nIter-1),1)
  call dcopy_(nInter,t_g,1,g(1,nIter),1)
  if (iPrint >= 99) then
    call RecPrt(' Newq: q ',' ',q,nInter,nIter+1)
    call RecPrt(' Newq: dq',' ',dq,nInter,nIter)
    call RecPrt(' Newq: g ',' ',g,nInter,nIter)
  end if

  call mma_deallocate(t_q)
  call mma_deallocate(t_g)
  call mma_deallocate(t_dq)
end if

! Estimate energy at the relaxed geometry

if (iand(iOptC,8) /= 8) then

  ! g(r)

  call dcopy_(nInter,g(1,nIter),1,Scrt1,1)
  call DScal_(nInter,-One,Scrt1,1)

  ! 1/2 H(r) (r -r)
  !            0

  call dGeMV_('N',nInter,nInter,Half,H,nInter,dq(1,nIter),1,One,Scrt1,1)
  E_Delta = DDot_(nInter,Scrt1,1,dq(1,nIter),1)
end if

do i=1,nInter
  q(i,nIter+1) = q(i,nIter)+dq(i,nIter)
end do

#ifdef _DEBUGPRINT_
write(Lu,*) ' E_Delta=',E_Delta
call RecPrt('Newq (Exit): q',' ',q,nInter,nIter+1)
call RecPrt('Newq (Exit): dq',' ',dq,nInter,nIter)
call RecPrt('Newq (Exit): g',' ',g,nInter,nIter)
#endif

return

end subroutine Newq
