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

subroutine Localize_LoProp_Drv(Ttot,Ttot_Inv,nBas,iCenter,iType,nBas1,nBas2,nSym,nBasMax,ipP,Restart)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nBas1, iCenter(nBas1), iType(nBas1), nBas2, nBasMax, ipP
real(kind=wp), intent(out) :: Ttot(nBas2), Ttot_Inv(nBas2)
logical(kind=iwp), intent(in) :: Restart
integer(kind=iwp) :: idum(1), iOffs, iOfft, iOpt0, iOpt1, ip_all_ints, ip_restart, ip_SSym, ip_Tmp, ipS, ipScr, iRc, iSyLbl, iSym, &
                     nElem, nInts, nInts_tot, nScr
character(len=8) :: Label
logical(kind=iwp) :: Found
#include "WrkSpc.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Get the overlap matrix

iOpt1 = 1
iOpt0 = 0
Label = 'Mltpl  0'
iRc = -1
iSyLbl = 1
if (Restart) then
  call Qpg_iArray('LoProp nInts',Found,nElem)
  call Allocate_iWork(ip_restart,nElem)
  call Get_iArray('LoProp nInts',iWork(ip_restart),nElem)
  nInts = iWork(ip_restart+0)-4
  call GetMem('Tmp','Allo','Real',ip_SSym,nInts+4)
  call Qpg_dArray('LoProp Integrals',Found,nInts_tot)
  if (.not. Found) then
    write(u6,*) 'LoProp Integrals not available on the RunFile.'
    call Abend()
  end if
  call Allocate_Work(ip_all_ints,nInts_Tot)
  call Get_dArray('LoProp Integrals',Work(ip_all_ints),nInts_tot)
  call dCopy_(iWork(ip_restart+0),Work(ip_all_ints),1,Work(ip_SSym),1)
  call Get_iArray('LoProp iSyLbl',iWork(ip_restart),nElem)
  iSyLbl = iWork(ip_restart+0)
  call Free_Work(ip_all_ints)
  call Free_iWork(ip_restart)
else
  call iRdOne(iRc,iOpt1,Label,1,idum,iSyLbl)
  if (iRc /= 0) then
    write(u6,*) 'Polar: error reading length of mu!'
    write(u6,*) 'Mu=',0
    call Abend()
  end if
  nInts = idum(1)
  call GetMem('Tmp','Allo','Real',ip_SSym,nInts+4)
  call RdOne(iRc,iOpt0,Label,1,Work(ip_SSym),iSyLbl)
  if (iRc /= 0) then
    write(u6,*) 'Polar: error reading mu!'
    write(u6,*) 'Mu=',0
    call Abend()
  end if
end if
#ifdef _DEBUGPRINT_
call RecPrt('SSym',' ',Work(ip_SSym),1,nInts)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call GetMem('SMatrix','Allo','Real',ip_Tmp,nBas2)
if (nSym == 1) then
  ipS = ip_Tmp
else
  call GetMem('SMatrix','Allo','Real',ipS,nBas1**2)
end if

iOfft = ip_SSym
iOffs = ip_Tmp
do iSym=1,nSym
  if (nBas(iSym) == 0) Go To 99

  ! Now I square the overlap matrix because it has been stored as a
  ! lower triangle

  call Square(Work(iOfft),Work(iOffs),1,nBas(iSym),nBas(iSym))

  iOfft = iOfft+nBas(iSym)*(nBas(iSym)+1)/2
  iOffs = iOffs+nBas(iSym)**2
99 continue
end do
call Free_Work(ip_SSym)

if (nSym /= 1) then

  ! Desymmetrize

  nScr = nBasMax*nBas1
  call Allocate_Work(ipScr,nScr)
  call FZero(Work(ipS),nBas1**2)
  call Desymmetrize(Work(ip_Tmp),nBas2,Work(ipScr),nScr,Work(ipS),nBas,nBas1,Work(ipP),nSym,iSyLbl)
  call Free_Work(ipScr)
  call Free_Work(ip_Tmp)

end if
#ifdef _DEBUGPRINT_
call RecPrt('Overlap matrix',' ',Work(ipS),nBas1,nBas1)
#endif

! Localize

call Localize_LoProp(Ttot,Ttot_Inv,nBas1,Work(ipS),iCenter,iType)
#ifdef _DEBUGPRINT_
call RecPrt('Ttot',' ',Ttot,nBas1,nBas1)
call RecPrt('Ttot_Inv',' ',Ttot_Inv,nBas1,nBas1)
call xSpot('Exit Localize_LoProp_Drv')
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_Work(ipS)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Localize_LoProp_Drv
