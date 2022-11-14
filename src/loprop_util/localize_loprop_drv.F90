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

subroutine Localize_LoProp_Drv(Ttot,Ttot_Inv,nBas,iCenter,iType,nBas1,nBas2,nSym,nBasMax,P,Restart)

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nBas1, iCenter(nBas1), iType(nBas1), nBas2, nBasMax
real(kind=wp), intent(out) :: Ttot(nBas1,nBas1), Ttot_Inv(nBas1,nBas1)
real(kind=wp), intent(in) :: P(*)
logical(kind=iwp), intent(in) :: Restart
integer(kind=iwp) :: iComp, idum(1), iOffs, iOfft, iOpt0, iOpt1, iRc, iSyLbl, iSym, nElem, nInts, nInts_tot, nScr
character(len=8) :: Label
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: irestart(:)
real(kind=wp), allocatable :: all_ints(:), S(:), Scr(:), SSym(:), Tmp(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Get the overlap matrix

iOpt1 = ibset(0,sOpSiz)
iOpt0 = 0
Label = 'Mltpl  0'
iRc = -1
iSyLbl = 1
if (Restart) then
  call Qpg_iArray('LoProp nInts',Found,nElem)
  call mma_allocate(irestart,nElem,label='irestart')
  call Get_iArray('LoProp nInts',irestart,nElem)
  nInts = irestart(1)-4
  call mma_allocate(SSym,nInts+4,label='Tmp')
  call Qpg_dArray('LoProp Integrals',Found,nInts_tot)
  if (.not. Found) then
    write(u6,*) 'LoProp Integrals not available on the RunFile.'
    call Abend()
  end if
  call mma_allocate(all_ints,nInts_Tot,label='all_ints')
  call Get_dArray('LoProp Integrals',all_ints,nInts_tot)
  SSym(1:nInts+4) = all_ints(1:nInts+4)
  call Get_iArray('LoProp iSyLbl',irestart,nElem)
  iSyLbl = irestart(1)
  call mma_deallocate(all_ints)
  call mma_deallocate(irestart)
else
  iComp = 1
  call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
  if (iRc /= 0) then
    write(u6,*) 'Polar: error reading length of mu!'
    write(u6,*) 'Mu=',0
    call Abend()
  end if
  nInts = idum(1)
  call mma_allocate(SSym,nInts+4,label='Tmp')
  call RdOne(iRc,iOpt0,Label,iComp,SSym,iSyLbl)
  if (iRc /= 0) then
    write(u6,*) 'Polar: error reading mu!'
    write(u6,*) 'Mu=',0
    call Abend()
  end if
end if
#ifdef _DEBUGPRINT_
call RecPrt('SSym',' ',SSym,1,nInts)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Tmp,nBas2,label='SMatrix')

iOfft = 1
iOffs = 1
do iSym=1,nSym
  if (nBas(iSym) == 0) cycle

  ! Now I square the overlap matrix because it has been stored as a
  ! lower triangle

  call Square(SSym(iOfft),Tmp(iOffs),1,nBas(iSym),nBas(iSym))

  iOfft = iOfft+nBas(iSym)*(nBas(iSym)+1)/2
  iOffs = iOffs+nBas(iSym)**2
end do
call mma_deallocate(SSym)

if (nSym == 1) then
  call move_alloc(Tmp,S)
else

  ! Desymmetrize

  call mma_allocate(S,nBas1**2)
  nScr = nBasMax*nBas1
  call mma_allocate(Scr,nScr,label='Scr')
  S(:) = Zero
  call Desymmetrize(Tmp,nBas2,Scr,nScr,S,nBas,nBas1,P,nSym,iSyLbl)
  call mma_deallocate(Scr)
  call mma_deallocate(Tmp)

end if
#ifdef _DEBUGPRINT_
call RecPrt('Overlap matrix',' ',S,nBas1,nBas1)
#endif

! Localize

call Localize_LoProp(Ttot,Ttot_Inv,nBas1,S,iCenter,iType)
#ifdef _DEBUGPRINT_
call RecPrt('Ttot',' ',Ttot,nBas1,nBas1)
call RecPrt('Ttot_Inv',' ',Ttot_Inv,nBas1,nBas1)
write(u6,*)
write(u6,*) 'Exit Localize_LoProp_Drv'
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(S)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Localize_LoProp_Drv
