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

subroutine EPotPoint(iPotte,nPick,ipPick,ipDPick,nEPP,ipT,ipTi,NucNr,nB,iAtom,jAtom,ip_Center)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPotte, nPick, ipPick, ipDPick, nEPP, ipT, ipTi, NucNr, nB, iAtom, jAtom, ip_Center
integer(kind=iwp) :: iB1, iB2, iC1, iC2, iComp, iDSq, iDTrans, iOpt, iPo, iPoint, iPP, iPSq, iPTr, irc, iSmLbl, iTEMP, kaunter, &
                     nB2, nB22, nDens
real(kind=wp) :: dEx
character(len=10) :: Label
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: D1ao(:)
#include "WrkSpc.fh"

! Loop through all points, pick out the relevant ones, obtain the
! partial expectation value from the appropriate basis and return.

nB2 = nB*(nB+1)/2
nB22 = nB**2
call GetMem('DSq','Allo','Real',iDSq,nB22)
call Qpg_darray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens,Label='D1ao')
else
  write(u6,*) 'EPotPoint: D1ao not found.'
  call abend()
end if
call Get_D1ao(D1ao,nDens)
call Dsq(D1ao,Work(iDSq),1,nB,nB)
call mma_deallocate(D1ao)
call GetMem('TEMP','Allo','Real',iTEMP,nB22)
call GetMem('DTrans','Allo','Real',iDTrans,nB22)

! Contravariant transformation of density matrix.

call DGEMM_('N','N',nB,nB,nB,One,Work(ipTi),nB,Work(iDSq),nB,Zero,Work(iTEMP),nB)
call DGEMM_('N','T',nB,nB,nB,One,Work(iTEMP),nB,Work(ipTi),nB,Zero,Work(iDTrans),nB)
call GetMem('Points','Allo','Real',iPP,nB2+4)
call GetMem('PointsSq','Allo','Real',iPSq,nB22)
call GetMem('PointsTr','Allo','Real',iPTr,nB22)
do iPoint=1,nPick
  iPo = iWork(ipPick+iPoint-1)
  write(Label,'(A3,I5)') 'EF0',iPo
  irc = -1
  iOpt = 0
  iSmLbl = 0
  iComp = 1
  call RdOne(irc,iOpt,Label,iComp,Work(iPP),iSmLbl)
  call Square(Work(iPP),Work(iPSq),1,nB,nB)

  ! Covariant transformation of the matrix for the potential in this particular point.

  call DGEMM_('T','N',nB,nB,nB,One,Work(ipT),nB,Work(iPSq),nB,Zero,Work(iTEMP),nB)
  call DGEMM_('N','N',nB,nB,nB,One,Work(iTEMP),nB,Work(ipT),nB,Zero,Work(iPTr),nB)
  dEx = Zero
  kaunter = 0

  ! The usual stuff to get the localized value.

  do iB1=1,nB
    do iB2=1,nB
      iC1 = iWork(ip_Center+iB1-1)
      iC2 = iWork(ip_Center+iB2-1)
      if (((iC1 == iAtom) .and. (iC2 == jAtom)) .or. ((iC1 == jAtom) .and. (iC2 == iAtom))) then
        dEx = dEx+Work(iDTrans+kaunter)*Work(iPTr+kaunter)
      end if
      kaunter = kaunter+1
    end do
  end do

  ! Accumulate in the return vector.

  if (iAtom == jAtom) then
    Work(iPotte+iPoint-1) = -dEx+real(NucNr,kind=wp)/Work(ipDPick+iPoint-1)
  else
    Work(iPotte+iPoint-1) = -dEx
  end if
end do

! Deallocate.

call GetMem('DSq','Free','Real',iDSq,nB22)
call GetMem('TEMP','Free','Real',iTEMP,nB22)
call GetMem('DTrans','Free','Real',iDTrans,nB22)
call GetMem('Points','Free','Real',iPP,nB2+4)
call GetMem('PointsSq','Free','Real',iPSq,nB22)
call GetMem('PointsTr','Free','Real',iPTr,nB22)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nEPP)

end subroutine EPotPoint
