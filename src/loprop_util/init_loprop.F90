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

subroutine Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,ipC,ipQ_Nuc,ip_ANr,ip_Type,ip_Center,nSize,nBas1,nBas2,nBasMax,ipP,ipPInv)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nSym, nBas(8), nOrb(8), nAtoms, ipC, ipQ_Nuc, ip_ANr, ip_Type, ip_Center, nSize, nBas1, nBas2, &
                                  nBasMax, ipP, ipPInv
real(kind=wp), intent(out) :: CoC(3)
integer(kind=iwp) :: i, iDum, ISING, iSym
real(kind=wp) :: DET
logical(kind=iwp) :: lOrb
integer(kind=iwp), parameter :: Occ = 1, Vir = 0
#include "WrkSpc.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call qpg_iarray('nOrb',lOrb,iDum)
if (lOrb) then
  call Get_iArray('nOrb',nOrb,nSym)
else
  call ICopy(nSym,nBas,1,nOrb,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nSize = 0
nBas1 = 0
nBas2 = 0
nBasMax = 0
do iSym=1,nSym
  nSize = nSize+nBas(iSym)*(nBas(iSym)+1)/2
  nBas1 = nBas1+nBas(iSym)
  nBas2 = nBas2+nBas(iSym)**2
  nBasMax = max(nBasMax,nBas(iSym))
end do
nSize = nSize+4
!                                                                      *
!***********************************************************************
!                                                                      *
! Center of Charge
call Get_dArray('Center of Charge',CoC,3)

! List coordinates of Coordinates
call Get_iScalar('LP_nCenter',nAtoms)
call Allocate_Work(ipC,3*nAtoms)
call Get_dArray('LP_Coor',Work(ipC),3*nAtoms)
call Allocate_Work(ipQ_Nuc,nAtoms)

! Effective charge at each center
call Get_dArray('LP_Q',Work(ipQ_Nuc),nAtoms)

! Atom number of each center
call Allocate_iWork(ip_ANr,nAtoms)
call Get_iArray('LP_A',iWork(ip_ANr),nAtoms)

! Pick up information of orbital type. Occ/Vir
call Allocate_iWork(ip_type,nbas1)
call Get_iArray('Orbital Type',iWork(ip_type),nBas1)
do i=ip_type,ip_type+nBas1-1
  if ((iWork(i) /= Occ) .and. (iWork(i) /= Vir)) then
    write(u6,*) 'Orbital type vector is corrupted!'
    call Abend()
  end if
end do

! Pick up index array of which center a basis function belong.
call Allocate_iWork(ip_center,nbas1)
call Get_iArray('Center Index',iWork(ip_center),nBas1)

#ifdef _DEBUGPRINT_
write(u6,*) '******* LoProp Debug Info *******'
call RecPrt('Coordinates',' ',Work(ipC),3,nAtoms)
call RecPrt('Charges',' ',Work(ipQ_Nuc),1,nAtoms)
write(u6,*) 'Atom Nr:',(iWork(i),i=ip_ANr,ip_ANr+nAtoms-1)
do iBas=1,nBas1
  if (iWork(ip_type+iBas-1) == Occ) then
    write(u6,'(A,I3,A)') 'Basis function ',iBas,' is occupied'
  else
    write(u6,'(A,I3,A)') 'Basis function ',iBas,' is virtual'
  end if
end do
write(u6,*)
do iBas=1,nBas1
  write(u6,'(A,I3,A,I3)') 'Basis function ',iBas,' belongs to center ',iWork(ip_center+iBas-1)
end do
write(u6,*) '*********************************'
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of symmetry we need the desymmetrization matrix

if (nSym == 1) Go To 99

call Allocate_Work(ipP,nbas1**2)
call Allocate_Work(ipPInv,nbas1**2)
call Get_dArray('SM',Work(ipP),nbas1**2)
#ifdef _DEBUGPRINT_
call RecPrt('SM',' ',Work(ipP),nbas1,nbas1)
#endif
call MINV(Work(ipP),Work(ipPInv),ISING,DET,nBas1)
#ifdef _DEBUGPRINT_
call RecPrt('SMInv',' ',Work(ipPInv),nbas1,nbas1)
#endif
call DGeTMi(Work(ipPInv),nbas1,nbas1)
!                                                                      *
!***********************************************************************
!                                                                      *
99 continue
return

end subroutine Init_LoProp
