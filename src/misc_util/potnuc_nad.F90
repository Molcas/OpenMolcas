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

subroutine PotNuc_nad(nSym,nAtoms,ReCharge,ZRE_nad)
!***********************************************************************
!                                                                      *
!     purpose: Computes NAD part of the Nuclear Repulsion Energy.      *
!              An array of reference charges (ReCharge) is used to     *
!              identify which atoms were alternating as ghosts         *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "angstr.fh"
integer nSym, nAtoms
real*8 ReCharge(nAtoms)
integer iGen(3), iCoSet(0:7,0:7), iStab(0:7), iOper(0:7)
real*8, allocatable :: Charge(:), Coor(:,:)

!----------------------------------------------------------------------*
! Prologue                                                             *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Read symm. oper per symm. species                                    *
!----------------------------------------------------------------------*
call Get_iArray('Symmetry operations',iOper,nSym)
!----------------------------------------------------------------------*
! Read atom Charges                                                    *
!----------------------------------------------------------------------*
call mma_allocate(Charge,8*nAtoms,Label='Charge')
call Get_dArray('Effective nuclear Charge',Charge,nAtoms)
!----------------------------------------------------------------------*
! Read coordinates of atoms                                            *
!----------------------------------------------------------------------*
call mma_allocate(Coor,3,8*nAtoms,Label='Coor')
call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
!----------------------------------------------------------------------*
! Apply the symmetry operations                                        *
!----------------------------------------------------------------------*
nGen = 0
if (nSym == 2) nGen = 1
if (nSym == 4) nGen = 2
if (nSym == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen >= 3) iGen(3) = iOper(4)

iAll_Atom = 0
MaxDCR = 0
iAll_Atom = nAtoms
do iAtom=1,nAtoms
  iChAtom = iChxyz(Coor(1:3,iAtom),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nSym/nStab
  Charge_ = Charge(iAtom)

  do iCo=1,nCoSet-1

    iAll_Atom = iAll_Atom+1
    Charge(iAll_Atom) = Charge_

    call OA(iCoSet(iCo,0),Coor(1:3,iAtom),Coor(1:3,iAll_Atom))

  end do

end do

!----------------------------------------------------------------------*
! Compute NAD part of the nuclear repulsion energy                     *
!----------------------------------------------------------------------*
ZRE_nad = 0.0d0

if (ReCharge(1) > 0.0d0) then

  do jAt=0,iAll_Atom-1
    pCharge = Charge(jAt+1)
    if (pCharge > 0.0d0) then
      do iAt=0,jAt-1  ! loop downwards
        kAt = mod(iAt+1,nAtoms)
        if (kAt == 0) kAt = nAtoms
        qCharge = ReCharge(kAt)
        if (qCharge > 0.0d0) then
          Xpq = Coor(1,iAt+1)-Coor(1,jAt+1)
          Ypq = Coor(2,iAt+1)-Coor(2,jAt+1)
          Zpq = Coor(3,iAt+1)-Coor(3,jAt+1)
          Rpq = sqrt(Xpq**2+Ypq**2+Zpq**2)
          pq_rep = pCharge*qCharge/Rpq
          ZRE_nad = ZRE_nad+pq_rep
        end if
      end do
    end if
  end do

else

  do jAt=0,iAll_Atom-1
    pCharge = Charge(jAt+1)
    if (pCharge > 0.0d0) then
      do iAt=jAt+1,iAll_Atom-1   ! loop upwards
        kAt = mod(iAt+1,nAtoms)
        if (kAt == 0) kAt = nAtoms
        qCharge = ReCharge(kAt)
        if (qCharge > 0.0d0) then
          Xpq = Coor(1,iAt+1)-Coor(1,jAt+1)
          Ypq = Coor(2,iAt+1)-Coor(2,jAt+1)
          Zpq = Coor(3,iAt+1)-Coor(3,jAt+1)
          Rpq = sqrt(Xpq**2+Ypq**2+Zpq**2)
          pq_rep = pCharge*qCharge/Rpq
          ZRE_nad = ZRE_nad+pq_rep
        end if
      end do
    end if
  end do

end if

#ifdef _DEBUGPRINT_
!----------------------------------------------------------------------*
! Print coordinates of the system  / ZRE_nad energy                    *
!----------------------------------------------------------------------*
write(6,*)
write(6,'(6X,A)') 'Atoms cartesian coordinates in Angstrom:'
write(6,'(6X,A)') '-----------------------------------------------'
write(6,'(6X,A)') 'No.  Charge A/B      X         Y         Z     '
write(6,'(6X,A)') '-----------------------------------------------'
do iAt=0,iAll_Atom-1
  kAt = mod(iAt+1,nAtoms)
  write(6,'(4X,I4,2X,F4.0,1X,F4.0,2X,3F10.5)') iAt+1,Charge(1+iAt),Recharge(kAt),Angstr*Coor(1:3,iAt+1)
end do
write(6,'(6X,A)') '-----------------------------------------------'
write(6,'(6X,A,F12.6)') 'Nuclear repulsion energy (NAD) =',ZRE_nad
write(6,*)
write(6,*)
#endif
!----------------------------------------------------------------------*
! Normal exit                                                          *
!----------------------------------------------------------------------*
call mma_deallocate(Coor)
call mma_deallocate(Charge)

return

end subroutine PotNuc_nad
