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
! Copyright (C) 2008, Bjorn O. Roos                                    *
!               2008, Valera Veryazov                                  *
!***********************************************************************

subroutine expbas(ireturn)
!***********************************************************************
!                                                                      *
!     Objective: Expand MOs to larger basis set                        *
!                                                                      *
!     B. O. Roos, University of Lund, April 2008.                      *
!                                                                      *
!***********************************************************************

use info_expbas_mod, only: EB_FileOrb
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "Molcas.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ib1, ib2, iErr, iLen, ind, ipCMO1, ipCMO2, ishift, ist1, ist2, iSym, Lu_, LuInpOrb, nb1, nb2, nDim1, &
                     nDim2, nSym1, nSym2, nTot1, nTot2, nBas1(mxsym), nBas2(mxsym), indt1(maxbfn), indt2(maxbfn), Indtype(56)

real(kind=wp) :: Occ1(maxbfn), Eorb1(maxbfn), Occ2(maxbfn), Eorb2(maxbfn)
character(len=LenIn8) :: Bas1(maxbfn), Bas2(maxbfn)
character(len=80) :: VecTit
character(len=512) :: FName
logical(kind=iwp) :: Exist_1, Exist_2, okay

!----------------------------------------------------------------------*
!     Read information from Runfile 1                                  *
!----------------------------------------------------------------------*
FName = 'RUNFIL1'
iLen = len_trim(FName)
call f_Inquire(FName(:iLen),Exist_1)
if (.not. Exist_1) then
  write(u6,*) 'Error finding file '//FName(:iLen)
  call Abend()
end if
call namerun(FName(:iLen))
call get_iScalar('nSym',nSym1)
call Get_iArray('nBas',nBas1,nSym1)
nDim1 = 0
nTot1 = 0
do iSym=1,nSym1
  nDim1 = nDim1+nBas1(iSym)
  nTot1 = nTot1+nBas1(iSym)**2
end do
call Get_cArray('Unique Basis Names',Bas1,(LENIN8)*nDim1)
!----------------------------------------------------------------------*
!     Read information from Runfile 2                                  *
!----------------------------------------------------------------------*
FName = 'RUNFIL2'
iLen = len_trim(FName)
call f_Inquire(FName(:iLen),Exist_2)
if (.not. Exist_2) then
  write(u6,*) 'Error finding file '//FName(:iLen)
  call Abend()
end if
call namerun(FName(:iLen))
call get_iScalar('nSym',nSym2)
call Get_iArray('nBas',nBas2,nSym2)
nDim2 = 0
ntot2 = 0
do iSym=1,nSym2
  nDim2 = nDim2+nBas2(iSym)
  ntot2 = ntot2+nBas2(iSym)**2
end do
call Get_cArray('Unique Basis Names',Bas2,(LENIN8)*nDim2)
!----------------------------------------------------------------------*
!     Read MO coefficients from a formatted vector file                *
!----------------------------------------------------------------------*
call GetMem('CMO1','Allo','Real',ipCMO1,nTot1)
call GetMem('CMO2','Allo','Real',ipCMO2,nTot2)
FName = EB_FileOrb
if (len_trim(FName) == 0) FName = 'INPORB'
iLen = len_trim(FName)
call f_Inquire(FName(:iLen),okay)
if (okay) then
  LuInpOrb = 50
  call RdVec(FName(:iLen),LuInpOrb,'COEI',nSym1,nBas1,nBas1,Work(ipCMO1),Occ1,Eorb1,indt1,VecTit,1,iErr)
else
  write(u6,*) 'RdCMO: Error finding MO file'
  call Abend()
end if

!----------------------------------------------------------------------*
!     Print and check input information                                *
!----------------------------------------------------------------------*
!write(u6,910) 'Start of option expand.'
write(u6,930) trim(Vectit)
write(u6,910) 'Information from input runfile'
write(u6,920) 'Number of symmetries',nSym1
write(u6,920) 'Number of basis functions',(nBas1(i),i=1,nSym1)
write(u6,910) 'Information from expanded basis set runfile'
write(u6,920) 'Number of symmetries',nSym2
write(u6,920) 'Number of basis functions',(nBas2(i),i=1,nSym2)
! Check for inconsistensies:
if (nSym1 /= nSym2) then
  write(u6,*) 'Symmetries are not equal. Stop here',nSym1,nSym2
  call Abend()
end if
do isym=1,nSym1
  if (nBas1(isym) > nBas2(isym)) then
    write(u6,*) 'Second basis set must be larger than first'
    write(u6,*) 'not fulfilled in sym',isym,'basis functions are',nBas1(isym),nBas2(isym)
    call Abend()
  end if
end do
!----------------------------------------------------------------------*
!     Build the new orbitals                                           *
!----------------------------------------------------------------------*
ist1 = 0
ist2 = 0
ib1 = 1
ib2 = 1
do isym=1,nsym1
  nb1 = nBas1(isym)
  nb2 = nBas2(isym)
  if (nb2 > 0) then
    call expandbas(Bas1(ib1),nb1,Bas2(ib2),nb2,Work(ist1+ipCMO1),Work(ist2+ipCMO2),occ1(ib1),eorb1(ib1),indt1(ib1),occ2(ib2), &
                   eorb2(ib2),indt2(ib2))
    ist1 = ist1+nb1**2
    ist2 = ist2+nb2**2
    ib1 = ib1+nb1
    ib2 = ib2+nb2
  end if
end do
!----------------------------------------------------------------------*
!     Write the new orbitals in to the file EXPORB                     *
!----------------------------------------------------------------------*
! First resort indt to standard
do i=1,56
  Indtype(i) = 0
end do
ind = 0
ishift = 0
do isym=1,nSym2
  nb2 = nBas2(isym)
  if (nb2 /= 0) then
    do ib2=1,nb2
      ind = ind+1
      Indtype(ishift+indt2(ind)) = Indtype(ishift+indt2(ind))+1
    end do
  end if
  ishift = ishift+7
end do

VecTit = 'Basis set expanded orbital file EXPORB'
Lu_ = 60
call WRVEC('EXPORB',LU_,'COEI',nSym2,nBas2,nBas2,Work(ipCMO2),occ2,eorb2,Indtype,VecTit)
write(u6,*) 'New orbitals have been built in file EXPORB'
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*
call GetMem('CMO1','Free','Real',ipCMO1,nTot2)
call GetMem('CMO2','Free','Real',ipCMO2,nTot2)

ireturn = 0

return

910 format(/1x,a)
920 format(1x,a30,8i5)
930 format(/1x,'Header on input orbitals file:'/a)

end subroutine expbas
