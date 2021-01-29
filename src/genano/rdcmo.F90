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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine RdCmo()

use Genano_globals, only: kSet, nSym, nBas, kRfSet, isUHF, wSet, Ssym, Cmo, Occ, Cmo2, Occ2, Eps, lftdeg, rydgen, LenIn, BasName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), parameter :: log1e3=log(1.0e3_wp)
integer(kind=iwp) :: i, iDummy(1), iErr, indx, iOrb, irc, iSym, iSymOne, iWFtype, Lu_, Lu_One, nCmo, nDim
real(kind=wp) :: Dummy(1), eps0, eta
character(len=6) :: OneInt, NatOrb, RunFile
character(len=72) :: line

if (isUHF == 1) then
  call dCopy_(size(Cmo2),Cmo2,1,Cmo,1)
  call dCopy_(size(Occ2),Occ2,1,Occ,1)
  return
end if
!-----------------------------------------------------------------------
iSymOne = 1
kSet = kSet+1
OneInt = 'ONE'
NatOrb = 'NAT'
RunFile = 'RUN'
write(OneInt(4:6),'(i3.3)') kSet
write(NatOrb(4:6),'(i3.3)') kSet
write(RunFile(4:6),'(i3.3)') kSet
call NameRun(RunFile)
!-----------------------------------------------------------------------
write(u6,*)
write(u6,*) '--------------------------------------------------'
write(u6,*)
write(u6,'(a,i3,a,f7.3)') ' Adding density matrix',kSet,' with weight',wSet(kSet)
write(u6,*)
write(u6,*) 'Reading one-el. file: ',OneInt
Lu_One = 2
call OpnOne(irc,0,OneInt,Lu_One)
call get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
nDim = 0
nCmo = 0
do iSym=1,nSym
  nDim = nDim+nBas(iSym)
  nCmo = nCmo+nBas(iSym)**2
end do
if (allocated(BasName)) call mma_deallocate(BasName)
call mma_allocate(BasName,nDim,label='BaName')
call Get_cArray('Unique Basis Names',BasName,len(BasName)*nDim)
call ClsOne(irc,0)
write(u6,'(a,i5)') ' nSym:',nSym
write(u6,'(a,8i5)') ' nBas:',(nBas(i),i=1,nSym)
!write(u6,'(a,1x,a)') (BasName(1,i),BasName(2,i),i=1,nDim)
!-----------------------------------------------------------------------
! hack to fix Roland's inconsistent labels
do i=1,nDim
  if (BasName(i)(LenIn+3:LenIn+3) == 's') BasName(i)(LenIn+1:LenIn+3) = '01s'
  if (BasName(i)(LenIn+3:LenIn+3) == 'p') BasName(i)(LenIn+1:LenIn+3) = '02p'
  if (BasName(i)(LenIn+3:LenIn+3) == 'd') BasName(i)(LenIn+1:LenIn+3) = '03d'
  if (BasName(i)(LenIn+3:LenIn+3) == 'f') BasName(i)(LenIn+1:LenIn+3) = '04f'
  if (BasName(i)(LenIn+3:LenIn+3) == 'g') BasName(i)(LenIn+1:LenIn+3) = '05g'
  if (BasName(i)(LenIn+3:LenIn+3) == 'h') BasName(i)(LenIn+1:LenIn+3) = '06h'
  if (BasName(i)(LenIn+3:LenIn+3) == 'i') BasName(i)(LenIn+1:LenIn+3) = '07i'
  if (BasName(i)(LenIn+3:LenIn+3) == 'k') BasName(i)(LenIn+1:LenIn+3) = '08k'
end do
!write(u6,'(a,1x,a)') (BasName(1,i),BasName(2,i),i=1,nDim)
!-----------------------------------------------------------------------
if (kSet == 1) then
  call Init_GenANO()
else
  call Check_genano()
end if
!-----------------------------------------------------------------------
if (kSet == kRfSet) then
  Lu_One = 2
  call OpnOne(irc,0,OneInt,Lu_One)
  if (allocated(Cmo)) call mma_deallocate(Cmo)
  call mma_allocate(Cmo,nCmo,label='Cmo')
  call RdOne(irc,6,'Mltpl  0',1,Cmo,iSymOne)
  call CpOvlp(Cmo,Ssym)
  call ClsOne(irc,0)
end if
!-----------------------------------------------------------------------
write(u6,*)
write(u6,*) 'Reading orbital file: ',NatOrb
Lu_ = 17
call chk_vec_UHF(NatOrb,Lu_,isUHF)
if (allocated(Occ)) call mma_deallocate(Occ)
call mma_allocate(Occ,nDim,label='Occ')
if (allocated(Eps)) call mma_deallocate(Eps)
call mma_allocate(Eps,nDim,label='Eps')
Eps(:) = Zero
if (isUHF == 1) then
  if (allocated(Cmo2)) call mma_deallocate(Cmo2)
  call mma_allocate(Cmo2,nCmo,label='Cmo2')
  if (allocated(Occ2)) call mma_deallocate(Occ2)
  call mma_allocate(Occ2,nDim,label='Occ2')
  call RdVec_(NatOrb,Lu_,'CO',1,nSym,nBas,nBas,Cmo,Cmo2,Occ,Occ2,Dummy,Dummy,iDummy,line,0,iErr,iWFtype)
  write(u6,'(a)') '***'
  write(u6,'(a)') '*** rdcmo: fix reading of eps for uhf!!!'
  write(u6,'(a)') '***'
else
  call RdVec(NatOrb,Lu_,'COE',nSym,nBas,nBas,Cmo,Occ,Eps,iDummy,line,0,iErr)
end if
write(u6,*) 'Orbital set: ',trim(line)
!-----------------------------------------------------------------------
if (lftdeg) then
  indx = 1
  do iSym=1,nSym
    do iOrb=1,nBas(iSym)
      Occ(indx) = (1.001_wp/iOrb)*Occ(indx)
      indx = indx+1
    end do
  end do
end if
!-----------------------------------------------------------------------
if (rydgen) then
  call RdVec(NatOrb,Lu_,'COE',nSym,nBas,nBas,Cmo,Occ,Eps,iDummy,line,0,iErr)
  eps0 = huge(eps0)
  indx = 1
  do iSym=1,nSym
    do iOrb=1,nBas(iSym)
      if (Occ(indx) < 1.0e-2_wp) then
        eps0 = min(eps0,Eps(indx))
      end if
      indx = indx+1
    end do
  end do
  !write(u6,'(a,f12.6)') 'eps0',eps0
  indx = 1
  do iSym=1,nSym
    do iOrb=1,nBas(iSym)
      if (Occ(indx) > 1.0e-2_wp) then
        Occ(indx) = Zero
      else if (Eps(indx) < Zero) then
        eta = exp(log1e3*(Eps(indx)/eps0-One))
        !write(u6,'(a,2f12.6)') 'eps/eta',eps(indx),eta
        Occ(indx) = eta
      else
        Occ(indx) = Zero
      end if
      indx = indx+1
    end do
  end do
end if
!-----------------------------------------------------------------------
call NameRun('#Pop')

return

end subroutine RdCmo
