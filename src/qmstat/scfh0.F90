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

! Read all MO-transformed integrals and construct zeroth Fock-matrix.
subroutine ScfH0(nBas)

use qmstat_global, only: AddExt, ExtLabel, HHmat, iCompExt, iOrb, iPrint, MxSymQ, nExtAddOns, ScalExt, SupM, V1
use Index_Functions, only: iTri, nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas(MxSymQ)
#include "Molcas.fh"
#include "tratoc.fh"
integer(kind=iwp) :: i, iDisk, iExt, ij, ik, il, iLu1, iLu2, iopt, irc, iSmLbl, iSup, iToc(64), j, jk, jl, k, kaunter, kl, l, &
                     llmax, Lu_One, nBasM(MxSymQ), nBTri, nBuf1, nBuf2, nDelM(MxSymQ), nFroM(MxSymQ), nMAX, nOrbM(MxSymQ), nSize, &
                     nSymM
real(kind=wp) :: Ecor
character(len=LenIn8) :: NameM(maxbfn)
character(len=10) :: firstind
real(kind=wp), allocatable :: AOx(:), Buff(:), Fine(:,:), MOx(:), SqAO(:,:), TEMP(:,:)
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"

! Wilkommen.

write(u6,*)
write(u6,*)
write(u6,*) 'Reading MO-transformed integrals. Zeroth hamiltonian constructed.'

! Numbers and files.

nSize = nTri_Elem(iOrb(1))
call mma_allocate(SupM,nSize,nSize,label='SUPER')
iLu1 = 56
iLu2 = 58
call DaName(iLu1,'TRAONE')
call DaName(iLu2,'TRAINT')
iDisk = 0
! This is special utility to read header of TRAONE.
call Wr_Motra_Info(iLu1,2,iDisk,iToc,64,Ecor,nSymM,nBasM,nOrbM,nFroM,nDelM,MxSymQ,NameM,LenIn8*maxbfn)

! One checks.

if (nBasM(1) /= nBas(1)) then
  write(u6,*)
  write(u6,*) '  ERROR! Conflict between one-electron file and MO-transformed one-electron file.'
  write(u6,*) '         nBas=',nBas(1),' MO-nBas=',nBasM(1)
  call Quit(_RC_GENERAL_ERROR_)
end if

! Read one-electron matrix elements.

iDisk = iToc(2)
call mma_allocate(HHmat,nSize,label='HHmat')
call dDaFile(iLu1,2,HHmat,nSize,iDisk)
call DaClos(iLu1)

! Add external perturbation if requested.

if (AddExt) then
  write(u6,*) '    -- Adding external perturbation.'
  nBTri = nTri_Elem(nBas(1))
  Lu_One = IsFreeUnit(49)
  iopt = 0
  call OpnOne(irc,iopt,'ONEINT',Lu_One)
  call mma_allocate(AOx,nBTri,label='AOExt')
  call mma_allocate(TEMP,iOrb(1),nBas(1),label='TEMP')
  call mma_allocate(Fine,iOrb(1),iOrb(1),label='Final')
  call mma_allocate(SqAO,nBas(1),nBas(1),label='Squared')
  call mma_allocate(MOx,nSize,label='MOExt')
  do iExt=1,nExtAddOns
    irc = -1
    iopt = ibset(ibset(0,sNoOri),sNoNuc)
    iSmLbl = 0
    call RdOne(irc,iopt,ExtLabel(iExt),iCompExt(iExt),AOx,iSmLbl)
    AOx(:) = AOx*ScalExt(iExt)
    if (irc /= 0) then
      write(u6,*)
      write(u6,*) 'ERROR when reading ',ExtLabel(iExt),'.'
      write(u6,*) 'Have Seward computed this integral?'
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call Square(AOx,SqAO,1,nBas(1),nBas(1))
    call Dgemm_('T','N',iOrb(1),nBas(1),nBas(1),One,V1,nBas(1),SqAO,nBas(1),Zero,TEMP,iOrb(1))
    call Dgemm_('N','N',iOrb(1),iOrb(1),nBas(1),One,TEMP,iOrb(1),V1,nBas(1),Zero,Fine,iOrb(1))
    call SqToTri_Q(Fine,MOx,iOrb(1))
    HHmat(:) = HHmat+MOx
  end do
  call mma_deallocate(AOx)
  call mma_deallocate(TEMP)
  call mma_deallocate(Fine)
  call mma_deallocate(SqAO)
  call mma_deallocate(MOx)
  call ClsOne(irc,Lu_One)
end if

! Now to the two-electron matrix elements.

iDisk = 0
call iDaFile(iLu2,2,iTraToc,nTraToc,iDisk)
iDisk = iTraToc(1)
iSup = 0
! Ooohhhh, lets get crude! Read ALL (yes, you read right) integrals and then order them.
nBuf1 = nTri_Elem(nOrbM(1))
nBuf2 = nTri_Elem(nBuf1)

! Let's check if this construct is possible. If not advise user what to do.

call mma_maxDBLE(nMAX)
if (nMAX < (nBuf2+nBuf1**2)) then
  write(u6,*)
  write(u6,*) '  Too many MO-transformed two-electron integrals from Motra. Do you need all?'
  write(u6,*) '  If not, then use the DELEte keyword in Motra to remove the superfluous ones.'
  call Quit(_RC_GENERAL_ERROR_)
end if

! Proceed!

call mma_allocate(Buff,nBuf2,label='Buffer')
call mma_allocate(TEMP,nBuf1,nBuf1,label='Temporary')
call dDaFile(iLu2,2,Buff,nBuf2,iDisk)
do i=1,nBuf1
  do j=i,nBuf1
    iSup = iSup+1
    if ((i <= nSize) .and. (j <= nSize)) then
      TEMP(j,i) = Buff(iSup)
      TEMP(i,j) = Buff(iSup)
    end if
  end do
end do

! and see to that right numbers get in the right place.

do i=1,iOrb(1)
  do j=1,i
    do k=1,i
      llmax = k
      if (i == k) llmax = j
      do l=1,llmax
        ij = iTri(i,j)
        ik = iTri(i,k)
        il = iTri(i,l)
        jk = iTri(j,k)
        jl = iTri(j,l)
        kl = iTri(k,l)
        SupM(kl,ij) = TEMP(kl,ij)-(TEMP(jl,ik)+TEMP(jk,il))*Quart
        SupM(ij,kl) = SupM(kl,ij)
        SupM(jl,ik) = TEMP(jl,ik)-(TEMP(kl,ij)+TEMP(jk,il))*Quart
        SupM(ik,jl) = SupM(jl,ik)
        SupM(jk,il) = TEMP(jk,il)-(TEMP(jl,ik)+TEMP(kl,ij))*Quart
        SupM(il,jk) = SupM(jk,il)
      end do
    end do
  end do
end do
call mma_deallocate(Buff)
call mma_deallocate(TEMP)
call DaClos(iLu2)

! Serious amount of printing!

if (iPrint >= 35) then
  write(u6,*)
  write(u6,*) 'The Super Matrix in all its divine g(l)ory:'
  kaunter = 0
  do i=1,iOrb(1)
    do j=1,i
      write(firstind,'(I3,A,I3)') i,',',j
      kaunter = kaunter+1
      call TriPrt(firstind,' ',SupM(:,kaunter),iOrb(1))
    end do
  end do
  write(u6,*) 'Super Matrix End.'
end if
write(u6,*) '...Done!'

! The end.

return

end subroutine ScfH0
