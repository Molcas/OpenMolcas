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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "qm1.fh"
#include "numbers.fh"
#include "tratoc.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension nBasM(MxSym), nOrbM(MxSym), nDelM(MxSym), nFroM(MxSym)
dimension nBas(MxSym)
dimension iToc(64)
parameter(lenin8=6+8)
parameter(maxbfn=10000)
character NameM*(lenin8*maxbfn), firstind*10

! Wilkommen.

write(6,*)
write(6,*)
write(6,*) 'Reading MO-transformed integrals. Zeroth hamiltonian constructed.'

! Numbers and files.

nSize = iOrb(1)*(iOrb(1)+1)/2
call GetMem('FockM','Allo','Real',iPointF,nSize)
call GetMem('SUPER','Allo','Real',iSupM,nSize**2)
iLu1 = 56
iLu2 = 58
call DaName(iLu1,'TRAONE')
call DaName(iLu2,'TRAINT')
iDisk = 0
! This is special utility to read header of TRAONE.
! Last argument depends on mxorb in Molcas.fh.
call Wr_Motra_Info(iLu1,2,iDisk,iToc,64,Ecor,nSymM,nBasM,nOrbM,nFroM,nDelM,MxSym,NameM,lenin8*maxbfn)

! One checks.

if (nBasM(1) /= nBas(1)) then
  write(6,*)
  write(6,*) '  ERROR! Conflict between one-electron file and MO-transformed one-electron file.'
  write(6,*) '         nBas=',nBas(1),' MO-nBas=',nBasM(1)
  call Quit(_RC_GENERAL_ERROR_)
end if

! Read one-electron matrix elements.

iDisk = iToc(2)
call dDaFile(iLu1,2,Work(iPointF),nSize,iDisk)
call dcopy_(nSize,Work(iPointF),iONE,HHmat,iONE)
call GetMem('FockM','Free','Real',iPointF,nSize)
call DaClos(iLu1)

! Add external perturbation if requested.

if (AddExt) then
  write(6,*) '    -- Adding external perturbation.'
  nBTri = nBas(1)*(nBas(1)+1)/2
  Lu_One = 49
  Lu_One = IsFreeUnit(Lu_One)
  call OpnOne(irc,0,'ONEINT',Lu_One)
  call GetMem('AOExt','Allo','Real',ipAOx,nBTri+4)
  call GetMem('TEMP','Allo','Real',iTEMP,nBas(1)*iOrb(1))
  call GetMem('Final','Allo','Real',iFine,iOrb(1)**2)
  call GetMem('Squared','Allo','Real',iSqAO,nBas(1)**2)
  call GetMem('MOExt','Allo','Real',ipMOx,nSize)
  do iExt=1,nExtAddOns
    irc = -1
    iopt = 0
    iSmLbl = 0
    call RdOne(irc,iopt,ExtLabel(iExt),iCompExt(iExt),Work(ipAOx),iSmLbl)
    call DScal_(nBTri,ScalExt(iExt),Work(ipAOx),iONE)
    if (irc /= 0) then
      write(6,*)
      write(6,*) 'ERROR when reading ',ExtLabel(iExt),'.'
      write(6,*) 'Have Seward computed this integral?'
      call Quit(_RC_IO_ERROR_READ_)
    end if
    call Square(Work(ipAOx),Work(iSqAO),iONE,nBas(1),nBas(1))
    call Dgemm_('T','N',iOrb(1),nBas(1),nBas(1),ONE,Work(iV1),nBas(1),Work(iSqAO),nBas(1),ZERO,Work(iTEMP),iOrb(1))
    call Dgemm_('N','N',iOrb(1),iOrb(1),nBas(1),ONE,Work(iTEMP),iOrb(1),Work(iV1),nBas(1),ZERO,Work(iFine),iOrb(1))
    call SqToTri_Q(Work(iFine),Work(ipMOx),iOrb(1))
    call DaxPy_(nSize,ONE,Work(ipMOx),iONE,HHmat,iONE)
  end do
  call GetMem('AOExt','Free','Real',ipAOx,nBTri+4)
  call GetMem('TEMP','Free','Real',iTEMP,nBas(1)*iOrb(1))
  call GetMem('Final','Free','Real',iFine,iOrb(1)**2)
  call GetMem('Squared','Free','Real',iSqAO,nBas(1)**2)
  call GetMem('MOExt','Free','Real',ipMOx,nSize)
  call ClsOne(irc,Lu_One)
end if

! Now to the two-electron matrix elements.

iDisk = 0
call iDaFile(iLu2,2,iTraToc,nTraToc,iDisk)
iDisk = iTraToc(1)
iSup = 0
! Ooohhhh, lets get crude! Read ALL (yes, you read right) integrals and then order them.
nBuf1 = nOrbM(1)*(nOrbM(1)+1)/2
nBuf2 = nBuf1*(nBuf1+1)/2

! Let's check if this construct is possible. If not advise user what to do.

call GetMem('MAX','Max','Real',iDum,nMAX)
if (nMAX < (nBuf2+nBuf1**2)) then
  write(6,*)
  write(6,*) '  Too many MO-transformed two-electron integrals from Motra. Do you need all?'
  write(6,*) '  If not, then use the DELEte keyword in Motra to remove the superfluous ones.'
  call Quit(_RC_GENERAL_ERROR_)
end if

! Proceed!

call GetMem('Buffer','Allo','Real',iBuff,nBuf2)
call GetMem('Temporary','Allo','Real',iTEMP,nBuf1**2)
call dDaFile(iLu2,2,Work(iBuff),nBuf2,iDisk)
do i=1,nBuf1
  do j=i,nBuf1
    iSup = iSup+1
    if ((i <= nSize) .and. (j <= nSize)) then
      Work(iTEMP+(i-1)*nBuf1+j-1) = Work(iBuff+iSup-1)
      Work(iTEMP+(j-1)*nBuf1+i-1) = Work(iBuff+iSup-1)
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
        ij = ipair_qmstat(i,j)
        ik = ipair_qmstat(i,k)
        il = ipair_qmstat(i,l)
        jk = ipair_qmstat(j,k)
        jl = ipair_qmstat(j,l)
        kl = ipair_qmstat(k,l)
        Work(iSupM+nSize*(ij-1)+kl-1) = Work(iTEMP+(ij-1)*nBuf1+kl-1)-(Work(iTEMP+(ik-1)*nBuf1+jl-1)+ &
                                        Work(iTEMP+(il-1)*nBuf1+jk-1))/4
        Work(iSupM+nSize*(kl-1)+ij-1) = Work(iSupM+nSize*(ij-1)+kl-1)
        Work(iSupM+nSize*(ik-1)+jl-1) = Work(iTEMP+(ik-1)*nBuf1+jl-1)-(Work(iTEMP+(ij-1)*nBuf1+kl-1)+ &
                                        Work(iTEMP+(il-1)*nBuf1+jk-1))/4
        Work(iSupM+nSize*(jl-1)+ik-1) = Work(iSupM+nSize*(ik-1)+jl-1)
        Work(iSupM+nSize*(il-1)+jk-1) = Work(iTEMP+(il-1)*nBuf1+jk-1)-(Work(iTEMP+(ik-1)*nBuf1+jl-1)+ &
                                        Work(iTEMP+(ij-1)*nBuf1+kl-1))/4
        Work(iSupM+nSize*(jk-1)+il-1) = Work(iSupM+nSize*(il-1)+jk-1)
      end do
    end do
  end do
end do
call GetMem('Buffer','Free','Real',iBuff,nBuf2)
call GetMem('Temporary','Free','Real',iTEMP,nBuf1**2)
call DaClos(iLu2)

! Serious amount of printing!

if (iPrint >= 35) then
  write(6,*)
  write(6,*) 'The Super Matrix in all its divine g(l)ory:'
  kaunter = 0
  do i=1,iOrb(1)
    do j=1,i
      write(firstind,'(I3,A,I3)') i,',',j
      call TriPrt(firstind,' ',Work(iSupM+nSize*kaunter),iOrb(1))
      kaunter = kaunter+1
    end do
  end do
  write(6,*) 'Super Matrix End.'
end if
write(6,*) '...Done!'

! The end.

return

end subroutine ScfH0
