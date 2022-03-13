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

subroutine RasRasTrans(nB,nStatePrim,iEig2,iPrint)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "files_qmstat.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension iTocBig(MxStOT)
character*30 OutLine


!Guten Tag.

write(6,*) '     ----- Transform from non-orthogonal RASSCF states to orthogonal RASSI states.'

! Set zeros and decide if transformation is at all possible.

LuIn = 66
kaunt = 0
iDisk = 0
call DaName(LuIn,RassiM)
nTriSP = nStatePrim*(nStatePrim+1)/2
call iDaFile(LuIn,2,iTocBig,nTriSP,iDisk)
nSize = nB*(nB+1)/2
nTriS = nState*(nState+1)/2
nSizeBig = nSize*nTriS
nSizeBigPrim = nSize*nTriSP
call GetMem('HOWMUCH','Max','Real',ipMAX,MEMMAX)

! This means that we do not have memory enough for TDM in contracted
! form. Then there is no use to proceed at all.

if (MEMMAX <= nSizeBig) then
  write(6,*)
  write(6,*) 'The transition density matrix is too big to put in memory!'
  write(6,*) 'Either,'
  write(6,*) '       (1) increase MOLCAS_MEM,'
  write(6,*) '       (2) contract number of states further.'
  call Quit(_RC_GENERAL_ERROR_)

! Here we go if there is enough memory for an in core transformation.

else if (MEMMAX >= (nSizeBig+nSizeBigPrim+nTriSP+nTriS+nStatePrim**2+nState*nStatePrim+nState**2)) then
  call GetMem('ALLES','Allo','Real',iBigT,nSizeBig)
  call GetMem('ALLESin','Allo','Real',iBigV,nSizeBigPrim)
  call GetMem('Int1','Allo','Real',iInt1,nTriSP)
  call GetMem('Int2','Allo','Real',iInt2,nTriS)
  call GetMem('Square1','Allo','Real',iSnt1,nStatePrim**2)
  call GetMem('Square2','Allo','Real',iSnt2,nState*nStatePrim)
  call GetMem('Square3','Allo','Real',iSnt3,nState**2)
  call dcopy_(nSizeBig,[ZERO],iZERO,Work(iBigT),iONE)
  kaunt = 0
  do i=1,nStatePrim
    do j=1,i
      kaunt = kaunt+1
      iDisk = iTocBig(kaunt)
      call dDaFile(LuIn,2,Work(iBigV+(kaunt-1)*nSize),nSize,iDisk)
    end do
  end do

  ! A lot of printing of TDM if requested.

  if (iPrint >= 25) then
    kaunt = 0
    do i=1,nStatePrim
      do j=1,i
        write(OutLine,'(A,I3,I3)') 'TDM, Piece ',i,j
        call TriPrt(OutLine,' ',Work(iBigV+kaunt),nB)
        kaunt = kaunt+nSize
      end do
    end do
  end if

  ! Proceed with transformation.

  kaunt = 0
  do,iBas = 1,nB
  do,jBas = 1,iBas
  call dcopy_(nTriSP,Work(iBigV+kaunt),nSize,Work(iInt1),iONE)
  call Square(Work(iInt1),Work(iSnt1),iONE,nStatePrim,nStatePrim)
  call Dgemm_('T','N',nState,nStatePrim,nStatePrim,ONE,Work(iEig2),nStatePrim,Work(iSnt1),nStatePrim,ZERO,Work(iSnt2),nState)
  call Dgemm_('N','N',nState,nState,nStatePrim,ONE,Work(iSnt2),nState,Work(iEig2),nStatePrim,ZERO,Work(iSnt3),nState)
  call SqToTri_Q(Work(iSnt3),Work(iInt2),nState)
  call dcopy_(nTriS,Work(iInt2),iONE,Work(iBigT+kaunt),nSize)
  kaunt = kaunt+1
  end do
  end do
  call GetMem('ALLESin','Free','Real',iBigV,nSizeBigPrim)
  call GetMem('Int1','Free','Real',iInt1,nTriSP)
  call GetMem('Int2','Free','Real',iInt2,nTriS)
  call GetMem('Square1','Free','Real',iSnt1,nStatePrim**2)
  call GetMem('Square2','Free','Real',iSnt2,nState*nStatePrim)
  call GetMem('Square3','Free','Real',iSnt3,nState**2)

  ! Here we go if both TDM's can not be put in memory. Might be a bit
  ! slow due to its nested nature with repeated IO.

  else
  call GetMem('ALLES','Allo','Real',iBigT,nSizeBig)
  call GetMem('AOGamma','Allo','Real',ipAOG,nSize)
  call dcopy_(nSizeBig,[ZERO],iZERO,Work(iBigT),iONE)
  do iiS=1,nStatePrim
    do jjS=1,nStatePrim
      if (iiS <= jjS) then
        indypop = jjS*(jjS+1)/2-jjS+iiS
      else
        indypop = iiS*(iiS+1)/2-iiS+jjS
      end if
      iDisk = iTocBig(indypop)
      call dDaFile(LuIn,2,Work(ipAOG),nSize,iDisk)
      kaunter = 0
      do iB=1,nB
        do jB=1,iB
          do iS=1,nState
            do jS=1,iS
              index = (iS*(iS-1)/2+jS-1)*nSize
              index = index+kaunter
              Work(iBigT+index) = Work(iBigT+index)+Work(iEig2+iiS-1+(iS-1)*nStatePrim)*Work(iEig2+jjS-1+(jS-1)*nStatePrim)* &
                                  Work(ipAOG+kaunter)
            end do
          end do
          kaunter = kaunter+1
        end do
      end do
    end do
  end do
  call GetMem('AOGamma','Free','Real',ipAOG,nSize)
end if

! Deallocations and finish up.

call GetMem('RedEigV1','Free','Real',iEig2,nStatePrim**2)
call DaClos(LuIn)

return

end subroutine RasRasTrans
