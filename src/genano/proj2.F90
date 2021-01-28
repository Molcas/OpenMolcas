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
! This routine projects out specified vectors from the density matrix. *
! The procedure is to read the projection orbitals and perform a       *
! sequence of rank-1 updated to zero out the contributions <c|D|c>.    *
! This procedure does not necessarily produce a complete projection    *
! of the density matrix, and routine 'proj1' is preferable.            *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine Proj2

use Genano_globals, only: MxLqn, nPrim, iSymBk, tDsym, Ssym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: err, i, iB, iBlk, ij, iLqn, indS, indT, iO, iShell, j, jB, jk, kB, kl, lB, nB, nO, Lu
real(kind=wp) :: cnorm, eval, t
real(kind=wp), allocatable :: ProjOrb(:), TmpDens(:), TmpOvlp(:)
integer(kind=iwp), external :: IsFreeUnit

Lu = IsFreeUnit(17)
call molcas_open(Lu,'PROJ')
iBlk = 0
!--- Loop over l quantum number ---*
do iLqn=0,MxLqn
  !--- Read projection orbitals ---*
  read(Lu,*,iostat=err) nB,nO
  if (err /= 0) exit
  if (nB*nO <= 0) cycle
  if (nB /= nPrim(iLqn)) then
    write(u6,*) 'Project: inconsistency in number of functions'
    call Abend()
  end if
  call mma_allocate(ProjOrb,nB*nB,label='ProjOrb')
  call mma_allocate(TmpDens,nB*nB,label='TmpDens')
  call mma_allocate(TmpDens,nB*nB,label='TmpOvlp')
  do iB=1,nB
    read(Lu,*) (ProjOrb(iB+nB*(iO-1)),iO=1,nO)
  end do
  !write(u6,'(a,i5)') ' Projection orbitals',iLqn
  !do iB=1,nB
  !  write(u6,'(1x,6f12.6)') (ProjOrb(iB+nB*(iO-1)),iO=1,nO)
  !end do
  !--- Loop over m quantum number ---*
  do iShell=-iLqn,iLqn
    iBlk = iBlk+1
    !write(u6,'(a,2i5)') ' Density matrix block',iLqn,iShell
    !call Triprt(' ','(6f12.6)',tDsym(iSymBk(iBlk)),nPrim(iLqn))
    !write(u6,'(a,2i5)') ' Overlap matrix block',iLqn,iShell
    !call Triprt(' ','(6f12.6)',Ssym(iSymBk(iBlk)),nPrim(iLqn))
    !--- Copy S and D into square form ---*
    do i=1,nB
      do j=1,nB
        indT = min(i,j)+max(i,j)*(max(i,j)-1)/2
        indS = i+nB*(j-1)
        TmpDens(indS) = tDsym(iSymBk(iBlk)-1+indT)
        TmpOvlp(indS) = Ssym(iSymBk(iBlk)-1+indT)
        !write(u6,'(1x,4i5,2f12.6)') i,j,IndT,IndS,TmpDens(indS),TmpOvlp(indS)
      end do
    end do
    !--- Project ---*
    do iO=1,nO
      eval = Zero
      do iB=1,nB
        do jB=1,nB
          do kB=1,nB
            do lB=1,nB
              ij = iB+nB*(jB-1)
              jk = jB+nB*(kB-1)
              kl = kB+nB*(lB-1)
              t = ProjOrb(iB+nB*(iO-1))*ProjOrb(lB+nB*(iO-1))
              t = t*TmpOvlp(ij)
              t = t*TmpDens(jk)
              t = t*TmpOvlp(kl)
              eval = eval+t
            end do
          end do
        end do
      end do
      cnorm = Zero
      do iB=1,nB
        do jB=1,nB
          ij = iB+nB*(jB-1)
          t = ProjOrb(iB+nB*(iO-1))*ProjOrb(jB+nB*(iO-1))
          t = t*TmpOvlp(ij)
          cnorm = cnorm+t
        end do
      end do
      !write(u6,'(a,f12.6)') ' Eigenvalue ',eval
      !write(u6,'(a,f12.6)') ' Norm       ',cnorm
      eval = eval/cnorm/cnorm
      do iB=1,nB
        do jB=1,nB
          t = ProjOrb(iB+nB*(iO-1))*ProjOrb(jB+nB*(iO-1))
          TmpDens(iB+nB*(jB-1)) = TmpDens(iB+nB*(jB-1))-eval*t
        end do
      end do
    end do
    !--- Copy back D ---*
    do i=1,nB
      do j=1,nB
        indT = min(i,j)+max(i,j)*(max(i,j)-1)/2
        indS = i+nB*(j-1)
        tDsym(iSymBk(iBlk)-1+indT) = TmpDens(indS)
      end do
    end do
    call mma_deallocate(ProjOrb)
    call mma_deallocate(TmpDens)
    call mma_deallocate(TmpOvlp)
    !write(u6,'(a,2i5)') ' Density matrix block',iLqn,iShell
    !call Triprt(' ','(6f12.6)',tDsym(iSymBk(iBlk)),nPrim(iLqn))
  end do
end do
!--- ---*
close(Lu)

return

end subroutine Proj2
