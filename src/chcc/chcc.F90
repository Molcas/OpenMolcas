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

subroutine chcc(ireturn)
! docasny drajver reorder procesu

use chcc_global, only: conv, deallocate_arrays, generkey, maxiter, nc, printkey, restkey, TCpu, TCpu_l, TCpu0, TWall, TWall_l, &
                       TWall0
use Para_Info, only: MyRank, nProcs
#ifdef _MOLCAS_MPP_
use chcc_global, only: NChLoc
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: iter, Jal1, Jal2, LunAux, maxdim, maxspace, NchBlk, NChHere, NvGrp, NvSGrp, wrksize
real(kind=wp) :: e1new, e1old, e2new, e2old, e2os, escf
real(kind=wp), allocatable :: wrk(:)

!mp

#ifdef _MOLCAS_MPP_
if (is_real_par()) then
  write(u6,*) ' Parallel run'
else
  write(u6,*) ' Serial run'
end if
!mp write(u6,*) ' MyRank, Nprocs',MyRank,Nprocs
#else
write(u6,*) ' Serial run'
!mp write(u6,*) ' MyRank, Nprocs',MyRank,Nprocs
#endif
!mp
! vynuluj hodiny
call CWTime(TCpu,TWall)
TCpu0 = TCpu
TWall0 = TWall
TCpu_l = TCpu
TCpu_l = TCpu
TWall_l = TWall
TWall_l = TWall
!mp

!@@ real*8, allocatable :: wrk(:)
!   real*8 wrk(1:40000000)
!   real*8 wrk(1:26593281)
!   wrksize=26593281
!   wrksize=40000000

!mp ##########################################################
!0  info o cholesky vektoroch

call frankie_drv_fake(NChHere)
write(u6,'(A,i9,A,i4)') ' Number of Cholesky vectors ',NChHere,' on node ',myRank

#ifdef _MOLCAS_MPP_

NChLoc(:) = 0

NChLoc(MyRank) = NChHere
call gaigop(NChLoc(0),NProcs,'+')

jal2 = sum(NChLoc)

nc = jal2
#else
nc = NChHere
#endif
!mp ##########################################################
!
! Get the maximum available memory

call mma_maxDBLE(maxspace)
maxspace = maxspace-8
write(u6,'(A,i13,A,f9.1,A,f5.1,A)') ' Max Size              : ',maxspace,' in r*8 Words,',real(maxspace*8,kind=wp)/1024**2,' Mb,', &
                                    real(maxspace,kind=wp)*8/1024**3,' Gb'

!1 Nacitanie vstupu (Docasne) + time Delay

call IniReord(NvGrp,NvSGrp,NchBlk,LunAux)
!mp call TimeDelay(LunAux) ! temporarily disabled (MP)

! Decide on automatic segmentation generation
if (NvGrp == 0) then
  call autoSegmentation(Nprocs,maxspace,Jal1,Jal2,NvGrp,NvSGrp,NchBlk,wrksize,maxdim)
else
  call checkMem(NvGrp,NvSGrp,NchBlk,Jal1,Jal2,wrksize,maxdim)

  if (wrksize > maxspace) then
    write(u6,*) ' Not Enough Memory! Increase large and/or small segmentation',real(wrksize,kind=wp)/real(maxspace,kind=wp)
    call abend()
  end if
end if

! write Large segmentation to RunFile
call Put_iScalar('CHCCLarge',NvGrp)

! ---------------------------------------------------------
!
! Transformation of Local L(ao) -> L(mo)
! CD1tmp file is created

call frankie_drv(NChHere)
if (printkey >= 10) write(u6,*) ' After Frankie',myRank,NChHere

! ---------------------------------------------------------

call mma_allocate(wrk,wrksize,label='CCSD')
write(u6,'(A,i13,A,f9.1,A,f5.1,A)') ' Real Allocated Memory : ',wrksize,' in r*8 Words,',real(wrksize,kind=wp)*8/1024**2,' Mb,', &
                                    real(wrksize,kind=wp)*8/1024**2,' Gb'
wrk(:) = Zero

!3 Priprava integralov

call Reord_chcc(wrk,wrksize,NvGrp,NvSGrp,NchBlk,LunAux)
if (generkey == 1) then
  !mp if (printkey >= 10) then ! uvidime ...
  write(u6,*) ' Generation of integrals (Reord_chcc) done'
  write(u6,*)
  !mp end if
else
  write(u6,*) ' Generation of integrals (Reord_chcc) skipped, only basic'
  write(u6,*)
end if
!mp
call CWTime(TCpu,TWall)

if (printkey > 1) then
  write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
  write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
  write(u6,*)
  write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
  write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
  write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
  write(u6,*)
end if
TCpu_l = TCpu
TWall_l = TWall
!mp

!4 Restart, ak treba

if (restkey == 1) then
  ! read T1o,niter,E1old,E2old (T2 are in T2files)
  call GetRest(wrk,wrksize,LunAux,iter,E1old,E2old)
else
  !Bug v originale chyba
  ! set T1o=0
  call VanishT1(wrk,wrksize)

  !@@
  !write(u6,*) ' Pred Chck'
  !call MakeChckData(wrk,wrksize,LunAux)
  !call SaveChckData(LunAux)
  !call GetChckData(LunAux)
  !write(u6,*) ' Chck  done'
  !@@

  !5 iteracny cyklus

  e1old = Zero
  e2old = Zero
  iter = 1

  !mp

  write(u6,*)
  write(u6,*) '------------------------'
  write(u6,*) 'Starting CCSD iterations'
  write(u6,*) '------------------------'
  write(u6,*)

  write(u6,*) '                  CCSD Energy      Difference'
  write(u6,*)
  !mp!
  call xflush(u6)
  !mp!

  call CWTime(TCpu,TWall)

  if (printkey > 1) then
    write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
    write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
    write(u6,*)
    write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
    write(u6,*)
  end if
  TCpu_l = TCpu
  TWall_l = TWall
  !mp
end if

do

  call o3v3ctl(wrk,wrksize,NvGrp,LunAux)
  if (printkey > 1) write(u6,*) ' o3v3 done'
  !mp
  call CWTime(TCpu,TWall)
  if (printkey > 1) then
    write(u6,*)
    write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
    write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
    write(u6,*)
    write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
    write(u6,*)
  end if
  TCpu_l = TCpu
  TWall_l = TWall
  !mp
  call o2v4ctl(wrk,wrksize,NvGrp,NvSGrp,LunAux)
  if (printkey > 1) write(u6,*) ' o2v4 done'
  !mp
  call CWTime(TCpu,TWall)
  if (printkey > 1) then
    write(u6,*)
    write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
    write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
    write(u6,*)
    write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
    write(u6,*)
  end if
  TCpu_l = TCpu
  TWall_l = TWall
  !mp
  call summary(wrk,wrksize,NvGrp,LunAux,maxdim,e1new,e2new,e2os)
  if (printkey > 1) write(u6,*) ' summary done'
  !mp
  call CWTime(TCpu,TWall)
  if (printkey > 1) then
    write(u6,*)
    write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
    write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
    write(u6,*)
    write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
    write(u6,*)
  end if
  TCpu_l = TCpu
  TWall_l = TWall
  !mp
  call SaveRest(wrk,wrksize,LunAux,(iter+1),E1new,E2new)

  !mp write(u6,91) ' Iteration :',iter,e1new,e2new,e1old+e2old-e1new-e2new
  !mp 91 format(a12,1x,i3,1x,3(f15.12,1x))

  if (iter == 1) then
    write(u6,91) ' Iteration :',iter,e2new
  else
    write(u6,93) ' Iteration :',iter,e2new,e1old+e2old-e1new-e2new
  end if

  call xflush(u6)

  if ((abs(e1old+e2old-e1new-e2new) <= conv) .or. (iter >= maxiter)) exit
  e1old = e1new
  e2old = e2new
  iter = iter+1

end do

write(u6,*)
write(u6,*) ' Final CCSD energy decomposition'
write(u6,92) ' E1 CCSD energy :',e1new
write(u6,92) ' E2 CCSD energy :',e2new
write(u6,92) ' E2 CCSD ss     :',e2new-e2os
write(u6,92) ' E2 CCSD os     :',e2os
write(u6,*)
!mp for MOLCAS verify
call Get_dScalar('SCF energy',escf)
call Add_Info('CHCCene',[e2new],1,6)
call Add_Info('E_CHCC',[e2new],1,6)
call Add_Info('E_HYPE',[e2new+escf],1,6)
! for NUMERICAL_GRADIENTS
call Put_cArray('Relax Method','CHCC    ',8)
call Store_Energies(1,e2new+escf,1)
!mp

!@@ deallocate(wrk)
call mma_deallocate(wrk)
call deallocate_arrays()

ireturn = 0

return

91 format(a12,1x,i3,1x,f15.12,1x)
92 format(a17,1x,f15.12)
93 format(a12,1x,i3,1x,2(f15.12,1x))

end subroutine chcc
