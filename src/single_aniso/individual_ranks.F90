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

subroutine individual_ranks(nDIMCF,BC,BS,Hinit,LJ,iprint)

use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, u6

implicit none
#include "stdalloc.fh"
integer, intent(in) :: nDIMCF, iprint
real(kind=8), intent(in) :: BC(nDIMcf,0:nDIMcf)
real(kind=8), intent(in) :: BS(nDIMcf,0:nDIMcf)
complex(kind=8), intent(in) :: Hinit(nDIMcf,nDIMcf)
character(len=1), intent(in) :: LJ
! local variables:
integer :: i, j, k, n, m, ik, jEnd, nfields, info, ikmax, ip, ir, iq
integer, allocatable :: rankKQ(:) ! nDIMCF*(2*nDIMCF+1) )
integer, allocatable :: projKQ(:) ! nDIMCF*(2*nDIMCF+1) )
real(kind=8) :: Tnrm, TnrmKQ, wt
real(kind=8) :: RnrmKQ(nDIMcf,-nDIMcf:nDIMcf)
real(kind=8), allocatable :: ListKQ(:) ! nDIMCF*(2*nDIMCF+1) )
real(kind=8), allocatable :: Rnrm(:) !nDIMCF)
real(kind=8), allocatable :: Snrm(:) !nDIMCF)
real(kind=8), allocatable :: Wk(:,:) !nDIMCF,nDIMCF)
real(kind=8), allocatable :: Ws(:,:) !nDIMCF,nDIMCF)
real(kind=8), allocatable :: Winit(:) !nDIMCF)
real(kind=8), external :: dznrm2_
complex(kind=8) :: zf, redME
complex(kind=8), allocatable :: O(:,:)  !nDIMCF,nDIMCF)! real ITO
complex(kind=8), allocatable :: W(:,:)  !nDIMCF,nDIMCF)! imag ITO
complex(kind=8), allocatable :: HCF(:,:,:) !nDIMCF,nDIMCF,nDIMCF)
complex(kind=8), allocatable :: HCFS(:,:,:)!nDIMCF,nDIMCF,nDIMCF)
complex(kind=8), allocatable :: Zk(:,:,:)  !nDIMCF,nDIMCF,nDIMCF)
complex(kind=8), allocatable :: Zs(:,:,:)  !nDIMCF,nDIMCF,nDIMCF)
complex(kind=8), allocatable :: Zinit(:,:) !nDIMCF,nDIMCF)
complex(kind=8), allocatable :: HKQ(:,:), Zkq(:,:)
real(kind=8) :: dznrm2
external :: dznrm2
character(len=16) :: field(8)
character(len=6) :: iprog

!-----------------------------------------------------------------------
call mma_allocate(Rnrm,nDIMCF,'Rnrm')
call mma_allocate(Snrm,nDIMCF,'Snrm')
call mma_allocate(Wk,nDIMCF,nDIMCF,'Wk')
call mma_allocate(Ws,nDIMCF,nDIMCF,'Ws')
call mma_allocate(Winit,nDIMCF,'Winit')
call mma_allocate(ListKQ,nDIMCF*(2*nDIMCF+1),'ListKQ')
call mma_allocate(rankKQ,nDIMCF*(2*nDIMCF+1),'rankKQ')
call mma_allocate(projKQ,nDIMCF*(2*nDIMCF+1),'projKQ')

call mma_allocate(O,nDIMCF,nDIMCF,'O')
call mma_allocate(W,nDIMCF,nDIMCF,'W')
call mma_allocate(HKQ,nDIMCF,nDIMCF,'HKQ')
call mma_allocate(ZKQ,nDIMCF,nDIMCF,'ZKQ')
call mma_allocate(HCF,nDIMCF,nDIMCF,nDIMCF,'HCF')
call mma_allocate(HCFS,nDIMCF,nDIMCF,nDIMCF,'HCFS')
call mma_allocate(Zk,nDIMCF,nDIMCF,nDIMCF,'Zk')
call mma_allocate(Zs,nDIMCF,nDIMCF,nDIMCF,'Zs')
call mma_allocate(Zinit,nDIMCF,nDIMCF,'Zinit')

call dcopy_(nDIMCF,[Zero],0,Rnrm,1)
call dcopy_(nDIMCF,[Zero],0,Snrm,1)
Tnrm = Zero
call dcopy_(nDIMCF*nDIMCF,[Zero],0,Wk,1)
call dcopy_(nDIMCF*nDIMCF,[Zero],0,Ws,1)
call dcopy_(nDIMCF,[Zero],0,Winit,1)

call zcopy_(nDIMCF*nDIMCF,[cZero],0,Zinit,1)
call zcopy_(nDIMCF*nDIMCF*nDIMCF,[cZero],0,HCF,1)
call zcopy_(nDIMCF*nDIMCF*nDIMCF,[cZero],0,HCFS,1)
call zcopy_(nDIMCF*nDIMCF*nDIMCF,[cZero],0,Zk,1)
call zcopy_(nDIMCF*nDIMCF*nDIMCF,[cZero],0,Zs,1)
!-----------------------------------------------------------------------
! re-construct the  initial CF matrix:
do N=2,nDIMcf-1,2
  do M=0,N
    call Liviu_ESO(nDIMcf,N,M,O,W,redME)
    if (M == 0) then
      zf = BC(N,0)*cOne
      call zaxpy_(nDIMcf*nDIMcf,zf,O,1,HCF(N,1:nDIMcf,1:nDIMcf),1)
    else
      zf = BC(N,M)*cOne
      call zaxpy_(nDIMcf*nDIMcf,zf,O,1,HCF(N,1:nDIMcf,1:nDIMcf),1)
      zf = BS(N,M)*cOne
      call zaxpy_(nDIMcf*nDIMcf,zf,W,1,HCF(N,1:nDIMcf,1:nDIMcf),1)
    end if
  end do
end do

do k=2,nDIMCF-1,2
  do ik=2,k
    ! add the ranks as follows:
    ! O2= O2
    ! O4= O2+O4
    ! O6= O2+O4+O6
    ! O8= O2+O4+O6+O8
    ! etc...
    ! compute the cumulative CF matrix
    do i=1,nDIMcf
      do j=1,nDIMcf
        HCFS(k,i,j) = HCFS(k,i,j)+HCF(ik,i,j)
      end do
    end do
  end do
end do

do N=2,nDIMCF-1,2
  Rnrm(N) = dznrm2(nDIMcf*nDIMcf,HCF(N,:,:),1)
  Tnrm = Tnrm+Rnrm(N)
  do i=2,N,2
    Snrm(N) = Snrm(N)+Rnrm(i)
  end do
end do
! compute the CF spinting of individual weight operators:
do k=2,nDIMcf-1,2
  call DIAG_C2(HCF(k,:,:),nDIMcf,info,Wk(k,:),Zk(k,:,:))
  call DIAG_C2(HCFS(k,:,:),nDIMcf,info,Ws(k,:),Zs(k,:,:))
end do
! set the initial energies as to the sum of all contributions:
call DIAG_C2(Hinit,nDIMcf,info,Winit,Zinit)

!-----------------------------------------------------------------------
! individual parameter contribution:
TnrmKQ = Zero
RnrmKQ(:,:) = Zero
ListKQ(:) = Zero
projKQ(:) = 0
rankKQ(:) = 0
ik = 0
do N=2,nDIMcf-1,2
  do M=-N,N
    ik = ik+1
    ! generate ITO operators:
    call Liviu_ITO(nDIMcf,N,abs(M),O,W,redME)
    ! generate HCF for each parameter rank and projection:
    call zcopy_(nDIMcf*nDIMcf,[cZero],0,HKQ,1)
    if (M < 0) then
      zf = BS(N,abs(M))*cOne
      call zaxpy_(nDIMcf*nDIMcf,zf,W,1,HKQ,1)
    else if (M >= 0) then
      zf = BC(N,abs(M))*cOne
      call zaxpy_(nDIMcf*nDIMcf,zf,O,1,HKQ,1)
    end if
    ! find the rank value of each NM operator
    ! for further estimate the effect of the corresponding CF parameter
    RnrmKQ(N,M) = dznrm2_(nDIMcf*nDIMcf,HKQ(:,:),1)
    TnrmKQ = TnrmKQ+RnrmKQ(N,M)

    ! indexing lists
    ListKQ(ik) = RnrmKQ(N,M)
    rankKQ(ik) = N
    projKQ(ik) = M
  end do
end do
ikmax = ik
! make an ordered list of CF parameters
! following their descending effect:
call sort_KQ(ikmax,ListKQ,rankKQ,projKQ,2)
!--------- below we just print the data --------------------------------
if (LJ == 'L') then
  iprog = 'RASSCF'
else
  iprog = ' RASSI'
end if

if (iprint >= 4) then
  write(u6,'(A,F20.12)') 'Tnrm=',Tnrm
  do k=1,nDIMCF-1
    write(u6,'(2(A,i2,A,F20.12,2x))') 'Rnrm(',k,')=',Rnrm(k),'Snrm(',k,')=',Snrm(k)
  end do
end if

! print the cumulative and individual rank weight:
write(u6,'(/)')
write(u6,'(A)') 'CUMULATIVE WEIGHT OF INDIVIDUAL-RANK OPERATORS ON THE CRYSTAL FIELD SPLITTING:'
if ((nDIMCF-1) >= 2) write(u6,'(2x,A,F10.6,A)') 'O2 :------------------------------------------: ',(Snrm(2)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 4) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 :-------------------------------------: ',(Snrm(4)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 6) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 :--------------------------------: ',(Snrm(6)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 8) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 :---------------------------: ',(Snrm(8)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 10) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 :---------------------: ',(Snrm(10)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 12) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 + O12 :---------------: ',(Snrm(12)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 14) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 + O12 + O14 :---------: ',(Snrm(14)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 16) write(u6,'(2x,A,F10.6,A)') 'O2 + O4 + O6 + O8 + O10 + O12 + O14 + O16 :---: ',(Snrm(16)/Tnrm)*100.0_wp,' %.'

!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(A)') 'ENERGY SPLITTING INDUCED BY CUMMULATIVE INDIVIDUAL-RANK OPERATORS (in cm-1).'
field(1) = '      O2       |'
field(2) = '     O2+O4     |'
field(3) = '   O2+O4+O6    |'
field(4) = '  O2+O4+O6+O8  |'
field(5) = ' O2+O4+...+O10 |'
field(6) = ' O2+O4+...+O12 |'
field(7) = ' O2+O4+...+O14 |'
field(8) = ' O2+O4+...+O16 |'
nfields = int((nDIMCF-1)/2)
do J=1,nfields,4
  jEnd = min(nfields,J+3)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)

  if (mod(nDIMcf,2) == 1) then
    write(u6,'(3A,I2,3x,A,5x,A,4x,A,10(5x,A,6x,A))') '  ',LJ,' =',(nDIMcf-1)/2,'|',iprog,'|',('ONLY','|',i=j,jEnd)
  else
    write(u6,'(3A,I3,A,5x,A,4x,A,10(5x,A,6x,A))') ' ',LJ,' =',(nDIMcf-1),'/2 |',iprog,'|',('ONLY','|',i=j,jEnd)
  end if

  write(u6,'(10x,10A)') '|     INITIAL   |',(field(i),i=j,jEnd)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
  do i=1,nDIMcf
    write(u6,'(1x,A,1x,i2,2x,a,12(f14.8,1x,a))') 'w.f.',i,'|',(Winit(i)-Winit(1)),'|',((Ws(k,i)-Ws(k,1)),'|',k=2*j,2*jEnd,2)
  end do
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
end do ! j

write(u6,'(/)')
write(u6,'(A)') 'WEIGHT OF INDIVIDUAL-RANK OPERATORS ON THE CRYSTAL FIELD SPLITTING:'
if ((nDIMCF-1) >= 2) write(u6,'(2x,A,F10.6,A)') 'O2  :-----------------------------------------: ',(Rnrm(2)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 4) write(u6,'(2x,A,F10.6,A)') 'O4  :-----------------------------------------: ',(Rnrm(4)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 6) write(u6,'(2x,A,F10.6,A)') 'O6  :-----------------------------------------: ',(Rnrm(6)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 8) write(u6,'(2x,A,F10.6,A)') 'O8  :-----------------------------------------: ',(Rnrm(8)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 10) write(u6,'(2x,A,F10.6,A)') 'O10 :-----------------------------------------: ',(Rnrm(10)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 12) write(u6,'(2x,A,F10.6,A)') 'O12 :-----------------------------------------: ',(Rnrm(12)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 14) write(u6,'(2x,A,F10.6,A)') 'O14 :-----------------------------------------: ',(Rnrm(14)/Tnrm)*100.0_wp,' %.'
if ((nDIMCF-1) >= 16) write(u6,'(2x,A,F10.6,A)') 'O16 :-----------------------------------------: ',(Rnrm(16)/Tnrm)*100.0_wp,' %.'

!-----------------------------------------------------------------------
write(u6,'(/)')
write(u6,'(A)') 'ENERGY SPLITTING INDUCED BY INDIVIDUAL-RANK OPERATORS (in cm-1).'
field(1) = '      O2       |'
field(2) = '      O4       |'
field(3) = '      O6       |'
field(4) = '      O8       |'
field(5) = '      O10      |'
field(6) = '      O12      |'
field(7) = '      O14      |'
field(8) = '      O16      |'
nfields = int((nDIMCF-1)/2)

do J=1,nfields,4
  jEnd = min(nfields,J+3)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)

  if (mod(nDIMcf,2) == 1) then
    write(u6,'(3A, I2,3x,A,5x,A,4x,A,10(5x,A,6x,A))') '  ',LJ,' =',(nDIMcf-1)/2,'|',iprog,'|',('ONLY','|',i=j,jEnd)
  else
    write(u6,'(3A,I3,a,5x,A,4x,A,10(5x,A,6x,A))') ' ',LJ,' =',(nDIMcf-1),'/2 |',iprog,'|',('ONLY','|',i=j,jEnd)
  end if

  write(u6,'(10x,10A)') '|     INITIAL   |',(field(i),i=j,jEnd)
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
  do i=1,nDIMcf
    write(u6,'(1x,A,1x,i2,2x,a,12(f14.8,1x,a))') 'w.f.',i,'|',(Winit(i)-Winit(1)),'|',((Wk(k,i)-Wk(k,1)),'|',k=2*j,2*jEnd,2)
  end do
  write(u6,'(16A)') '----------|',('---------------|',i=j,jEnd+1)
end do ! j

!-----------------------------------------------------------------------
write(u6,'(/)')
write(u6,'(A)') 'WEIGHT OF INDIVIDUAL CRYSTAL FIELD PARAMETERS ON THE CRYSTAL FIELD SPLITTING: (in descending order):'
write(u6,'(A)') 'CFP are given in ITO used in J. Chem. Phys. 137, 064112 (2012).'

write(u6,'(100A)') ('-',i=1,55),'|'
write(u6,'(A)') '  k |  q  |         B(k,q)        |    Weight (in %)   |'
write(u6,'(A)') '----|-----|-----------------------|--------------------|'
do ik=1,ikmax
  ip = projKQ(ik)
  iq = abs(projKQ(ik))
  ir = rankKQ(ik)
  wt = 100.0_wp*ListKQ(ik)/TnrmKQ

  if (projKQ(ik) >= 0) then
    write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),(ES22.14,1x,A),F19.14,1x,A)') ir,'|',ip,'|',BC(ir,iq),'|',wt,'|'
  else
    write(u6,'((1x,I2,1x,A),(1x,I3,1x,A),(ES22.14,1x,A),F19.14,1x,A)') ir,'|',ip,'|',BS(ir,iq),'|',wt,'|'
  end if
end do
write(u6,'(A)') '----|-----|-----------------------|--------------------|'

!-----------------------------------------------------------------------
call mma_deallocate(Rnrm)
call mma_deallocate(Snrm)
call mma_deallocate(Wk)
call mma_deallocate(Ws)
call mma_deallocate(Winit)

call mma_deallocate(O)
call mma_deallocate(W)
call mma_deallocate(HKQ)
call mma_deallocate(HCF)
call mma_deallocate(HCFS)
call mma_deallocate(Zk)
call mma_deallocate(Zkq)
call mma_deallocate(Zs)
call mma_deallocate(Zinit)
call mma_deallocate(ListKQ)
call mma_deallocate(rankKQ)
call mma_deallocate(projKQ)

return

end subroutine individual_ranks
