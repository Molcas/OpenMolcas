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

subroutine Cho_VecTransp(Vec,Jin,Jfi,iSym,iRed,iPass)

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs
use Cholesky, only: Cho_AdrVec, Cho_Real_Par, iiBstR, iiBstR_G, iL2G, IndRed, InfVec_G, IndRed, LuCho_G, LuPri, MaxVec, myNumCho, &
                    nnBstR, nnBstR_G
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: u6
#endif
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Vec(*)
integer(kind=iwp), intent(in) :: Jin, Jfi, iSym, iRed, iPass
character(len=*), parameter :: SecNam = 'Cho_VecTransp'
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
integer(kind=iwp) :: g_a, i, i1, iAdr, iCount, iNode, iOpt, irc, iRSL, iStart, iVec, iVec1, j, j1, Jin0, jRed, jv, jVec, LastV, &
                     lTot, MxRSL, MyEnd, myStart, nProcs_eff, nRS_g, nRS_l, nV, nVR
logical(kind=iwp) :: ok
integer(kind=iwp), allocatable :: iAdrLG(:,:), iVecR(:), Map(:), MapRS2RS(:), nRSL(:)
real(kind=wp), allocatable :: VecR(:,:)
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LocDbg = _DBG_
logical(kind=iwp), external :: ga_create_irreg, ga_destroy
#ifndef _GA_
!VVP:2014 DGA is here
#include "WrkSpc.fh"
integer(kind=iwp) :: iGAL, nelm
integer(kind=iwp), external :: ga_local_woff
logical(kind=iwp), external :: ga_create_local
#endif

if (.not. Cho_Real_Par) then
  if (LocDbg) then
    write(u6,'(A,A,A)') 'Illegal call to ',SecNam,':'
    write(u6,*) 'Should only be called in parallel, but Cho_Real_Par = ',Cho_Real_Par
  end if
  call Cho_Quit('Illegal call to '//SecNam,103)
end if

if (iRed == 2) then
  jRed = 3
else if (iRed == 3) then
  jRed = 2
else
  call Cho_Quit('iRed must be 2 or 3 in '//SecNam,104)
end if
nRS_l = nnBstR(iSym,iRed)   ! local  red set dimension
nRS_g = nnBstR_G(iSym,iRed) ! global red set dimension
nV = Jfi-Jin+1
nVR = 0

call mma_allocate(iVecR,nV,Label='iVecR')
call cho_p_distrib_vec(Jin,Jfi,iVecR,nVR)
call mma_allocate(VecR,nRS_g,nVR+1,Label='VecR')

call mma_allocate(nRSL,nProcs,Label='nRSL')
nRSL(:) = 0
nRSL(1+MyRank) = nRS_l  ! MyRank starts from 0
call Cho_GAIGOP(nRSL,nProcs,'+')

MxRSL = nRSL(1)
do i=2,nProcs
  MxRSL = max(MxRSL,nRSL(i))
end do
call mma_allocate(iAdrLG,MxRSL,nProcs,Label='iAdrLG')

call mma_allocate(Map,nProcs,Label='Map')
nProcs_eff = 0
iStart = 1
myStart = 0
do i=1,nProcs
  if (nRSL(i) > 0) then
    nProcs_eff = nProcs_eff+1
    Map(nProcs_eff) = iStart
    if ((i-1) == myRank) myStart = iStart
    iStart = iStart+nRSL(i)
  end if
end do

if (LocDbg) then
  write(LuPri,*)
  write(LuPri,*) SecNam,': debug info.'
  write(LuPri,*) '#nodes: ',nProcs,'  myRank: ',myRank
  write(LuPri,*) '#contributing nodes: ',nProcs_eff
  write(LuPri,*) 'Symmetry block: ',iSym
  write(LuPri,*) 'On this node:'
  write(LuPri,*) 'Vector dimension : ',nRS_l
  write(LuPri,*) 'Number of vectors: ',nV,' (',Jin,'-',Jfi,')'
  write(LuPri,*) 'Global vector dimension : ',nRS_g
  write(LuPri,*) 'Number of global vectors: ',nVR
  write(Lupri,*) 'MAP:'
  write(LuPri,*) (map(i),i=1,nProcs_eff)
end if
!VVP:2014 Local rather than Global
#ifdef _GA_
ok = ga_create_irreg(mt_dbl,nRS_g,nV,'Ga_Vec',Map,nProcs_eff,1,1,g_a)
#else
ok = ga_create_local(mt_dbl,nRS_g,nV,'Ga_Vec',g_a)
#endif
if (.not. ok) call Cho_Quit(SecNam//': ga_create_irreg error',101)

if (nRS_l > 0) then
  myEnd = myStart+nRS_l-1
# ifdef _GA_
  call ga_put(g_a,myStart,myEnd,1,nV,Vec,nRS_l)
# else
  !VVP:2014 the minimal latency and scalable putC call
  call ga_putc(g_a,myStart,myEnd,1,nV,Vec,nRS_l)
# endif
end if
#ifndef _GA_
nelm = nRS_g*nV
iGAL = ga_local_woff(g_a)
call Cho_GAdGOP(Work(iGAL),nelm,'+')
#else
call GASync()
#endif
Jin0 = Jin-1
do i=1,nVR
  jv = iVecR(i)-Jin0
# ifdef _GA_
  call ga_get(g_a,1,nRS_g,jv,jv,VecR(:,i),nRS_g)
# else
  !VVP:2014 the minimal latency and scalable getC call
  call ga_getc(g_a,1,nRS_g,jv,jv,VecR(:,i),nRS_g)
# endif
end do

ok = ga_destroy(g_a)
if (.not. ok) call Cho_Quit(SecNam//': ga_destroy error',101)

! write the reordered vec on disk

call Cho_P_IndxSwp()
irc = -1
call Cho_X_RSCopy(irc,1,jRed)
if (irc /= 0) call Cho_Quit(SecNam//': Non-zero return code from Cho_X_RSCopy',104)

call mma_allocate(MapRS2RS,nnBstR(iSym,1),Label='MapRS2RS')
call Cho_RS2RS(mapRS2RS,size(mapRS2RS),jRed,iRed,iPass,iSym)
call Cho_P_IndxSwp()

iAdrLG(:,:) = 0
do i=1,nRS_l
  i1 = IndRed(iiBstR(iSym,iRed)+i,iRed) ! addr in local rs1
  j1 = iL2G(i1) ! addr in global rs1
  j = mapRS2RS(j1-iiBstR_G(iSym,1)) ! addr in glob. rs
  iAdrLG(i,myRank+1) = j
end do
call Cho_GAIGOP(iAdrLG,size(iAdrLG),'+')

call mma_deallocate(MapRS2RS)

if (LocDbg) then
  iCount = sum(nRSL(1:nProcs))
  if (iCount /= nRS_g) call Cho_Quit('nRSL error in '//SecNam,104)
end if

do j=1,nVR
  VecR(:,nVr+1) = VecR(:,j)
  iCount = 0
  do iNode=1,nProcs
    do iRSL=1,nRSL(iNode)
      VecR(iAdrLG(iRSL,iNode),j) = VecR(1+iCount,nVR+1)
      iCount = iCount+1
    end do
  end do
end do

if (CHO_ADRVEC /= 1) then ! only WA files!!
  call Cho_Quit('CHO_ADRVEC error in '//SecNam,102)
else ! write to disk and update InfVec_G(*,3,iSym)
  iVec1 = myNumCho(iSym)+1
  lTot = nRS_g*nVR
  if (lTot > 0) then
    iOpt = 1
    iAdr = InfVec_G(iVec1,3,iSym)
    call dDAfile(LuCho_G(iSym),iOpt,VecR,lTot,iAdr)
  end if
  do iVec=1,nVR
    jVec = iVec1+iVec-1
    if (jVec < MaxVec) InfVec_G(jVec+1,3,iSym) = InfVec_G(jVec,3,iSym)+nRS_g
  end do
  LastV = myNumCho(iSym)+nVR
  if (LastV > MaxVec) call Cho_Quit('Max. number of vectors exceeded in '//SecNam,104)
end if
myNumCho(iSym) = myNumCho(iSym)+nVR

! deallocations

call mma_deallocate(Map)
call mma_deallocate(iAdrLG)
call mma_deallocate(nRSL)
call mma_deallocate(VecR)
call mma_deallocate(iVecR)
#else

#include "macros.fh"
unused_var(Vec(1))
unused_var(Jin)
unused_var(Jfi)
unused_var(iSym)
unused_var(iRed)
unused_var(iPass)
call Cho_Quit(SecNam//' should never be called in serial installation',103)

#endif

end subroutine Cho_VecTransp
