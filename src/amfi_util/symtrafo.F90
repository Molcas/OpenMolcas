!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SymTrafo(LUPROP,lOper,nComp,nBas,nIrrep,Label,MolWgh)
!bs Purpose: combine SO-integrals from amfi to symmetry-adapted
!bs integrals on one file AOPROPER_MF_SYM

use AMFI_global, only: Lmax, MxCart
use index_functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LUPROP, nComp, lOper(nComp), nIrrep, nBas(0:nIrrep-1), MolWgh
character(len=8), intent(in) :: Label
#include "Molcas.fh"
integer(kind=iwp) :: I, iBas, icc, iCent, icentprev, icoeff, iComp, idummy(8), iIrrep, ijSO, ilcentprev, imcentprev, indx, indexi, &
                     indexj, iOff, iOff2, iOpt, iorb, ipSCR, iRC, irun, isame, iSmLbl, iSO, iSO_a, iSO_r, istatus, isymunit, &
                     iunit, j1, j12, j2, jcent, jcentprev, jlcentprev, jmcentprev, jrun, jsame, jSO, jSO_r, lauf, laufalt, &
                     length3, length3_tot, LenInt, LenTot, LLhigh, Lrun, Mrun, ncontcent(0:Lmax), not_defined, nSOs, numbofcent, &
                     numboffunct, numboffunct3, numbofsym
real(kind=wp) :: coeff, Sgn
!BS character(len=20) :: filename
character(len=8) :: xa2(4)
character(len=3) :: send
logical(kind=iwp) :: EX
integer(kind=iwp), allocatable :: C(:), ifirstLM(:,:,:), ip(:), iSO_info(:,:), Lcent(:), Lhighcent(:), Lval(:), Mcent(:), mval(:), &
                                  nadpt(:), ncent(:), nphase(:,:), numballcart(:)
real(kind=wp), allocatable :: AMFI_Int(:,:), Scr(:,:), SOInt(:)
integer(kind=iwp), external :: iPntSO, isfreeunit, n2Tri

! These variables are just placeholders for reading
#include "macros.fh"
unused_var(idummy)
unused_var(xa2)

!#######################################################################
call mma_allocate(ifirstLM,[0,Lmax],[-Lmax,Lmax],[1,MxAtom],label='ifirstLM')

! read information from SYMINFO
isymunit = isfreeunit(58)
call f_inquire('SYMINFO',EX)
if (.not. EX) call SysAbendMsg('systrafo','SYMINFO not present','Sorry')
call molcas_open(isymunit,'SYMINFO')
rewind(isymunit)
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) 'Symmetry adapation of the SO-integrals'
#endif
read(isymunit,*)
read(isymunit,*)
read(isymunit,*)
numboffunct = 0
send = ''
do while (send /= 'END')
  numboffunct = numboffunct+1
  read(isymunit,'(A3)') send
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'there are totally ',numboffunct,' functions'
#endif
if (numboffunct > MxOrb) call SysAbendMsg('symtrafo','increase MxOrb in Molcas.fh',' ')
rewind isymunit
read(isymunit,*)
read(isymunit,*)
numbofcent = 0
call mma_allocate(ncent,numboffunct,label='ncent')
call mma_allocate(Lval,numboffunct,label='Lval')
call mma_allocate(mval,numboffunct,label='mval')
call mma_allocate(nadpt,numboffunct,label='nadpt')
call mma_allocate(nphase,8,numboffunct,label='nphase')
do irun=1,numboffunct
  read(isymunit,*) indx,ncent(irun),lval(irun),mval(irun),nadpt(irun),(nphase(I,irun),I=1,nadpt(irun))
  numbofcent = max(numbofcent,ncent(irun))
  if (indx /= irun) call SysAbendMsg('symtrafo','weird numbering on SYMINFO',' ')
end do
close(iSymUnit)
#ifdef _DEBUGPRINT_
write(u6,*) 'number of unique centres',numbofcent
#endif

! clean up arrays for new integrals
numboffunct3 = (numboffunct*numboffunct+numboffunct)/2

call mma_allocate(AMFI_Int,numboffunct3,3,label='AMFI_Int')
AMFI_Int(:,:) = Zero

nSOs = 0
do iIrrep=0,nIrrep-1
  nSOs = nSOs+nBas(iIrrep)
end do
call mma_allocate(iSO_info,2,nSOs,label='iSO_info')
iSO_a = 0
do iIrrep=0,nIrrep-1
  iSO_r = 0
  do iBas=1,nBas(iIrrep)
    iSO_a = iSO_a+1
    iSO_r = iSO_r+1
    iSO_info(1,iSO_a) = iIrrep
    iSO_info(2,iSO_a) = iSO_r
  end do
end do

! loop over unique centres to read integrals and information

iunit = LUPROP
call mma_allocate(Scr,numboffunct3,3,label='Scr')
Scr(:,:) = Zero
ipSCR = 1
length3_tot = 0

call mma_allocate(C,numbofcent,label='C')
call mma_allocate(numballcart,numbofcent,label='numballcart')
call mma_allocate(Lcent,MxCart,label='Lcent')
call mma_allocate(Mcent,MxCart,label='Mcent')
call mma_allocate(Lhighcent,MxAtom,label='Lhighcent')

! In a MPI run not all atomic block will be available in
! all processes. Make up so we know later if a particular
! atom is present.

C(:) = -99
do jcent=1,numbofcent

# ifdef _DEBUGPRINT_
  write(u6,*) 'read integrals and info for centre ',jcent
# endif

  ! Note that when running in parallel this list is incomplete.
  ! Hence, we process the centers which each process host.

  read(iunit,iostat=istatus) iCent
  if (istatus < 0) exit
  read(iunit) xa2,numbofsym,(idummy(I),i=1,numbofsym),numballcart(icent),(Lcent(i),I=1,numballcart(icent)), &
              (mcent(i),I=1,numballcart(icent)),Lhighcent(icent),(ncontcent(I),I=0,Lhighcent(icent))
# ifdef _DEBUGPRINT_
  write(u6,*) numballcart(icent),'functions on centre ',icent
# endif
  length3 = iTri(numballcart(icent),numballcart(icent))
  C(iCent) = ipSCR
  read(iunit) (Scr(i,1),i=ipSCR,ipSCR+length3-1)
  read(iunit) xa2
  read(iunit) (Scr(i,2),i=ipSCR,ipSCR+length3-1)
  read(iunit) xa2
  read(iunit) (Scr(i,3),i=ipSCR,ipSCR+length3-1)
  ipScr = ipScr+length3
  length3_tot = length3_tot+length3
  !ulf
  ! check if any L-value is missing
  LLhigh = Lhighcent(icent)
  do i=1,Lhighcent(icent)
    if (ncontcent(I) == 0) LLhigh = LLhigh-1
  end do
  Lhighcent(icent) = LLhigh
  !bs determize where the first function of a special type is..
  not_defined = iTri(numboffunct,numboffunct)+1
  do Lrun=0,Lhighcent(icent)
    ifirstLM(Lrun,-Lrun:Lrun,icent) = not_defined
  end do
  do iorb=1,numballcart(icent)
    Lrun = Lcent(iorb)
    Mrun = Mcent(iorb)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iorb,Lrun,mrun',iorb,Lrun,mrun
#   endif
    ifirstLM(Lrun,Mrun,icent) = min(iorb,ifirstLM(Lrun,Mrun,icent))
  end do

  !bs determined..
  !bs check if all of them were found
  do Lrun=0,Lhighcent(icent)
    do Mrun=-Lrun,Lrun
      if (ifirstLM(Lrun,Mrun,icent) == not_defined) then
        write(u6,*) 'problems for centre,L,M ',icent,Lrun,Mrun
        call SysAbendMsg('symtrafo','problems with L- and M-values',' ')
      end if
    end do
  end do
end do    !end of loop over centres
call mma_deallocate(numballcart)
call mma_deallocate(Lcent)
call mma_deallocate(Mcent)
call mma_deallocate(Lhighcent)
#ifdef _DEBUGPRINT_
write(u6,*) 'length3_tot=',length3_tot
call RecPrt('SCR(1,1)',' ',Scr(1,1),1,length3_tot)
call RecPrt('SCR(1,2)',' ',Scr(1,2),1,length3_tot)
call RecPrt('SCR(1,3)',' ',Scr(1,3),1,length3_tot)
do iCent=1,numbofcent
  write(u6,*) C(iCent)
end do
#endif
! If this process does not have any blocks of integrals proceed
! directly to the distribution step.
if (Length3_tot /= 0) then

  !bs Finally the transformation!!!!

  icentprev = 0
  jcentprev = 0
  ilcentprev = -1
  jlcentprev = -1
  imcentprev = 20
  jmcentprev = 20
  isame = 1
  jsame = 1
  lauf = 0
  do irun=1,numboffunct
    ! Skip if center corresponding to this basis function is
    ! not available at this process.
    if ((ncent(irun) == icentprev) .and. (ilcentprev == lval(irun)) .and. (imcentprev == mval(irun))) then
      isame = isame+1
    else
      isame = 1
      icentprev = ncent(irun)
      ilcentprev = lval(irun)
      imcentprev = mval(irun)
    end if
    do jrun=1,irun
      lauf = lauf+1
      if ((ncent(jrun) == jcentprev) .and. (jlcentprev == lval(jrun)) .and. (jmcentprev == mval(jrun))) then
        jsame = jsame+1
      else
        jsame = 1
        jcentprev = ncent(jrun)
        jlcentprev = lval(jrun)
        jmcentprev = mval(jrun)
      end if
      !bs check for same centers
      if (ncent(irun) == ncent(jrun)) then
        if ((lval(irun) == lval(jrun)) .and. (lval(irun) > 0)) then
          if (abs(abs(mval(irun))-abs(mval(jrun))) <= 1) then

            !bs the only cases where non-zero integrals occur
            if (nadpt(irun) == 1) then
              coeff = One
            else
              icoeff = 0
              do icc=1,nadpt(irun)
                icoeff = icoeff+nphase(icc,irun)*nphase(icc,jrun)
              end do
              coeff = real(icoeff,kind=wp)
              if (MolWgh == 2) then
                coeff = coeff/real(nadpt(irun),kind=wp)
              else
                coeff = coeff/real(nadpt(irun)*nadpt(irun),kind=wp)
              end if
            end if
            !bs determine indices of atomic integrals
            indexi = ifirstLM(lval(irun),mval(irun),ncent(irun))+isame-1
            indexj = ifirstLM(lval(irun),mval(jrun),ncent(irun))+jsame-1
            laufalt = iTri(indexi,indexj)

            if (C(nCent(iRun)) /= -99) then
              ipSCR = C(ncent(irun))-1+laufalt
              ! DebugDebug
              !write(u6,*) 'laufalt=',laufalt
              !write(u6,*) 'ip''s:',ipSCR
              !write(u6,*) Scr(ipSCR,1),Scr(ipSCR,2),Scr(ipSCR,3)
              ! DebugDebug
              Sgn = One
              if (indexi > indexj) Sgn = -Sgn
              AMFI_Int(lauf,1) = Sgn*coeff*Scr(ipScr,1)
              AMFI_Int(lauf,2) = Sgn*coeff*Scr(ipScr,2)
              AMFI_Int(lauf,3) = Sgn*coeff*Scr(ipScr,3)
            end if

          end if
        end if
      end if
    end do ! jrun
  end do ! irun
end if
call mma_deallocate(C)
call mma_deallocate(ncent)
call mma_deallocate(Lval)
call mma_deallocate(mval)
call mma_deallocate(nadpt)
call mma_deallocate(nphase)

! Allocate memory for symmetry adapted one electron integrals.
! Will just store the unique elements, i.e. low triangular blocks
! and lower triangular elements in the diagonal blocks.

call mma_allocate(ip,nComp,label='ip')
LenTot = 0
do iComp=1,nComp
  ip(iComp) = 1+LenTot
  LenInt = n2Tri(lOper(iComp))
  LenTot = LenTot+LenInt+4
end do
call mma_allocate(SOInt,LenTot,label='SOInt')
SOInt(:) = Zero

! This test is not valid for parallel execution, since here
! we have only incomplete lists.
!if (lauf /= numboffunct3) call SysAbendMsg('symtrafo','error in numbering ',' ')
do iComp=1,nComp
  do iSO=1,numboffunct
    j1 = iSO_info(1,iSO)
    iSO_r = iSO_info(2,iSO)
    do jSO=1,iSO
      j2 = iSO_info(1,jSO)
      jSO_r = iSO_info(2,jSO)
      j12 = ieor(j1,j2)
      if (.not. btest(lOper(iComp),j12)) cycle

      iOff = iPntSO(max(j1,j2),min(j1,j2),lOper(iComp),nBas)
      iOff = iOff+ip(iComp)-1

      iOff2 = iTri(iSO,jSO)

      if (j1 == j2) then
        ijSO = iTri(iSO_r,jSO_r)
      else
        ijSO = (jSO_r-1)*nBas(j1)+iSO_r
      end if
      SOInt(iOff+ijSO) = -AMFI_Int(iOff2,iComp)

    end do
  end do

  ! Write out integrals to ONEINT for this specific component of the operator.

  iOpt = 0
  iRC = -1
  iSmLbl = lOper(iComp)
  call GADSum(SOInt(ip(iComp)),n2Tri(iSmLbl))
  call WrOne(iRC,iOpt,Label,iComp,SOInt(ip(iComp)),iSmLbl)
  if (iRC /= 0) call SysAbendMsg('symtrafo','     Error in subroutine ONEEL ','     Abend in subroutine WrOne')

end do

#ifdef _DEBUGPRINT_
call PrMtrx(Label,lOper,nComp,ip,SOInt)
#endif

call mma_deallocate(ifirstLM)
call mma_deallocate(iSO_info)
call mma_deallocate(Scr)
call mma_deallocate(AMFI_Int)
call mma_deallocate(ip)
call mma_deallocate(SOInt)
!BS write(u6,*) 'Symmetry transformation successfully done'

return

end subroutine SymTrafo
