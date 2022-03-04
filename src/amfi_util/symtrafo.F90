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

subroutine SymTrafo(LUPROP,ip,lOper,nComp,nBas,nIrrep,Label,MolWgh,SOInt,LenTot)
!bs   Purpose: combine SO-integrals from amfi to symmetry-adapted
!bs   integrals on one file AOPROPER_MF_SYM

implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: AMFI_Int(:,:), Scr(:,:)
integer, allocatable :: iSO_info(:,:)
parameter(maxorbs=MxOrb)
parameter(maxcent=MxAtom)
real*8 SOInt(LenTot)
character*8 ya, za, xa2
character*3 end
!BS character*20 filename
logical EX
!BS namelist /SYMTRA/ none
dimension ya(4), za(4)
dimension xa2(4)
dimension ncent(maxorbs), Lval(maxorbs), mval(maxorbs), nadpt(maxorbs), nphase(8,maxorbs), idummy(8), Lhighcent(maxcent), &
          Lcent(MxCart), Mcent(MxCart), ncontcent(0:Lmax), numballcart(maxcent)
allocatable ifirstLM(:,:,:)
integer ip(nComp), nBas(0:nIrrep-1), lOper(nComp), ipC(MxAtom)
character Label*8
!Statement function
IPNT(I,J) = (max(i,j)*max(i,j)-max(i,j))/2+min(i,j)

!#######################################################################
end = '   '  ! added due to cray warnings. B.S. 04/10/04
ya(1) = '********'
za(1) = '********'
ya(2) = '        '
Za(2) = '        '
ya(3) = 'ANTISYMM'
Za(3) = 'ANTISYMM'
ya(4) = 'Y1SPNORB'
ZA(4) = 'Z1SPNORB'
call mma_allocate(ifirstLM,[0,Lmax],[-Lmax,Lmax],[1,maxcent],label='ifirstLM')

! read information from SYMINFO
isymunit = isfreeunit(58)
call f_inquire('SYMINFO',EX)
if (.not. EX) call SysAbendMsg('systrafo','SYMINFO not present','Sorry')
call molcas_open(isymunit,'SYMINFO')
rewind(isymunit)
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(6,*) 'Symmetry adapation of the SO-integrals'
#endif
read(isymunit,*)
read(isymunit,*)
read(isymunit,*)
numboffunct = 0
do while (end /= 'END')
  numboffunct = numboffunct+1
  read(isymunit,'(A3)') end
end do
#ifdef _DEBUGPRINT_
write(6,*) 'there are totally ',numboffunct,' functions'
#endif
if (numboffunct > maxorbs) call SysAbendMsg('symtrafo','increase maxorbs in symtrafo',' ')
rewind isymunit
read(isymunit,*)
read(isymunit,*)
numbofcent = 0
do irun=1,numboffunct
  read(isymunit,*) index,ncent(irun),lval(irun),mval(irun),nadpt(irun),(nphase(I,irun),I=1,nadpt(irun))
  numbofcent = max(numbofcent,ncent(irun))
  if (index /= irun) call SysAbendMsg('symtrafo','weird numbering  on SYMINFO',' ')
end do
close(iSymUnit)
#ifdef _DEBUGPRINT_
write(6,*) 'number of unique centres',numbofcent
#endif

! clean up arrays for new integrals
numboffunct3 = (numboffunct*numboffunct+numboffunct)/2

call mma_allocate(AMFI_Int,numboffunct3,3,Label='AMFI_Int')
AMFI_Int(:,:) = Zero

nSOs = 0
do iIrrep=0,nIrrep-1
  nSOs = nSOs+nBas(iIrrep)
end do
call mma_allocate(iSO_info,2,nSOs,Label='iSO_info')
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
call mma_allocate(Scr,numboffunct3,3,Label='Scr')
Scr(:,:) = Zero
ipSCR = 1
length3_tot = 0

! In a MPI run not all atomic block will be available in
! all processes. Make up so we know later if a particular
! atom is present.

ipC(1:numbofcent) = -99
do jcent=1,numbofcent

# ifdef _DEBUGPRINT_
  write(6,*) 'read integrals and info for centre ',jcent
# endif

  ! Note that when running in parallel this list is incomplete.
  ! Hence, we process the centers which each process host.

  read(iunit,end=199) iCent
  read(iunit) xa2,numbofsym,(idummy(I),i=1,numbofsym),numballcart(icent),(Lcent(i),I=1,numballcart(icent)), &
              (mcent(i),I=1,numballcart(icent)),Lhighcent(icent),(ncontcent(I),I=0,Lhighcent(icent))
# ifdef _DEBUGPRINT_
  write(6,*) numballcart(icent),'functions on centre ',icent
# endif
  length3 = ipnt(numballcart(icent),numballcart(icent))
  ipC(iCent) = ipSCR
  read(iunit) (Scr(i,1),i=ipSCR,ipSCR+length3-1)
  read(iunit) Ya
  read(iunit) (Scr(i,2),i=ipSCR,ipSCR+length3-1)
  read(iunit) Za
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
  not_defined = ipnt(numboffunct,numboffunct)+1
  do Lrun=0,Lhighcent(icent)
    do Mrun=-Lrun,Lrun
      ifirstLM(Lrun,Mrun,icent) = not_defined
    end do
  end do
  do iorb=1,numballcart(icent)
    Lrun = Lcent(iorb)
    Mrun = Mcent(iorb)
#   ifdef _DEBUGPRINT_
    write(6,*) 'iorb,Lrun,mrun',iorb,Lrun,mrun
#   endif
    ifirstLM(Lrun,Mrun,icent) = min(iorb,ifirstLM(Lrun,Mrun,icent))
  end do

  !bs determined..
  !bs check if all of them were found
  do Lrun=0,Lhighcent(icent)
    do Mrun=-Lrun,Lrun
      if (ifirstLM(Lrun,Mrun,icent) == not_defined) then
        write(6,*) 'problems for centre,L,M ',icent,Lrun,Mrun
        call SysAbendMsg('symtrafo','problems with L- and M-values',' ')
      end if
    end do
  end do
end do    !end of loop over centres
199 continue
#ifdef _DEBUGPRINT_
write(6,*) 'length3_tot=',length3_tot
call RecPrt('SCR(1,1)',' ',Scr(1,1),1,length3_tot)
call RecPrt('SCR(1,2)',' ',Scr(1,2),1,length3_tot)
call RecPrt('SCR(1,3)',' ',Scr(1,3),1,length3_tot)
do iCent=1,numbofcent
  write(6,*) ipC(iCent)
end do
#endif
! If this process does not have any blocks of integrals proceed
! directly to the distribution step.
if (Length3_tot == 0) Go To 299

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
            coeff = 1d0
          else
            icoeff = 0
            do icc=1,nadpt(irun)
              icoeff = icoeff+nphase(icc,irun)*nphase(icc,jrun)
            end do
            coeff = dble(icoeff)
            if (MolWgh == 2) then
              coeff = coeff/dble(nadpt(irun))
            else
              coeff = coeff/dble(nadpt(irun)*nadpt(irun))
            end if
          end if
          !bs determine indices of atomic integrals
          indexi = ifirstLM(lval(irun),mval(irun),ncent(irun))+isame-1
          indexj = ifirstLM(lval(irun),mval(jrun),ncent(irun))+jsame-1
          laufalt = ipnt(indexi,indexj)

          if (ipC(nCent(iRun)) /= -99) then
            ipSCR = ipC(ncent(irun))-1+laufalt
            ! DebugDebug
            !write(6,*) 'laufalt=',laufalt
            !write(6,*) 'ip''s:',ipSCR
            !write(6,*) Scr(ipSCR,1),Scr(ipSCR,2),Scr(ipSCR,3)
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
299 continue
! This test is not valid for parallel execusion, since here
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
      if (iand(lOper(iComp),2**j12) == 0) Go To 99

      iOff = iPntSO(max(j1,j2),min(j1,j2),lOper(iComp),nBas)
      iOff = iOff+ip(iComp)-1

      iOff2 = iPnt(iSO,jSO)

      tmp = -AMFI_Int(iOff2,iComp)
      if (j1 == j2) then
        ijSO = iPnt(iSO_r,jSO_r)
      else
        ijSO = (jSO_r-1)*nBas(j1)+iSO_r
      end if
      SOInt(iOff+ijSO) = tmp

99    continue
    end do
  end do

  ! Write out integrals to ONEINT for this specific component of the operator.

  iOpt = 0
  iRC = -1
  iSmLbl = lOper(iComp)
  call GADSum(SOInt(ip(iComp)),n2Tri(iSmLbl))
  call WrOne(iRC,iOpt,Label,iComp,SOInt(ip(iComp)),iSmLbl)
  if (iRC /= 0) then
    call SysAbendMsg('symtrafo','     Error in subroutine ONEEL ','     Abend in subroutine WrOne')
  end if

end do

#ifdef _DEBUGPRINT_
call PrMtrx(Label,lOper,nComp,ip,SOInt)
#endif

call mma_deallocate(ifirstLM)
call mma_deallocate(iSO_info)
call mma_deallocate(Scr)
call mma_deallocate(AMFI_Int)
!BS write(6,*) 'Symmetry transformation successfully done'

return
#ifdef _WARNING_WORKAROUND_
if (.false.) then
  call Unused_integer_array(idummy)
  call Unused_character(ya)
  call Unused_character(za)
  call Unused_character(xa2)
end if
#endif

end subroutine SymTrafo
