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

subroutine Get_Density_Matrix(ip_D,nBas1,nBas2,nBasMax,nBas,nSym,ipP,UserDen,PrintDen,SubtractDen,SubScale,Q_Nuc,nAtoms,iPert, &
                              Restart,Utility,TDensity,nStateI,nStateF)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
! OBSERVE! MxState has to be the same as MxStat in cntrl.fh in src/rassi.
parameter(MxState=200,MxStOT=MxState*(MxState+1)/2)
integer nBas(nSym)
character*16 Label, filename
dimension Q_Nuc(nAtoms)
dimension iToc(MxStOT)
logical UserDen, PrintDen, SubtractDen, Exist, Restart, Utility
logical TDensity, ok1, ok2, Found
character*8 Method
real*8, allocatable :: DTmp(:), DSym(:)

write(Label,'(A,I1)') 'LoProp Dens ',iPert
if (Restart) then

  ! Retrieve density matrix from the runfile

  call qpg_dArray(Label,Exist,nDens)
  if ((.not. Exist) .or. (nDens == 0)) then
    call SysAbendMsg('get_density_matrix','Could not locate:',Label)
  end if
  call Allocate_Work(ip_D,nDens)
  call Get_dArray(Label,Work(ip_D),nDens)

  nSize = nDens

else if (nSym == 1) then

  nSize = nBas(1)*(nBas(1)+1)/2

  if (UserDen) then
    ! Additions made by A.Ohrn to enable LoProp to read in a density matrix
    ! provided by the user. Generalized by P. Soderhjelm so that also
    ! perturbed density matrices may be read for polarizability calculation
    Lu_Ud = 56
    Lu_Ud = IsFreeUnit(Lu_Ud)
    if (ipert == 0) then
      filename = 'USERDEN'
    else
      write(filename,'(A7,I1)') 'USERDEN',ipert
    end if
    call OpnFl(filename,Lu_Ud,Exist)
    if (.not. Exist) then
      write(6,*)
      write(6,*) ' Unable to locate user density matrix.'
      call Abend()
    end if
    call GetMem('UserDen','Allo','Real',ipUser,nSize)
    read(Lu_Ud,*) (Work(ipUser+k),k=0,nSize-1)
    call Put_D1ao(Work(ipUser),nSize)
    call GetMem('UserDen','Free','Real',ipUser,nSize)
    close(Lu_Ud)
  end if
  ! Addition of A.Ohrn, to collect transition density matrix from ToFile.
  if (TDensity) then
    LuIn = 57
    LuIn = IsFreeUnit(LuIn)
    call DaName(LuIn,'TOFILE')
    iDisk = 0
    ! Table-of-contents
    call iDaFile(LuIn,2,iToc,MxStOT,iDisk)
    ! Allocation of density
    call GetMem('TDMden','Allo','Real',iTDMden,nSize)
    ! Loop to 'suck-out' the relevant matrix from the ToFile.
    nStateM = max(nStateI,nStateF)
    kaunter = 0
    do iS1=1,nStateM
      do iS2=1,iS1
        kaunter = kaunter+1
        iDisk = iToc(kaunter)
        call dDaFile(LuIn,2,Work(iTDMden),nSize,iDisk)
        ok1 = (iS1 == nStateI) .and. (iS2 == nStateF)
        ok2 = (iS1 == nStateF) .and. (iS2 == nStateI)
        if (ok1 .or. ok2) then
          call Put_D1ao(Work(iTDMden),nSize)
        end if
      end do
    end do
    call GetMem('TDMden','Free','Real',iTDMden,nSize)
    call DaClos(LuIn)
  end if
  ! End addition A.Ohrn.

  if (SubtractDen) then
    ! Additions made by P.Soderhjelm to enable LoProp to read in a density matrix
    ! provided by the user and subtract that from the current one (which could
    ! in fact be provided by the Userdens keyword. After subtraction the matrix
    ! is scaled by a constant SubScale. Nuclear charges are set to zero.
    Lu_Ud = 56
    Lu_Ud = IsFreeUnit(Lu_Ud)
    call OpnFl('SUBDEN',Lu_Ud,Exist)
    if (.not. Exist) then
      write(6,*)
      write(6,*) ' Unable to locate density matrix to subtract.'
      call Abend()
    end if
    call GetMem('UserDen','Allo','Real',ipUser,nSize)
    read(Lu_Ud,*) (Work(ipUser+k),k=0,nSize-1)
    call mma_allocate(DTmp,nSize)
    call Get_D1ao(Dtmp,nSize)
    do k=0,nSize-1
      Dtmp(1+k) = SubScale*(Dtmp(1+k)-Work(ipUser+k))
    end do
    call Put_D1ao(Dtmp,nSize)
    call mma_deallocate(DTmp)
    call GetMem('UserDen','Free','Real',ipUser,nSize)
    close(Lu_Ud)
    do k=1,nAtoms
      Q_Nuc(k) = 0.0d0
    end do
  end if
  ! End addition P.Soderhjelm.

  ! Addition by J.Bostrom to check if loprop should pick
  ! the variational density instead.
  call Get_cArray('Relax Method',Method,8)
  iMp2Prpt = 0
  if (Method == 'MBPT2   ') then
    call Get_iScalar('mp2prpt',iMp2Prpt)
  end if
  call getmem('D','Allo','Real',ip_D,nSize)
  if (iMp2Prpt /= 0) then
    call Get_D1ao_var(Work(ip_D),nSize)
  else
  ! End Addition J.Bostrom (well, the End If too ofc.)
    call Get_D1ao(Work(ip_D),nSize)
  end if
  nDens = nSize
# ifdef _DEBUGPRINT_
  call RecPrt('D',' ',Work(ip_D),1,nDens)
# endif
  if (PrintDen) then
    ! Addition made by P.Soderhjelm to enable LoProp to print the density matrix
    ! to a text file for later use.
    Lu_Ud = 56
    Lu_Ud = IsFreeUnit(Lu_Ud)
    call OpnFl('PRDEN',Lu_Ud,Exist)
    write(Lu_Ud,'(10d25.16)') (Work(ip_D+k),k=0,nDens-1)
    close(Lu_Ud)
  end if

else

  call Allocate_Work(ip_D,nBas1*(nBas1+1)/2)
  call Allocate_Work(ip_D_sq,nBas1**2)
  call Allocate_Work(ip_Tmp,nBas2)

  call Qpg_darray('D1ao',Found,nDens)
  if (Found .and. (nDens /= 0)) then
    call mma_allocate(DSym,nDens,Label='DSym')
    call Get_D1ao(DSym,nDens)
  else
    write(6,*) 'Get_density_matrix: not found.'
    call Abend()
  end if
  iSyLbl = 1
# ifdef _DEBUGPRINT_
  call RecPrt('DSym',' ',DSym,1,nDens)
# endif
  iOfft = 1
  iOffs = ip_Tmp
  do iSym=1,nSym
    if (nBas(iSym) == 0) Go To 99
#   ifdef _DEBUGPRINT_
    call TriPrt('DSym',' ',DSym(iOfft),nBas(iSym))
#   endif
    call Square(DSym(iOfft),Work(iOffs),1,nBas(iSym),nBas(iSym))
    call DScal_(nBas(iSym)**2,Half,Work(iOffs),1)
    call DScal_(nBas(iSym),Two,Work(iOffs),nBas(iSym)+1)
#   ifdef _DEBUGPRINT_
    call RecPrt('DSym',' ',Work(iOffs),nBas(iSym),nBas(iSym))
#   endif
    iOfft = iOfft+nBas(iSym)*(nBas(iSym)+1)/2
    iOffs = iOffs+nBas(iSym)**2
99  continue
  end do
  call mma_deallocate(DSym)

  nScr = nBasMax*nBas1
  call Allocate_Work(ipScr,nScr)
  call Desymmetrize(Work(ip_Tmp),nBas2,Work(ipScr),nScr,Work(ip_D_sq),nBas,nBas1,Work(ipP),nSym,iSyLbl)
  call Free_Work(ipScr)
  call Free_Work(ip_Tmp)

  call Triangularize(Work(ip_D_sq),Work(ip_D),nBas1,.true.)
  call Free_Work(ip_D_sq)
end if
#ifdef _DEBUGPRINT_
call TriPrt('Density Matrix',' ',Work(ip_D),nBas1)
#endif

! Copy Density Matrix to the runfile

if ((.not. Restart) .and. (.not. Utility)) then
  !call put_dArray(Label,Work(ip_D),nDens)
  call put_dArray(Label,Work(ip_D),nBas1*(nBas1+1)/2)
end if

return

end subroutine Get_Density_Matrix
