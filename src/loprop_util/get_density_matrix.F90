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

subroutine Get_Density_Matrix(D,nBas1,nBas2,nBasMax,nBas,nSym,P,UserDen,PrintDen,SubtractDen,SubScale,Q_Nuc,nAtoms,iPert,Restart, &
                              Utility,TDensity,nStateI,nStateF)

use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
type(Alloc1DArray_Type), intent(out) :: D
integer(kind=iwp), intent(in) :: nBas1, nBas2, nBasMax, nSym, nBas(nSym), nAtoms, iPert, nStateI, nStateF
real(kind=wp), intent(in) :: P(*), SubScale
logical(kind=iwp), intent(in) :: UserDen, PrintDen, SubtractDen, Restart, Utility, TDensity
real(kind=wp), intent(inout) :: Q_Nuc(nAtoms)
integer(kind=iwp) :: iDisk, iMp2Prpt, iOffs, iOfft, iS1, iS2, iSyLbl, iSym, k, kaunter, Lu_Ud, LuIn, nDens, nScr, nSize, nStateM
character(len=16) :: Label, filename
character(len=8) :: Method
logical(kind=iwp) :: Exists, Found, ok1, ok2
integer(kind=iwp), allocatable :: iToc(:)
real(kind=wp), allocatable :: D_sq(:,:), DTmp(:), DSym(:), Scr(:), TDMden(:), Tmp(:), User(:)
! OBSERVE! MxState has to be the same as MxStat in cntrl.fh in src/rassi.
!          (but there's no such variable)
integer(kind=iwp), parameter :: MxState = 200, MxStOT = MxState*(MxState+1)/2
integer(kind=iwp), external :: IsFreeUnit

write(Label,'(A,I1)') 'LoProp Dens ',iPert
if (Restart) then

  ! Retrieve density matrix from the runfile

  call qpg_dArray(Label,Exists,nDens)
  if ((.not. Exists) .or. (nDens == 0)) then
    call SysAbendMsg('get_density_matrix','Could not locate:',Label)
  end if
  call mma_allocate(D%A,nDens,label='D')
  call Get_dArray(Label,D%A,nDens)

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
    call OpnFl(filename,Lu_Ud,Exists)
    if (.not. Exists) then
      write(u6,*)
      write(u6,*) ' Unable to locate user density matrix.'
      call Abend()
    end if
    call mma_allocate(User,nSize,label='UserDen')
    read(Lu_Ud,*) User(:)
    call Put_dArray('D1ao',User,nSize)
    call mma_deallocate(User)
    close(Lu_Ud)
  end if
  ! Addition of A.Ohrn, to collect transition density matrix from ToFile.
  if (TDensity) then
    LuIn = 57
    LuIn = IsFreeUnit(LuIn)
    call DaName(LuIn,'TOFILE')
    iDisk = 0
    ! Table-of-contents
    call mma_allocate(iToc,MxStOT,label='iToc')
    call iDaFile(LuIn,2,iToc,MxStOT,iDisk)
    ! Allocation of density
    call mma_allocate(TDMden,nSize,label='TDMden')
    ! Loop to 'suck-out' the relevant matrix from the ToFile.
    nStateM = max(nStateI,nStateF)
    kaunter = 0
    do iS1=1,nStateM
      do iS2=1,iS1
        kaunter = kaunter+1
        iDisk = iToc(kaunter)
        call dDaFile(LuIn,2,TDMden,nSize,iDisk)
        ok1 = (iS1 == nStateI) .and. (iS2 == nStateF)
        ok2 = (iS1 == nStateF) .and. (iS2 == nStateI)
        if (ok1 .or. ok2) then
          call Put_dArray('D1ao',TDMden,nSize)
        end if
      end do
    end do
    call mma_deallocate(iToc)
    call mma_deallocate(TDMden)
    call DaClos(LuIn)
  end if
  ! End addition A.Ohrn.

  ! Additions made by P.Soderhjelm to enable LoProp to read in a density matrix
  ! provided by the user and subtract that from the current one (which could
  ! in fact be provided by the Userdens keyword. After subtraction the matrix
  ! is scaled by a constant SubScale. Nuclear charges are set to zero.
  if (SubtractDen) then
    Lu_Ud = 56
    Lu_Ud = IsFreeUnit(Lu_Ud)
    call OpnFl('SUBDEN',Lu_Ud,Exists)
    if (.not. Exists) then
      write(u6,*)
      write(u6,*) ' Unable to locate density matrix to subtract.'
      call Abend()
    end if
    call mma_allocate(User,nSize,label='UserDen')
    read(Lu_Ud,*) User(:)
    call mma_allocate(DTmp,nSize)
    call Get_dArray_chk('D1ao',Dtmp,nSize)
    Dtmp(:) = SubScale*(Dtmp(:)-User(:))
    call Put_dArray('D1ao',Dtmp,nSize)
    call mma_deallocate(DTmp)
    call mma_deallocate(User)
    close(Lu_Ud)
    do k=1,nAtoms
      Q_Nuc(k) = Zero
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
  call mma_allocate(D%A,nSize,label='D')
  if (iMp2Prpt /= 0) then
    call Get_D1ao_var(D%A,nSize)
  else
  ! End Addition J.Bostrom (well, the End If too ofc.)
    call Get_dArray_chk('D1ao',D%A,nSize)
  end if
  nDens = nSize
# ifdef _DEBUGPRINT_
  call RecPrt('D',' ',D%A,1,nDens)
# endif
  if (PrintDen) then
    ! Addition made by P.Soderhjelm to enable LoProp to print the density matrix
    ! to a text file for later use.
    Lu_Ud = 56
    Lu_Ud = IsFreeUnit(Lu_Ud)
    call OpnFl('PRDEN',Lu_Ud,Exists)
    write(Lu_Ud,'(10F25.16)') D%A(:)
    close(Lu_Ud)
  end if

else

  call mma_allocate(D%A,nBas1*(nBas1+1)/2,label='D')
  call mma_allocate(D_sq,nBas1,nBas1,label='D_sq')
  call mma_allocate(Tmp,nBas2,label='Tmp')

  call Qpg_darray('D1ao',Found,nDens)
  if (Found .and. (nDens /= 0)) then
    call mma_allocate(DSym,nDens,Label='DSym')
    call Get_dArray_chk('D1ao',DSym,nDens)
  else
    write(u6,*) 'Get_density_matrix: not found.'
    call Abend()
  end if
  iSyLbl = 1
# ifdef _DEBUGPRINT_
  call RecPrt('DSym',' ',DSym,1,nDens)
# endif
  iOfft = 1
  iOffs = 1
  do iSym=1,nSym
    if (nBas(iSym) == 0) cycle
#   ifdef _DEBUGPRINT_
    call TriPrt('DSym',' ',DSym(iOfft),nBas(iSym))
#   endif
    call Square(DSym(iOfft),Tmp(iOffs),1,nBas(iSym),nBas(iSym))
    call DScal_(nBas(iSym)**2,Half,Tmp(iOffs),1)
    call DScal_(nBas(iSym),Two,Tmp(iOffs),nBas(iSym)+1)
#   ifdef _DEBUGPRINT_
    call RecPrt('DSym',' ',Tmp(iOffs),nBas(iSym),nBas(iSym))
#   endif
    iOfft = iOfft+nBas(iSym)*(nBas(iSym)+1)/2
    iOffs = iOffs+nBas(iSym)**2
  end do
  call mma_deallocate(DSym)

  nScr = nBasMax*nBas1
  call mma_allocate(Scr,nScr,label='Scr')
  call Desymmetrize(Tmp,nBas2,Scr,nScr,D_sq,nBas,nBas1,P,nSym,iSyLbl)
  call mma_deallocate(Scr)
  call mma_deallocate(Tmp)

  call Triangularize(D_sq,D%A,nBas1,.true.)
  call mma_deallocate(D_sq)
end if
#ifdef _DEBUGPRINT_
call TriPrt('Density Matrix',' ',D%A,nBas1)
#endif

! Copy Density Matrix to the runfile

if ((.not. Restart) .and. (.not. Utility)) then
  !call put_dArray(Label,D%A,nDens)
  call put_dArray(Label,D%A,nBas1*(nBas1+1)/2)
end if

return

end subroutine Get_Density_Matrix
