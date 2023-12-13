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

subroutine Bond_List(nq,nsAtom,iIter,nIter,Cx,Process,Valu,nB,qLbl,fconst,rMult,LuIC,Indq,Proc_dB,iTabBonds,nBonds,iTabAI,mAtoms, &
                     mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: ANr, AtomLbl, Fragments_Bond, iOptC, jStab, Magic_Bond, nStab, vdW_Bond
use ddvdt, only: A_StrH, aAV, alpha_vdW, f_Const_Min, r_ref_vdW, rAV, rkr, rkr_vdW
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nsAtom, iIter, nIter, nB, LuIC, nBonds, iTabBonds(3,nBonds), mAtoms, iTabAI(2,mAtoms), nB_Tot, &
                                 ndB_Tot
integer(kind=iwp), intent(inout) :: nq, Indq(3,nB), mB_Tot, mdB_Tot, iBM(nB_Tot), idBM(2,ndB_Tot), mqB(nB)
real(kind=wp), intent(in) :: Cx(3,nsAtom,nIter)
logical(kind=iwp), intent(in) :: Process, Proc_dB
real(kind=wp), intent(inout) :: Valu(nB,nIter), fconst(nB), rMult(nB), BM(nB_Tot), dBM(ndB_Tot)
character(len=14), intent(inout) :: qLbl(nB)
#include "Molcas.fh"
integer(kind=iwp), parameter :: mB = 2*3
integer(kind=iwp) :: iAtom, iAtom_, iBond, iBondType, iCase, iDeg, iDCR(2), iDCRR(0:7), iE1, iE2, iF1, iF2, Ind(2), iRow, &
                     iStabM(0:7), jAtom, jAtom_, jRow, kDCRR, Lambda, nCent, nDCRR, nqB, nStabM
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: A(3,2), Alpha, Deg, f_Const, Grad(mB), Hess(mB**2), r0, Rab, RabCov, Rij2, Val
logical(kind=iwp) :: Help
character(len=14) :: Label
character(len=LenIn4) :: Lbls(2)
integer(kind=iwp), parameter :: iChOp(0:7) = [1,1,1,2,1,2,2,3]
character(len=*), parameter :: ChOp(0:7) = ['E  ','X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ']
real(kind=wp), external :: CovRad
logical(kind=iwp), external :: R_Stab_A

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (nBonds < 1) return

nqB = 0
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ---> Enter Bonds.'
write(u6,*)
write(u6,*) 'Process=',Process
call RecPrt('CX',' ',CX,3*nsAtom,nIter)
write(u6,'(20(1X,A))') (AtomLbl(i),i=1,nsAtom)
write(u6,*)
write(u6,*) ' iTabAI'
write(u6,*)
do iAtom=1,mAtoms
  write(u6,*) iTabAI(1,iAtom),iTabAI(2,iAtom)
end do
#endif

! Loop over bonds

nCent = 2
do iBond=1,nBonds
  iBondType = iTabBonds(3,iBond)

  ! We will only incorpotate covalent and fragment bonds

  if (iBondType == vdW_Bond) cycle   ! vdW bonds
  if (iBondType > Magic_Bond) cycle  ! magic bonds

  do iCase=1,2

    if (iCase == 1) then
      iAtom_ = iTabBonds(1,iBond)
      jAtom_ = iTabBonds(2,iBond)
    else
      iAtom_ = iTabBonds(2,iBond)
      jAtom_ = iTabBonds(1,iBond)
    end if
    iAtom = iTabAI(1,iAtom_)
    jAtom = iTabAI(1,jAtom_)

    iDCR(1) = iTabAI(2,iAtom_)
    iDCR(2) = iTabAI(2,jAtom_)
    if (jAtom > iAtom) cycle
    if (iDCR(1) /= iOper(0)) cycle
    if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(2) /= iOper(0))) cycle
    iRow = max(ANr(iAtom),1)
    jRow = max(ANr(jAtom),1)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iAtom,jAtom=',iAtom,jAtom
#   endif
    Help = (iRow > 3) .or. (jRow > 3)
    Ind(1) = iAtom
    Ind(2) = jAtom
    A(:,1) = Cx(:,iAtom,iIter)
    write(Label,'(A,I2,A,I2,A)') 'B(',iAtom,',',jAtom,')'

#   ifdef _DEBUGPRINT_
    call RecPrt('A',' ',Cx(:,iAtom,iIter),1,3)
    call RecPrt('B',' ',Cx(:,jAtom,iIter),1,3)
#   endif

    ! Form double coset representatives

    call DCR(Lambda,jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCRR,nDCRR)
    kDCRR = iDCR(2)

#   ifdef _DEBUGPRINT_
    write(u6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
    write(u6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
    write(u6,'(10A)') 'R={',(ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
    write(u6,'(2A)') 'R=',ChOp(iDCR(2))
#   endif

    call OA(iDCR(2),Cx(:,jAtom,iIter),A(:,2))

    ! Compute the stabilizer of A & R(B), this is done in two ways.
    !
    ! A=/=B, the stabilizer is formed as the intersection of
    !        the stabilizers of A and B.
    !
    ! A=B, the stabilizer is formed as union of U and R(U)

    if (iAtom /= jAtom) then
      call Inter(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iStabM,nStabM)
    else
      call Union(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),kDCRR,iStabM,nStabM)
    end if
#   ifdef _DEBUGPRINT_
    write(u6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#   endif

    ! Now evaluate the degeneracy of the bond.

    iDeg = nIrrep/nStabM
    Deg = sqrt(real(iDeg,kind=wp))
#   ifdef _DEBUGPRINT_
    write(u6,*) ' nIrrep,nStabM=',nIrrep,nStabM
#   endif

    nq = nq+1
    if (.not. Process) mB_Tot = mB_Tot+mB
    if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2

    nqB = nqB+1
    iF1 = 1
    call NxtWrd(AtomLbl(iAtom),iF1,iE1)
    Lbls(1) = AtomLbl(iAtom)(iF1:iE1)
    iF2 = 1
    call NxtWrd(AtomLbl(jAtom),iF2,iE2)
    Lbls(2) = AtomLbl(jAtom)(iF2:iE2)
    if (kDCRR /= 0) then
      Lbls(2)(iE2+1:iE2+1) = '('
      Lbls(2)(iE2+2:iE2+1+iChOp(kDCRR)) = ChOp(kDCRR)(1:iChOp(kDCRR))
      Lbls(2)(iE2+2+iChOp(kDCRR):iE2+2+iChOp(kDCRR)) = ')'
      call NxtWrd(Lbls(2),iF2,iE2)
    end if
    write(LuIC,'(A,I3.3,4A)') 'b',nqB,' = Bond ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,I3.3,4A)') 'b',nqB,' = Bond ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2)
#   endif
    Label = ' '
    write(Label,'(A,I3.3)') 'b',nqB
    if (.not. Proc_dB) Hess(:) = Zero
    call Strtch(A,nCent,Val,Grad,.false.,'        ',Hess,Proc_dB)

    if (Process) then

      Indq(1,nq) = 1
      Indq(2,nq) = (jAtom-1)*nsAtom+iAtom
      Indq(3,nq) = kDCRR+1

      Rij2 = (A(1,1)-A(1,2))**2+(A(2,1)-A(2,2))**2+(A(3,1)-A(3,2))**2
      Rab = sqrt(Rij2)
      if (Help) then
        RabCov = CovRad(ANr(iAtom))+CovRad(ANr(jAtom))
        !if ((iRow == 1) .and. (jRow == 1)) then
        !  ! Bond a la Fischer & Almlof
        f_Const = A_StrH(1)*exp(-A_StrH(2)*(Rab-RabCov))
        !else
        !  ij = iTri(iRow+1,jRow+1)-1
        !  f_Const = A_Str/(Rab-B_Str(ij))**3
        !end if
      else
        if (btest(iOptC,11) .and. (iBondType == 1)) then
          r0 = r_ref_vdW(iRow,jRow)
          f_Const = rkr_vdW*exp(-Alpha_vdW*(Rab-r0)**2)
        else
          r0 = rAV(iRow,jRow)
          Alpha = aAV(iRow,jRow)
          f_Const = rkr*exp(Alpha*(r0**2-rij2))
        end if
      end if

      f_Const = max(f_Const,f_Const_Min)
      if (iBondType == Fragments_Bond) f_Const = f_Const*1.0e3_wp
      fconst(nq) = sqrt(f_Const)
      rMult(nq) = Deg

      Valu(nq,iIter) = Val
      qLbl(nq) = Label

      ! Project the gradient vector

      call ProjSym(nCent,Ind,A,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,mqB,nB,nq,rMult(nq))

    end if

  end do     ! iCase

end do       ! iBond

return

end subroutine Bond_List
