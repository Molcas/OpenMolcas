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

subroutine Bond_List(nq,nsAtom,iIter,nIter,Cx,Process,value,nB,qLbl,fconst,rMult,LuIC,Indq,Proc_dB,iTabBonds,nBonds,iTabAI,mAtoms, &
                     mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: jStab, nStab, AtomLbl, ANr
use Slapaf_Parameters, only: iOptC

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
parameter(mB=2*3)
real*8 Cx(3,nsAtom,nIter), A(3,2), Grad(mB), Hess(mB**2), fconst(nB), value(nB,nIter), rMult(nB), BM(nB_Tot), dBM(ndB_Tot)
integer iDCRR(0:7), iStabM(0:7), Ind(2), iDCR(2), iChOp(0:7), Indq(3,nB), iTabBonds(3,nBonds), iTabAI(2,mAtoms), iBM(nB_Tot), &
        idBM(2,ndB_Tot), mqB(nB)
logical Process, Proc_dB, Help, R_Stab_A
character*14 Label, qLbl(nB)
character*3 ChOp(0:7)
character*(LenIn4) Lbls(2)
#include "bondtypes.fh"
#define _FMIN_
#define _VDW_
#include "ddvdt.fh"
#define _SCHLEGEL_
#include "ddvdt_bond.fh"
data ChOp/'E  ','X ','Y ','XY ','Z  ','XZ ','YZ ','XYZ'/
data iChOp/1,1,1,2,1,2,2,3/

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (nBonds < 1) return

nqB = 0
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) ' ---> Enter Bonds.'
write(6,*)
write(6,*) 'Process=',Process
call RecPrt('CX',' ',CX,3*nsAtom,nIter)
write(6,'(20(1X,A))') (AtomLbl(i),i=1,nsAtom)
write(6,*)
write(6,*) ' iTabAI'
write(6,*)
do iAtom=1,mAtoms
  write(6,*) iTabAI(1,iAtom),iTabAI(2,iAtom)
end do
#endif

! Loop over bonds

nCent = 2
do iBond=1,nBonds
  iBondType = iTabBonds(3,iBond)

  ! We will only incorpotate covalent and fragment bonds

  if (iBondType == vdW_Bond) Go To 1   ! vdW bonds
  if (iBondType > Magic_Bond) Go To 1  ! magic bonds

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
    if (jAtom > iAtom) Go To 2
    if (iDCR(1) /= iOper(0)) Go To 2
    if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(2) /= iOper(0))) Go To 2
    iRow = max(ANr(iAtom),1)
    jRow = max(ANr(jAtom),1)
#   ifdef _DEBUGPRINT_
    write(6,*) 'iAtom,jAtom=',iAtom,jAtom
#   endif
    Help = (iRow > 3) .or. (jRow > 3)
    Ind(1) = iAtom
    Ind(2) = jAtom
    call dcopy_(3,Cx(1,iAtom,iIter),1,A,1)
    write(Label,'(A,I2,A,I2,A)') 'B(',iAtom,',',jAtom,')'

#   ifdef _DEBUGPRINT_
    call RecPrt('A',' ',Cx(1,iAtom,iIter),1,3)
    call RecPrt('B',' ',Cx(1,jAtom,iIter),1,3)
#   endif


    ! Form double coset representatives

    call DCR(Lambda,jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCRR,nDCRR)
    kDCRR = iDCR(2)

#   ifdef _DEBUGPRINT_
    write(6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
    write(6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
    write(6,'(10A)') 'R={',(ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
    write(6,'(2A)') 'R=',ChOp(iDCR(2))
#   endif

    call OA(iDCR(2),Cx(1:3,jAtom,iIter),A(1:3,2))

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
    write(6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#   endif

    ! Now evaluate the degeneracy of the bond.

    iDeg = nIrrep/nStabM
    Deg = sqrt(dble(iDeg))
#   ifdef _DEBUGPRINT_
    write(6,*) ' nIrrep,nStabM=',nIrrep,nStabM
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
    write(6,'(A,I3.3,4A)') 'b',nqB,' = Bond ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2)
#   endif
    Label = ' '
    write(Label,'(A,I3.3)') 'b',nqB
    if (.not. Proc_dB) call FZero(Hess,36)
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
        !  ij = Max(iRow,jRow)*(Max(iRow,jRow)+1)/2+Min(iRow,jRow)
        !  f_Const = A_Str/(Rab-B_Str(ij))**3
        !end if
      else
        if ((iand(iOptC,2048) == 2048) .and. (iBondType == 1)) then
          r0 = r_ref_vdW(iRow,jRow)
          f_Const = rkr_vdW*exp(-Alpha_vdW*(Rab-r0)**2)
        else
          r0 = rAV(iRow,jRow)
          Alpha = aAV(iRow,jRow)
          f_Const = rkr*exp(Alpha*(r0**2-rij2))
        end if
      end if

      f_Const = max(f_Const,f_Const_Min)
      if (iBondType == Fragments_Bond) f_Const = f_Const*1.d3
      fconst(nq) = sqrt(f_Const)
      rMult(nq) = Deg

      value(nq,iIter) = Val
      qLbl(nq) = Label

      ! Project the gradient vector

      call ProjSym(nCent,Ind,A,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,mqB,nB,nq,rMult(nq))

    end if

2   continue
  end do     ! iCase

1 continue
end do       ! iBond

return

end subroutine Bond_List
