************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Compute_Xhole_Int(nBasLop,nSym,ipSqMom,Func,nSize)
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "nq_info.fh"
#include "status.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: D1ao(:)
      Dimension nBasLop(nSym)
!     Logical Do_Gamma, Do_Grad, On_Top, Do_Tau, Do_MO, Do_TwoEl
      Logical DSCF
      Character*4 DFTFOCK, KSDFT
      Integer nSize
      Real*8, Allocatable:: CMO(:)

*
*-- Check symmetry
*
      If(nSym.ne.1) then
        Write(6,*)
        Write(6,*)' You should not run LoProp with symmetry!'
        Call Abend()
      Endif
*
*-- Set a lot of numbers and labels. See below for more help.
*
      nB=nBasLop(1)     !-- number of basis functions
      nTri=nB*(nB+1)/2  !-- obvious!
      Func=Zero         !-- initialize
!IFG: uncomment when DrvNQ call fixed (see below)
!     Dens=Zero         !-- initialize
!     nFckDim=1         !-- for open-shell the fock-matrix has
!                       !   two versions, but this implementation
!                       !   is not for UHF (but for CASSCF)
!     nD=1              !-- much like nFckDim
!     Do_Gamma=.false.  !-- optional stuff in DFT-theory; for our
!                       !   purpose, just scheiss.
!     Do_Grad=.false.   !-- see above.
!     On_Top=.false.    !-- see above.
!     Do_Tau=.false.    !-- see above.
!     Do_MO=.true.      !-- in our integration kernel, we need MOs.
!     Do_TwoEl=.false.  !-- scheiss in our functional.

      Write(DFTFOCK,'(A)')'XHOL'  !-- Tell the routines that we wish
                                  !   to run a xhole calculation.
      Write(KSDFT,'(A)')'LDA'     !-- Just to get the routines to
                                  !   take the fastest path, no other
                                  !   meaning.
      Call mma_allocate(D1ao,nTri)
      nDens=nTri
      Call Get_D1ao(D1ao,nDens)!-- The density matrix.
      Functional_type=LDA_type    !-- Number from nq_info.fh.
      EThr=1.0d-9
      Call Put_dScalar('EThr',EThr)!-- A threshold for energy accuracy
                                   !   in DFT/SCF. Dummy.
      Call Seward_Init             !-- Initialize a lot of shit.
      nDiff=0                      !-- like above
      DSCF=.false.                 !-- like above
      Call GetInf(DSCF,nDiff)      !-- like above
      Call SetUp_iSD               !-- like above
      Call Get_iScalar('nSym',mIrrep)    !-- Stupid number in nq_info
                                         !   needed to fool do_mo
                                         !   stuff, since they normally
                                         !   are only activated by
                                         !   CASDFT
      Call Get_iArray('nBas',mBas,mIrrep)!-- Like above

*
*-- Then we need orbital density dipoles.
*
      nCMO=nB**2
      Call mma_allocate(CMO,nCMO,Label='CMO')
      nOrb=INT(sqrt(dble(nCMO)))
      Call Get_CMO(CMO,nCMO)
      nB=mBas(0)
      Call GetMem('MultSq','Allo','Real',iMultSq,nB**2)
      Call GetMem('TEMP','Allo','Real',iTEMP,nB*nOrb)
      Call GetMem('MultiKulti','Allo','Real',iMult1,nOrb**2+4)
      Call GetMem('OrbDipsX','Allo','Real',ip_OrbDip(1),nOrb*(nOrb+1)/2)
      Call GetMem('OrbDipsY','Allo','Real',ip_OrbDip(2),nOrb*(nOrb+1)/2)
      Call GetMem('OrbDipsZ','Allo','Real',ip_OrbDip(3),nOrb*(nOrb+1)/2)
      irc=-1
      Lu_One=49
      Lu_One=IsFreeUnit(Lu_One)
      Call OpnOne(irc,0,'ONEINT',Lu_One)
      If(irc.ne.0) then
        Write(6,*)
        Write(6,*)'ERROR! Could not open one-electron integral file.'
        Call Abend()
      Endif
      Do i=1,3
        irc=-1
        iOpt=1
        iSmLbl=1
        Call iRdOne(irc,iOpt,'Mltpl  1',i,nSize,iSmLbl)
        irc=-1
        iOpt=0
        iSmLbl=0
        Call RdOne(irc,iOpt,'Mltpl  1',i,Work(iMult1),iSmLbl)
        Call Square(Work(iMult1),Work(iMultSq),1,nB,nB)
      Call DGEMM_('T','N',nOrb,nB,nB,
     &            One,CMO,nB,
     &                Work(iMultSq),nB,
     &            Zero,Work(iTEMP),nOrb)
      Call DGEMM_('N','N',nOrb,nOrb,nB,
     &            One,Work(iTEMP),nOrb,
     &                CMO,nB,
     &                Zero,Work(iMult1),nOrb)
        kaunt1=0
        kaunt2=0
        Do iO1=1,nOrb
          Do iO2=1,nOrb
            If(iO1.ge.iO2) then
              Work(ip_OrbDip(i)+kaunt1)=Work(iMult1+kaunt2)
              kaunt1=kaunt1+1
            Endif
            kaunt2=kaunt2+1
          Enddo
        Enddo
      Enddo
      Call ClsOne(irc,Lu_One)
      Call GetMem('MultSq','Free','Real',iMultSq,nB**2)
      Call GetMem('TEMP','Free','Real',iTEMP,nB*nOrb)
      Call GetMem('MultiKulti','Free','Real',iMult1,nOrb**2+4)
      Call mma_deallocate(CMO)

*
*-- Anders' little helper on DrvNQ:
*     Argument (1): The name on the integration kernel. Is a routine
*                   in src/nq_util/ directory.
*              (2): Through some "fooling" of the integration routines
*                   this argument will upon return contain the matrix
*                   elements.
*              (3): nFckDim is described above. Just one.
*              (4): Here the functional "energy" comes, in other words
*                   the molecular expectation value
*              (5): Just shit.
*              (6): The one-electron density matrix.
*              (7): Dimension on all one-electron matrices.
*              (8): Like (3).
*              (9-17): Total scheiss from our perspective.
*              (18): Label to signal that we will need MOs.
*              (19): We use no two-electron stuff, hence false.
*              (20): Usually, this tells what type of Fock-matrix
*                    we use, but now we "rob" this variable for our
*                    purpose and let the value XHOL signify that we
*                    are computing the xhole dipole.
      Call GetMem('X-Dipole elements','Allo','Real',ip_MatEl,nTri)
!IFG: The call to DrvNQ below is completely messed up, please fix
      Call WarningMessage(2,'There is surely a bug here!')
!     Call DrvNQ(Do_XHoleDip,Work(ip_MatEl),nFckDim,Func,Dens
!    &          ,D1ao,nTri,nD,Do_Gamma,Do_Grad,Dummy
!    &          ,iDummy,Dummy,Dummy,iDummy,On_Top,Do_Tau,Do_MO
!    &          ,Do_TwoEl,DFTFOCK)
*      FFF=ddot_(nDens,D1ao,1,Work(ip_MatEl),1)
*      write(6,*)'YYY:',nDens,FFF,Func,ip_MatEl
*
*-- Put the second-moments in square form.
*
      Call GetMem('2MomSq','Allo','Real',ipSqMom,nB**2)
      Call Square(Work(ip_MatEl),Work(ipSqMom),1,nB,nB)

*
*-- Deallocate
*
      Call mma_deallocate(D1ao)
      Call Free_iSD()
      Call GetMem('X-Dipole elements','Free','Real',ip_MatEl,nTri)
      Call GetMem('OrbDipsX','Free','Real',ip_OrbDip(1),nOrb*(nOrb+1)/2)
      Call GetMem('OrbDipsY','Free','Real',ip_OrbDip(2),nOrb*(nOrb+1)/2)
      Call GetMem('OrbDipsZ','Free','Real',ip_OrbDip(3),nOrb*(nOrb+1)/2)
      Call Free_HerRW()

      Return
      End
