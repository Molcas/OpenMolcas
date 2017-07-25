      SubRoutine Drv_Fck_DMET(Label,ip,lOper,nComp,CCoor,
     &                   nOrdOp,rNuc,rHrmt,iChO,
     &                   DMET_f,nBfn)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
#include "warnings.fh"
      Real*8 DMET_f(nBfn,nBfn)
      Real*8, Dimension(:), Allocatable :: Array
      Character Label*8
      Real*8 CCoor(3,nComp), rNuc(nComp)
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 112
      iPrint = nPrint(iRout)
      Call qEnter('Drv_Fck')
      If (iPrint.ge.19) Then
         Write (6,*) ' In OneEl: Label', Label
         Write (6,*) ' In OneEl: nComp'
         Write (6,'(1X,8I5)') nComp
         Write (6,*) ' In OneEl: lOper'
         Write (6,'(1X,8I5)') lOper
         Write (6,*) ' In OneEl: n2Tri'
         Do iComp = 1, nComp
            ip(iComp) = n2Tri(lOper(iComp))
         End Do
         Write (6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
         Call RecPrt(' CCoor',' ',CCoor,3,nComp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
      nIC = 0
      llOper = 0
      Do iComp = 1, nComp
         llOper = iOr(llOper,lOper(iComp))
         Do iIrrep = 0, nIrrep-1
            If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0) nIC = nIC + 1
         End Do
      End Do
      If (iPrint.ge.20) Write (6,*) ' nIC =',nIC
      If (nIC.eq.0) Go To 999
      Call SOS(iStabO,nStabO,llOper)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for symmetry adapted one electron integrals.
*     Will just store the unique elements, i.e. low triangular blocks
*     and lower triangular elements in the diagonal blocks.
*
      Call ICopy(nComp,-1,0,ip,1)
      LenTot=0
      Do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      Call mma_allocate(Array,LenTot,label='Array')
      ip(1)=1
      Call DCopy_(LenTot,Zero,0,Array(ip(1)),1)
      iadr=ip(1)
      do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         ip(icomp)=iadr
         iadr=iadr+LenInt+4
*        Copy center of operator to work area.
         Call DCopy_(3,Ccoor(1,iComp),1,Array(ip(iComp)+LenInt),1)
*        Copy nuclear contribution to work area.
         Array(ip(iComp)+LenInt+3) = rNuc(iComp)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Compute all SO integrals for all components of the operator.

        Do iBfn = 1, nBfn
           Do jBfn = 1, iBfn
                 ijBfn = iBfn*(iBfn-1)/2 + jBfn - 1 + ip(1)
                 Array(ijBfn)=DMET_f(iBfn,jBfn)
           End Do
        End Do
*                                                                      *
************************************************************************
*                                                                      *
*                    P O S T P R O C E S S I N G                       *
*                                                                      *
************************************************************************
*                                                                      *
      write(6,*) 'Fock Matrix'
      Call PrMtrx(Label,lOper,nComp,ip,Array)
      If (iPrint.ge.10) Call PrMtrx(Label,lOper,nComp,ip,Array)
*                                                                      *
************************************************************************
*                                                                      *
*---- Write integrals to disc.
*
      mpp_state=1
      Do iComp = 1, nComp
         iSmLbl = lOper(iComp)
*                                                                      *
************************************************************************
*                                                                      *
*------- Write integrals to disc
*
         iOpt = 0
         iRC = -1
         write(6,*) "before if"
         If (Label(1:3).eq.'PAM')
     &      Write(Label,'(A5,I3.3)') 'PAM  ',iPAMcount
         Write(6,*) ' oneel *',Label,'*'

         Call WrOne(iRC,iOpt,Label,iComp,Array(ip(iComp)),iSmLbl)
         write(6,*) 'after wrone'

         If (Label(1:3).eq.'PAM')
     &      Call WrOne(iRC,iOpt,Label,1,Array(ip(iComp)),iSmLbl)
         iPAMcount=iPAMcount+1

         If (iRC.ne.0) then
            Call qTrace
            Write(6,*) ' *** Error in subroutine ONEEL ***'
            Write(6,*) '     Abend in subroutine WrOne'
            Call Quit(_RC_IO_ERROR_WRITE_)
         End If
      End Do  ! iComp
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory for integral
*
      Call mma_deallocate(Array)
*                                                                      *
************************************************************************
*                                                                      *
 999  Continue
      Call qExit('Drv_Fck_DMET')
      Return
      End
