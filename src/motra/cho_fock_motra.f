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
      SUBROUTINE Cho_Fock_MoTra(nSym,nBas,nFro,
     &                DLT,DSQ,FLT,nFLT,FSQ,ExFac)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym)
      Real*8  DLT(*), DSQ(*), FLT(nFLT), FSQ(*)
      Real*8  ExFac
#include "WrkSpc.fh"

      Integer  ip_of_Work
      External ip_of_Work

      Integer NScreen
      Real*8  dmpk, dFKmat
      Integer nDen, nXorb(8), nChMOs(8)


*****************************************************************
*  CALCULATE AND RETURN FMAT DUE TO FROZEN ORBITALS ONLY
*****************************************************************


      NScreen=10
      dmpK=1.0d-1
      dFKmat=0.0d0
      nDen=1
      Call IZero(nXorb,nSym)

* --- Initialize Cholesky information

      ChFracMem=0.0d0
      CALL CHO_X_INIT(irc,ChFracMem)
      if (irc.ne.0) then
         write(6,*)'Cho_Fock_Motra: Cho_X_Init returns error code ',irc
         Call AbEnd()
      endif

      ipDSQ = ip_of_Work(DSQ(1))  ! not needed on exit

      MOdim=0
      Do i=1,nSym
         MOdim=MOdim+nBas(i)**2
      End Do
      CALL GETMEM('choMOs','allo','real',ipMOs,MOdim)

      ipd = ipDSQ
      ipV = ipMOs
      Do i=1,nSym
          if(nBas(i).gt.0)then
            Ymax=0.0d0
            do ja=1,nBas(i)
               jaa=ipd-1+nBas(i)*(ja-1)+ja
               Ymax=Max(Ymax,Work(jaa))
            end do
            Thr = 1.0d-8*Ymax
            CALL CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),
     &                     NumV,Thr,irc)
            if(irc.ne.0)then
              write(6,*)'Cho_Fock_Motra: CD_incore returns rc ',irc
              Call AbEnd()
            endif

            nChMOs(i)= NumV

            if ( NumV .ne. nFro(i) ) then
             write(6,'(a,a,i6,a,i6,a,i6,a,i6,a,i6)')
     &       'Warning! Cho_Fock_Motra: nr of Frozen orbitals from the ',
     &       'decomposition of the density matrix is ',numV,
     &       ' in symm. ',i, '; Expected value = ',nFro(i),
     &       '; Max diagonal of the density in symm. ',i,
     &       ' is equal to ',Ymax
            endif

          else

            nChMOs(i)= 0

          endif

          ipd = ipd + nBas(i)**2
          ipV = ipV + nBas(i)**2
      End Do

      ipDLT = ip_of_Work(DLT(1))
      ipFLT = ip_of_Work(FLT(1))
      ipFSQ = ip_of_Work(FSQ(1))  ! not needed on exit

      CALL CHO_LK_SCF(irc,nDen,[ipFLT],[ipFSQ],nXorb,nFro,
     &                [ipMOs],[ipDLT],0.5d0*ExFac,NScreen,dmpk,dFKmat)
      if (irc.ne.0) then
         write(6,*)'Cho_Fock_Motra: Cho_LK_scf returns error code ',irc
         Call AbEnd()
      endif

      Call GADSUM(FLT,nFLT)

      CALL GETMEM('choMOs','free','real',ipMOs,MOdim)

* --- Finalize Cholesky information

      CALL CHO_X_FINAL(irc)
      if (irc.ne.0) then
         write(6,*)'Cho_Fock_Motra: Cho_X_Final returns error code ',irc
         write(6,*)'Try recovery -- continue.'
      endif

      RETURN
      END
