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
      SUBROUTINE full2red(XLT,Xab)
      use definitions, only: iwp, wp
      use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nBas, nnBstR,
     &                    nSym
      Implicit None
      real(kind=wp), intent(in):: XLT(*)
      real(kind=wp), intent(out):: Xab(*)

      Integer(kind=iwp) ISLT(8)
      Integer(kind=iwp), External:: cho_isao
      Integer(kind=iwp) iLoc,jSym,IS,ISYM,NB
      Integer(kind=iwp) jRab,kRab,iRab,iag,ibg,iSyma,ias,ibs,iab,kfrom

* Select table column for use with caspt2:
      iLoc=3
* jSym=1 always: Used for density matrices.
      jSym = 1
c Offsets to symmetry block in the LT matrix
      IS=0
      DO ISYM=1,NSYM
       ISLT(ISYM)=IS
       NB=NBAS(ISYM)
       IS=IS+(NB*(NB+1))/2
      END DO


      Do jRab=1,nnBstR(jSym,iLoc)
         kRab = iiBstr(jSym,iLoc) + jRab
         iRab = IndRed(kRab,iLoc)
         iag   = iRS2F(1,iRab)
         ibg   = iRS2F(2,iRab)
         iSyma = cho_isao(iag)
         ias   = iag - ibas(iSyma)
         ibs   = ibg - ibas(iSyma)
         if(ias.ge.ibs) then
          iab=(ias*(ias-1))/2+ibs
         else
          iab=(ibs*(ibs-1))/2+ias
         end if
         kfrom = isLT(iSyma) + iab
         Xab(jRab) = Xab(jRab)+XLT(kfrom)
      End Do

      End SUBROUTINE full2red

      SUBROUTINE red2full(XLT,Xab)
      use definitions, only: iwp, wp
      use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nBas, nnBstR,
     &                    nSym
      Implicit None
      real(kind=wp), intent(out):: XLT(*)
      real(kind=wp), intent(in):: Xab(*)

      Integer(kind=iwp) ISLT(8)
      Integer(kind=iwp), External:: cho_isao
      Integer(kind=iwp) iLoc,jSym,IS,ISYM,NB
      Integer(kind=iwp) jRab,kRab,iRab,iag,ibg,iSyma,ias,ibs,iab,kto

* Select table column for use with caspt2:
      iLoc=3
* jSym=1 always: Used for density matrices.
      jSym = 1
c Offsets to symmetry block in the LT matrix
      IS=0
      DO ISYM=1,NSYM
       ISLT(ISYM)=IS
       NB=NBAS(ISYM)
       IS=IS+(NB*(NB+1))/2
      END DO

      Do jRab=1,nnBstR(jSym,iLoc)
         kRab = iiBstr(jSym,iLoc) + jRab
         iRab = IndRed(kRab,iLoc)
         iag   = iRS2F(1,iRab)
         ibg   = iRS2F(2,iRab)
         iSyma = cho_isao(iag)
         ias   = iag - ibas(iSyma)
         ibs   = ibg - ibas(iSyma)
         if(ias.ge.ibs) then
          iab=(ias*(ias-1))/2+ibs
         else
          iab=(ibs*(ibs-1))/2+ias
         end if
         kto = isLT(iSyma) + iab
         XLT(kto) = XLT(kto)+Xab(jRab)
      End Do

      End SUBROUTINE red2full
