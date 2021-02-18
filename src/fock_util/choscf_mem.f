************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE CHOSCF_MEM(nSym,nBas,iUHF,DoExchange,ipNocc,
     &                      ALGO,REORD,MinMem,lOff1)

*****************************************************************
*  Author : F. Aquilante
*
*  Purpose:
*            Returns an array MinMem that contains the minimum
*            amount of memory required for each Cholesky
*            vector for the building of the frozen AO-Fock matrix
*
*            The first vector is stored in a larger memory block
*            in some cases in order to re-use the allocated memory
*            for computing the exchange intermediate.
*            The value of LOFF1 specifies the amount of memory
*            reserved for the 1st vector
*
*            Note that this is a specialized code and can be used
*                 only in connection with the truth table defined
*                 in the calling routine
******************************************************************
      use ChoArr, only: nDimRS
      Implicit Real*8 (a-h,o-z)
      Integer nSym,nBas(nSym),MinMem(nSym),iUHF,ALGO
      Integer Moccmx(nSym),Mabmx(nSym),MxBas(nSym)
      Logical REORD,xToDo,DoExchange(*)
      Integer ipNocc(*)

#include "WrkSpc.fh"

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
**************************************************



      If (iUHF .eq. 0) Then
         nDen = 1
         xToDo=DoExChange(1)
      Else
         nDen = 3
         xToDo=DoExChange(2)
      End If

C============================
      lOff1=0
      Do i=1,nDen
         Do j=1,nSym
            lOff1=Max(lOff1,iWork(ipNocc(i)+j-1))
         End Do
      End Do

      do j=1,nSym
         Moccmx(j)=0
         do i=1,nDen
            Moccmx(j)=Max(Moccmx(j),iWork(ipNocc(i)+j-1))
         end do
      end do

C Max dimension of a read symmetry block
        Nmax=0
        Do iSym=1,nSym
           if(NBAS(iSym).gt.Nmax .and. Moccmx(isym).ne.0 )then
           Nmax = NBAS(iSym)
           endif
        End Do

        lOff1 = Nmax*lOff1

C============================
      Do jSym=1,nSym

C Total length of a vector
       Mabmx(jSym)=0
       MxBas(jSym)=0
       Nab=0
       NSab=0
       Do ksym=1,nSym
          iSymp=MulD2h(ksym,jSym)
          If (iSymp.gt.ksym .and. (Moccmx(iSymp).ne.0 .or.
     &                             Moccmx(ksym).ne.0) )  Then
             Nab = Nab + nBas(ksym)*nBas(iSymp)
             Mabmx(jSym)= Max(Mabmx(jSym),Max(nBas(ksym)*Moccmx(iSymp),
     &                        nBas(iSymp)*Moccmx(ksym)))
             MxBas(jSym)= Max(MxBas(jSym),nBas(ksym)*nBas(iSymp))
          Else
               If (iSymp.eq.ksym) Then ! all are needed for the Coulomb
                  Nab = Nab + nBas(ksym)*(nBas(ksym)+1)/2
                  NSab= NSab + nBas(ksym)**2
               End If
          End If
       End Do

C ====================================================================
C === The following memory management is bound to the truth table  ===
C === assigned in the calling routine
C=====================================================================

        IF(.not.xToDo) THEN

         MinMem(jSym) = Nab + 1  ! to store 1 vector of L(rJ),s + V(J)
         if(.not.REORD) then
           MinMem(jSym)=Nab+nDimRS(jSym,1)
c           ! L + read 1 vect reduced set1
         end if

        ELSE ! Memory for "off-diagonal" exchange only (jSym.ne.1)

             if (REORD) then
               MinMem(jSym) = 2*Nab ! 1 vector of L(rJ),s + W(rJ),s
             else
                if(ALGO.eq.2)then
                  MinMem(jSym) = Nab + Max(nDimRS(jSym,1),Mabmx(jSym))
                else
                  MinMem(jSym) = Nab + Max(nDimRS(jSym,1),MxBas(jSym))
                endif
             endif

        END IF

C-----
        IF(xToDo.and.jSym.eq.1) THEN  ! jSym=1 is a special case

           If(nSym.eq.1) then

               if (ALGO.eq.2) then
                   if (lOff1.gt.Nab) then
                      MinMem(jSym) = lOff1 + NSab
                   else
                      MinMem(jSym) = Nab + NSab
                      lOff1 = Nab
                   endif
               else

                   MinMem(jSym) = 2*NSab
                   lOff1 = NSab

               endif

           Else  ! nSym > 1   jSym=1

               if (ALGO.eq.2) then
                   if (lOff1.gt.nBas(1)*(nBas(1)+1)/2) then
                      MinMem(jSym) = Nab - nBas(1)*(nBas(1)+1)/2
     &                             + lOff1 + Nmax**2
                   else
                      MinMem(jSym) = Nab + Nmax**2
                      lOff1 = nBas(1)*(nBas(1)+1)/2
                   endif
               else

                   MinMem(jSym) = 2*(Nmax**2)
     &                          + (Nab-nBas(1)*(nBas(1)+1)/2)
                   lOff1 = Nmax**2

               endif

           EndIf

        ENDIF   ! JSYM=1  special case

C ====================================================================

      End Do  ! loop over jSym

      Return
      END
