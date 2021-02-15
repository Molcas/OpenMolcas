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
* Copyright (C) 2005, Giovanni Ghigo                                   *
************************************************************************
*  Cho_GenE
*
*> @brief
*>   Routine for the generation of the ``A,B`` block of exchange integrals (in symmetries \p iSymA, \p iSymB)
*>   for occupied MO \p iI, \p iJ in symmetries \p iSymI, \p iSymJ
*> @author Giovanni Ghigo
*>
*> @details
*> The routine generates the ``A,B`` block of integrals gathering 9
*> sub-blocks. These are combination of inactive, active, and
*> secondary ``A,B`` MO.
*>
*> @param[in] iSymI,iSymJ,iSymA,iSymB Symmetry block of the two-electrons integrals
*> @param[in] NumV                    Number of Cholesky vectors to transform in the current batch
*> @param[in] iAddEx                  Memory pointer of the ``A,B`` integrals block
*> @param[in] LenEx                   Length of the ``A,B`` integrals block
************************************************************************
      Subroutine Cho_GenE(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                      iAddEx,LenEx)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           Jannuary-June 2005                                         *
*----------------------------------------------------------------------*
* This is the routine that really generates the A,B block of exchange  *
* integrals (in symmetries iSymA,iSymB) for occupied MO iI,iJ in       *
* symmetries iSymI,iSymJ. The 3 x 3 A,B block is built gathering 9     *
* sub-blocks. These are combination of inactive, active, and secondary *
* A,B MO                                                               *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Integer iAddSB(3,3)

      Integer LenA(3), LenB(3)

      Do iSB_A = 1, 3
        Do iSB_B = 1, 3
          iAddSB(iSB_A,iSB_B)= ip_Dummy  ! Mem Address of the SubBlocks
        EndDo
      EndDo

* --- GENERATION of SubBlocks
      If (SubBlocks(1,1)) Call MkExSB11(iAddSB(1,1),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
      If (SubBlocks(1,2)) Call MkExSB12(iAddSB(1,2),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
      If (SubBlocks(1,3)) Call MkExSB13(iAddSB(1,3),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
      If (SubBlocks(2,1)) Call MkExSB21(iAddSB(2,1),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV, iAddSB(1,2))
      If (SubBlocks(2,2)) Call MkExSB22(iAddSB(2,2),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
      If (SubBlocks(2,3)) Call MkExSB23(iAddSB(2,3),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
      If (SubBlocks(3,1)) Call MkExSB31(iAddSB(3,1),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV, iAddSB(1,3))
      If (SubBlocks(3,2)) Call MkExSB32(iAddSB(3,2),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV, iAddSB(2,3))
      If (SubBlocks(3,3)) Call MkExSB33(iAddSB(3,3),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
* --- END GENERATION of SubBlocks

* --- GATERING of SubBlocks
      iAddExSB = iAddEx
      IF (DoTCVA) THEN

       LenB(1) = nIsh(iSymB)
       LenB(2) = nAsh(iSymB)
       LenB(3) = nSsh(iSymB)
       LenA(1) = nIsh(iSymA)
       LenA(2) = nAsh(iSymA)
       LenA(3) = nSsh(iSymA)
       If (iSymA.NE.iSymB) then

         Do iSB_A = 1, 3
           Do iA = 1,LenA(iSB_A)
             Do iSB_B = 1, 3
               If (LenB(iSB_B)==0) Cycle                  ! SB(1,iSB_A)

               iAddSBi = iAddSB(iSB_A,1) + LenB(iSB_B) * (iA-1)
               call daxpy_(LenB(iSB_B),1.0d0,Work(iAddSBi),1,
     &                      Work(iAddExSB),1)
               iAddExSB = iAddExSB + LenB(iSB_B)
             EndDo
           EndDo  ! iA
         EndDo ! iSB_A

       else

         Do iSB_B = 1, 3
           Do iB = 1,LenB(iSB_B)
             Do iSB_A = 1, 3
               If (LenA(iSB_A)==0) Cycle                   ! SB(1,iSB_B)

               iAddSBi = iAddSB(iSB_B,iSB_A) + LenA(iSB_A) * (iB-1)
               call daxpy_(LenA(iSB_A),1.0d0,Work(iAddSBi),1,
     &                     Work(iAddExSB),1)
               iAddExSB = iAddExSB + LenA(iSB_A)
             EndDo  ! iSB_A

           EndDo  ! iB
         EndDo ! iSB_B

       EndIf

      ELSE

       If (iSymA.NE.iSymB) then
         iLenBb = nSsh(iSymB)
         If (iLenBb.GT.0) then                          ! SB(3,3)
           iLenA = nSsh(iSymA)
           Do iA = 1,iLenA
             iAddSBi = iAddSB(3,3) + iLenBb * (iA-1)
             call daxpy_(iLenBb,1.0d0,Work(iAddSBi),1,Work(iAddExSB),1)
             iAddExSB = iAddExSB + iLenBb
           EndDo  ! iA
         EndIf
       else
         iLenAa = nSsh(iSymA)
         If (iLenAa.GT.0) then                          ! SB(3,3)
           iLenB = nSsh(iSymB)
           Do iB = 1,iLenB
             iAddSBi = iAddSB(3,3) + iLenAa * (iB-1)
             call daxpy_(iLenAa,1.0d0,Work(iAddSBi),1,Work(iAddExSB),1)
             iAddExSB = iAddExSB + iLenAa
           EndDo  ! iB
         EndIf
       EndIf

      ENDIF
* --- END GATERING of SubBlocks

      iDum=0
      Do iSB_A = 1, 3
        Do iSB_B = 1, 3
          If (iAddSB(iSB_A,iSB_B).NE.ip_Dummy)
     &      Call GetMem('SB','Free','Real',
     &                  iAddSB(iSB_A,iSB_B), iDum)
        EndDo
      EndDo

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LenEx)
      End
