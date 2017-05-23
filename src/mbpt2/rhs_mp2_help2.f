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

      subroutine rhs_mp2_help2(iSymA, iSymB, iSymI, iSymJ)
*
      Implicit Real*8 (a-h,o-z)
*
*defining One etc.
#include "real.fh"
#include "trafo.fh"
#include "WrkSpc.fh"
#include "files_mbpt2.fh"
#include "mp2grad.fh"
#include "corbinf.fh"
*
      iVirVir(i,j,k) = nFro(k) + nOcc(k) +  j-1 +
     &                (nOrb(k) + nDel(k))*
     &                (i + nFro(k) + nOcc(k) -1)
      iOccOcc(i,j,k) = j-1 + (nOrb(k)+nDel(k))
     &                               * (i -1)
      iMp2Lagr(i,j,k) = ip_Mp2Lagr(k) +
     &                         j-1 + (nOcc(k)+nFro(k))*(i-1)
*--------------------------------------------------
*      write(6,*) 'Starting subroutine with symmetries: '
*      write(6,*) 'I,J,A,B', I, J, iA, iB
*      write(6,*) 'NewSym',NewSym
*      write(6,*) 'IpoDens', IpoDens
*--------------------------------------------------
*
*     Reads the <..|AB> - integrals for the L3-term
*
      nB = nExt(iSymB) + nDel(iSymB)
      Do iA = 1, nExt(iSymA) + nDel(iSymA)
*        Decide i the AB-matrix is symmetric, if so
*        only loop over lower triangle
         if(iSymA.eq.iSymB) nB = iA
         Do iB = 1, nB
*
            fac_ab = One
            if((iA.ne.iB).and.(iSymA.eq.iSymB)) fac_ab = Two
*           If I is not equal J there is a symmetry J I B A
*           that is identical and not stored.
            if(iSymA.ne.iSymB) fac_ab = Two
*
            Call Exch(iSymI,iSymA,iSymJ,iSymB,
     &                iA+nOcc(iSymA)+nFro(iSymA),
     &                iB+nOcc(iSymB)+nFro(iSymB),
     &                Work(ipInt1), Work(ipScr1))
            If (iSymA.ne.iSymB) then
               Call Exch(iSymJ,iSymA,iSymI,iSymB,
     &                   iA+nOcc(iSymA)+nFro(iSymA),
     &                   iB+nOcc(iSymB)+nFro(iSymB),
     &                   Work(ipInt2),Work(ipScr1))
            End If
*            write(6,*) 'iOrbA',iA , 'iOrbB',iB
*            Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
*     &                nOrb(iSymA)+nDel(iSymA),
*     &                nOrb(iSymB)+nDel(iSymA))
*            If (iSymA.Ne.iSymB) Then
*               Call RecPrt('Int2:','(8F10.6)',Work(ipInt2),
*     &                     nOrb(iSymB)+nDel(iSymB),
*     &                     nOrb(iSymA)+nDel(iSymA))
*            End If
*
****************************************************************************************
*
*                  Lagrangian term 3
*
***************************************************************************************
*
*
*
            If((iSymA.Eq.iSymI).And.(iSymB.Eq.iSymJ)) Then
               If((iA.Eq.iB).And.(iSymA.Eq.iSymB)) Then
                  Fac = Half
               Else
                  Fac = One
               Endif
               iIC =(nFro(iSymB) + nOcc(iSymB)) *
     &              (nOrb(iSymA) + nDel(iSymA))
               iCI = nFro(iSymB) + nOcc(iSymB)
               Do iC = 1, nExt(iSymB) + nDel(iSymB)
                  Do iI = 1, nOcc(iSymA) + nFro(iSymA)
                     xaibc = Work(ipInt1 + iIC + iI-1 + (iC-1)
     &                     * (nOrb(iSymA)+nDel(iSymA)))
                     If(iSymA.Eq.iSymB) Then
                        xacbi = Work(ipInt1 + iCI +  iC-1 +
     &                          (iI-1)*(nOrb(iSymB)+nDel(iSymB)))
                     Else
                        xacbi = Work(ipInt2 + iCI + iC-1 +
     &                          (iI-1)*(nOrb(iSymB)+nDel(iSymB)))
                     End If
                    Work(iMp2Lagr(iA,iI,iSymA)) =
     &                    Work(iMp2Lagr(iA,iI,iSymA)) -
     &                        Fac*Work(ip_Density(iSymB) +
     &                                 iVirVir(iB,iC,iSymB))
     &                              * (Two*xaibc - xacbi)
                     If(iSymA.Eq.iSymB) Then
                       Work(iMp2Lagr(iB,iI,iSymA)) =
     &                    Work(iMp2Lagr(iB,iI,iSymA)) -
     &                        Fac*Work(ip_Density(iSymB) +
     &                                 iVirVir(iA,iC,iSymB))
     &                              * (Two*xacbi - xaibc)
                   End If
*                        Write(6,*) 'AIBC', iA, iI,iB,iC
*                        Write(6,*) 'AIBC2',iB,iI,iA, iC
*                        Write(6,*) 'icab', xaibc
*                        Write(6,*) 'ibac', xacbi
*                        Write(6,*) 'Dens', Work(ip_Density(iSymB) +
*     &                                 iVirVir(iB,iC,iSymB))
*                        Write(6,*) 'Dens2',Work(ip_Density(iSymB) +
*     &                                 iVirVir(iA,iC,iSymB))
*                        Write(6,*) 'A', (Two*xaibc - xacbi)
*                        Write(6,*) 'A2', (Two*xacbi - xaibc)
*                        Write(6,*) 'Bidrag', Work(ip_Density(iSymB)+
*     &                                 iVirVir(iB,iC,iSymB))
*     &                              * (Two*xaibc - xacbi)
*                        Write(6,*) 'Bidrag2',Work(ip_Density(iSymB) +
*     &                                 iVirVir(iA,iC,iSymB))
*     &                              * (Two*xacbi - xaibc)

*
                  EndDo
               EndDo
*
*           And again we need another loop for the symmetries that isnt written
*           explicitly. all references to A and B need to change too.
*
               If(iSymA.ne.iSymB) then
                  iIC = (nFro(iSymA) + nOcc(iSymA)) *
     &                  (nOrb(iSymB) + nDel(iSymB))
                  iCI = nFro(iSymA) + nOcc(iSymA)
                  Do iC = 1, nExt(iSymA) + nDel(iSymA)
                     Do iI = 1, nFro(iSymB) + nOcc(iSymB)
                        xaibc = Work(ipInt1 + iCI + iC-1 +
     &                           (iI-1)*(nOrb(iSymA)+nDel(iSymA)))
                        xacbi = Work(ipInt2 + iIC + iI-1 +
     &                           (iC-1)*(nOrb(iSymB)+nDel(iSymB)))
                        Work(iMp2Lagr(iB,iI,iSymB)) =
     &                    Work(iMp2Lagr(iB,iI,iSymB)) -
     &                        Fac*Work(ip_Density(iSymA) +
     &                                 iVirVir(iA,iC,iSymA))
     &                              * (Two*xaibc - xacbi)
                     EndDo
                  EndDo
               EndIf
            Endif
*

         EndDo
      EndDo
*     Reads one block of <ij|..> at the time.
*
      nJ = nFro(iSymJ) + nOcc(iSymJ)
      Do iI = 1, nFro(iSymI) + nOcc(iSymI)
*        Decide if the IJ-block i symmetric and
*        if so only read a triangle
         if(iSymI.eq.iSymJ) nJ = iI
******** Check if the IJ-matrix is symmetric, if so
***      only loop over lower triangle.
*
*****************************************
         Do iJ = 1, nJ
*           If we are only looping over a triangular matrix
*           we need to count each offdiagonal element twice.
            fac_ij = One
            If((iI.ne.iJ).and.(iSymI.eq.iSymJ)) fac_ij = Two
*           If I is not equal J there is a symmetry J I B A
*           that is identical and not stored.
            If(iSymI.Ne.iSymJ) fac_ij = Two
*           Copy one AB-block of integrals to the workspace
            Call Exch(iSymA,iSymI,iSymB,iSymJ,iI,
     &                iJ,Work(ipInt1),
     &                Work(ipScr1))
            if (iSymI.ne.iSymJ) then
               Call Exch(iSymB,iSymI,iSymA,iSymJ,
     &                   iI,iJ,
     &                   Work(ipInt2),Work(ipScr1))
            endif
*            Write(6,*)
*            Write(6,*) ' *  i,j = ',iI,iJ
*            Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),nOrb(iSymA),
*     &                nOrb(iSymB))
*            if (iSymI.ne.iSymJ) then
*            Call RecPrt('Int2:','(8F10.6)',Work(ipInt2),nOrb(iSymB),
*     &                nOrb(iSymA))
*            endif


****************************************************************************************
*
*                  Lagrangian term 4
*
***************************************************************************************
*
*
*
            If((iSymI.eq.iSymA).and.(iSymJ.eq.iSymB)) then
               If((iI.eq.iJ).and.(iSymI.eq.iSymJ)) Then
                  Fac = Half
               Else
                  Fac = One
               Endif
*
*
*
               iKA = (nFro(iSymJ) + nOcc(iSymJ)) *
     &               (nOrb(iSymI) + nDel(iSymI))
               iAK = nFro(iSymJ) + nOcc(iSymJ)
               Do iK = 1, nFro(iSymI) + nOcc(iSymI)
                  Do iA = 1, nExt(iSymJ) + nDel(iSymJ)
                     xikja = Work(ipInt1 + iKA + iK-1 +
     &                            (iA-1)*(nOrb(iSymI)+nDel(iSymI)))
                     if(iSymI.eq.iSymJ) then
                        xiajk = Work(ipInt1 + iAK + iA-1 +
     &                               (iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                     else
                        xiajk = Work(ipInt2 + iAK + iA-1 +
     &                               (iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                     endif
                     Work(iMp2Lagr(iA,iJ,iSymJ)) =
     &                    Work(iMp2Lagr(iA,iJ,iSymJ)) -
     &                        Fac*Work(ip_Density(iSymI) +
     &                                 iOccOcc(iI,iK,iSymI))
     &                              * (Two*xikja - xiajk)
                     If (iSymI.eq.iSymJ) then
                        Work(iMp2Lagr(iA,iI,iSymJ)) =
     &                    Work(iMp2Lagr(iA,iI,iSymJ)) -
     &                        Fac*Work(ip_Density(iSymI) +
     &                                 iOccOcc(iJ,iK,iSymI))
     &                              * (Two*xiajk - xikja)
                     EndIf
*                     write(6,*) 'xikja', xikja
*                     write(6,*) 'xiajk', xiajk
                  EndDo
               EndDo

               If(iSymI.ne.iSymJ) then
                  iKA =(nFro(iSymI) + nOcc(iSymI)) *
     &                 (nOrb(iSymJ) + nDel(iSymJ))
                  iAK = nFro(iSymI) + nOcc(iSymI)
                  Do iK = 1, nFro(iSymJ) + nOcc(iSymJ)
                     Do iA = 1, nExt(iSymI) + nDel(iSymI)
                        xikja =Work(ipInt1 + iAK + iA-1 +
     &                            (iK-1)*(nOrb(iSymI)+nDel(iSymI)))
                        xiajk =Work(ipInt2 + iKA + iK-1 +
     &                               (iA-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                        Work(iMp2Lagr(iA,iI,iSymI)) =
     &                    Work(iMp2Lagr(iA,iI,iSymI)) -
     &                        Fac*Work(ip_Density(iSymJ) +
     &                                 iOccOcc(iJ,iK,iSymJ))
     &                              * (Two*xikja - xiajk)
                     EndDo
                  EndDo
               Endif
            Endif
         EndDo
      EndDo

      Return
      End
