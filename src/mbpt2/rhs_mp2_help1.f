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
      subroutine rhs_mp2_help1(iSymA, iSymB, iSymI, iSymJ)
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
*
*defining One etc.
#include "real.fh"
#include "trafo.fh"
#include "WrkSpc.fh"
#include "files_mbpt2.fh"
#include "mp2grad.fh"
#include "corbinf.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      iVirVir(i,j,k) = nFro(k) + nOcc(k) +  j-1 +
     &                (nOrb(k) + nDel(k))*
     &                (i + nFro(k) + nOcc(k) -1)
      iOccOcc(i,j,k) = nFro(k) + j-1 + (nOrb(k)+nDel(k))
     &                               * (i + nFro(k)-1)
      iVirDel(i,j,k) = nFro(k) + nOcc(k) + j-1 +
     &                (nOrb(k) + nDel(k)) *
     &                (i + nFro(k) + nOcc(k) + nExt(k) - 1)
      iOccFro(i,j,k) = j-1 + (nOrb(k) +nDel(k))
     &                               * (i + nFro(k) -1)
      iVirOcc(i,j,k) = j-1 + nFro(k) + nOcc(k) + (nOrb(k) + nDel(k))
     &                                *(i-1)
      iMp2Lagr(i,j,k) = ip_Mp2Lagr(k) +
     &                           j-1 + (nOcc(k)+nFro(k))*(i-1)
*                                                                      *
************************************************************************
*                                                                      *
*    First we read in blocks of all integrals of the form <I.|J.>
*    If the I and J is of the same symmetry we only need to
*    consider blocks with OrbJ < OrbI
      nJ = nOcc(iSymJ)
      Do iI = 1, nOcc(iSymI)
*        Check if we need to consider all combinations of orbitals
         if(iSymI.eq.iSymJ) nJ = iI
         Do iJ = 1, nJ
*            Write(6,*) 'iSymI, iSymJ', iSymI,iSymJ
*            Write(6,*) 'iI, iJ', iI, iJ
*           If only blocks with OrbJ < OrbI is used we need
*           to multiply the result with two for offdiag terms.
            fac_ij = One
            if((iI.ne.iJ).and.(iSymI.eq.iSymJ)) fac_ij = Two
*           If I is not equal J there is a symmetry J I B A
*           that is identical and not stored.
            if(iSymI.ne.iSymJ) fac_ij = Two
*           Needed?
*
*           Copy one AB-block of integrals to the workspace
            Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &                iI+nFro(iSymI),iJ+nFro(iSymJ),
     &                Work(ipInt1),Work(ipScr1))
            If(iSymI.ne.iSymJ) then
               Call Exch(iSymB,iSymI,iSymA,iSymJ,
     &                   iI+nFro(iSymI),iJ+nFro(iSymJ),
     &                   Work(ipInt2),Work(ipScr1))
            End If
*            Write(6,*) 'Rhs_mp2_help1'
*            Write(6,*) ' *  i,j = ',iI, iJ
*            Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
*     &                  nOrb(iSymA)+nDel(iSymA),
*     &                  nOrb(iSymB)+nDel(iSymB))
*            If(iSymI.ne.iSymJ) Then
*               Call RecPrt('Int2:','(8F10.6)',Work(ipInt2),
*     &                     nOrb(iSymB)+nDel(iSymB),
*     &                     nOrb(iSymA)+nDel(iSymA))
*            End If
*            write(6,*) ''
*
*           Calculates the EMP2-energy as a test
            iAB = nFro(iSymA) + nOcc(iSymA)  +
     &           (nFro(iSymB) + nOcc(iSymB))
     &          *(nOrb(iSymA) + nDel(iSymA))
            iBA = nFro(iSymB) + nOcc(iSymB)  +
     &           (nFro(iSymA) + nOcc(iSymA))
     &          *(nOrb(iSymB) + nDel(iSymB))
*
            Do iA = 1, nExt(iSymA)
            nB = nExt(iSymB)
            if(iSymA.eq.iSymB) nB = iA
               Do iB = 1, nB
                   EDiff = One/(work(mAdOcc(iSymI)+iI-1)
     &                       + work(mAdOcc(iSymJ)+iJ-1)
     &                       - work(mAdVir(iSymA)+iA-1)
     &                       - work(mAdVir(iSymB)+iB-1) )

                  xiajb=work(ipInt1 + iAB + iA-1 + (iB-1)*
     &                       (nOrb(iSymA)+nDel(iSymA)))
                  If(iSymI.eq.iSymJ) then
                     xibja=work(ipInt1 + iAB +
     &                          iB-1 + (iA-1)*(nOrb(iSymB)+nDel(iSymB)))
                  Else
                     xibja=work(ipInt2 + iBA +
     &                          iB-1 + (iA-1)*(nOrb(iSymB)+nDel(iSymB)))
                  End If
                  EMP2 = EMP2 + Ediff*Fac_ij * (xiajb*(Two*xiajb-xibja))
                  VECL2=VECL2 + Ediff**2*Fac_ij*
     &                  (xiajb*(Two*xiajb-xibja))
                  If(((iA.ne.iB).and.(iSymA.eq.iSymB))
     &                .or. (iSymI.ne.iSymJ)) Then
                     EMP2= EMP2+ Ediff*Fac_ij*(xibja*(Two*xibja-xiajb))
                     VECL2=VECL2 + Ediff**2*Fac_ij*
     &                     (xibja*(Two*xibja-xiajb))
                  End If
               End Do           !iOrbB
            End Do              !iOrbA
*                                                                      *
************************************************************************
*                                                                      *
*           P Virtual-Virtual
*                                                                      *
************************************************************************
*                                                                      *
*           Calculate the virtual-virtual perturbation to the
*           P-matrix
*           In this and all other things we calculate from now on
*           we need to do an extra loop over symmetries where J is
*           not equal to I since only symmetries where I > J is saved.
*           In the extra loops we do some special handling to get the
*           results from the symmetries with I < J.
*
               iAC = nFro(iSymA) + nOcc(iSymA)  +
     &              (nFro(iSymB) + nOcc(iSymB))
     &             *(nOrb(iSymA) + nDel(iSymA))
               iCA = nFro(iSymB) + nOcc(iSymB)  +
     &              (nFro(iSymA) + nOcc(iSymA))
     &             *(nOrb(iSymB) + nDel(iSymB))
               iBC = iAC
               iCB = iCA
               Do iC = 1, nExt(iSymB)
                  Do iA = 1, nExt(iSymA)
                     xiajc=work(ipInt1 + iAC + iA-1 +
     &                             (iC-1)*(nOrb(iSymA) + nDel(iSymA)))
                     If(ISymI.Eq.iSymJ) then
                        xicja=work(ipInt1 + iCA + iC-1 +
     &                       (iA-1)*(nOrb(iSymB)+nDel(iSymB)))
                     Else
                        xicja=work(ipInt2 + iCA + iC-1 +
     &                       (iA-1)*(nOrb(iSymB)+nDel(iSymB)))
                     End If
                     EDiffac = work(mAdOcc(iSymI)+iI-1)
     &                       + work(mAdOcc(iSymJ)+iJ-1)
     &                       - work(mAdVir(iSymA)+iA-1)
     &                       - work(mAdVir(iSymB)+iC-1)
                     T_ij = (Two*xiajc-xicja)/EDiffac
                     T_ji = (Two*xicja-xiajc)/EDiffac
                     Do iB = 1, nExt(iSymA)
                        xibjc=work(ipInt1 + iBC + iB-1 +
     &                       (iC-1)*(nOrb(iSymA) + nDel(iSymA)))
                        If(ISymI.Eq.iSymJ) then
                           xicjb=work(ipInt1 + iCB + iC-1 +
     &                          (iB-1)*(nOrb(iSymB) +nDel(iSymB)))
                        Else
                           xicjb=work(ipInt2 + iCB + iC-1 +
     &                          (iB-1)*(nOrb(iSymB) +nDel(iSymB)))
                        endif
                        EDiffbc = work(mAdOcc(iSymI)+iI-1)
     &                          + work(mAdOcc(iSymJ)+iJ-1)
     &                          - work(mAdVir(iSymB)+iC-1)
     &                          - work(mAdVir(iSymA)+iB-1)
                        Work(ip_Density(iSymA) +
     &                       iVirVir(iA,iB,iSymA)) =
     &                       Work(ip_Density(iSymA) +
     &                       iVirVir(iA,iB,iSymA)) +
     &                       Fac_ij * ((xibjc*T_ij) +
     &                       (xicjb*T_ji)) /EDiffbc
                         Work(ip_WDensity(iSymA) +
     &                       iVirVir(iB,iA,iSymA)) =
     &                       Work(ip_WDensity(iSymA) +
     &                       iVirVir(iB,iA,iSymA)) -
     &                       Fac_ij * ((xibjc*T_ij) +
     &                       (xicjb*T_ji))
                     End Do
*                     Write(6,*) 'IJABC',iI,iJ,iA,iB,iC
*                     Write(6,*) 'Symm', iSymI,iSymJ,iSymA,iSymB,iSymC
*                        write(6,*) 'EDiffBC', EDiffbc
*                        Write(6,*) 'EDiffAC', EDiffac
*                        write(6,*) 'xicja', xicja
*                        write(6,*) 'xicjb', xicjb
*                        write(6,*) 'xiajc', xiajc
*                        write(6,*) 'xibjc', xibjc
                     Do iBDel = 1, nDel(iSymA)
                        xibjc=work(ipInt1 + iBC +
     &                        iBDel-1 + nExt(iSymA)+
     &                        (iC-1)*(nOrb(iSymA) + nDel(iSymA)))
                        If(ISymI.Eq.iSymJ) then
                           xicjb=work(ipInt1 + iBC + iC-1 +
     &                          (iBDel + nExt(iSymA)-1)
     &                          *(nOrb(iSymB) +nDel(iSymB)))
                        Else
                           xicjb=work(ipInt2 + iCB + iC-1 +
     &                          (iBDel + nExt(iSymA)-1)
     &                          *(nOrb(iSymB) +nDel(iSymB)))

                        End If
                        EDiffbc = Work(madVir(iSymA)+iA-1)
     &                             - Work(mAdDel(iSymA)+iBDel-1)
                        Work(ip_Density(iSymA) +
     &                       iVirDel(iBDel,iA,iSymA)) =
     &                       Work(ip_Density(iSymA) +
     &                       iVirDel(iBDel,iA,iSymA)) +
     &                       Fac_ij *
     &                       ((xibjc*T_ij) +
     &                       (xicjb*T_ji)) /EDiffbc
                        Work(ip_WDensity(iSymA) +
     &                       iVirDel(iBDel,iA,iSymA)) =
     &                       Work(ip_WDensity(iSymA) +
     &                       iVirDel(iBDel,iA,iSymA)) -
     &                       Fac_ij *
     &                       ((xibjc*T_ij) +
     &                       (xicjb*T_ji))
*--------------------------------------------------------
*                        Write(6,*) 'index Dens', iVirDel(iBDel,iA,iSymA)
*                        Write(6,*) 'xicjb', xicjb
*                        Write(6,*) 'xibjc', xibjc
*                        Write(6,*) 'EDiffbc', EDiffbc
*                        Write(6,*) 'E_a', Work(madVir(iSymA)+iA-1)
*                        Write(6,*) 'E_B', Work(mAdDel(iSymA)+iBDel-1)
*                        Write(6,*) 'Tij', T_ij
*                        Write(6,*) 'Tji', T_ji
*                        Write(6,*) 'iB, iC', iBDel+nExt(iSymA), iC
*                        Write(6,*) 'iI, iJ', iI, iJ
*-------------------------------------
                        End Do
*                         work(indexW) = work(indexW) -
*     &                     Fac_ij*EDiffbc*(Two*(xijac*xijbc+
*     &                     xijca*xijcb) - xijac*xijcb-xijca*xijbc)
*-----------------------------------------------------------------------
*                        write (6,*) 'The Denom is', EDenom
*                        write(6,*) 'xicja', xicja
*                        write(6,*) 'xicjb', xicjb
*                        write(6,*) 'xiajc', xiajc
*                        write(6,*) 'xibjc', xibjc
*                        write(6,*) 'this is with the unmirrored'
*------------------------------------------------------------------------
*
                  enddo         !iOrbB
               enddo            !iOrbA
*
*
*              Here is the extra loop.
               if (iSymA.ne.iSymB) then
               iAC = nFro(iSymB) + nOcc(iSymB)  +
     &              (nFro(iSymA) + nOcc(iSymA))
     &             *(nOrb(iSymB) + nDel(iSymB))
               iCA = nFro(iSymA) + nOcc(iSymA) +
     &              (nFro(iSymB) + nOcc(iSymB))
     &             *(nOrb(iSymA) + nDel(iSymA))
               iBC = iAC
               iCB = iCA
               Do iC = 1, nExt(iSymA)
                  Do iA = 1, nExt(iSymB)
                     EDiffac = work(mAdOcc(iSymI)+iI-1)
     &                       + work(mAdOcc(iSymJ)+iJ-1)
     &                       - work(mAdVir(iSymB)+iA-1)
     &                       - work(mAdVir(iSymA)+iC-1)
                     xiajc=work(ipInt1 + iCA + iC-1 +
     &                             (iA-1)*(nOrb(iSymA) + nDel(iSymA)))
                     xicja=work(ipInt2 + iAC + iA-1 +
     &                             (iC-1)*(nOrb(iSymB) + nDel(iSymB)))
                     T_ij = (Two*xiajc-xicja)/EDiffac
                     T_ji = (Two*xicja-xiajc)/EDiffac
                     Do iB = 1,nExt(iSymB)
                        EDiffbc = work(mAdOcc(iSymI)+iI-1)
     &                          + work(mAdOcc(iSymJ)+iJ-1)
     &                          - work(mAdVir(iSymA)+iC-1)
     &                          - work(mAdVir(iSymB)+iB-1)
                        xibjc=work(ipInt1 + iCB + iC-1 +
     &                             (iB-1)*(nOrb(iSymA) + nDel(iSymA)))
                        xicjb=work(ipInt2 + iBC + iB-1 +
     &                             (iC-1)*(nOrb(iSymB) + nDel(iSymB)))
                        Work(ip_Density(iSymB) + iVirVir(iA,iB,iSymB)) =
     &                  Work(ip_Density(iSymB) + iVirVir(iA,iB,iSymB)) +
     &                                 Fac_ij  *
     &                                ((xibjc*T_ij) +
     &                                 (xicjb*T_ji)) / EDiffbc
                        Work(ip_WDensity(iSymB) + iVirVir(iB,iA,iSymB))=
     &                  Work(ip_WDensity(iSymB) + iVirVir(iB,iA,iSymB))-
     &                                 Fac_ij  *
     &                                ((xibjc*T_ij) +
     &                                 (xicjb*T_ji))

                     End Do
                     Do iBDel = 1, nDel(iSymB)
                        xibjc=work(ipInt1 + iBC +
     &                             iBDel-1 + nExt(iSymB) +
     &                        (iC-1)
     &                       *(nOrb(iSymB) + nDel(iSymB)))
                        xicjb=work(ipInt2 + iCB + iC-1 +
     &                             (iBdel + nExt(iSymB)-1)*
     &                             (nOrb(iSymA) + nDel(iSymA)))
                        EDiffbc = Work(madVir(iSymB)+iA-1)
     &                          - Work(mAdDel(iSymB)+iBDel-1)
                        Work(ip_Density(iSymA) +
     &                       iVirDel(iBDel,iA,iSymB)) =
     &                       Work(ip_Density(iSymB) +
     &                       iVirDel(iBDel,iA,iSymB)) +
     &                       Fac_ij *
     &                       ((xibjc*T_ij) +
     &                       (xicjb*T_ji)) /EDiffbc
                        Work(ip_WDensity(iSymA) +
     &                       iVirDel(iBDel,iA,iSymB)) =
     &                       Work(ip_WDensity(iSymB) +
     &                       iVirDel(iBDel,iA,iSymB)) -
     &                       Fac_ij *
     &                       ((xibjc*T_ij) +
     &                       (xicjb*T_ji))
                     End Do
*-----------------------------------------------------------------------
*                        write (6,*) 'The Denom is', EDenom
*                        write(6,*) 'xicja', xicja
*                        write(6,*) 'xicjb', xicjb
*                        write(6,*) 'xiajc', xiajc
*                        write(6,*) 'xibjc', xibjc
*                        write(6,*) 'this is with the mirrored symmetry'
*------------------------------------------------------------------------


                  enddo         !iOrbB
               enddo            !iOrbA
            endif
*                                                                      *
************************************************************************
*                                                                      *
*           Lagrangian term 2
*                                                                      *
************************************************************************
*                                                                      *
            iAB = nFro(iSymA) + nOcc(iSymA) +
     &           (nFro(iSymB) + nOcc(iSymB)) *
     &           (nOrb(iSymA) + nDel(iSymA))
            iBA = nFro(iSymB) + nOcc(iSymB) +
     &           (nFro(iSymA) + nOcc(iSymA)) *
     &           (nOrb(iSymB) + nDel(iSymB))
            iKB = (nFro(iSymB) + nOcc(iSymB))
     &            *(nOrb(iSymA) + nDel(iSymA))
            iBK = nFro(iSymB) + nOcc(iSymB)
            Do iA = 1, nExt(iSymA)
               Do iB = 1,nExt(iSymB)
                  EDenom =(work(mAdOcc(iSymI)+iI-1)
     &                   + work(mAdOcc(iSymJ)+iJ-1)
     &                   - work(mAdVir(iSymA)+iA-1)
     &                   - work(mAdVir(iSymB)+iB-1) )
                  xiajb=work(ipInt1 + iAB + iA-1 +
     &                  (iB-1)*(nOrb(iSymA) + nDel(iSymA)))
                  If(isymI.Eq.iSymJ) Then
                     xibja=work(ipInt1 + iBA + iB-1 +
     &                          (iA-1)*(nOrb(iSymB)+nDel(iSymB)))
                  Else
                     xibja=work(ipInt2 + iBA + iB-1 +
     &                          (iA-1)*(nOrb(iSymB)+nDel(iSymB)))
                  End If
                  Tij= (Two*xiajb-xibja)/EDenom
                  Tji= (Two*xibja-xiajb)/EDenom
                  Do iK = 1 , nFro(iSymA) + nOcc(iSymA)
                     xikjb=work(ipInt1 + iKB + iK-1 +
     &                     (iB-1)*(nOrb(iSymA) + nDel(iSymA)))
                     If(isymI.Eq.iSymJ) Then
                        xibjk=work(ipInt1 + iBK + iB-1 +
     &                             (iK-1)*(nOrb(iSymB)+nDel(iSymB)))
                     Else
                        xibjk=work(ipInt2 + iBK + iB-1 +
     &                             (iK-1)*(nOrb(iSymB)+nDel(iSymB)))
                     Endif
                     Work(iMp2Lagr(iA,iK,iSymA)) =
     &                 Work(iMp2Lagr(iA,iK,iSymA)) + Fac_ij*
     &                               (Tij*Xikjb+Tji*Xibjk)
                     Work(ip_WDensity(iSymA) + iVirOcc(iK,iA,iSymA))=
     &                 Work(ip_WDensity(iSymA) + iVirOcc(iK,iA,iSymA))-
     &                     2.0d0*Fac_ij*(Tij*Xikjb+Tji*Xibjk)
*-------------------------------------------------------
*                     write(6,*) 'The Denom is', EDenom
*                     write(6,*) 'Tij  ', Tij
*                     write(6,*) 'Tji  ', Tji
*                     write(6,*) 'xiajb', xiajb
*                     write(6,*) 'xiajb', xiajb
*                     write(6,*) 'xibja', xibja
*                     write(6,*) 'xikjb', xikjb
*                     write(6,*) 'xibjk', xibjk
*-----------------------------------------------------------------------
                  End Do         !OrbK
               End Do            !OrbA
            End Do               !OrbB
*
*     In the case that the symmetries A and B are not the same
*     only the case where A is greater than B is written on disk
*     so we must do an extra loop for where B is greater than A
*
            If(iSymA.ne.iSymB) Then
               iAB = nFro(iSymB) + nOcc(iSymB) +
     &              (nFro(iSymA) + nOcc(iSymA))
     &             *(nOrb(iSymB) + nDel(iSymB))
               iBA = nFro(iSymA) + nOcc(iSymA) +
     &              (nFro(iSymB) + nOcc(iSymB))
     &             *(nOrb(iSymA) + nDel(iSymB))
               iKB = (nFro(iSymA) +nOcc(iSymA))
     &             *(nOrb(iSymB)+nDel(iSymB))
               iBK = nFro(iSymA) + nOcc(iSymA)
               Do iA = 1, nExt(iSymB)
                  Do iB = 1, nExt(iSymA)
                     EDenom = (work(mAdOcc(iSymI)+iI-1)
     &                       + work(mAdOcc(iSymJ)+iJ-1)
     &                       - work(mAdVir(iSymB)+iA-1)
     &                       - work(mAdVir(iSymA)+iB-1) )
                     xiajb=work(ipInt1 + iBA + iB-1 +
     &                             (iA-1)*(nOrb(iSymA)+nDel(iSymA)))
                     xibja=work(ipInt2 + iAB + iA-1 +
     &                             (iB-1)*(nOrb(iSymB)+nDel(iSymB)))
                     Tij= (Two*xiajb-xibja)/EDenom
                     Tji= (Two*xibja-xiajb)/EDenom
                     Do iK = 1 , nFro(iSymB) + nOcc(iSymB)
                        xikjb=work(ipInt1 + iBK + iB-1 +
     &                             (iK-1)*(nOrb(iSymA)+nDel(iSymA)))

                        xibjk=work(ipInt2 + iKB + iK-1 +
     &                             (iB-1)*(nOrb(iSymB)+nDel(iSymB)))
*
                        Work(iMp2Lagr(iA,iK,iSymB)) =
     &                    Work(iMp2Lagr(iA,iK,iSymB)) + Fac_ij*
     &                                 (Tij*Xikjb+Tji*Xibjk)
                        Work(ip_WDensity(iSymB) + iVirOcc(iK,iA,iSymB))=
     &                  Work(ip_WDensity(iSymB) + iVirOcc(iK,iA,iSymB))-
     &                     2.0d0*Fac_ij*(Tij*Xikjb+Tji*Xibjk)
*-----------------------------------------------------------------------
*                        write (6,*) 'The Denom is', EDenom
*                        write(6,*) 'Tij  ', Tij
*                        write(6,*) 'Tji  ', Tji
*                        write(6,*) 'xijab', xijab
*                        write(6,*) 'xijba', xijba
*                        write(6,*) 'xijkb', xijkb
*                        write(6,*) 'xijbk', xijbk
*                        write(6,*) 'this is with the mirrored symmetry'
*------------------------------------------------------------------------
                     EndDo      !OrbK
                  EndDo         !OrbA
               EndDo            !OrbB
            EndIf
*
*           When we are done with a block we free the memory again.
*
         End Do                  !OrbI
      End Do
*
******************************************************************************************
*
*             Terms that need <A.|B.>-integrals
*
******************************************************************************************
*
*
*
      nB = nExt(iSymB)
      Do iA = 1, nExt(iSymA)
*        Decide i the AB-matrix is symmetric, if so
*        only loop over lower triangle
         If(iSymA.Eq.iSymB) nB = iA

         Do iB = 1, nB
*            Write(6,*) 'iSymA, iSymB', iSymA, iSymB
*            Write(6,*) 'iA, iB', iA, iB
            fac_ab = One
            if((iA.Ne.iB).and.(iSymA.eq.iSymB)) fac_ab = Two
*           If I is not equal J there is a symmetry J I B A
*           that is identical and not stored.
            if(iSymA.Ne.iSymB) fac_ab = Two
*
            Call Exch(iSymI,iSymA,iSymJ,iSymB,
     &                iA + nOcc(iSymA)+nFro(iSymA),
     &                iB + nOcc(iSymB)+nFro(iSymB),
     &                Work(ipInt1), Work(ipScr1))
            If(iSymA.Ne.iSymB) Then
               Call Exch(iSymJ,iSymA,iSymI,iSymB,
     &                   iA + nOcc(iSymA)+nFro(iSymA),
     &                   iB + nOcc(iSymB)+nFro(iSymB),
     &                   Work(ipInt2),Work(ipScr1))
            End If
*            Write(6,*)
*            Write(6,*) ' *  A,B = ',iA,iB
*            Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
*     &                  nOrb(iSymI)+nDel(iSymI),
*     &                  nOrb(iSymJ)+nDel(iSymJ))
*            If (iSymA.Ne.iSymB) Then
*               Call RecPrt('Int2:','(8F10.6)',Work(ipInt2),
*     &                     nOrb(iSymJ) + nDel(iSymJ),
*     &                     nOrb(iSymI) + nDel(iSymI))
*            End If
*            Write(6,*)

*********************************************************************
*
*                   Calculate Pij
*
*********************************************************************
*
*           We start by taking i and j in symI and k in symJ
            iIK= nFro(iSymI) +
     &           nFro(iSymJ) * (nOrb(iSymI)+nDel(iSymI))
            iKI= nFro(iSymJ) +
     &           nFro(iSymI) * (nOrb(iSymJ)+nDel(iSymJ))
            iJK= iIK
            iKJ= iKI
            Do iK = 1, nOcc(iSymJ)
               Do iI=1, nOcc(iSymI)
                   EDiffik = work(mAdOcc(iSymI) + iI-1)
     &                     + work(mAdOcc(iSymJ)+ iK-1)
     &                     - work(mAdVir(iSymA)+iA-1)
     &                     - work(mAdVir(iSymB)+iB-1)
                   xiakb=Work(ipInt1 +iIK + iI-1 + (iK-1)*
     &                          (nOrb(iSymI) + nDel(iSymI)))
                   If(iSymA.Eq.iSymB) then
                      xkaib=Work(ipInt1+iKI +iK-1 +(iI-1)*
     &                             (nOrb(iSymJ) + nDel(iSymJ)))
                   Else
                      xkaib=Work(ipInt2+iKI+iK-1 +(iI-1)*
     &                             (nOrb(iSymJ) + nDel(iSymJ)))
                   End If
                   T_ab = (Two*xiakb-xkaib)/EDiffik
                   T_ba = (Two*xkaib-xiakb)/EDiffik
                   Do iJ= 1, nOcc(iSymI)
*                    Calculating the denominator

                     EDiffjk = work(mAdOcc(iSymI)+iJ-1)
     &                       + work(mAdOcc(iSymJ)+iK-1)
     &                       - work(mAdVir(iSymA)+iA-1)
     &                       - work(mAdVir(iSymB)+iB-1)

                     xjakb=Work(ipInt1 +iJK + iJ-1 + (iK-1)*
     &                          (nOrb(iSymI) + nDel(iSymI)))
                     if(iSymA.Eq.iSymB) then
                        xkajb=Work(ipInt1+iKJ +iK-1 +(iJ-1)*
     &                             (nOrb(iSymJ) + nDel(iSymJ)))
                     else
                        xkajb=Work(ipInt2+iKJ+iK-1 +(iJ-1)*
     &                             (nOrb(iSymJ) + nDel(iSymJ)))
                     endif
*---------------------------------------------------------------------
*                     write (6,*) 'The Denom is', EDenom
*                     write(6,*) 'xikab', xikab
*                     write(6,*) 'xjkab', xjkab
*                     write(6,*) 'xkiab', xkiab
*                     write(6,*) 'xkjab', xkjab
*-----------------------------------------------------------------------
                     Work(ip_Density(iSymI) + iOccOcc(iI,iJ,iSymI)) =
     &                  Work(ip_Density(iSymI) + iOccOcc(iI,iJ,iSymI)) -
     &                                         Fac_ab *
     &                                  (xjakb*T_ab +
     &                                   xkajb*T_ba)/EDiffjk
                     Work(ip_WDensity(iSymI) + iOccOcc(iI,iJ,iSymI)) =
     &                  Work(ip_WDensity(iSymI) + iOccOcc(iI,iJ,iSymI))-
     &                                         Fac_ab *
     &                                  (xjakb*T_ab +
     &                                   xkajb*T_ba)
*------------------------------------------------------------------------
                  End Do
                  Do iJFroz = 1, nFro(iSymI)
                     xjakb=Work(ipInt1 + iJFroz-1 +
     &                    (iK+nFro(iSymJ)-1)*
     &                    (nOrb(iSymI) + nDel(iSymI)))
                     if(iSymI.Eq.iSymJ) then
                        xkajb=Work(ipInt1 +iK-1
     &                       + nFro(iSymJ) +(iJFroz-1)*
     &                       (nOrb(iSymJ) + nDel(iSymJ)))
                     else
                        xkajb=Work(ipInt2+
     &                       iK-1 + nFro(iSymJ) +(iJFroz-1)*
     &                       (nOrb(iSymJ) + nDel(iSymJ)))
                     endif
                     EDiffjk = Work(madOcc(iSymI)+iI-1)
     &                       - Work(mAdFro(iSymI)+iJFroz-1)
                        Work(ip_Density(iSymI) +
     &                       iOccFro(iI,iJFroz,iSymI)) =
     &                       Work(ip_Density(iSymI) +
     &                       iOccFro(iI,iJFroz,iSymI)) +
     &                       Fac_ab *
     &                       ((xjakb*T_ab) +
     &                       (xkajb*T_ba)) /EDiffjk
                        Work(ip_WDensity(iSymI) +
     &                       iOccFro(iI,iJFroz,iSymI)) =
     &                       Work(ip_WDensity(iSymI) +
     &                       iOccFro(iI,iJFroz,iSymI)) -
     &                       Fac_ab *
     &                       ((xjakb*T_ab) +
     &                       (xkajb*T_ba))
                  End Do
               EndDo
            EndDo
*           And then we do the same thing with i and j in symJ and k in symI
            if(iSymI.ne.iSymJ) then
               iIK = nFro(iSymJ) +
     &               nFro(iSymI) * (nOrb(iSymJ)+nDel(iSymJ))
               iKI = nFro(iSymI) +
     &               nFro(iSymJ) * (nOrb(iSymI)+nDel(iSymI))
               iJK = iIK
               iKJ = iKI
                Do iK = 1, nOcc(iSymI)
                  Do iI = 1, nOcc(iSymJ)
                     EDiffik = work(mAdOcc(iSymJ)+iI-1)
     &                       + work(mAdOcc(iSymI)+iK-1)
     &                       - work(mAdVir(iSymA)+iA-1)
     &                       - work(mAdVir(iSymB)+iB-1)
                     xiakb = Work(ipInt1 + iKI +  iK-1 +
     &                               (iI-1)*(nOrb(iSymI)+nDel(iSymI)))
                     xkaib = Work(ipInt2 + iIK + iI-1 +
     &                               (iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                     T_ab = (Two*xiakb-xkaib)/EDiffik
                     T_ba = (Two*xkaib-xiakb)/EDiffik
                     Do iJ = 1, nOcc(iSymJ)
*                       Calculating the denominator

                        EDiffjk = work(mAdOcc(iSymJ)+iJ-1)
     &                          + work(mAdOcc(iSymI)+iK-1)
     &                          - work(mAdVir(iSymA)+iA-1)
     &                          - work(mAdVir(iSymB)+iB-1)
*
                        xjakb = Work(ipInt1 + iKJ + iK-1 +
     &                               (iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
                        xkajb = Work(ipInt2 + iJK + iJ-1 +
     &                               (iK-1)*(nOrb(iSymJ)+nDel(iSymJ)))
*
                        Work(ip_Density(iSymJ) + iOccOcc(iI,iJ,iSymJ)) =
     &                  Work(ip_Density(iSymJ) + iOccOcc(iI,iJ,iSymJ)) -
     &                               Fac_ab *
     &                              (xjakb*T_ab +
     &                               xkajb*T_ba) /EDiffjk
                        Work(ip_WDensity(iSymJ)+iOccOcc(iI,iJ,iSymJ)) =
     &                   Work(ip_WDensity(iSymJ)+iOccOcc(iI,iJ,iSymJ)) -
     &                               Fac_ab *
     &                              (xjakb*T_ab +
     &                               xkajb*T_ba)


*--------------------------------------------------------------------
*                     write (6,*) 'The Denom is', EDenom
*                     write(6,*) 'xikab', xikab
*                     write(6,*) 'xjkab', xjkab
*                     write(6,*) 'xkiab', xkiab
*                     write(6,*) 'xkjab', xkjab
*--------------------------------------------------------------------
                     EndDo
                     Do iJFroz = 1, nFro(iSymJ)
                        xjakb=Work(ipInt1 + iK+nFro(iSymI)-1 +
     &                            (iJFroz-1)*
     &                            (nOrb(iSymI) + nDel(iSymI)))
                        xkajb=Work(ipInt2+
     &                             iJFroz-1 +(iK + nFro(iSymI)-1)*
     &                            (nOrb(iSymJ) + nDel(iSymJ)))
                        EDiffjk = Work(madOcc(iSymJ)+iI-1)
     &                          - Work(mAdFro(iSymJ)+iJFroz-1)
                        Work(ip_Density(iSymJ) +
     &                       iOccFro(iI,iJFroz,iSymJ)) =
     &                       Work(ip_Density(iSymJ) +
     &                       iOccFro(iI,iJFroz,iSymJ)) +
     &                       Fac_ab *
     &                       ((xjakb*T_ab) +
     &                       (xkajb*T_ba)) /EDiffjk
                        Work(ip_WDensity(iSymJ) +
     &                       iOccFro(iI,iJFroz,iSymJ)) =
     &                       Work(ip_WDensity(iSymJ) +
     &                       iOccFro(iI,iJFroz,iSymJ)) -
     &                       Fac_ab *
     &                       ((xjakb*T_ab) +
     &                       (xkajb*T_ba))
                     End Do
                  EndDo
               EndDo
            EndIf

***************************************************************************
*
*                  Lagrangian term 1
*
***************************************************************************
*
            iCJ = nFro(iSymI) + nOcc(iSymI) +
     &            nFro(iSymJ)*
     &           (nOrb(iSymI)+nDel(iSymI))
            iJC = nFro(iSymJ) +
     &           (nFro(iSymI) + nOcc(iSymI)) *
     &           (nOrb(iSymJ) + nDel(iSymJ))
            iIJ = nFro(iSymI) +
     &            nFro(iSymJ) *
     &           (nOrb(iSymI) + nDel(iSymI))
            IJI = nFro(iSymJ) +
     &            nFro(iSymI) *
     &           (nOrb(iSymJ) + nDel(iSymJ))
            Do iC = 1, nExt(iSymI) + nDel(iSymI)
               Do iI = 1, nOcc(iSymI)
                  Do iJ = 1, nOcc(iSymJ)
                     EDenom =  One/(work(mAdOcc(iSymI)+iI-1)
     &                    + work(mAdOcc(iSymJ)+iJ-1)
     &                    - work(mAdVir(iSymA)+iA-1)
     &                    - work(mAdVir(iSymB)+iB-1) )
                     xaibj = Work(ipInt1 + iIJ + iI-1 +
     &                            (iJ-1)*(nOrb(iSymI)+nDel(iSymI)))
                     xacbj = Work(ipInt1 + iCJ + iC-1 +
     &                            (iJ-1) * (nOrb(iSymI)+nDel(iSymI)))
                     if(iSymA.eq.iSymB) then
                        xbiaj = Work(ipInt1 + iJI + iJ-1
     &                        + (iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                        xbcaj = Work(ipInt1 + iJC + iJ-1
     &                        + (iC-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                     else
                        xbiaj = Work(ipInt2 + iJI + iJ-1
     &                        + (iI-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                        xbcaj = Work(ipInt2 + iJC + iJ-1
     &                        + (iC-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                     endif
                     Work(iMp2Lagr(iC,iI+nFro(iSymI),iSymI)) =
     &                    Work(iMp2Lagr(iC,iI+nFro(iSymI),iSymI)) +
     &                              EDenom * Fac_ab *(-Two*(xaibj*xacbj
     &                              + xbiaj*xbcaj) + xaibj*xbcaj +
     &                                xbiaj*xacbj)

*-------------------------------------------------------------------
*                     write(*,*) 'EDenom', EDenom
*                     write(*,*) 'Fac_ab', Fac_ab
*                     write(*,*) 'xaibj', xaibj
*                     write(*,*) 'xacbj', xacbj
*                     write(*,*) 'xbiaj', xbiaj
*                     write(*,*) 'xbcaj', xbcaj
*-------------------------------------------------------------------
                  EndDo
               EndDo
            EndDo
*           And now we switch symmetries and go again
            If(iSymA.Ne.iSymB) Then
               iCJ = nFro(iSymJ) + nOcc(iSymJ) +
     &               nFro(iSymI)*
     &              (nOrb(iSymJ)+nDel(iSymJ))
                 iJC = nFro(iSymI) +
     &                (nFro(iSymJ) + nOcc(iSymJ)) *
     &                (nOrb(iSymI) + nDel(iSymI))
                 iIJ = nFro(iSymJ) +
     &                 nFro(iSymI) *
     &                (nOrb(iSymJ) + nDel(iSymJ))
                 IJI = nFro(iSymI) +
     &                 nFro(iSymJ) *
     &                (nOrb(iSymI) + nDel(iSymI))
               Do iC = 1, nExt(iSymJ) + nDel(iSymJ)
                  Do iI = 1, nOcc(iSymJ)
                     Do iJ = 1, nOcc(iSymI)
                        EDenom =  One/(work(mAdOcc(iSymJ)+iI-1)
     &                          + work(mAdOcc(iSymI)+iJ-1)
     &                          - work(mAdVir(iSymA)+iA-1)
     &                          - work(mAdVir(iSymB)+iB-1) )
                        xaibj = Work(ipInt1 + iJI + iJ-1
     &                        + (iI-1)*(nOrb(iSymI)+nDel(iSymI)))
                        xacbj = Work(ipInt1 + iJC + iJ-1
     &                        + (iC-1)*(nOrb(iSymI)+nDel(iSymI)))
                        xajbi = Work(ipInt2 + iIJ + iI-1
     &                        + (iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                        xajbc = Work(ipInt2 + iCJ + iC-1
     &                        + (iJ-1)*(nOrb(iSymJ)+nDel(iSymJ)))
                        Work(iMp2Lagr(iC,iI+nFro(iSymJ),iSymJ)) =
     &                    Work(iMp2Lagr(iC,iI+nFro(iSymJ),iSymJ)) +
     &                              EDenom * Fac_ab *(-Two*(xaibj*xacbj
     &                              + xajbi*xajbc) + xaibj*xajbc +
     &                              xajbi*xacbj)

*-----------------------------------------------------------------
*                        Write(*,*) 'EDenom', EDenom
*                        write(*,*) 'Fac_ab', Fac_ab
*                        write(6,*) 'xijab', xijab
*                        write(6,*) 'xcjab', xcjab
*                        write(6,*) 'xjiab', xjiab
*                        write(6,*) 'xjcab', xjcab
*-------------------------------------------------------------------
                  End Do
               End Do
            End Do
         End If
*
*
*
         EndDo
      EndDo

      Return
      End
