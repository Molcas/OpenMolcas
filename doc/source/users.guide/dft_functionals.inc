DFT functionals:
::::::::::::::::

Below is listed the keywords for the DFT functionals currently implemented in the package.

.. class:: keywordlist

:kword:`LSDA`, :kword:`LDA`, :kword:`SVWN`
  Vosko, Wilk, and Nusair 1980 correlation functional fitting the
  RPA solution to the uniform electron gas :cite:`Vosko:80` (functional III in the paper).

:kword:`LSDA5`, :kword:`LDA5`, :kword:`SVWN5`
  Functional V from the VWN80 paper :cite:`Vosko:80` which fits the Ceperley--Alder solution to the uniform electron gas.

:kword:`HFB`
  Becke's 1988 exchange functional which includes
  the Slater exchange along with corrections involving the gradient of the density :cite:`Becke:88`.

:kword:`HFS`
  :math:`\rho^{4/3}` with the theoretical coefficient of 2/3 also known as Local Spin Density exchange
  :cite:`Hohenberg:64,Kohn:65,Slater:74`.

:kword:`HFB86`
  Becke's 1986 two-parameter exchange functional which includes
  the Slater exchange along with corrections involving the gradient of the density :cite:`B86,B07`.

:kword:`HFO`
  Handy's stand-alone OPTX exchange functional :cite:`OPTX`.

  .. :kword:`XALPHA`
       :math:`\rho^{4/3}` with the empirical coefficient of 0.7 usually used when this exchange functional
       is used without correlation functional :cite:`Hohenberg:64,Kohn:65,Slater:74`.

:kword:`BLYP`
  Becke's 1988 exchange functional which includes
  the Slater exchange along with corrections involving the gradient of the density :cite:`Becke:88`.
  Correlation functional of Lee, Yang, and Parr, which includes both local and non-local terms :cite:`Lee:88,Miehlich:89`.

:kword:`BPBE`
  Becke's 1988 exchange functional which includes
  the Slater exchange along with corrections involving the gradient of the density :cite:`Becke:88` ,
  combined with the GGA correlation functional by Perdew, Burke and Ernzerhof :cite:`PBE:96`.

:kword:`B3LYP`
  .. compound::

    Becke's 3 parameter functional :cite:`Becke:93` with the form

    .. math:: A E_{\text{X}}^{\text{Slater}} + (1-A) E_{\text{X}}^{\text{HF}} + B \Delta E_{\text{X}}^{\text{Becke}} + E_{\text{C}}^{\text{VWN}} + C \Delta E_{\text{C}}^{\text{non-local}}

    where the non-local correlation functional is the LYP functional and the VWN is functional III (not functional V).
    The constants :math:`A`, :math:`B`, :math:`C` are those determined by Becke by fitting to the G1 molecule set, namely
    :math:`A`\=0.80, :math:`B`\=0.72, :math:`C`\=0.81.

:kword:`B3LYP5`
  .. compound::

    Becke's 3 parameter functional :cite:`Becke:93` with the form

    .. math:: A E_{\text{X}}^{\text{Slater}} + (1-A) E_{\text{X}}^{\text{HF}} + B \Delta E_{\text{X}}^{\text{Becke}} + E_{\text{C}}^{\text{VWN}} + C \Delta E_{\text{C}}^{\text{non-local}}

    where the non-local correlation functional is the LYP functional and the VWN is functional V.
    The constants :math:`A`, :math:`B`, :math:`C` are those determined by Becke by fitting to the G1 molecule set, namely
    :math:`A`\=0.80, :math:`B`\=0.72, :math:`C`\=0.81.

:kword:`B2PLYP`
  Grimme's double-hybrid density functional :cite:`B2PLYP` based on Becke's 1988 exchange and LYP correlation GGA functionals with the form

  .. math:: A E_{\text{X}}^{\text{Slater}} + (1-A) E_{\text{X}}^{\text{HF}} + A \Delta E_{\text{X}}^{\text{Becke}} + C E_{\text{C}}^{\text{LYP}} + (1-C) E_{\text{C}}^{\text{PT2}}

  The constants :math:`A`, and :math:`C` are :math:`A`\=0.47, :math:`C`\=0.73. The SCF program computes only the DFT part of the B2PLYP energy.
  In order to get the PT2 term, one has to run the MBPT2 program on converged B2PLYP orbitals, and scale the MP2
  correlation energy by the factor :math:`(1-C)`\=0.27.

:kword:`B86LYP`
  Becke's 1986 exchange :cite:`B86,B07` functional combined with the LYP correlation :cite:`Lee:88,Miehlich:89`.

:kword:`BWig`
  Becke's 1988 GGA exchange functional combined with the local Wigner correlation functional :cite:`GillWig`.

:kword:`GLYP`
  Gill's 1996 GGA exchange functional :cite:`G96` combined with the combined with the LYP correlation :cite:`Lee:88,Miehlich:89`.

:kword:`OLYP`
  Handy's OPTX exchange functional :cite:`OPTX` combined with the LYP correlation :cite:`Lee:88,Miehlich:89`.

:kword:`OPBE`
  Handy's OPTX exchange functional :cite:`OPTX` combined with the PBE correlation:cite:`PBE:96`.

:kword:`O3LYP`
  A hybrid density functional based on the OPTX exchange :cite:`O3LYP` , with the form

  .. math:: A E_{\text{X}}^{\text{HF}} + B E_{\text{X}}^{\text{Slater}} + C \Delta E_{\text{X}}^{\text{OPTX}} + 0.19 E_{\text{C}}^{\text{VWN}} + 0.81 \Delta E_{\text{C}}^{\text{LYP}}

  The constants :math:`A`, :math:`B`, :math:`C` are as provided in Ref. :cite:`O3LYP`: :math:`A`\=0.1161, :math:`B`\=0.9262, :math:`C`\=0.8133.

:kword:`KT3`
  The exchange-correlation functional by Keal and Tozer, 2004 :cite:`KT2,KT3`.

:kword:`TLYP`
  .. compound::

    .. math:: E_{\text{X}}^{\text{HF}} + E_{\text{C}}^{\text{non-local}}

    where the non-local correlation functional is the LYP functional.

:kword:`PBE`
  The Perdew, Burke, Ernzerhof GGA functional 1996 :cite:`PBE:96`.

:kword:`PBE0`
  The Perdew, Ernzerhof, Burke non-empirical hybrid functional 1996 :cite:`PBE0:96`.

:kword:`PBEsol`
  The Perdew et al. 2008 modification of PBE for solids.

  .. :cite:`PBEsol`

:kword:`RGE2`
  The regularized gradient approximation (RGE2) exchange functional by Ruzsinszky, Csonka, and Scuseria, 2009 that contains
  higher-power s terms in the exchange functional, as compared to the PBEsol. It is coupled with the PBEsol correlation :cite:`RGE2`.

:kword:`PTCA`
  The correlation functional by Tognetti, Cortona, and Adamo combined with the PBE exchange :cite:`TCA`.

:kword:`SSB`
  The exchange functional SSB-sw by Swart, Solà, and Bickelhaupt, 2008 :cite:`SSB` that switches between the OPTX exchange
  for small values of the reduced density gradient and the PBE exchange for the large ones. It is combined with the PBE correlation functional.

:kword:`M06`
  The M06 functional of the Minnesota 2006 class of functionals by Zhao and Truhlar :cite:`ZT1,ZT2,ZT3,ZT4`.

:kword:`M06L`
  The M06-L functional of the Minnesota 2006 class of functionals by Zhao and Truhlar :cite:`ZT1,ZT2,ZT3,ZT4`.

:kword:`M06HF`
  The M06-HF functional of the Minnesota 2006 class of functionals by Zhao and Truhlar :cite:`ZT1,ZT2,ZT3,ZT4`.

:kword:`M062X`
  The M06-2X functional of the Minnesota 2006 class of functionals by Zhao and Truhlar :cite:`ZT1,ZT2,ZT3,ZT4`.
