/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/

window.MathJax = {
  // For MathJax 2
  "HTML-CSS": {
    scale: 90,
    preferredFont: "STIX",
    webFont: "STIX-Web",
  },
  SVG: {
    scale: 90,
    font: "STIX-Web",
  },
  // For MathJax 3
  chtml: {
    scale: 1,
    mtextInheritFont: true,
    matchFontHeight: true,
  },
  svg: {
    scale: 1,
    mtextInheritFont: true,
    matchFontHeight: true,
  },
  tex: {
    extensions: ["AMSmath.js","AMSsymbols.js","mhchem.js"],
    macros: {
      mat: ["\\boldsymbol{\#1}", 1],
      sign: ["\\operatorname{sign}", 0],
      Tr: ["\\operatorname{Tr}", 0],
      abs: ["\\operatorname{abs}", 0],
      bra: ["\\left<\#1\\right|", 1],
      ket: ["\\left|\#1\\right>", 1],
      braket: ["\\left<\#1\\middle|\#2\\right>", 2],
      braopket: ["\\left<\#1\\middle|\#2\\middle|\#3\\right>", 3],
    },
  },
};

// For MathJax 2
window.MathJax.TeX = window.MathJax.tex;
window.MathJax.TeX.Macros = window.MathJax.tex.macros;
