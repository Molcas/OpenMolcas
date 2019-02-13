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

var maxHeight = "10em";
var minWinHeight = 400;

// Escapes special characters and returns a valid jQuery selector
function jqSelector(str) {
  return str.replace(/([;&,\.\+\*\~':"\!\^#$%@\[\]\(\)=>\|])/g, '\\$1');
}

function expand(item) {
  var sel = getSelection().toString();
  if (!sel) { /* Do not trigger if there is selected text */
    item.css("max-height",item.data("natHeight")+20); /* Add an offset to get rid of the horizontal scrollbar) */
    item.off("click");
    item.click(function() {
      collapse($(this));
    });
  };
}

function collapse(item) {
  var sel = getSelection().toString();
  if (!sel) { /* Do not trigger if there is selected text */
    item.css("max-height",maxHeight);
    item.off("click");
    item.click(function() {
      expand($(this));
    });
  };
}

function fixednav(){
  /* Hide duplicate navigation links
     (do this before measuring the space */
  $("div.footer div[role='navigation']").each(function(i) {
    $(this).css("display","none");
  });
  $("div.footer div.source_link").each(function(i) {
    $(this).css("margin-top","0");
  });

  /* Compute available space */
  var header = 0;
  $("div.header-wrapper").each(function(i) {
    header = header + $(this).outerHeight();
  });
  var footer = 0;
  $("div.footer-wrapper").each(function(i) {
    footer = footer + $(this).outerHeight();
  });
  var available = $(window).height() - header - footer;

  /* Save the original bottom padding of the sidebar */
  $("div.sidebar").each(function(i) {
    var sb = $(this).attr("style");
    if (typeof sb == typeof undefined || sb == false) {
      sbpad = $(this).css("padding-bottom");
    }
  });

  /* If the window is big enough, change styles for fixed navigation */
  if (available > minWinHeight) {
    /* Fixed header and footer */
    $("div.header-wrapper").each(function(i) {
      $(this).css("position","fixed");
      $(this).css("top","0");
      $(this).css("width","100%");
      $(this).css("z-index","1");
    });
    $("div.footer-wrapper").each(function(i) {
      $(this).css("position","fixed");
      $(this).css("bottom","0");
      $(this).css("width","100%");
      $(this).css("z-index","1");
    });
    /* Fixed sidebar */
    $("div.sidebar").each(function(i) {
      $(this).css("position","fixed");
      $(this).css("top",header);
      $(this).css("height",available);
      $(this).css("overflow","auto");
      $(this).css("margin","0");
      $(this).css("border-bottom","none");
      /* Apply the padding as margin (see below) instead) */
      $(this).css("padding-bottom",0);
      $(this).css("z-index","1");
    });
    /* Workaround for the bottom padding, which is not properly applied */
    $("div.sidebar div[role='search']").each(function(i) {
      $(this).css("margin-bottom",sbpad);
    });
    /* Add margins to the main content */
    $("div.content-wrapper").each(function(i) {
      $(this).css("margin-top",header);
      $(this).css("margin-bottom",footer);
    });
    /* Hide duplicate navigation links */
    $("div.footer div[role='navigation']").each(function(i) {
      $(this).css("display","none");
    });
    /* Make title in header visible */
    $("div.header div.currenttitle").each(function(i) {
      $(this).css("display","block");
    });
    /* adjust scroll position */
    hash = "#" + jqSelector(window.location.hash.substring(1));
    if (hash.length > 1) {
      $("html, body").animate({
        scrollTop: $(hash).offset().top - header
      }, 500);
    };

  /* Otherwise, remove all added styles */
  } else {
    $("div.header-wrapper").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.footer-wrapper").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.sidebar").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.sidebar div[role='search']").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.content-wrapper").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.footer div[role='navigation']").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.footer div.source_link").each(function(i) {
      $(this).removeAttr("style");
    });
    $("div.header div.currenttitle").each(function(i) {
      $(this).removeAttr("style");
    });
  };
};

$(document).ready(function() {
  $("div.document div.highlight").each(function(i) {
    $(this).data("natHeight",$(this).height()); /* Store natural height */
    $(this).css("max-height",maxHeight);
    var dif = $(this).data("natHeight") - $(this).height();
    if (dif > 60) { /* Do not assing events if it is not collapsed */
      $(this).css("transition-duration",Math.min(Math.max(0.5,dif/1000),3) + "s");
      $(this).click(function() {
        expand($(this));
      });
    } else { /* Do not collapse if the difference is not large enough */
      $(this).css("max-height","none");
    };
  });
  /* cludge for clicking the same link twice */
  $("a").click(function() {
    location.hash == $(this).attr("href") ? fixednav() : null;
  });
});

window.onload = function() {
  try {
    MathJax.Hub.Queue(function() {
      fixednav();
    });
  }
  catch(err) {
    fixednav();
  }
};

window.onhashchange = function() {
  fixednav();
};

$(window).resize(function() {
  fixednav();
});

