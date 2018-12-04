#!/usr/bin/perl -wc

# TEST!!!
# Package for molcas style

# (C)reated by Galex.
# Last-modified: 10-10-2003
# (C)left ValeSoft 


=for REMOVED
# internal package for utility functions
package local::Molcas_style::utils;

sub Replace # modifies '$_'!
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;

 my ($left,$right,$from,$to)=@_;
 my ($from_len,$to_len)=(length($from),length($to));
 my $was_replaced=0;
 
 pos()=0;
 while (/\Q$left\E/g) {
  while (m/(\Q$from\E|\Q$right\E)/g) {
   my $p=pos();
   # warn "p=$p";
   if ($1 eq $from) {
    # warn "Found!";
    pos()=$p-$from_len;
    s/\G\Q$from\E/$to/;
    pos()=$p-$from_len+$to_len;
    $was_replaced=1;
   } else {
    last;
   }
  }
 }

 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return $was_replaced;
}
=cut

package hidden; # small package for testing purposes

$Debug=0;

package main; # the main part of the module

warn "\n\nMolcas.perl Loading, pwd is ".`pwd`."\nNow" if $hidden::Debug;


# Hash describing \molcas-like commands (w/o args)
my %molcas_cmds=qw
(
 molcas    MOLCAS
 molcasi   MOLCAS-1
 molcasii  MOLCAS-2
 molcasiii MOLCAS-3
 molcasiv  MOLCAS-4
 molcasv   MOLCAS-5
 molcasvi  MOLCAS-6
 molcasvii MOLCAS-7
 molcasviii MOLCAS-8
 MolcasWWW http://www.molcas.org/
);

# sometimes in future, they will be read from outside:
my %molcas_const=qw
(
 molcasversion 8.1
 molcasversionnodots 81
 MolcasRoot /usr/local/lib/molcas81
);

=for COMMENT (galexv)
 Below is ugly eval()ed-code. Yes, I know, it can be done
 in a little bit cleaner way with closure-generators and
 modifications of the perl symbol table, like this:
 
 *{"do_cmd_".$cmd_name}=generate_noarg_sub($constant,$opentag,$closetag);
  and
 *{"do_cmd_".$cmd_name}=generate_onearg_sub($opentag,$closetag);
 
 but do we really need it?
 Besides, it can lead us to a mess with global, local- and my-variables,
 which probably doesn't cost slight decrease in compilation time...
  
=cut

# Generate 'do_cmd_molcas'-like subs to handle \molcas-like cmds
foreach (keys %molcas_cmds)
{
 my $code=<<"~~END_OF_CODE~~";
sub do_cmd_$_
{ 
 return "<I>$molcas_cmds{$_}</I>".\$_[0];
}
~~END_OF_CODE~~

 print "Generated code:\n$code\n" if ($hidden::Debug);
 eval($code);
}

# Generate 0-argument subs to handle constant-like cmds
foreach (keys %molcas_const)
{
 my $code=<<"~~END_OF_CODE~~";
sub do_cmd_$_
{ 
 return "$molcas_const{$_}".\$_[0];
}
~~END_OF_CODE~~

 print "Generated code:\n$code\n" if ($hidden::Debug);
 eval($code);
}

# command depending on other ones
sub do_cmd_molcasthis
{
 my ($macro_body,$value)=($new_command{'molcasversion'},'???');
 if (defined($macro_body)) {
  # \molcasversion definition was seen by latex2html
  my ($nargs,$opt);
  ($nargs,$value,$opt)=split /:!:/,$body; # parse definition
  
  if ($nargs!=0 || $opt ne '}') {
   warn "\\molcasthis: wrong \\molcasversion definition --- taking as '$value'\n";
  }
 } elsif (defined(&do_cmd_molcasversion)) { # is it defined as do_cmd_...?
  # use \molcasversion from style (.perl) file
  $value=do_cmd_molcasversion();
 } else {
  warn "\\molcasthis: \\molcasversion not defined --- taking as '$value'\n";
 }  
 return "<I>MOLCAS-$value</I>".$_[0];
}

# Command that makes visible space
sub do_cmd_visiblespace
{
#VV return '&nbsp;<img src="visible_space.png" alt=" ">&nbsp;'.$_[0];
 return '&nbsp;&nbsp;'.$_[0];
}

  

# one-argument command tag definitions
my %one_arg_cmd=
(
# 'keyword'=>['<font size=+1>','</font>'],
 'keyword'=>['<b>','</b>'],
 'program'=>['<tt><b>','</b></tt>'],
 'command'=>['<b><i>','</i></b>'],
 'file'=>['<tt><i><font size="+1">','</font></i></tt>'],
 'ftncode'=>['<tt>','</tt>']
);
 
# generate \keyword{..}-like command handlers
foreach (keys %one_arg_cmd)
{
 my ($otag,$ctag)=@{$one_arg_cmd{$_}};
 my $code=<<"~~END_OF_CODE~~";
sub do_cmd_$_
{
 local \$_=shift;
 my \$arg=(s/\$next_pair_pr_rx// || s/\$next_pair_rx//)?\$2:&missing_braces;
 return '$otag'.\$arg.'$ctag'.\$_;
}
~~END_OF_CODE~~
 print "Generated code:\n$code\n" if ($hidden::Debug);
 eval($code);
}


# \namelist has to be coded separately due to \obeyspaces:

sub do_cmd_namelist
{
 local $_=shift;
 my $arg=(s/$next_pair_pr_rx// || s/$next_pair_rx//)?$2:&missing_braces;
 $arg=~s/ /\&nbsp;/g;     # \obeyspaces
 return '<table border=1><tr><td>'.$arg.'</table>'.$_;
}

# Sub to copy 'visible space' to current (document) directory.
sub copy_visible_space_img
{
 my $img_f=$MOLCAS_STYLE->{VISIBLE_SPACE_IMG};
 if (!defined($img_f)) {
  warn "No image for \"visible space\" given in configuration,\n".
       "so none is copied.\n.";
  return;
 }
 L2hos->Copy($img_f,".${dd}visible_space.png");
}

# and call it as part of initialization process
copy_visible_space_img();


=for COMMENT:
Preprocessor hook, to process command symbols
in our verbatim-like environments specially.
=cut

=for  REMOVED
sub preprocess_alltt # modifies '$_' !
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 
 # process all command symbols in verbatim-like environments
 for my $env (qw(inputlisting sourcelisting)) {
  # quote math symbols
  local::Molcas_style::utils::Replace("\\begin{$env}",
                                      "\\end{$env}",
	  			      '$','\\$');

  # remove all unquoted braces (FIXME: quick ugly hack!!!)
  for my $b (qw({ })) {
  # replace quoted braces by special symbol:
   local::Molcas_style::utils::Replace("\\begin{$env}",
                                       "\\end{$env}",
	  			       "\\$b","\x01");
   # replace braces by empty string:
   local::Molcas_style::utils::Replace("\\begin{$env}",
                                       "\\end{$env}",
	  			       $b,'');
   # replace special symbol back by quoted braces:
   local::Molcas_style::utils::Replace("\\begin{$env}",
                                       "\\end{$env}",
	  			       "\x01","\\$b",);
				       
  }

  # remove trailing spaces
  while (
   local::Molcas_style::utils::Replace("\\begin{$env}",
                                       "\\end{$env}",
	  			       " \n","\n")
  ) {
   warn "Trailing spaces removed\n" if $hidden::Debug;
  };				       

  # remove empty lines
  while (
   local::Molcas_style::utils::Replace("\\begin{$env}",
                                       "\\end{$env}",
	  			       "\n\n","\n")
  ) {
   warn "Empty lines removed\n" if $hidden::Debug;
  };				       
 }

 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
}
=cut

# Another variant of preprocessing.

# What environments should be preprocessed?
my %PreprocessEnv=(
 inputlisting  => \&quote_specials,
 sourcelisting => \&quote_specials,
);
 
sub preprocess_alltt # modifies '$_' !
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 
 # for environments which need preprocessing...
 for $env (keys %PreprocessEnv) {
  my $sub=$PreprocessEnv{$env};
  s/(^[^%\n]*\\begin\{$env\}\n?)(.*?)(^[^%\n]*\\end\{$env\})/$1.$sub->($2).$3/msge; # nesting won't work!!
 }
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
}
 
sub quote_specials
{
 my $in=shift;
 $in=~s/[%{}\$\\\n]/'<<molcas_sty_quoted>><<'.ord($&).'>>'/ge;
 return $in;
}

sub restore_specials # modifies argument!
{
 $_[0]=~s/<<molcas_sty_quoted>><<(\d+)>>/chr($1)/ge;
}

# Helper for verbatim-like environments
sub verbatim_like_helper # ($input_text --> $processed_text)
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;

 # The following code is taken (with modifications)
 # from 'verbatim.perl' style file
 local($closures,$reopens) = &preserve_open_tags; # ensure proper tag nesting
 $verbatim{++$global{'verbatim_counter'}} = $_[0]; # save content verbatim
 
 my ($tag_open,$tag_close) = ('<TT>','</TT>');
# if ($USING_STYLES) {
#     $env_id .= ' CLASS="verbatim"' unless ($env_id =~ /(^|\s)CLASS\s*\=/i);
#     $verb_pre =~ s/>/ $env_id>/;
# }
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return $closures."<BR>\n". # close all tags
        $tag_open. # open our tag
	$verbatim_mark.'verbatim'.$global{'verbatim_counter'}.'#'. # set verbatim mark
        $tag_close. # close our tag
	$reopens; # reopen all tags, closed just above
 
}


# process argument as needed by {inputlisting} and {sourcelisting}
sub process_molcas_verb_environment # modifies argument!
{
 # expand \symbol to quoted-specials
 $_[0]=~s/\\symbol\{(\d+)\}/<<molcas_sty_quoted>><<$1>>/g;
 $_[0]=~s/\\symbol\{('[0-7])\}/'<<molcas_sty_quoted>><<'.oct($1).'>>'/ge;
 $_[0]=~s/\\symbol\{("[0-9A-Fa-f])\}/'<<molcas_sty_quoted>><<'.hex($1).'>>'/ge;
 
# # treat some commands specially:
# s/\\hspace\*?\{\d*cm\}/\\quad\\quad/g;

 # remove unescaped braces and percents 
 # FIXME: percented lines should be removed completely!
 $_[0]=~s/(^|[^\\])[{}%]/$1/g; 

 # unescape braces and percents
 $_[0]=~s/\\([{}%])/$1/g;

 # expand simple 0-argument commands (will break on other ones)
 $_[0]=~s/\\(\w+)/defined(&{"do_cmd_$1"}) && &{"do_cmd_$1"}('')/eg;

 # expand quoted-specials:
 restore_specials($_[0]); 
}
 

sub do_env_inputlisting
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 my $visible_space_tag=do_cmd_visiblespace();
 
 restore_specials($in);
 process_molcas_verb_environment($in);
 
 $in=~s/ +(\r?\n)/$1/g; # remove trailing spaces
 $in=~s/(\r?\n)/<BR>$1/g;   # \obeylines
 $in=~s/\ {2,}/\ /g; # reduce spaces
 $in=~s/ /$visible_space_tag/g; # show spaces
 
 warn("\nSubroutine ",(caller(0))[3]," is about to exit") if $hidden::Debug;
 return verbatim_like_helper($in); # process as verbatim-like
}

sub do_env_sourcelisting
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 
 restore_specials($in);
 process_molcas_verb_environment($in);

 $in=~s/(\r?\n)/<BR>$1/g;   # \obeylines
 
 warn("\nSubroutine ",(caller(0))[3]," is about to exit") if $hidden::Debug;
 return verbatim_like_helper($in); # process as verbatim-like
}


# Helper for list-like environments,
# derived from 'do_env_description'

sub molcas_list_helper 
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 local($_, $bullet) = @_;
 &protect_useritems($_);
 $_ = &translate_environments($_); # attempt to fix

 s/\n?$item_description_rx\s*($labels_rx8)?\s*/"\n<TR>\n<TD>". 
     (($9)? "<A NAME=\"$9\">$1<\/A>" : $1 ) ."<TD>"/egm;
 s/\n?\\item\b\s*([^$letters\\]|)\s*/\n<TR>\n<TD>$bullet<TD>$1/gm;
 
 s/^\s+//m;
 s/\n$//s;
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return $_;
}


my $table_attr='BORDER=1 CELLSPACING=0 CELLPADDING=4';

# List-like environments

sub do_env_programlist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 $in=molcas_list_helper($in,'<i>program</i>');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD width="20%"><b>Program</b><TD><b>Purpose</b>
$in
</TABLE>
~~END_OF_TABLE~~
}

sub do_env_filelist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 $in=molcas_list_helper($in,'<i>file</i>');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD width="20%"><b>File</b><TD><b>Contents</b>
$in
</TABLE>
~~END_OF_TABLE~~
}

sub do_env_commandlist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 $in=molcas_list_helper($in,'<i>command</i>');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD width="20%"><b>Command</b><TD><b>Purpose</b>
$in
</TABLE>
~~END_OF_TABLE~~
}

sub do_env_variablelist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 $in=molcas_list_helper($in,'<i>&nbsp;</i>');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD width="20%"><b>Variable</b><TD><b>Purpose</b>
$in
</TABLE>
~~END_OF_TABLE~~
}

sub do_env_keywordlist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 my $in=shift;
 $in=molcas_list_helper($in,'&nbsp;');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD width="20%"><b>Keyword</b><TD><b>Meaning</b>
$in
</TABLE>
~~END_OF_TABLE~~
}

# Environment with 1 arg

sub do_env_syntaxlist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 local $_=shift;
 my $arg=(s/$next_pair_pr_rx// || s/$next_pair_rx//)?$2:&missing_braces;
 
 $_=molcas_list_helper($_,'---');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD><b>Syntax:</b><TD>$arg
$_
</TABLE>
~~END_OF_TABLE~~
}
 
# Environment with 3 args

sub do_env_whatnotlist
{
 warn("\nSubroutine ",(caller(0))[3]," called") if $hidden::Debug;
 local $_=shift;
 my $arg1=(s/$next_pair_pr_rx// || s/$next_pair_rx//)?$2:&missing_braces;
 my $arg2=(s/$next_pair_pr_rx// || s/$next_pair_rx//)?$2:&missing_braces;
 my $arg3=(s/$next_pair_pr_rx// || s/$next_pair_rx//)?$2:&missing_braces;
 
 $_=molcas_list_helper($_,'&nbsp;');
 warn("\nSubroutine ",(caller(0))[3]," exited") if $hidden::Debug;
 return <<"~~END_OF_TABLE~~";
<TABLE $table_attr> 
<TR><TD>$arg2<TD>$arg3
$_
</TABLE>
~~END_OF_TABLE~~
}


# Declarations (\keyfont-like)

# Tag definitions:
my %font_decl=
(
 'inputfont'=>['<TT>','</TT>'],
 'sourcefont'=>['<TT>','</TT>'],
 'keyfont'=>['<font size=+1>','</font>'],
 'cmndfont'=>['<TT>','</TT>'],
 'prgmfont'=>['<TT>','</TT>'],
 'filefont'=>['<TT>','</TT>'],
 'varfont'=>['<TT>','</TT>'],
 'ftnfont'=>['<TT>','</TT>']
);
 
# generate \keyfont-like command handlers
foreach (keys %font_decl)
{
 my ($otag,$ctag)=@{$font_decl{$_}};
 my $code=<<"~~END_OF_CODE~~";
sub do_cmd_$_
{
 return "$otag".\$_[0]."$ctag";
}
~~END_OF_CODE~~
 print "Generated code:\n$code\n" if ($hidden::Debug);
 eval($code);
}
 

=for REMOVED
&process_commands_wrap_deferred(<<"~~IN_TEX~~");
inputlisting
~~IN_TEX~~
=cut


=for COMMENT
&ignore_commands(<<"==IGNORED==");
molcas
==IGNORED==
=cut

&process_commands_in_tex(<<"~~IN_TEX~~");
logo # []
flowchart# {}
myincludegraphics# {}
~~IN_TEX~~

1;
