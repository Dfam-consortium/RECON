#!/usr/local/bin/perl
use strict;
##
## Prototype script to explore RECON results for a given family
##   ./graphFamilyElementRelationships.pl <directory_containing_RECON_results>
##                                        <family_id_to_visualize>
##
## e.g.
##   ./graphFamilyElementRelationships.pl . 3 > recon_family_3.html
##
##
my $reconDir  = $ARGV[ 0 ];
my $family_ID = $ARGV[ 1 ];

$reconDir = "." if ( ! -d $reconDir );

# Double check that we can see the summary, ele_redef_res and 
# edge_redef_res subdirectories
die "Could not locate RECON summary directory: $reconDir/summary\n"
    if ( ! ( -d "$reconDir/summary" && -d "$reconDir/ele_redef_res" &&
             -d "$reconDir/edge_redef_res" ) );

#
# Parse summary/eles file to obtain family to element mapping
#
open IN, "<$reconDir/summary/eles" or die;
my %elements      = ();
my %eleInFamilies = ();
while ( <IN> )
{
  next if ( /^#/ );
  s/^\s+//;
  my @flds = split( /\s+/ );
  $eleInFamilies{ $flds[ 1 ] } = $flds[ 0 ];

  #  fam    ele   dir  sequence    start     end
  #       1      1  1       gi|1     1519     1709
  #       1  11727  1     gi|596    29279    29564
  #       3  11454 -1     gi|588     7468     7549
  if ( /^\s*($family_ID\s+.*)/ )
  {
    $elements{ $flds[ 1 ] } = {};
  }
}
close IN;


#
# parse all elements files for the query family
#
my $maxImagesInASingleElement = 0;
my $maxElementGraphLen        = 0;
foreach my $element ( keys( %elements ) )
{
  # Which is the most authoritative?
  #open IN, "<$reconDir/edge_redef_res/e$element"
  #open IN, "<$reconDir/ele_def_res/e$element"
  open IN, "<$reconDir/ele_redef_res/e$element"
      or die "Could not open e$element\n";
  my $eleSeq;
  my $eleStart;
  my $eleEnd;
  my $img_no;
  my $eleFileData = "";
  while ( <IN> )
  {
    # Example from ele_redef_res
    #
    #    frag gi|10 184 317
    #    img_no 3
    #    edge_no 3
    #    msp 13550 p 254 79.0 -1 17 gi|10 184 264 11454 gi|588 7472 7554
    #    msp 13552 s 260 81.7 -1 17 gi|10 191 250 11607 gi|592 12260 12320
    #    msp 13554 p 287 86.7 1 17 gi|10 220 317 3838 gi|256 1587 1688
    #    edge 24 p 1 86 17 3838
    #    edge 25 P -1 79 17 11454
    #    edge 26 s 0 0 17 11607
    #
    # Save the entire file in the record
    $elements{$element}->{'file_data'} .= $_;
    if ( /^frag\s+(\S+)\s+(\d+)\s+(\d+)/ )
    {
      $eleSeq                            = $1;
      $eleStart                          = $2;
      $eleEnd                            = $3;
      $elements{$element}->{'ele_seq'}   = $eleSeq;
      $elements{$element}->{'ele_start'} = $eleStart;
      $elements{$element}->{'ele_end'}   = $eleEnd;
    } elsif ( /^img_no\s+(\d+)/ )
    {
      if ( $1 > $maxImagesInASingleElement )
      {
        $maxImagesInASingleElement = $1;
      }
    } elsif ( /^msp/ )
    {
      # The element numbers are in the msp line in col 7/11
      my @flds   = split();
      my $msp_id = $flds[ 1 ];
      if ( $msp_id & 1 )
      {
        # Odd ids decrement to even so that we can match MSPs
        # ....msps always come in pairs
        $msp_id--;
      }
      $elements{$element}->{$msp_id}->{'line'}       = $_;
      $elements{$element}->{$msp_id}->{'edge_type'}  = $flds[ 2 ];
      $elements{$element}->{$msp_id}->{'msp_orient'} = $flds[ 5 ];
      my $otherElementID = -1;
      if ( $flds[ 6 ] == $element )
      {
        # First coordinate is covering this element
        $elements{$element}->{$msp_id}->{'msp_start'} = $flds[ 8 ];
        $elements{$element}->{$msp_id}->{'msp_end'}   = $flds[ 9 ];
        $otherElementID                               = $flds[ 10 ];
      } else
      {
        # Second coordinate is covering this element
        $elements{$element}->{$msp_id}->{'msp_start'} = $flds[ 12 ];
        $elements{$element}->{$msp_id}->{'msp_end'}   = $flds[ 13 ];
        $otherElementID                               = $flds[ 6 ];
      }
      if ( !defined $elements{$otherElementID} )
      {
        $elements{$element}->{$msp_id}->{'paired'} = 0;
      } else
      {
        $elements{$element}->{$msp_id}->{'paired'} = 1;
      }
      if ( defined $eleInFamilies{$otherElementID} )
      {
        $elements{$element}->{$msp_id}->{'family_link'} =
            $eleInFamilies{$otherElementID};
      } else
      {
        $elements{$element}->{$msp_id}->{'family_link'} = "";
      }

      if ( !defined $elements{$element}->{'min_msp_start'}
           || $elements{$element}->{'min_msp_start'} >
           $elements{$element}->{$msp_id}->{'msp_start'} )
      {
        $elements{$element}->{'min_msp_start'} =
            $elements{$element}->{$msp_id}->{'msp_start'};
      }
      if ( !defined $elements{$element}->{'max_msp_end'}
           || $elements{$element}->{'max_msp_end'} <
           $elements{$element}->{$msp_id}->{'msp_end'} )
      {
        $elements{$element}->{'max_msp_end'} =
            $elements{$element}->{$msp_id}->{'msp_end'};
      }
      my $ele_len =
          $elements{$element}->{'max_msp_end'} -
          $elements{$element}->{'min_msp_start'} + 1;
      if ( $ele_len > $maxElementGraphLen )
      {
        $maxElementGraphLen = $ele_len;
      }
    }
  }
  close IN;
}


##
## Prepare the HTML output
##
print "<html>\n";
print "<h2>RECON Family $family_ID</h2>\n";
# Head stanza
print <<EOH;
  <head>
    <title>RECON Family Element Relationships</title>
    <script>
            var prevClass = "";
            var prevColor = "yellow";
     
            function highlightPair(event) {
              if (prevClass != "") {
                var rects = document.getElementsByClassName(prevClass);
                for (i = 0; i < rects.length; ++i) {
                  var element = rects[i];
                  element.setAttribute("stroke", prevColor);
                }
              }
              var myclass = event.target.getAttribute("class");
              var rects = document.getElementsByClassName(myclass);
              for (i = 0; i < rects.length; ++i) {
                var element = rects[i];
                prevColor = element.getAttribute("stroke");
                element.setAttribute("stroke", "blue");
                var aEles = element.getElementsByTagName("animate");
                if (aEles.length > 0) {
                  aEles[0].beginElement();
                }
              }
              prevClass = myclass;
            }
    </script>
  </head>
EOH

# visual constants
my $graphWidth = 1600;
my $annotHeight = 4;
my $annotSpacer = 3;
my $fontSize = "1.0em";

# Dense
#my $graphWidth = 200;
#my $annotHeight = 1;
#my $annotSpacer = 2;
#my $fontSize = "0.5em";

# So that all graphs represent the same scale
my $graphPixelsPerBP = $graphWidth / $maxElementGraphLen;
foreach my $element ( keys( %elements ) )
{
  my $ele = $elements{$element};
  my $rectID      = 1;
  my $curY        = 5;
  my $maxWidth    = $ele->{'max_msp_end'} - $ele->{'min_msp_start'} + 1;
  my $graphWidth  = $maxWidth * $graphPixelsPerBP;
  my $graphHeight = ( $annotHeight + $annotSpacer ) * scalar( keys( %{$ele} ) );
  my $eleSeq      = $ele->{'ele_seq'};
  my $eleStart    = sprintf( "%0d",
                          ( $ele->{'ele_start'} - $ele->{'min_msp_start'} ) *
                              $graphPixelsPerBP );
  my $eleEnd = sprintf( "%0d",
                        ( $ele->{'ele_end'} - $ele->{'min_msp_start'} ) *
                            $graphPixelsPerBP );
  my $svgHeight  = $graphHeight + 10;
  my $textHeight = $svgHeight - 2;

  print "<svg width=\"$graphWidth\" height=\"$svgHeight\" style=\"border: 1px solid black;\">\n";
  print "  <g>\n";
  # Draw the "frag" boundaries for this element as a dashed vertical grey lines
  print "  <line x1=\"$eleStart\" y1=\"0\" x2=\"$eleStart\" y2=\"$graphHeight\" stroke=\"grey\" stroke-width=\"0.5\" stroke-dasharray=\"5\"/>\n";
  print "  <line x1=\"$eleEnd\" y1=\"0\" x2=\"$eleEnd\" y2=\"$graphHeight\" stroke=\"grey\" stroke-width=\"0.5\" stroke-dasharray=\"5\"/>\n";
  # Label the plot
  print "  <text font-size=\"$fontSize\" x=\"0\" y=\"$textHeight\">ele-$element $eleSeq:$ele->{'ele_start'}-$ele->{'ele_end'} [" . ($ele->{'ele_end'}-$ele->{'ele_start'}+1) . "]</text>\n";

  # Sort the alignments (aka "images") on start positions 
  foreach my $mspID (
    sort {
      $ele->{$a}->{'msp_start'} <=> $ele->{$b}->{'msp_start'}
    } grep
    {
      /^\d+$/
    } keys( %{$ele} )
      )
  {
    # Not sure why this was necessary.  Are there e# files with
    # msps without a numerical id?
    next if ( $mspID !~ /^\d+$/ );

    my $msp        = $ele->{$mspID};
    my $relStart   = $msp->{'msp_start'} - $ele->{'min_msp_start'};
    my $relEnd     = $msp->{'msp_end'} - $ele->{'min_msp_start'};
    my $graphStart = sprintf( "%0d", $relStart * $graphPixelsPerBP );
    my $graphEnd   = sprintf( "%0d", $relEnd * $graphPixelsPerBP );
    my $graphLen   = $graphEnd - $graphStart + 1;

    my $additionalAttr = "";
    # secondary edge
    if ( $msp->{'edge_type'} eq 's' )
    {
      $additionalAttr .= " stroke-opacity=\"0.5\"";
    }
    if ( $msp->{'family_link'} eq "" )
    {
      $additionalAttr .= " stroke-dasharray=\"10\"";
    }

    if ( $msp->{'paired'} )
    {
      # Draw inter-family alignment
      print
"    <line id=\"rect-$rectID\" class=\"msp-$mspID\" x1=\"$graphStart\" y1=\"$curY\" x2=\"$graphEnd\" y2=\"$curY\" stroke=\"green\" stroke-width=\"$annotHeight\" onClick=\"highlightPair(event)\" $additionalAttr>\n";
    } else
    {
      # Draw intra-family alignment
      print
"    <line id=\"rect-$rectID\" class=\"msp-$mspID\" x1=\"$graphStart\" y1=\"$curY\" x2=\"$graphEnd\" y2=\"$curY\" stroke=\"red\" stroke-width=\"$annotHeight\" onClick=\"highlightPair(event)\" $additionalAttr>\n";
    }
    # Add hover text indicating the msp details an the family it links to
    print "      <title>$msp->{'line'}: family-$msp->{'family_link'}</title>\n";
    print
"      <animate attributeType=\"CSS\" attributeName=\"opacity\" values=\"0.75;0.5;0.25;0;0.25;0.5;0.75;1\" begin=\"indefinite\" dur=\"0.15s\" repeatCount=\"3\" /></line>\n";
    $rectID++;
    $curY += $annotHeight + $annotSpacer;
  }
  print " </g>\n";
  print "</svg>\n";
}

## Fiddle for playing around with interactivity
## https://jsfiddle.net/9830qscr/13/

