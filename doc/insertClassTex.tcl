
proc glob-r {{dir .} args} {
    set res {}
    foreach i [lsort [glob -nocomplain -dir $dir *]] {
        if {[file isdirectory $i]} {
            eval [list lappend res] [eval [linsert $args 0 glob-r $i]]
        } else {
            if {[llength $args]} {
                foreach arg $args {
                    if {[string match $arg $i]} {
                        lappend res $i
                        break
                    }
                }
            } else {
                lappend res $i
            }
        }
    }
    return $res
}

proc insertTexData {texFile texGeneratedPDF htmlFile} {
	set fp [open $htmlFile "r"];
	set contents [read $fp];
	close $fp;
	set fp [open $htmlFile "w"]; #no append

	set lines [split $contents "\n"];

	set texGeneratedPDF [file tail $texGeneratedPDF];

	set next 0;
	foreach line $lines {
		# puts "line: $line";
		if {$next == 0 && [regexp {^<div class="dynheader">$} $line] != 0} {
			puts $fp "
			<script type=\"text/x-mathjax-config\">
			MathJax.Hub.Config({
			  tex2jax: {
			    inlineMath: \[\['$','$'\]\],
			    processEscapes: true
			  }
			});
			</script>"
			puts $fp "<script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>"
			#now insert tex
			puts $fp "<div class=\"insertedTex\">";
			
			set texFP [open $texFile "r"];
			set texContents [read $texFP];
			set texContents [split $texContents "\n"];
			# set texContents [string map {\\ \\\\} $texContents];
			close $texFP;

			foreach texLine $texContents {
				if {[regexp {^(.*)\\emph\{(.*)\}(.*)$} $texLine -> first emph second] != 0} {
					set texLine "${first} <i>${emph}</i> ${second}";
					# puts $fp "${first} <i>${emph}</i> ${second}";
				}
				if {[regexp {^(.*)\\url\{(.*)\}(.*)$} $texLine -> first url second] != 0} {
					set texLine "${first} <a href=\"${url}\">${url}</a> ${second}";
					# puts $fp "${first} <a href=\"${url}\">${url}</a> ${second}";
				}
				if {[regexp {^(.*)\\small\{(.*)\}(.*)$} $texLine -> first small second] != 0} {
					set texLine "${first} <small>${small}</small> ${second}";
				}
				if {[regexp {^(.*)\\(.*)\{(.*)\}(.*)$} $texLine -> first command brackets second] != 0} {
					set texLine "${first} \$\\${command}{${brackets}}\$ ${second}";
					# puts $fp "${first} \$\\${command}{${brackets}}\$ ${second}";
				}
				if {[regexp {^[^\\]*\{(.*)\}(.*)$} $texLine -> first brackets second] != 0} {
					set texLine "${first} \${${brackets}}\$ ${second}";
					# puts $fp "${first} \${${brackets}}\$ ${second}";
				}
		
				puts $fp $texLine;
			}

			# puts $fp $texContents;
			# puts $fp "Blah blah tex";
			puts $fp "<a href=\"${texGeneratedPDF}\"> More... (pdf)</a>"
			puts $fp "</div><!-- tex prompt -->";			

			set next 1;
		}
		
		puts $fp $line;
	}

	close $fp;
}


proc parseClassFiles {texDir insertDir} {
	set texFiles {};
	lappend texFiles {*}[glob-r $texDir "*.tex"];
	
	set partTex {};
	foreach tex $texFiles {
		if {[regexp {.*-doc\.tex} $tex] != 0} {
			continue;
		}
		if {[regexp {.*-rest\.tex} $tex] != 0} {
			continue;
		}

		lappend partTex $tex;
	}

	set pdfFiles {};
	lappend pdfFiles {*}[glob-r $texDir "*.pdf"];

	# foreach tex $partTex {
	# 	puts "part tex: $tex";
	# }
	# foreach pdf $pdfFiles {
	# 	puts "pdf: $pdf";
	# }

	set texMap {};
	foreach texF $partTex {
		set texTail [file root [file tail $texF]];
		foreach pdf $pdfFiles {
			if {[string match -nocase "*/${texTail}-doc.pdf" $pdf]} {
				dict set texMap $texF $pdf;
			}
		}
	}

	set htmlFiles [glob -nocomplain -directory $insertDir "*.html"];
	set htmlFile "";
	foreach texF $partTex {
		set texTail [file root [file tail $texF]];
		if {![dict exists $texMap $texF]} {
			puts "Warning: could not handle tex file: $texF";
			continue;
		}

		foreach htmlF $htmlFiles {
			if {[string match -nocase "*/class${texTail}.html" $htmlF] ||
					[string match -nocase "*/union${texTail}.html" $htmlF]} {
				# puts "matched $texF and $htmlF";
				set pdfF [dict get $texMap $texF];
				set htmlDir [file dirname $htmlF];
				set newPDF [file join $htmlDir [file tail $pdfF]];
				file copy -force -- $pdfF $newPDF;
				# puts "copying from $pdfF to $newPDF";
				insertTexData $texF $newPDF $htmlF;
				break;
			}
		}
	}


}



if {$::argc > 1} {

	lassign $::argv texDir insertDir;
	puts "texDir: $texDir"
	puts "insertDir: $insertDir";
	parseClassFiles $texDir $insertDir;
	# insertInheritsDiagram $texDir $insertDir 
	# fixHierarchyList $texDir $insertDir;
	# insertSimpleSVG $newSVG $htmlFile;
} else {
	puts "USAGE: tclsh insertSimpleDiagram.tcl INSERT_SVG_FILE HTML_FILE"
}