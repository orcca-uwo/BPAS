

proc insertSimpleSVG {simpleSVGFile htmlFile} {
	
	set fp [open $htmlFile "r"];
	set contents [read $fp];
	close $fp;
	set fp [open $htmlFile "w"]; #no append

	set lines [split $contents "\n"];

	set next 0;
	foreach line $lines {
		# puts "line: $line";
		if {[regexp {.*class=\"dynheader\".*} $line] != 0} {
			#insert the simplified diagram
			set next 1;
		}
		if {$next == 1 && [regexp {Inheritance diagram for (.*):</div>$} $line -> className]} {
			puts $fp "Simplified semantic inheritance diagram for ${className}:</div>";
			puts $fp "<div class=\"dyncontent\">"
			puts $fp "<div class=\"center\"><div class=\"zoom\"><iframe scrolling=\"no\" frameborder=\"0\" src=\"${simpleSVGFile}\" width=\"100%\" height=\"600\"><p><b>This browser is not able to show SVG: try Firefox, Chrome, Safari, or Opera instead.</b></p></iframe></div>"
			puts $fp "</div>"
			puts $fp "</div>"

			#insert javascript
			puts $fp "<script type=\"text/javascript\">\n
				function toggleDiagramShow (id) {
					var rows =\$('div.dyncontent.' + id);
					var img = \$('div.dynheader.'+id+' img');
					var src = \$(img).attr('src');
					if (rows.filter(':first').is(':visible')===true) {
					    rows.css('display','none');
					    \$(img).attr('src',src.substring(0,src.length-8)+'closed.png');
					} else {
					    rows.css('display','block');
					    \$(img).attr('src',src.substring(0,src.length-10)+'open.png');
					}
				}</script>"

			set jsClassName [string map {" " "_" "&lt;" "_" "&gt;" "_" "," "_"} $className];
			puts $fp "<div class=\"dynheader ${jsClassName}\" onclick=\"javascript:toggleDiagramShow('${jsClassName}')\"><img src=\"closed.png\" alt=\"-\"/>&#160;Full inheritance diagram for ${className}:</div>";
			set next 2;
			continue;
		}

		if {$next == 2 && [regexp {^<div class="dyncontent">$} $line]} {
			puts $fp "<div class=\"inherit dyncontent ${jsClassName}\">";
			set next 0;
			continue;
		}

		puts $fp $line;
	}

	close $fp;
}


proc parseClassFiles {simpleDir insertDir} {
	set classFiles {};
	lappend classFiles {*}[glob -directory $simpleDir "class*.html"];
	lappend classFiles {*}[glob -directory $simpleDir "union*.html"];

	set noSVG {};
	set hasSVG {};
	foreach classFile $classFiles {
		set fname [file tail $classFile];
		# puts $fname
		set fp [open $classFile "r"];
		set lines [split [read $fp] "\n"];

		set found 0;
		foreach line $lines {
			if {[regexp {.*"(.*\.svg)".*} $line -> svgname]} {
				set fullSVGName [file join $simpleDir $svgname];
				set svgNameRoot [file root $svgname];
				set svgname "${svgNameRoot}-simple.svg";
				set newSVGName [file join $insertDir $svgname]
				file copy -force $fullSVGName $newSVGName;
				puts "copied $fullSVGName to $newSVGName"
				insertSimpleSVG $svgname [file join $insertDir $fname];
				set found 1;
				break;
			}
		}

		if {$found} {
			lappend hasSVG $fname;
		} else {
			lappend noSVG $fname;
		}
	}

	# puts "noSVG: "
	# foreach fname $noSVG {
	# 	puts $fname;
	# }
	# puts "hasSVG: " 
	# foreach fname $hasSVG {
	# 	puts $fname;
	# }
}

proc insertInheritsDiagram {simpleDir insertDir} {

	set simpFile [file join $simpleDir "inherits.html"];
	set insertFile [file join $insertDir "inherits.html"];
	set inheritDiagrams [glob -directory $simpleDir "inherit*.svg"];
	foreach diagram $inheritDiagrams {
		set diaRoot [file root [file tail $diagram]];
		set diaRoot "${diaRoot}-simple.svg";
		set destFile [file join $insertDir $diaRoot];
		file copy -force $diagram $destFile;
		puts "copied $diagram to $destFile"
	}

	set sfp [open $simpFile "r"];
	set ifp [open $insertFile "r"];
	set simpLines [split [read $sfp] "\n"];
	set insertLines [split [read $ifp] "\n"];

	close $sfp;
	close $ifp;

	set fp [open $insertFile "w"];

	set next 0;
	while {1} {
		set insertLines [lassign $insertLines insertLine];
		if {$next} {
			break; #throw away next line
		}
		if {[regexp {.*Go to the textual class hierarchy.*} $insertLine] != 0} {
			set next 1;
			continue;
		}

		puts $fp $insertLine;
	}

	set next 0;
	while {1} {
		set simpLines [lassign $simpLines insertLine];
		if {[regexp {.*Go to the textual class hierarchy.*} $insertLine] != 0} {
			puts $fp $insertLine;
			set next 1;
			continue;
		}
		if {$next && [regexp {(.*src=)(.*)\.svg(.*)} $insertLine -> begin svgname end]} {
			set svgname "${svgname}-simple.svg";
			puts $fp "${begin}${svgname}${end}";
			continue;
		}
		if {$next && [regexp {</table>} $insertLine] != 0} {
			puts $fp $insertLine;
			break;
		}
		
		if {$next} {
			puts $fp $insertLine
		}
	}


	puts $fp "<script type=\"text/javascript\">\n
				function toggleHierarchyShow (id) {
					var rows =\$('table.toggletable.' + id);
					var img = \$('div.dynheader.'+id+' img');
					var src = \$(img).attr('src');
					if (rows.filter(':first').is(':visible')===true) {
					    rows.css('display','none');
					    \$(img).attr('src',src.substring(0,src.length-8)+'closed.png');
					} else {
					    rows.css('display','block');
					    \$(img).attr('src',src.substring(0,src.length-10)+'open.png');
					}
				}</script>"
	puts $fp "<div class=\"dynheader full_class_hierarchy\" onclick=\"javascript:toggleHierarchyShow('full_class_hierarchy')\"><img src=\"closed.png\" alt=\"-\"/>&#160;<b>Full Class Hierarchy:</b></div>"
	puts $fp "<table class=\"toggletable full_class_hierarchy\" border=\"0\" cellspacing=\"10\" cellpadding=\"0\" style=\"display: none;\">";

	set next 1;
	foreach line $insertLines {
		if {$next && [regexp {</table>} $line] != 0} {
			puts $fp "</table></div>";
			set next 0;
			continue;
		}
		puts $fp $line;
	}

	close $fp;

}

proc fixHierarchyList {simpleDir insertDir} {
	set simpFile [file join $simpleDir "hierarchy.html"];
	set insertFile [file join $insertDir "hierarchy.html"];
	file copy -force $simpFile $insertFile;
}


if {$::argc > 1} {

	lassign $::argv simpleDir insertDir;
	puts "simpleDir: $simpleDir"
	puts "insertDir: $insertDir";
	parseClassFiles $simpleDir $insertDir
	insertInheritsDiagram $simpleDir $insertDir 
	fixHierarchyList $simpleDir $insertDir;
	# insertSimpleSVG $newSVG $htmlFile;
} else {
	puts "USAGE: tclsh insertSimpleDiagram.tcl INSERT_SVG_FILE HTML_FILE"
}