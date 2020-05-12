#!/bin/bash

# NAME:         extract_pdf_results
# VERSION:      0.1
# AUTHOR:       (c) 2014 Glutanimate
# DESCRIPTION:  Extracts PDF pages that contain supplied string and concatenates them to a new file.
# FEATURES:     
# DEPENDENCIES: pdfgrep pdftk
#               ➥install on Ubuntu/Debian with sudo apt-get install pdfgrep pdftk
#
# LICENSE:      GNU GPLv3 (http://www.gnu.de/documents/gpl-3.0.en.html)
#
# NOTICE:       THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
#               EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES 
#               PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR 
#               IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
#               AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND 
#               PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE,
#               YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#
#               IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY 
#               COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS 
#               PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, 
#               INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE 
#               THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED 
#               INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE 
#               PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER 
#               PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
#
# USAGE:        extract_pdf_results <string> <pdffile>

STRING="$1"
FILE="$2"
FILENAME="${FILE##*/})"
BASENAME="${FILENAME%.*}"
DIRNAME="${FILE%/*}"

echo "Processing $FILE..."

## find pages that contain string, remove duplicates, convert newlines to spaces

echo "Looking for $STRING..."

PAGES="$(pdfgrep -n "$STRING" "$FILE" | cut -f1 -d ":" | uniq | tr '\n' ' ')"

echo "Matching pages:
$PAGES"

## extract pages to new file in original directory

echo "Extracting result pages..."

pdftk "$FILE" cat $PAGES output "${DIRNAME}/${BASENAME}_pages_with_${STRING}.pdf"

echo "Done."

