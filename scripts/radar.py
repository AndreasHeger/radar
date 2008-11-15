####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: WrapperRadar.py,v 1.1.1.1 2002/07/02 10:46:58 heger Exp $
##
##
####
####
USAGE = """radar.py [OPTIONS] filenames

run radar on one or more fasta-formatted files.
"""

import re, sys, os, string, tempfile, optparse, shutil, types, subprocess

rx_repeat = re.compile("^\s*(\d+)\-\s*(\d+)\s+\(\s*(\S+)\/\s*(\S+)\)\s+(\S+)")

##----------------------------------------------------------------------------------------
DEBUG = 0

import pyradar

class RadarRepeat:
    def __init__( self, file):
        """Build members from a file."""

        ## read up to header
        while 1:
            line = file.readline()
            if not line: break
            
            if line[0:3] == "No.": break

        ## read header
        line = file.readline()
        data = string.split( line, "|")

        self.mNRepeats = string.atoi(string.strip(data[0]))
        self.mScore    = string.atof(string.strip(data[1]))
        self.mLength   = string.atoi(string.strip(data[2]))
        self.mDiagonal = string.atoi(string.strip(data[3]))
        self.mBW_from  = string.atoi(string.strip(data[4]))
        self.mBW_to    = string.atoi(string.strip(data[5]))
        self.mLevel    = string.atoi(string.strip(data[6]))

        ## read single repeats
        file.readline()

        self.mRepeatUnits = []
        while 1:
            line = file.readline()
            if not line: break
            
            result = rx_repeat.search( line )

            if not result: break
            
            self.mRepeatUnits.append( result.groups())
        
    def GetBWFrom(self):
        return self.mBW_from
    def GetBWTo(self):
        return self.mBW_to
    def GetNRepeats(self):
        return self.mNRepeats
    def GetDiagonal(self):
        return self.mDiagonal
    def GetLength(self):
        return self.mLength
    def GetScore(self):
        return self.mScore
    def GetLevel(self):
        return self.mLevel
    def GetRepeatUnits( self ):
        return self.mRepeatUnits
    
    def Print( self ):
        print  string.join ( 
            map( str, (self.mNRepeats, self.mScore, self.mLength, self.mLevel)) ,
            "\t" )
        print self.mRepeatUnits
        


    ##------------------------------------------------------------------------            
    def BuildResult( self, file ):
        """parse the output file and return a tuple of repeats.
        """

        repeats = []
        
        while 1:
            line = file.readline()
            if not line: break
            
            if line[:-1] == "repeatfinder finished without problems":
                break
            elif line[:-1] == "Not enough dots given":
                break
            elif line[:-1] == "No repeats found":                               # Version 15a
                break
            else:
                try:
                    repeats.append( RadarRepeat( file ))
                except ValueError:
                    break

        return tuple(repeats)

##------------------------------------------------------------
def multiple_fasta_iterator (filenames, regex_identifier = None ):
    """iterate over multiple fasta files."""
    
    identifier = None
        
    for filename in filenames:

        if filename[-3:] == ".gz":
            infile = gzip.open( filename, "r" )
        elif filename == "-":
            infile = sys.stdin
        else:
            infile = open( filename, "r")

        for line in infile:

            if line[0] == "#":  continue

            if line[0] == ">" :

                if identifier:
                    yield identifier, "".join(fragments)
                    
                if regex_identifier:
                    try:
                        identifier = re.search(regex_identifier, line[1:-1]).groups()[0]
                    except AttributeError:
                        raise "could not parse identifier from line %s - check the input" % line[1:-1]
                else:
                    identifier = re.split("\s", line[1:-1])[0]
                fragments = []
            else:
                fragments.append (re.sub( "\s", "", line.strip() ) )

        if filename != "-": infile.close()

    if identifier:
        yield identifier, "".join(fragments)
    raise StopIteration

def runLFasta( filename_sequence, filename_output, lfasta_options ):
    
    retcode = subprocess.call( "lfasta %s %s %s > %s " % \
                               (filename_sequence, filename_sequence, lfasta_options, filename_output),
                               stderr=subprocess.PIPE,
                               shell=True)
    if retcode != 0:
        raise IOError( "error while running lfasta: code = %i" % retcode )    

def main():

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "-a", "--filename-fasta", dest="filename_fasta", type="string", action="append",
                       help="filename with fasta sequence [default=%default]")

    parser.add_option("-v", "--verbose", dest="loglevel", type="int",
                      help="loglevel [%default]. The higher, the more output." )

    parser.add_option("-L", "--log", dest="stdlog", type="string",
                      help="file with logging information [default = stdout].",
                      metavar = "FILE"  )
    
    parser.add_option("-E", "--error", dest="stderr", type="string",
                      help="file with error information [default = stderr].",
                      metavar = "FILE"  )
    
    parser.add_option("-S", "--stdout", dest="stdout", type="string",
                      help="file where output is to go [default = stdout].",
                      metavar = "FILE"  )
        
    parser.set_defaults(
                        stderr = sys.stderr,
                        stdout = sys.stdout, 
                        stdlog = sys.stdout,
                        filename_fasta = [],
                        loglevel = 0,
                        tmpdir = "/tmp",
                        regex_identifier = "^(\S+)" )

    (options, args) = parser.parse_args()
             
    if options.stdout != sys.stdout: 
        options.stdout = open(options.stdout, "w")
    if options.stderr != sys.stderr:
        options.stderr = open(options.stderr, "w")
    if options.stdlog != sys.stdout:
        options.stdlog = open(options.stdlog, "a") 
 
    options.filename_fasta += args
    
    if len(options.filename_fasta) == 0:
        options.filename_fasta.append("-")
    
    pyradar.setLogLevel( options.loglevel )
    
    iterator = multiple_fasta_iterator( options.filename_fasta,
                                        regex_identifier = options.regex_identifier )
    
    for identifier, sequence in iterator:
    
        options.stdout.write(">%s\n" % identifier )
        tmpdir = tempfile.mkdtemp( prefix = options.tmpdir + "/" )
        filename_sequence = tmpdir + "/seq"
        filename_lfasta1 = tmpdir + "/lfasta1"
        filename_lfasta2 = tmpdir + "/lfasta2"            
        outfile = open(filename_sequence, "w")
        outfile.write(sequence)
        outfile.close()
        
        runLFasta( filename_sequence, filename_lfasta1, "")
        runLFasta( filename_sequence, filename_lfasta2, "-f -12 -g -2 -s 250" )
        
        pyradar.run_from_files( filename_sequence, 
                                filename_sequence,
                                filename_lfasta1, 
                                filename_lfasta2 )
        
        shutil.rmtree( tmpdir )
    
if __name__ == '__main__':
    main()
    






