#!/usr/bin/env python

import sys


def construct(*args, **kwargs):
    start = int(kwargs.get('start'))
    end = int(kwargs.get('end'))
    
    for i in range(start, end):
        PID = "P"+str(i)
        makeProjectBarcode(PID)
    
def makeProjectBarcode(PID):
    print "^XA" #start of label
    print "^DFFORMAT^FS" #download and store format, name of format, end of field data (FS = field stop)
    print "^LH0,0" # label home position (label home = LH)
    print "^FO400,20^AFN,60,20^FN1^FS" #AF = assign font F, field number 1 (FN1), print text at position field origin (FO) rel. to home
    print "^FO120,5^BCN,70,N,N^FN2^FS" #BC=barcode 128, field number 2, Normal orientation, height 70, no interpreation line. 
    print "^XZ" #end format
    
    for i in range (1,6):
        PlateID="P"+str(i)
        plateCode=PID+PlateID
        print "^XA" #start of label format
        print "^XFFORMAT^FS" #label home posision
        print "^FN1^FD"+plateCode+"^FS" #this is readable
        print "^FN2^FD"+plateCode+"^FS" #this is the barcode
        print "^XZ"

def getArgs():
    from optparse import OptionParser
    ''' Options '''
    parser = OptionParser(
        description="""Tool for constructing barcode labels for NGI Genomics Projects""",
        usage='-s <start project ID> -e <end project id>',
        version="%prog V0.02 sverker.lundin@scilifelab.se")
    parser.add_option('-s', '--start', type=int,
        help='the starting project ID (numeric, e.g. 123)')
    parser.add_option('-e', '--end', type=int,
        help='the last project ID (numerig, e.g. 234)')
    return parser

def main():
    parser = getArgs()
    (options, args) = parser.parse_args()
    if not (options.start):
        print >> sys.stderr, 'Usage: %s %s' % \
            (parser.get_prog_name(), parser.usage)
        
    if not (options.end):
        options.end = options.start +1
    if options.start >= options.end:
        print >> sys.stderr, 'end value has to be > start value'
        sys.exit()
    try:
        construct(start=options.start, end=options.end)
    except KeyboardInterrupt:
        print >> sys.stderr, 'Interupted!'
        quit()
   
if __name__ == '__main__':
    main()