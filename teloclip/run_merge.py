#!/usr/bin/env python

from __future__ import print_function
from teloclip import __version__
from teloclip import log
import wrappers
import teloclip
import argparse
import sys

# Select version-compatible user input function
if not sys.version_info[0] > 2:
    def input(s):
        return raw_input(s)

def mainArgs():
    parser = argparse.ArgumentParser(description='xxxx.',prog='teloclip-merge')
    # Input options
    parser.add_argument('samfile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('--refIdx',type=str,required=True,help='Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`')
    parser.add_argument('--ref',type=str,required=True,help='Path to reference fasta.')
    # Output options
    parser.add_argument('--out',type=str,required=True,help='Output updated contigs to file.')
    parser.add_argument('--temp',type=str,default=None,help='Path to temp directory.')
    # Settings
    parser.add_argument('--minClip',type=int,default=1,help='Require clip to extend past ref contig end by at least N bases.')
    parser.add_argument('--maxBreak',type=int,default=50,help='Tolerate max N unaligned bases at contig ends.')
    parser.add_argument('--readType',default=None,choices=[None,'PB','ONT'],help='Set platform type for aligned long-reads. PB = PacBio, ONT = Oxford Nanopore Technologies.')
    # Third party apps
    parser.add_argument('--minimap2',type=str,default='minimap2',help='Path to minimap2.')
    parser.add_argument('--racon',type=str,default='racon',help='Path to racon.')
    parser.add_argument('--miniasm',type=str,default='miniasm',help='Path to miniasm.')
    # Version info
    parser.add_argument('--version', action='version',version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()
    return args



# Make temp fasta
# manageTemp(record=(name,seq), tempPath=filePath, scrub=False)
# Del fasta
# manageTemp(tempPath=filePath, scrub=True)


def clipdown(args,ID=None,alnList=None,contig=None,temp=None,mode=None):
    # Note: add flags for having already run error correction + 
    ''' Interactive contig extension'''
    log('Processing contig: %s' % str(ID))
    pad = len(str(len(alnList)))
    contigLen = len(contig)
    idx = 0
    if mode == "R":
        for aln in alnList: # i.e [(alnStart,alnEnd,rightClipLen,readSeq,readname)]
            idx += 1
            gap = contigLen - aln[1] # alnEnd
            overhang = aln[3][-aln[2]:]
            maplen = aln[1] - aln[0]
            teloclip.writeClip(idx,pad,gap,overhang,maplen)
    elif mode == "L":
        log("Reporting left end overhangs as reverse complement. \n")
        for aln in alnList: # i.e #[(alnStart,alnEnd,leftClipLen,readSeq,readname)]
            idx += 1
            gap = aln[0] - 1 # alnStart
            overhang = teloclip.revComp(aln[3][:aln[2]])
            maplen = aln[1] - aln[0]
            teloclip.writeClip(idx,pad,gap,overhang,maplen)
    # Prompt user to accept a read, skip contig end, or continue
    Q1 = input("Enter read ID to accept alignment [or RETURN to continue]: ")
    if Q1:
        targetIdx = int(Q1.strip().lstrip(0)) - 1
        if mode == "R":
            log("Selected alignment %s: %s" % (str(targetIdx + 1), alnList[targetIdx][3][-aln[2]:]))
        elif mode == "L":
            # Edit to report rc L end clip:
            # log("Selected alignment %s: %s" % (str(targetIdx + 1), alnList[targetIdx][3][-aln[2]:]))
            pass
        if input("Is this correct? [Y/N]: ") == 'Y':
            # When an alignment is accepted for clipdown:
            # Select extraction protocol for L or R contig end.
            # Check for terminal gaps (may need to trim from raw contig)
            # Extract clipped OH from aligned read
            # Append to contig
            # Update header string with L/R mod tag "TC:Left:5C,E400" = left end clipped 5 bases, extended by 400 bases.
            # return modified contig seq and header line
            # return (header, seq)
            pass
    else:
        # If rejected all candidate alignments and running int A or CA mode, end here.
        if 'A' in mode:
            return None
        else:
            Q2 = input("Type 'S' to skip this contig end or RETURN to continue]: ")
            if Q2 == 'S':
                return None
            else:
                # Prompt user to exclude suspect alignments
                Q3 = input("To exclude alignments enter IDs as space-delimited list [or RETURN to continue]: ")
                if Q3:
                    # Parse idx IDs, split on whitespace, scrub left zeros
                    # -1 to match list index
                    exclude = [int(x.lstrip('0')) - 1 for x in Q3.split()]
                    keepAln = teloclip.filterList(alnList, exclude)
                else: 
                    keepAln = alnList
                if not alnList:
                    log("No remaining alignments. Moving on.")
                    return None
                # If running in C mode (alignments already from corrected reads) then only present option 'A'
                if mode == 'C':
                    log("There are %s remaining alignments. Aligned reads have already been error corrected. \n\
                        You may choose from the following options:\n\
                        - Assemble raw reads and align de novo contigs to the reference contig [A]" % str(len(keepAln)))
                else:
                    log("There are %s remaining alignments. You may choose from the following options: \n\
                        - Error-correct reads and re-align to the reference contig [C] \n\
                        - Assemble raw reads and align de novo contigs to the reference contig [A] \n\
                        - Error-correct reads before running the assembly step [CA] \n" % str(len(keepAln)))
                Q4 = input("Select an option: ").strip()
                if Q4 == 'C' or 'CA':
                    # Write out accepted reads
                    # Write contig
                    # Align reads to self `minimap2 -x ava-pb  reads.fq reads.fq > ovlp.paf`
                    # Correct with Racon
                    # Store corrected reads
                    if Q4 == 'C':
                        # Align corrected reads back to single contig
                        # Import new samfile
                        # Call clipdown() with 'C' flag
                        # clipdown(args,ID=None,alnList=None,contig=None,temp=None,mode='C')
                        # return selected alignment
                        pass
                    elif Q4 == 'CA':
                        # Run miniasm
                        # Align any assembled contigs to current ref contig
                        # Call clipdown() with 'CA' flag
                        # clipdown(args,ID=None,alnList=None,contig=None,temp=None,mode='CA')
                        # return selected alignment
                        pass
                elif Q4 == 'A':
                    # Write out accepted reads
                    # Write contig
                    # Run miniasm
                    # Align any assembled contigs to current ref contig
                    # Call clipdown() with 'A' flag
                    # clipdown(args,ID=None,alnList=None,contig=None,temp=None,mode='A')
                    # return selected alignment
                    pass
                else: 
                    log("No valid option detected. Moving on.")
                    return None




def main():
    # Get cmd line args
    args = mainArgs()

    # Check for required programs.
    tools = [args.minimap2,args.racon,args.miniasm]
    missing_tools = []
    for tool in tools:
        missing_tools += teloclip.missing_tool(tool)
    if missing_tools:
        log('WARNING: Some tools required by teloclip-merge could not be found: ' +
              ', '.join(missing_tools))
        log('You may need to install them to use all features.')
    
    # Load reference genome to dict with format
    # {'name':(header, seq)}
    contigs = teloclip.fasta2dict(args.ref)
    # Load ref contigs lengths as dict
    contigInfo = teloclip.read_fai(args.refIdx)

    # Load alignments from samfile or stdin
    #{'contig':{"L":[(alnStart,alnEnd,leftClipLen,readSeq,readname)],"R":[(alnStart,alnEnd,rightClipLen,readSeq,readname)]}}
    alignments = teloclip.loadSam(samfile=args.samfile,contigs=contigInfo,maxBreak=args.maxBreak,minClip=args.minClip)

    # Create output and temp paths as required
    # Update to use args.out and args.temp
    tempDir = teloclip.dochecks(args)

    for name in contigs.keys():
        # Count alignments for each end of contig
        leftCount = len(alignments[name]['L'])
        rightCount = len(alignments[name]['R'])
        newSeqR = None
        newSeqL = None
        # Process right end alignments if found
        if rightCount:
            # Run interactive mode
            # If edits made return updated (header, seq) object
            newSeqR = clipdown(args,ID=name,alnList=alignments[name]['R'],contig=contigs[name],temp=tempDir,mode="R")
            if newSeqR:
                contigs[name] = newSeqR
        # Process left end alignments if found
        if leftCount:
            # Run interactive mode
            # If edits made return updated (header, seq) object
            newSeqL = clipdown(args,ID=name,alnList=alignments[name]['L'],contig=contigs[name],temp=tempDir,mode="L")
            if newSeqL:
                contigs[name] = newSeqL
        # Update counters
            # Left only
            # Right only
            # Both
            # Neither

    # Write updated contigs to outfile
    output = open(args.out, 'w')
    for name in contigs:
       header,seq = contigs[name]
       teloclip.writefasta(output,name,seq,header=header,length=80)
    output.close()

    # Report counters


