#!/usr/bin/env python

import cairo, argparse, sys

##################
#    ARGPARSE    #
##################

def get_args():
    parser = argparse.ArgumentParser(description="A script to take a FASTA file and a text file with binding motifs to visualize binding motifs on sequences")
    parser.add_argument("-f", "--fasta", help="FASTA file", type=str)
    parser.add_argument("-m", "--motifs", help="File of binding motifs", type=str)
    args = parser.parse_args()

    if not args.fasta.lower().endswith(".fasta", ".fa", ".fna"):
        sys.exit(f'Error: Input file {args.fasta} must be a FASTA file')

    return args


##################
#   CLASSES      #
##################

class Motif:
    def __init__(self, the_binding_motif:str, the_color:list, the_x_cord:int, the_y_cord:int):
        self.binding_motif = the_binding_motif
        self.color = the_color
        self.x_cord = the_x_cord,
        self.y_cord = the_y_cord
    
    def set_x(self, new_x_cord:int):
        self.x_cord = new_x_cord

    def set_y(self, new_y_cord:int):
        self.y_cord = new_y_cord

    def set_color(self, new_color:list):
        self.color = new_color

class NucleotideSequence:
    def __init__(self, the_start:int, the_end:int, the_y_cord:int, exon_or_intron:str):
        self.start = the_start
        self.end = the_end
        self.y_cord = the_y_cord
        self.type = exon_or_intron

    def set_start(self, new_start:int):
        self.start = new_start   

    def set_end(self, new_end:int):
        self.end = new_end

    def set_y(self, new_y_cord:int):
        self.y_cord = new_y_cord

    def set_intron_exon(self, intron_or_exon:str):
        self.type = intron_or_exon

class Header:
    def __init__(self, the_header_name:str, the_x_cord:int, the_y_cord:int):
        self.name = the_header_name
        self.x_cord = the_x_cord
        self.y_cord = the_y_cord

    def set_name(self, new_name:str):
        self.name = new_name

    def set_x(self, new_x_cord:int):
        self.x_cord = new_x_cord

    def set_y(self, new_y_cord:int):
        self.y_cord = new_y_cord


##################
#HELPER FUNCTIONS#
##################


##################
#   MAIN BODY    #
##################

def main():
    return 0

if __name__ == "__main__":
    main()