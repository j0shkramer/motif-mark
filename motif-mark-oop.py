#!/usr/bin/env python

import cairo, argparse, sys, random, re

##################
#    ARGPARSE    #
##################

def get_args():
    parser = argparse.ArgumentParser(description="A script to take a FASTA file and a text file with binding motifs to visualize binding motifs on sequences")
    parser.add_argument("-f", "--fasta", help="FASTA file", type=str)
    parser.add_argument("-m", "--motifs", help="File of binding motifs", type=str)
    args = parser.parse_args()

    if not args.fasta.lower().endswith((".fasta", ".fa", ".fna")):
        sys.exit(f'Error: Input file {args.fasta} must be a FASTA file')

    return args


##################
#   CLASSES      #
##################

class Motif:
    def __init__(self, the_binding_motif:str, the_color:list[float], the_x_cord:int, the_y_cord:int, the_header_idx:int):
        self.binding_motif = the_binding_motif
        self.color = the_color
        self.x_cord = the_x_cord
        self.y_cord = the_y_cord
        self.header_idx = the_header_idx
    
    def set_x(self, new_x_cord:int):
        self.x_cord = new_x_cord

    def set_y(self, new_y_cord:int):
        self.y_cord = new_y_cord

    def set_color(self, new_color:list):
        self.color = new_color

class NucleotideSequence:
    def __init__(self, the_start:int, the_end:int, the_y_cord:int, exon_or_intron:str, the_header_idx:int, the_seq:str):
        self.start = the_start
        self.end = the_end
        self.y_cord = the_y_cord
        self.type = exon_or_intron
        self.header_idx = the_header_idx
        self.seq = the_seq

    def set_start(self, new_start:int):
        self.start = new_start   

    def set_end(self, new_end:int):
        self.end = new_end

    def set_y(self, new_y_cord:int):
        self.y_cord = new_y_cord

    def set_intron_exon(self, intron_or_exon:str):
        self.type = intron_or_exon

class Header:
    def __init__(self, the_header_name:str, the_x_cord:int, the_y_cord:int, the_header_idx:int):
        self.name = the_header_name
        self.x_cord = the_x_cord
        self.y_cord = the_y_cord
        self.header_idx = the_header_idx

    def set_name(self, new_name:str):
        self.name = new_name

    def set_x(self, new_x_cord:int):
        self.x_cord = new_x_cord

    def set_y(self, new_y_cord:int):
        self.y_cord = new_y_cord


IUPAC_TO_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "[UT]",
    "U": "[UT]",
    "W": "[AT]",      
    "S": "[CG]",     
    "K": "[GT]",
    "M": "[AC]",
    "R": "[AG]",
    "Y": "[CT]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


##################
#HELPER FUNCTIONS#
##################

def read_motifs(motif_file:str):
    # Hold of the motifs we are looking for and their randomly assigned RGB value
    motif_colors: dict[str,list] = {}
    motif_list: list = []
    with open(motif_file, "r") as motifs:
        for line in motifs:
            line = line.strip()
            motif_list.append(line.upper())
            motif_colors[line.upper()] = [random.uniform(0,1), 
                                      random.uniform(0,1), 
                                      random.uniform(0,1)]
    return motif_colors, motif_list

def process_intron(intron_list:list, nucleotide_idx:int, header_ct:int, curr_seq:str):
    end = 100 + nucleotide_idx
    curr_intron = NucleotideSequence(end - len(curr_seq), end, 100, "Intron", header_ct, curr_seq)
    intron_list.append(curr_intron)
    return intron_list

def process_exon(exon_list:list, nucleotide_idx:int, header_ct:int, curr_seq:str):
    end = 100 + nucleotide_idx
    curr_exon = NucleotideSequence(end - len(curr_seq), end, 100, "Intron", header_ct, curr_seq)
    exon_list.append(curr_exon)
    return exon_list

def process_header(line:str, header_list: list, header_name:str, header_ct:int):
    header_name = line[1:]
    header_ct += 1
    curr_header = Header(header_name, 25, 75, header_ct)
    header_list.append(curr_header)
    return header_list

def process_motif(motifs_seen:list, read_ct:int, motif_colors:dict, motif_hits:list, motif:str):
    for hits in motif_hits:
        curr_motif = Motif(motif, motif_colors[motif], 100 + hits.start(), 100, read_ct)
        motifs_seen.append(curr_motif)
    return motifs_seen

def read_fasta(fasta_life:str, intron_list: list, exon_list: list, motifs_seen: list, header_list: list, motif_colors: list, motif_list: list, read_list: list):
    with open(fasta_life, "r") as fasta:
        curr_sequence:str = ""
        curr_nucelotide:str = ""
        nucleotide_idx: int = 0
        header_ct:int = 0
        curr_read = ""
        for line in fasta:
            line = line.strip()
            if line.startswith(">") and curr_sequence == "":
                header_list = process_header(line, header_list, line[1:], header_ct)
                header_ct += 1
                nucleotide_idx = 0
                curr_nucelotide = ""
                curr_sequence = ""
            elif line.startswith(">") and curr_sequence.islower():
                curr_read += curr_sequence
                header_list = process_header(line, header_list, line[1:], header_ct)
                intron_list = process_intron(intron_list, nucleotide_idx, header_ct, curr_sequence)
                read_list.append(curr_read.upper())
                header_ct += 1
                nucleotide_idx = 0
                curr_nucelotide = ""
                curr_sequence = ""
                curr_read = ""
            elif line.startswith(">") and curr_sequence.isupper():
                curr_read += curr_sequence
                header_list = process_header(line, header_list, line[1:], header_ct)
                exon_list = process_exon(exon_list, nucleotide_idx, header_ct, curr_sequence)
                read_list.append(curr_read.upper())
                header_ct += 1
                nucleotide_idx = 0
                curr_nucelotide = ""
                curr_sequence = ""
                curr_read = ""
            else:
                for nucleotide in line:
                    if curr_nucelotide != "":
                        if ((curr_nucelotide.isupper() and nucleotide.isupper()) or (curr_nucelotide.islower() and nucleotide.islower())):
                            curr_sequence += nucleotide
                            curr_nucelotide = nucleotide
                            nucleotide_idx += 1
                        else:
                            if curr_sequence.islower():
                                intron_list = process_intron(intron_list, nucleotide_idx, header_ct, curr_sequence)
                                curr_read += curr_sequence
                                curr_nucelotide = ""
                                curr_sequence = ""
                            else:
                                exon_list = process_exon(exon_list, nucleotide_idx, header_ct, curr_sequence)
                                curr_read += curr_sequence
                                curr_nucelotide = ""
                                curr_sequence = ""
                                
                    else:
                        curr_sequence += nucleotide
                        curr_nucelotide = nucleotide
                        nucleotide_idx += 1
        if curr_sequence.islower():
            curr_read += curr_sequence
            intron_list = process_intron(intron_list, nucleotide_idx, header_ct, curr_sequence)
            read_list.append(curr_read.upper())
        else:
            curr_read += curr_sequence
            exon_list = process_intron(exon_list, nucleotide_idx, header_ct, curr_sequence)
            read_list.append(curr_read.upper())
    read_ct:int = 0
    for read in read_list:
        read_ct += 1
        for motif in motif_list:
            regex_motif = ""
            for base in motif:
                regex_motif += IUPAC_TO_REGEX[base]
            regex_motif = re.compile(regex_motif)
            motif_hits = re.finditer(regex_motif, read)
            motifs_seen = process_motif(motifs_seen, read_ct, motif_colors, motif_hits, motif)
    return intron_list, exon_list, motifs_seen, header_list, read_list

def visualize_motifs(intron_list: list, exon_list: list, motifs_seen: list, header_list: list, filename:str, motif_colors:dict):
    reads = len(header_list)
    adjustment_idx: int = int(1000 / (reads - 1))
    figure_name = filename.replace(".", "_")
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,1200,1500)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1,1,1)
    ctx.rectangle(0,0,1200,1500)
    ctx.fill()
    ctx.stroke()

    for intron in intron_list:
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(10)
        ctx.move_to(intron.start, (intron.y_cord + (adjustment_idx * (intron.header_idx - 1))))
        ctx.line_to(intron.end, (intron.y_cord + (adjustment_idx * (intron.header_idx - 1))))
        ctx.stroke()

    for exon in exon_list:
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(30)
        ctx.move_to(exon.start, (exon.y_cord + (adjustment_idx * (exon.header_idx - 1))))
        ctx.line_to(exon.end, (exon.y_cord + (adjustment_idx * (exon.header_idx - 1))))
        ctx.stroke()

    for motif in motifs_seen:
        ctx.set_source_rgb(*motif.color)
        ctx.rectangle(motif.x_cord, 
                      (motif.y_cord + (adjustment_idx * (motif.header_idx - 1))) - 15,
                      len(motif.binding_motif),
                      30)
        ctx.fill()

    for header in header_list:
        ctx.set_source_rgb(0,0,0)
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(20)
        ctx.move_to(header.x_cord, header.y_cord + ((adjustment_idx * (header.header_idx - 1))))
        ctx.show_text(f'{header.name}')

    legend_adjustment_idx = 300 / len(motif_colors)
    legend_ct = 0

    for motif in sorted(motif_colors.keys()):
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(max(8, legend_adjustment_idx - 5))
        y = 1200 + (legend_ct * legend_adjustment_idx)
        x = 25
        ctx.set_source_rgb(*motif_colors[motif])
        ctx.rectangle(x, y, 25, max(8, legend_adjustment_idx - 5))              
        ctx.fill()
        label = str(motif)
        ext = ctx.text_extents(label)
        swatch_center_y = y + (max(8, legend_adjustment_idx - 5)) / 2
        text_y = swatch_center_y - (ext.height / 2 + ext.y_bearing)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(55, text_y)
        ctx.show_text(f'{motif}')
        legend_ct += 1

    ctx.set_source_rgb(0,0,0)
    ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(40)
    ctx.move_to(25, 40)
    ctx.show_text(f'Binding Motifs Identified in {filename}')

    surface.write_to_png(f'{figure_name}.png')



##################
#   MAIN BODY    #
##################

def main():
    args = get_args()
    intron_list: list[NucleotideSequence] = []
    exon_list: list[NucleotideSequence] = []
    motifs_seen: list[Motif] = []
    header_list: list[Header] = []
    read_list: list[str] = []
    motif_colors, motif_list = read_motifs(args.motifs)
    intron_list, exon_list, motifs_seen, header_list, read_list = read_fasta(args.fasta, intron_list, exon_list, motifs_seen, header_list, motif_colors, motif_list, read_list)
    t = visualize_motifs(intron_list, exon_list, motifs_seen, header_list, args.fasta, motif_colors)
    t

if __name__ == "__main__":
    main()