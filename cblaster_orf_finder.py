#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import csv
import subprocess
from Bio import SeqIO, Entrez

from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio import pairwise2
from collections import defaultdict
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import matplotlib.pyplot as plt
from C_blaster_extraction_custom import NCBICDSExtractor
from weblogo import LogoData, LogoOptions, LogoFormat, eps_formatter

# Configuration

Entrez.email = "friederike@biermann-erfurt.de"  # Replace with your email for NCBI Entrez queries

# Define a list of potential RBS patterns
rbs_patterns = ["AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA"]
alternative_start_codons = ["ATG", "GTG", "TTG"]
stop_codons = ["TAA", "TAG", "TGA"]



class ORFProcessor:
    def __init__(self, orfs):
        """
        Initializes the ORFProcessor class with a list of ORFs.
        
        Args:
            orfs (list): List of ORFs and their details (orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer, record_id).
        """
        self.orfs = orfs

    def save_orfs_to_csv(self, csv_output_file):
        """
        Saves ORF information to a CSV file.

        Args:
            csv_output_file (str): Path to the CSV file to write the ORFs.
        """
        with open(csv_output_file, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            
            # Write the header
            csv_writer.writerow(['Record ID', 'ORF Sequence', 'Start', 'End', 'Strand', 'Protein Sequence', 'RBS Type'])
            
            # Write ORF data into the CSV file
            for orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, _, record_id in self.orfs:
                csv_writer.writerow([
                    record_id,                # Record ID from GenBank
                    str(orf_seq),             # ORF nucleotide sequence
                    start_pos,                # Start position of ORF
                    end_pos,                  # End position of ORF
                    strand,                   # Strand ("+" or "-")
                    str(protein_seq),         # Protein sequence (translated ORF)
                    rbs_type                  # Ribosome Binding Site (RBS) type or GA-rich
                ])
        print(f"ORFs written to {csv_output_file}")

    def filter_orfs_by_length(self, target_orf, tolerance=0.3):
        """
        Filters ORFs based on length similarity (±30% tolerance).

        Args:
            target_orf (tuple): The target ORF (orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer, record_id).
            tolerance (float): Percentage tolerance for ORF length similarity (default is 30%).

        Returns:
            list: Filtered ORFs that fall within the length tolerance.
        """
        target_length = len(target_orf[1])
        length_range = (target_length * (1 - tolerance), target_length * (1 + tolerance))
        
        return [orf for orf in self.orfs if length_range[0] <= len(orf[1]) <= length_range[1] and orf[-1] != target_orf[-1]]

    def is_similar_orf(self, orf1, orf2):
        """
        Checks if two ORFs are similar based on sequence similarity (global or 6 AA motif).

        Args:
            orf1 (tuple): ORF tuple containing the protein sequence.
            orf2 (tuple): ORF tuple containing the protein sequence.

        Returns:
            bool: True if the ORFs are similar, False otherwise.
        """
        seq1 = orf1[1][:-1]
        seq2 = orf2[1][:-1]

        # Global sequence alignment for 30% similarity
        alignments = pairwise2.align.globalxx(seq1, seq2, score_only=True)
        print(alignments)
        similarity = alignments / min(len(seq1), len(seq2))
        
        if similarity >= 0.5:
            return True

        # Search for 6 AA motif with 80% similarity
        for i in range(len(seq1) - 5):
            motif1 = seq1[i:i + 6]
            for j in range(len(seq2) - 5):
                motif2 = seq2[j:j + 6]
                motif_similarity = pairwise2.align.globalxx(motif1, motif2, score_only=True) / 6
                
                if motif_similarity >= 0.8:
                    print(seq1, seq2, motif1, motif2, "similar")
                    return True

        return False

    def cluster_orfs(self):
        """
        Clusters ORFs based on length and sequence similarity.

        Returns:
            list: A list of ORF clusters.
        """
        clusters = []
        visited = set()

        for orf in self.orfs:
            if orf in visited:
                continue

            # Start a new cluster with the current ORF
            cluster = [orf]
            visited.add(orf)

            # Filter ORFs by length (±30% tolerance)
            length_filtered_orfs = self.filter_orfs_by_length(orf)

            # Add ORFs to the cluster based on sequence similarity
            for candidate_orf in length_filtered_orfs:
                if candidate_orf not in visited and self.is_similar_orf(orf, candidate_orf):
                    cluster.append(candidate_orf)
                    visited.add(candidate_orf)

            clusters.append(cluster)

        return clusters

    def align_orfs(self, cluster, output_file):
        """
        Aligns ORFs using MUSCLE and writes the alignment to a file.

        Args:
            cluster (list): List of ORFs in the cluster.
            output_file (str): Path to the output file for the alignment.

        Returns:
            str: Path to the aligned FASTA file.
        """
        # Write the ORFs to a temporary FASTA file
        with open("temp_cluster.fasta", "w") as fasta_handle:
            for i, orf in enumerate(cluster):
                fasta_handle.write(f">orf_{i}\n{orf[1]}\n")

        # Run MUSCLE alignment
        muscle_cmd = [
            "muscle", 
            "-in", "temp_cluster.fasta", 
            "-out", output_file
        ]
        
        try:
            # Execute MUSCLE command
            subprocess.check_call(muscle_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"Error running MUSCLE on {output_file}")
            return None

        # Cleanup temporary file
        os.remove("temp_cluster.fasta")
        
        return output_file

    def create_sequence_logo(self, alignment_file, output_logo_file):
        """
        Creates a sequence logo from an aligned FASTA file.

        Args:
            alignment_file (str): Path to the aligned FASTA file.
            output_logo_file (str): Path to the output sequence logo file.
        """
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Convert alignment to amino acid sequences (Bio.AlignIO can handle this directly)
        sequences = [str(record.seq) for record in alignment]
        
        # Prepare the data for WebLogo
        data = LogoData.from_seqs(sequences)
        options = LogoOptions()
        options.title = "Sequence Logo"
        
        format = LogoFormat(data, options)
        
        # Write the logo to the output file
        with open(output_logo_file, "wb") as logo_handle:
            eps_formatter(data, format, logo_handle)

        print(f"Sequence logo saved to {output_logo_file}")

    def process_clusters_and_logos(self, output_root):
        """
        Main function to cluster ORFs, align them, and create sequence logos.

        Clusters the ORFs, sorts clusters by size, aligns the largest 5 clusters, and generates sequence logos.
        """
        clusters = self.cluster_orfs()
        
        # Sort clusters by size and take the top 15 largest clusters
        sorted_clusters = sorted(clusters, key=len, reverse=True)[:15]

        for i, cluster in enumerate(sorted_clusters, 1):
            print(f"Processing cluster {i} with {len(cluster)} ORFs...")
            
            # Align the ORFs in the cluster
            aligned_fasta = self.align_orfs(cluster, f"{output_root}_aligned_cluster_{i}.fasta")
            
            # Generate a sequence logo from the aligned sequences
            #self.create_sequence_logo(aligned_fasta, f"{output_root}_sequence_logo_cluster_{i}.eps")
            
    def save_orfs_to_fasta(self, fasta_output_file):
        """
        Saves ORF information to a FASTA file.

        Args:
            orfs (list): List of ORFs and their details.
            fasta_output_file (str): Path to the FASTA file to write the ORFs.
        """
        fasta_records = []
        
        for i, (orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, _, record_id) in enumerate(self.orfs, 1):
            header = f"{record_id}_ORF_{i}_start_{start_pos}_end_{end_pos}_strand_{strand}_rbs_{rbs_type}"
            # Create a SeqRecord object for each ORF
            fasta_record = SeqRecord(Seq(str(orf_seq)), id=header, description="")
            fasta_records.append(fasta_record)

        # Write all ORF records to a FASTA file
        with open(fasta_output_file, "w") as output_handle:
            SeqIO.write(fasta_records, output_handle, "fasta")
        
        print(f"ORFs saved to {fasta_output_file}")

class ORFFinder:
    def __init__(self, gbk_file):
        """
        Initializes the ORFFinder with a GenBank file and sets up the output directory
        based on the filename of the GenBank file.

        Args:
            gbk_file (str): Path to the GenBank file.
        """
        self.gbk_file = gbk_file
        
        # Define RBS patterns, start codons, and stop codons
        self.rbs_patterns = ["AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA"]
        self.alternative_start_codons = ["ATG", "GTG", "TTG"]
        self.stop_codons = ["TAA", "TAG", "TGA"]

    def ga_content(self, seq):
        """
        Calculate G and A content percentage.
        """
        ga_count = seq.count('G') + seq.count('A')
        return ga_count / len(seq) if len(seq) > 0 else 0

    def search_orfs(self, sequence, is_reverse=False):
        """

        Search for ORFs given a sequence and RBS patterns.
        """
        found_orfs = []
        seq_length = len(sequence)
        list_start_positions = []
        list_end_positions = []
        
        # Search for RBS and ORFs in the forward direction
        for rbs_pattern in self.rbs_patterns:
            for rbs_start in range(seq_length - len(rbs_pattern)):
                if sequence[rbs_start : rbs_start + len(rbs_pattern)] == rbs_pattern:
                    for start_codon in range(rbs_start + len(rbs_pattern), min(rbs_start + len(rbs_pattern) + 20, seq_length - 3)):
                        if sequence[start_codon : start_codon + 3] in self.alternative_start_codons:
                            for end_codon in range(start_codon + 12, min(start_codon + 60, seq_length), 3):
                                if sequence[end_codon : end_codon + 3] in self.stop_codons:
                                    orf_seq = Seq(sequence[start_codon : end_codon + 3])
                                    protein_seq = orf_seq.translate()
                                    strand = "-" if is_reverse else "+"
                                    if is_reverse:
                                        start_pos = seq_length - end_codon - 3
                                        end_pos = seq_length - start_codon
                                    else:
                                        start_pos = start_codon
                                        end_pos = end_codon + 3
                                    if start_pos not in list_start_positions and end_pos not in list_end_positions:
                                        found_orfs.append((orf_seq, protein_seq, strand, start_pos, end_pos, rbs_pattern, rbs_start - start_codon))
                                        list_start_positions.append(start_pos)
                                        list_end_positions.append(end_pos)
                                        break

        # If no RBS pattern is found, search for GA-rich regions (>= 80% GA)
        if not found_orfs:
            for start_codon in range(seq_length - 20):
                if sequence[start_codon : start_codon + 3] in self.alternative_start_codons:
                    for buffer in range(0, 20):
                        upstream_region = sequence[max(0, start_codon - 8 - buffer):start_codon - buffer]
                        if len(upstream_region) == 8 and self.ga_content(upstream_region) >= 0.8:
                            for end_codon in range(start_codon + 12, min(start_codon + 60, seq_length), 3):
                                if sequence[end_codon : end_codon + 3] in self.stop_codons:
                                    orf_seq = Seq(sequence[start_codon : end_codon + 3])
                                    protein_seq = orf_seq.translate()
                                    strand = "-" if is_reverse else "+"
                                    if is_reverse:
                                        start_pos = seq_length - end_codon - 3
                                        end_pos = seq_length - start_codon
                                    else:
                                        start_pos = start_codon
                                        end_pos = end_codon + 3
                                    if start_pos not in list_start_positions and end_pos not in list_end_positions:
                                        found_orfs.append((orf_seq, protein_seq, strand, start_pos, end_pos, f"GA-rich", buffer))
                                        list_start_positions.append(start_pos)
                                        list_end_positions.append(end_pos)
                                        break
        return found_orfs

    def find_orfs_with_criteria(self, dna_sequence):
        """
        Finds ORFs with RBS and sequence criteria for both forward and reverse strands.
        """
        found_orfs = []
        orfs_forward = self.search_orfs(dna_sequence)
        reverse_complement = str(Seq(dna_sequence).reverse_complement())
        orfs_reverse = self.search_orfs(reverse_complement, is_reverse=True)
        found_orfs = orfs_forward + orfs_reverse
        return found_orfs


    def process_gbk_file(self):
        """
        Processes the GenBank file, finds ORFs, re-annotates already annotated genes
        with less than 100 amino acids, allows up to 20% overlap, and returns ORF details.
        Removes ORFs with identical sequences coming from the same record.
        """
        orfs_found = []
        updated_records = []
        
        for record in SeqIO.parse(self.gbk_file, "genbank"):
            dna_sequence = str(record.seq)
            existing_orfs_set = set()

            # Process existing CDS features and find new ORFs
            self.annotate_existing_short_cdss(record, orfs_found, existing_orfs_set)
            self.find_and_add_new_orfs(record, dna_sequence, orfs_found, existing_orfs_set)

            updated_records.append(record)

        # Remove duplicates based on sequence and record ID
        orfs_found = self.remove_duplicate_orfs(orfs_found)

        return orfs_found, updated_records


    def remove_duplicate_orfs(self, orfs_found):
        """
        Removes ORFs with identical sequences and that come from the same record.

        Args:
            orfs_found (list): List of tuples containing ORF details.
                            Each ORF tuple contains (orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer, record_id).
        
        Returns:
            list: A list of unique ORFs, with duplicates removed.
        """
        unique_orfs = []
        seen = set()  # Set to track seen (orf_seq, record_id) pairs

        for orf in orfs_found:
            orf_seq = orf[0]
            record_id = orf[7]  # The record ID is the last item in the tuple

            # Only add ORF if it has not been encountered yet
            print(orf_seq, record_id)
            if (orf_seq, record_id) not in seen:
                unique_orfs.append(orf)
                seen.add((orf_seq, record_id))  # Mark this (sequence, record) pair as seen

        return unique_orfs


    def annotate_existing_short_cdss(self, record, orfs_found, existing_orfs_set):
        """
        Re-annotate existing CDS features with less than 100 amino acids and add them to orfs_found.
        
        Args:
            record (SeqRecord): The GenBank record to process.
            orfs_found (list): List to store found ORFs (new and existing short CDS).
            existing_orfs_set (set): A set to track existing ORFs and prevent duplicates.
        """
        for feature in record.features:
            if feature.type == "CDS":
                protein_seq = self.get_protein_sequence(record, feature)
                
                # Check if the translated protein has less than 100 amino acids
                if len(protein_seq) < 100:
                    feature.qualifiers["note"] = "Re-annotated CDS, less than 100 AA"
                    orf_key = (int(feature.location.start), int(feature.location.end), feature.strand)

                    # Add to orfs_found if it hasn't been added yet
                    if orf_key not in existing_orfs_set:
                        orfs_found.append((
                            str(feature.extract(record.seq)),
                            protein_seq,
                            "+" if feature.strand == 1 else "-",
                            int(feature.location.start),
                            int(feature.location.end),
                            "Existing CDS < 100 AA",
                            "",  # buffer is left empty for existing features
                            record.id
                        ))
                        existing_orfs_set.add(orf_key)


    def get_protein_sequence(self, record, feature):
        """
        Extract or translate the protein sequence from a given CDS feature.
        
        Args:
            record (SeqRecord): The GenBank record containing the feature.
            feature (SeqFeature): The CDS feature to extract/translate the protein sequence.
        
        Returns:
            str: The translated protein sequence.
        """
        if "translation" in feature.qualifiers:
            return feature.qualifiers["translation"][0]
        else:
            protein_seq = str(feature.extract(record.seq).translate(to_stop=True))
            feature.qualifiers["translation"] = [protein_seq]
            return protein_seq


    def find_and_add_new_orfs(self, record, dna_sequence, orfs_found, existing_orfs_set):
        """
        Find new ORFs in the DNA sequence and add them if they don't overlap more than 20% with existing CDS.
        
        Args:
            record (SeqRecord): The GenBank record to update.
            dna_sequence (str): The DNA sequence of the record.
            orfs_found (list): List to store found ORFs (new and existing short CDS).
            existing_orfs_set (set): A set to track existing ORFs and prevent duplicates.
        """
        found_orfs = self.find_orfs_with_criteria(dna_sequence)

        for orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer in found_orfs:
            if not self.check_orf_overlap(record, start_pos, end_pos):
                self.add_new_orf(record, orfs_found, orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer)


    def check_orf_overlap(self, record, start_pos, end_pos):
        """
        Check if a new ORF overlaps more than 20% with any existing CDS feature in the record.
        
        Args:
            record (SeqRecord): The GenBank record to check for overlaps.
            start_pos (int): Start position of the new ORF.
            end_pos (int): End position of the new ORF.
        
        Returns:
            bool: True if the ORF overlaps more than 20%, False otherwise.
        """
        for feature in record.features:
            if feature.type == "CDS":
                overlap_start = max(start_pos, int(feature.location.start))
                overlap_end = min(end_pos, int(feature.location.end))
                overlap_length = max(0, overlap_end - overlap_start)
                
                existing_length = int(feature.location.end) - int(feature.location.start)
                orf_length = end_pos - start_pos

                overlap_percentage = overlap_length / orf_length
                
                if overlap_percentage > 0.20:
                    return True

        return False


    def add_new_orf(self, record, orfs_found, orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer):
        """
        Add a new ORF to the record and update the orfs_found list.
        
        Args:
            record (SeqRecord): The GenBank record to update.
            orfs_found (list): List to store found ORFs.
            orf_seq (str): DNA sequence of the ORF.
            protein_seq (str): Translated protein sequence.
            strand (str): The strand where the ORF is located ("+" or "-").
            start_pos (int): Start position of the ORF.
            end_pos (int): End position of the ORF.
            rbs_type (str): Ribosome binding site type.
            buffer (str): Buffer information (if applicable).
        """
        orf_feature = SeqFeature(
            FeatureLocation(start_pos, end_pos, strand=1 if strand == "+" else -1),
            type="CDS",
            qualifiers={"translation": str(protein_seq), "note": rbs_type}
        )
        record.features.append(orf_feature)
        orfs_found.append((orf_seq, protein_seq, strand, start_pos, end_pos, rbs_type, buffer, record.id))


    def run_orf_finder(self):
        """
        Runs the ORF finder and processes the GenBank file. This function will return the ORFs
        and the updated GenBank records for further processing.
        """
        print(f"Running ORF finder on {self.gbk_file}...")
        orfs, updated_records = self.process_gbk_file()
        print(f"ORF finding completed on {self.gbk_file}.")
        return orfs, updated_records

class CDSExtractor:
    def __init__(self, gbk_file):
        self.gbk_file = gbk_file

    def extract_cds(self, fasta_output, num_cds=4):
        """
        Extract up to `num_cds` CDS sequences, prioritizing removal from the ends, and including those that
        overlap with features that overlap with 'predicted_BGC_type'.
        Parameters:
        fasta_output (str): Path to the output FASTA file.
        num_cds (int): Number of CDS sequences to extract.
        """
        # Extract features that overlap with "predicted_BGC_type"
        overlapping_with_predicted_bgc = self._get_features_overlapping_predicted_bgc()

        # Extract CDS features that overlap with those features
        all_cds = self._get_cds_features(overlapping_with_predicted_bgc)

        # Split CDS features into overlapping and non-overlapping groups
        overlapping_cds, non_overlapping_cds = self._split_cds_by_overlap(all_cds)

        # Select CDS sequences, prioritizing middle ones and overlapping features
        selected_cds = self._select_cds(overlapping_cds, non_overlapping_cds, num_cds)

        # Write the final selection to the output file
        self._write_to_fasta(selected_cds, fasta_output)

    def _get_features_overlapping_predicted_bgc(self):
        """
        Parse the GenBank file and extract features that overlap with 'predicted_BGC_type' features,
        but exclude the 'predicted_BGC_type' features themselves.
        """
        overlapping_locations = []
        predicted_bgc_locations = []

        for record in SeqIO.parse(self.gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "predicted_BGC_type":
                    # Collect 'predicted_BGC_type' locations
                    start, end = feature.location.start, feature.location.end
                    predicted_bgc_locations.append((start, end))

        for record in SeqIO.parse(self.gbk_file, "genbank"):
            for feature in record.features:
                if feature.type != "predicted_BGC_type":
                    start, end = feature.location.start, feature.location.end
                    # Check if this feature overlaps with any 'predicted_BGC_type' region
                    overlaps_bgc = any(start <= bgc_end and end >= bgc_start for bgc_start, bgc_end in predicted_bgc_locations)
                    if overlaps_bgc:
                        overlapping_locations.append((start, end))

        return overlapping_locations

    def _get_cds_features(self, overlapping_with_predicted_bgc):
        """
        Parse the GenBank file and extract CDS features, marking those that overlap with features
        overlapping with 'predicted_BGC_type'.
        """
        all_cds = []
        for record in SeqIO.parse(self.gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and 'translation' in feature.qualifiers:
                    sequence = feature.qualifiers['translation'][0].strip()
                    if sequence:
                        gene_name = feature.qualifiers.get('gene', ['unknown_gene'])[0]
                        header = f">{record.id}_{gene_name}_{feature.location.start}-{feature.location.end}"
                        start, end = feature.location.start, feature.location.end
                        overlaps_with_overlapping = any(start <= overlap_end and end >= overlap_start for overlap_start, overlap_end in overlapping_with_predicted_bgc)
                        all_cds.append((header, sequence, overlaps_with_overlapping))
        return all_cds

    def _split_cds_by_overlap(self, all_cds):
        """
        Split the CDS features into two lists: those that overlap with features overlapping 'predicted_BGC_type'
        and those that don't.
        """
        overlapping_cds = [f"{header}\n{sequence}" for header, sequence, overlaps_with_overlapping in all_cds if overlaps_with_overlapping]
        non_overlapping_cds = [f"{header}\n{sequence}" for header, sequence, overlaps_with_overlapping in all_cds if not overlaps_with_overlapping]
        return overlapping_cds, non_overlapping_cds

    def _select_cds(self, overlapping_cds, non_overlapping_cds, num_cds):
        """
        Select up to `num_cds` CDS sequences, including all overlapping ones and filling the rest
        from the middle of the non-overlapping list.
        """
        if len(overlapping_cds) + len(non_overlapping_cds) <= num_cds:
            return overlapping_cds + non_overlapping_cds
        remaining_slots = num_cds - len(overlapping_cds)
        total_non_overlapping = len(non_overlapping_cds)

        if remaining_slots > 0:
            # Prioritize middle portion of non-overlapping CDS
            mid_start = total_non_overlapping // 4
            mid_end = total_non_overlapping - (total_non_overlapping // 4)
            middle_cds = non_overlapping_cds[mid_start:mid_end]
            selected_cds = overlapping_cds + middle_cds[:remaining_slots]
        else:
            # If no remaining slots, only use overlapping CDS
            selected_cds = overlapping_cds

        return selected_cds

    def _write_to_fasta(self, sequences, fasta_output):
        """
        Write the selected CDS sequences to a FASTA file.
        """
        if sequences:
            with open(fasta_output, "w") as fasta_file:
                fasta_file.write("\n".join(sequences))
            print(f"Up to {len(sequences)} CDS sequences written to {fasta_output}")


class CDHitReducer:
    def __init__(self, input_dir, output_fasta="combined.fasta", identity=0.97):
        """
        Initialize the CDHitReducer.
        
        Args:
            input_dir (str): Directory containing GenBank files.
            output_fasta (str): Path for the combined FASTA file.
            identity (float): Sequence identity threshold for CD-HIT (default 0.97).
        """
        self.input_dir = input_dir
        self.output_fasta = output_fasta
        self.identity = identity
        self.nr_fasta = os.path.splitext(output_fasta)[0] + "_nr.fasta"
        self.clstr_file = os.path.splitext(output_fasta)[0] + "_nr.fasta.clstr"

    def combine_sequences_to_fasta(self):
        """
        Combine all sequences from GenBank files into a single FASTA file.
        """
        seq_mapping = {}  # Mapping from sequence ID to (GenBank file, SeqRecord)
        with open(self.output_fasta, "w") as fasta_handle:
            for gbk_file in os.listdir(self.input_dir):
                if gbk_file.endswith(".gbk"):
                    gbk_path = os.path.join(self.input_dir, gbk_file)
                    for record in SeqIO.parse(gbk_path, "genbank"):
                        SeqIO.write(SeqRecord(record.seq, id=record.id, description=""), fasta_handle, "fasta")
                        seq_mapping[record.id] = (gbk_file, record)  # Store the sequence with its original GenBank file
                    print(f"Sequences from {gbk_file} added to {self.output_fasta}")
        return seq_mapping

    def run_cd_hit(self):
        """
        Run CD-HIT on the combined FASTA file to remove redundant sequences.
        """
        try:
            # Run CD-HIT
            subprocess.run(['cd-hit', '-i', self.output_fasta, '-o', self.nr_fasta, '-c', str(self.identity)], check=True)
            print(f"Non-redundant sequences saved to {self.nr_fasta}")

            # Check if the .clstr file was created
            if not os.path.exists(self.clstr_file):
                raise FileNotFoundError(f"CD-HIT did not generate the expected cluster file: {self.clstr_file}")

        except subprocess.CalledProcessError as e:
            print(f"Error running CD-HIT: {e}")
        except FileNotFoundError as fnf_error:
            print(fnf_error)

    def parse_cdhit_clusters(self):
        """
        Parse the .clstr file from CD-HIT to identify redundant sequences.
        Returns a dictionary mapping clusters to sequence IDs, with cluster 0 being the non-redundant ones.
        """
        clusters = {}
        cluster_id = None

        if not os.path.exists(self.clstr_file):
            raise FileNotFoundError(f"The CD-HIT cluster file '{self.clstr_file}' was not found. Please ensure CD-HIT ran successfully.")

        with open(self.clstr_file, 'r') as clstr:
            for line in clstr:
                if line.startswith('>Cluster'):
                    cluster_id = int(line.split()[1])  # Get the cluster number
                    clusters[cluster_id] = []
                elif line.strip():  # Process sequence lines
                    seq_id = line.split('>')[1].split('...')[0]  # Extract the sequence ID
                    clusters[cluster_id].append(seq_id)
        return clusters

    def map_and_remove_redundant_sequences(self, seq_mapping, clusters):
        """
        Map non-redundant sequences back to their original GenBank files and remove redundant sequences.
        Always keep the first sequence in each CD-HIT cluster.
        
        Args:
            seq_mapping (dict): Mapping from sequence ID to (GenBank file, SeqRecord).
            clusters (dict): CD-HIT cluster dictionary.
        """
        # Keep a record of sequences to retain (non-redundant)
        sequences_to_keep = set()

        # Always keep the first sequence in each cluster
        for cluster_id, seq_ids in clusters.items():
            sequences_to_keep.add(seq_ids[0])  # Keep the first sequence in each cluster
            print(f"Keeping {seq_ids[0]} from Cluster {cluster_id}")

        # Remove redundant sequences from each GenBank file
        for gbk_file in os.listdir(self.input_dir):
            if gbk_file.endswith(".gbk"):
                gbk_path = os.path.join(self.input_dir, gbk_file)
                original_records = list(SeqIO.parse(gbk_path, "genbank"))
                new_records = []

                # Keep only non-redundant sequences
                for record in original_records:
                    if record.id in sequences_to_keep:
                        new_records.append(record)

                # Overwrite the original GenBank file with non-redundant sequences or delete it if empty
                if new_records:
                    with open(gbk_path, "w") as gbk_handle:
                        SeqIO.write(new_records, gbk_handle, "genbank")
                    print(f"Updated {gbk_file} with {len(new_records)} non-redundant sequences.")
                else:
                    # If there are no non-redundant sequences, remove the GenBank file
                    os.remove(gbk_path)
                    print(f"Removed {gbk_file} as it contained no non-redundant sequences.")


    def remove_redundant_sequences(self):
        """
        Execute the full process of combining sequences, running CD-HIT, and removing redundant sequences.
        """
        # Step 1: Combine sequences into a single FASTA file
        seq_mapping = self.combine_sequences_to_fasta()
        
        # Step 2: Run CD-HIT to reduce sequence redundancy
        self.run_cd_hit()

        # Step 3: Parse the CD-HIT output (.clstr file) to find redundant sequences
        clusters = self.parse_cdhit_clusters()

        # Step 4: Map and remove redundant sequences from GenBank files
        self.map_and_remove_redundant_sequences(seq_mapping, clusters)



class CBLASTProcessor:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file

    def run_cblaster(self):
        """
        Run cblaster search using the input FASTA file.
        """
        try:
            subprocess.run(['cblaster', 'search', '-qf', self.fasta_file, '-s', f"{os.path.splitext(self.fasta_file)[0]}_session.json"], check=True) #"--rid", "M1004R3B016",, 
            print("cblaster search completed.")
        except subprocess.CalledProcessError as e:
            print(f"Error running cblaster search: {e}")
            raise

    def extract_clusters(self, output_dir, max_clusters = 50, upstream=3000, downstream=3000):
        """
        Extract gene clusters from the cblaster session.
        """
        try:
            extractor = NCBICDSExtractor(session_file=f"{os.path.splitext(self.fasta_file)[0]}_session.json", upstream=upstream, downstream=downstream)
            extractor.extract_clusters(output_dir=output_dir, prefix="cluster_", max_clusters=max_clusters)
        except subprocess.CalledProcessError as e:
            print(f"Error extracting clusters: {e}")
            raise


class ORFFinderDirectoryProcessor:
    def __init__(self, input_dir):
        self.input_dir = input_dir

    def run_orf_finder_on_directory(self):
        """
        Processes all GenBank files in a directory, finds ORFs, and returns the ORFs.
        """
        combined_orfs = []
        for gbk_file in os.listdir(self.input_dir):
            if gbk_file.endswith(".gbk"):
                gbk_path = os.path.join(self.input_dir, gbk_file)
                print(f"Running ORF finder on {gbk_path}...")
                orf_finder = ORFFinder(gbk_path)
                orfs, updated_records = orf_finder.process_gbk_file()
                combined_orfs.extend(orfs)
                self.save_updated_records(updated_records)
                print(f"ORF finding completed on {gbk_file}")
        return combined_orfs
    
    def save_updated_records(self, updated_records):
        """
        Save updated GenBank records to files using the locus as the filename.
        
        Args:
            updated_records (list): List of updated SeqRecord objects.
            output_dir (str): Directory to save the updated GenBank files.
        """
        # Ensure the output directory exists
        if not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)

        # Save each updated record as a GenBank file named after its locus (record.id)
        for record in updated_records:
            output_file = os.path.join(self.input_dir, f"{record.id}.gbk")
            with open(output_file, "w") as gbk_handle:
                SeqIO.write(record, gbk_handle, "genbank")
            print(f"Updated record saved to {output_file}")


def main():
    if len(sys.argv) != 2:
        print("Usage: script.py <input_gbk_file>")
        sys.exit(1)

    input_gbk_file = sys.argv[1]
    input_root = os.path.splitext(input_gbk_file)[0]
    output_fasta = f"{input_root}_cds_sequences.fasta"
    output_protein_fasta = f"{input_root}_orfs_protein.fasta"

    if not os.path.isfile(input_gbk_file):
        print(f"Error: {input_gbk_file} is not a valid file")
        sys.exit(1)

    # Step 1: Extract CDS from the input GenBank file and write to a FASTA file
    cds_extractor = CDSExtractor(input_gbk_file)
    cds_extractor.extract_cds(output_fasta)

    # Step 2: Run cblaster using the generated FASTA file
    cblaster_processor = CBLASTProcessor(output_fasta)
    cblaster_processor.run_cblaster()

    # Step 3: Extract gene clusters from the cblaster session file (assumed to be session.json)
    output_directory = input_root
    cblaster_processor.extract_clusters(output_directory, upstream=50000, downstream=50000)

    # Step 4: Remove redundant sequences in-place in the output directory
    cd_hit_reducer = CDHitReducer(output_directory, identity=0.97)
    cd_hit_reducer.remove_redundant_sequences()

    # Step 5: Run the ORF finder on the extracted gene clusters
    orf_finder_processor = ORFFinderDirectoryProcessor(output_directory)
    orfs = orf_finder_processor.run_orf_finder_on_directory()

    # Step 6: Process the found ORFs and generate CSV, FASTA, and sequence logos
    processor = ORFProcessor(orfs)
    processor.save_orfs_to_csv(f"{input_root}_orfs_output.csv")
    processor.save_orfs_to_fasta(f"{input_root}_orfs.fasta")

    # Step 7: Process the clusters and generate sequence logos
    processor.process_clusters_and_logos(input_root)

if __name__ == "__main__":
    main()

