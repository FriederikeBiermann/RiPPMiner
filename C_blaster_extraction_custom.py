#!/usr/bin/env python3

import logging
from pathlib import Path
from Bio import SeqIO, Entrez
from cblaster.classes import Session
from cblaster.extract import organism_matches, parse_organisms, parse_scaffolds

# Configure logging
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

Entrez.email = "friederike@biermann-erfurt.de"  # Set your email to use NCBI Entrez API

class NCBICDSExtractor:
    def __init__(self, session_file, upstream=3000, downstream=3000):
        """
        Initialize the extractor with the session file and upstream/downstream flanking regions.
        Args:
            session_file (str): Path to the cblaster session file.
            upstream (int): Size of the upstream region to extract.
            downstream (int): Size of the downstream region to extract.
        """
        self.session_file = session_file
        self.upstream = upstream
        self.downstream = downstream

    def extract_clusters(self, output_dir, prefix="", cluster_numbers=None, score_threshold=None, organisms=None, scaffolds=None, format_="genbank", max_clusters=50):
        """
        Extract clusters from the cblaster session file and write them to GenBank and FASTA files.
        """
        LOG.info(f"Loading session from: {self.session_file}")
        with open(self.session_file) as fp:
            session = Session.from_json(fp)

        LOG.info("Extracting clusters that match the filters")
        cluster_hierarchy = self.get_sorted_cluster_hierarchies(session, cluster_numbers, score_threshold, organisms, scaffolds, max_clusters)

        LOG.info(f"Extracted {len(cluster_hierarchy)} clusters.")
        if len(cluster_hierarchy) == 0:
            LOG.info("No clusters meet the filtering criteria. Exiting...")
            raise SystemExit(0)

        LOG.info("Writing GenBank and FASTA files")
        self.create_genbanks_from_clusters(cluster_hierarchy, output_dir, prefix, format_)
        LOG.info(f"Clusters written to {output_dir}")

    def fetch_gbk_file(self, accession, start, end):
        """
        Fetch the GenBank file from NCBI for the given accession and range, including flanking regions.
        Args:
            accession (str): GenBank accession number.
            start (int): Start coordinate for the sequence (1-based).
            end (int): End coordinate for the sequence (1-based).
        Returns:
            SeqRecord: GenBank record for the specified range.
        """
        # Adjust start and end to include upstream and downstream regions
        adjusted_start = max(1, start - self.upstream)  # Ensure we don't go below 1 (NCBI uses 1-based indexing)
        adjusted_end = end + self.downstream

        LOG.info(f"Fetching GenBank file for {accession} from {adjusted_start} to {adjusted_end}")
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text", seq_start=adjusted_start, seq_stop=adjusted_end)
        return SeqIO.read(handle, "genbank")

    def create_genbanks_from_clusters(self, cluster_hierarchy, output_dir, prefix, format_):
        """
        Create GenBank and FASTA files for each selected cluster.
        Args:
            cluster_hierarchy (list): Filtered list of clusters.
            output_dir (str): Directory to store the output files.
            prefix (str): Prefix for the output file names.
            format_ (str): File format for the GenBank output (e.g., 'genbank').
        """
        output_dir = Path(output_dir)
        if not output_dir.is_dir():
            output_dir.mkdir()

        for cluster, scaffold, organism in cluster_hierarchy:
            # Fetch the GenBank file for this cluster, including flanking regions
            record = self.fetch_gbk_file(scaffold.accession, cluster.intermediate_start + 1, cluster.intermediate_end)

            # Generate output file paths
            output_gbk = output_dir / f"{prefix}cluster{cluster.number}.gbk"
            output_fasta = output_dir / f"{prefix}cluster{cluster.number}.fasta"

            # Write GenBank output
            with open(output_gbk, "w") as gbk_handle:
                SeqIO.write(record, gbk_handle, format_)
                LOG.info(f"Created {output_gbk} for cluster {cluster.number}")

            # Write FASTA output
            self.write_fasta(record, output_fasta)

    def write_fasta(self, record, output_fasta):
        """
        Write the cluster sequence to a FASTA file.
        Args:
            record (SeqRecord): GenBank record containing the sequence.
            output_fasta (str): Path to the output FASTA file.
        """
        with open(output_fasta, "w") as fasta_handle:
            fasta_handle.write(f">{record.id} {record.description}\n{record.seq}\n")
        LOG.info(f"Created {output_fasta}")

    def get_sorted_cluster_hierarchies(self, 
        session,
        cluster_numbers=None,
        score_threshold=None,
        organisms=None,
        scaffolds=None,
        max_clusters=50,
    ):
        """
        Filter out selected clusters with their associated scaffold and organism.
        Args:
            session (Session): A session object.
            cluster_numbers (list): Numbers of clusters or a range of numbers (e.g., '1-5').
            score_threshold (float): Minimum score a cluster needs to have in order to be included.
            organisms (list): Regex patterns for organisms of which all clusters need to be extracted.
            scaffolds (list): Names of scaffolds of which all clusters need to be extracted.
            max_clusters (int): The maximum number of clusters extracted regardless of filters. Set to None to extract all clusters.
        Returns:
            List of tuples of the form (Cluster object, Scaffold object, organism_name).
        """
        # No filter options return all clusters
        selected_clusters = set()

        # Prepare the filters defined by the user
        if cluster_numbers:
            cluster_numbers = parse_numbers(cluster_numbers)
        if organisms:
            organisms = parse_organisms(organisms)
        if scaffolds:
            scaffolds = parse_scaffolds(scaffolds)

        # Apply the filters to the clusters
        for organism in session.organisms:
            if organisms and not organism_matches(organism.name, organisms):
                continue
            for scaffold in organism.scaffolds.values():
                if scaffolds and scaffold.accession not in scaffolds:
                    continue
                for cluster in scaffold.clusters:
                    if cluster_numbers and cluster.number not in cluster_numbers:
                        continue
                    if score_threshold and cluster.score < score_threshold:
                        continue
                    selected_clusters.add((cluster, scaffold, organism.full_name))

        # Sort clusters by score and return the top clusters
        return sorted(
            selected_clusters,
            key=lambda x: (x[0].score, -x[0].start, -x[0].end, x[1].accession),
            reverse=True,
        )[:max_clusters]

def cluster_in_range(start, end, cluster):
    """
    Checks if a cluster is within a given range.
    Args:
        start (int): Start of the range.
        end (int): End of the range.
        cluster (Cluster): cblaster cluster object.
    Returns:
        boolean: True if the cluster is within the start and end values.
    """
    if start is end is None:
        return True
    return start <= cluster.start and end >= cluster.end

def parse_numbers(cluster_numbers):
    """
    Parses cluster numbers from user input.
    Args:
        cluster_numbers (list): A list of numbers or ranges of numbers.
    Returns:
        list: A list of integer numbers.
    """
    chosen_cluster_numbers = set()
    for number in cluster_numbers:
        try:
            if "-" in number:
                start, stop = number.split("-")
                chosen_cluster_numbers.update(range(int(start), int(stop) + 1))
            else:
                chosen_cluster_numbers.add(int(number))
        except ValueError:
            LOG.warning(f"Cannot extract cluster '{number}': number is not a valid integer")
    return chosen_cluster_numbers

