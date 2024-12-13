#!/usr/bin/env python3
import argparse
import itertools
import pysam
import numpy as np
import tqdm.auto as tqdm
import logging

def parse_args():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge haplotype alignments from multiple sources')
    parser.add_argument('--alignment_source', required=True, type=str, action='append', help='Alignment sources (specify multiple times)')
    parser.add_argument('--output', required=True, type=str,help='The output realigned file')
    parser.add_argument('--min_mapping_quality', type=int, default=60, help='low mapping quality threshold')
    args = parser.parse_args()
    args.alignment_sources = args.alignment_source
    return args

class MappingKey:
     # keep read name, chromosome, read start, read end, mapping quality and if the alignment is supplementary
    def __init__(self, rec:pysam.AlignedSegment):
        self.name = rec.query_name
        self.chrom = rec.reference_name
        self.start = rec.query_alignment_start
        self.end = rec.query_alignment_end
        self.mapq = rec.mapping_quality
        self.is_supplementary = rec.is_supplementary
        self.sa_tag = None
        if rec.has_tag('SA'):
            self.sa_tag = rec.get_tag('SA')

    def poor_sa(self, mq_threshold) -> bool:
        if self.sa_tag is not None and 'decoy' in self.sa_tag:
            return True
        if self.is_supplementary and self.mapq < mq_threshold:
            return True
        return False 
    
def _poor_mq_has_sa(mappings: list, mq_threshold: int) -> bool:
    '''Check if there is a poor mapping quality read with supplementary alignment

    Parameters
    ----------
    mappings : list
        List of MappingKey objects

    Returns
    -------
    bool
        True if there is a poor mapping quality read with supplementary alignment, False otherwise
    '''
    has_sa = False
    for mapping in mappings:
        if mapping.is_supplementary:
            has_sa = True
            break
    poor_mq = False
    if has_sa:
        for mapping in mappings:
            if not mapping.is_supplementary and mapping.mapq < mq_threshold:
                poor_mq = True
                break
            
    return poor_mq

def _q_span(mappings: list) -> int:
    '''Calculate the query span of the mappings

    Parameters
    ----------
    mappings : list
        List of MappingKey objects

    Returns
    -------
    int
        Query span
    '''
    return max([x.end for x in mappings]) - min([x.start for x in mappings])

def compare_read_mappings(first: list, second: list, mq_threshold: int) -> int:
    '''Compare two lists of MappingKey objects

    Parameters
    ----------
    first : list
        List of MappingKey objects
    second : list
        List of MappingKey objects
    mq_threshold : int
        Mapping quality threshold

    Returns
    -------
    int
        1 if first > second, -1 if second > first, 0 if equal
    '''
    first_poor_mq = _poor_mq_has_sa(first, mq_threshold) or np.any([x.poor_sa(mq_threshold) for x in first])
    second_poor_mq = _poor_mq_has_sa(second, mq_threshold) or np.any([x.poor_sa(mq_threshold) for x in second])
    if not first_poor_mq and second_poor_mq:
        return 1
    if first_poor_mq and not second_poor_mq:
        return -1   
    if first_poor_mq and second_poor_mq:
        return 0
    first_q_span = _q_span(first)
    second_q_span = _q_span(second)
    if first_q_span > second_q_span:
        return 1
    if first_q_span < second_q_span:
        return -1
    return 0
    
def find_best_mapping_index(mappings: list, mq_threshold: int) -> int:
    '''Find the best mapping index

    Parameters
    ----------
    mappings : list
        List of MappingKey objects
    mq_threshold : int
        Mapping quality threshold

    Returns
    -------
    int
        Best mapping index
    '''
    best_index = 0
    for i in range(1, len(mappings)):
        comp = compare_read_mappings(mappings[best_index], mappings[i], mq_threshold)
        if comp == -1:
            best_index = i
    return best_index

def run():
    args = parse_args()
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(__name__)    
    logger.info(f"Merging alignments: start")
    logger.info(f"Alignment sources: {args.alignment_sources}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Min mapping quality: {args.min_mapping_quality}")
    logger.info(f"Reading alignment files and creating mapping keys:start")
    mappings = [sorted([ MappingKey(rec) for rec in pysam.AlignmentFile(x, "rb") ], key=lambda x: x.name) for x in args.alignment_sources]
    logger.info(f"Reading alignment files and creating mapping keys:done")
    mappings = [ {k: list(v) for k, v in itertools.groupby(x, key=lambda x: x.name)} for x in mappings]
    best_mappings = {}
    logger.info(f"Choosing best alignment per read: start")    
    for k in mappings[0].keys():
        best_mappings[k] = find_best_mapping_index([x[k] for x in mappings], args.min_mapping_quality)
    logger.info(f"Choosing best alignment per read: end")    


    # Create a new BAM file for output
    counters = np.zeros(len(args.alignment_sources))
    logger.info("Reading from the sources and writing into the output:start")
    with pysam.AlignmentFile(args.output, "wb", template=pysam.AlignmentFile(args.alignment_sources[0])) as output:
        for i, r in enumerate(args.alignment_sources):
            logger.info(f"Reading alignments from {r}")
            input = pysam.AlignmentFile(r, "rb")
            for rec in tqdm.tqdm(input):
                if best_mappings[rec.query_name] == i:
                    output.write(rec)
                    counters[i] += 1
            logger.info(f"Wrote {counters[i]} alignments from {r}")
    output.close()

if __name__ == "__main__":
    run()                       