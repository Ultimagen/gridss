library(StructuralVariantAnnotation, quietly=TRUE)
library(VariantAnnotation, quietly=TRUE)
library(GenomicRanges, quietly=TRUE)
library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE)
library(parallel, quietly = TRUE)
library(pbapply, quietly = TRUE)
library(stringr, quietly=TRUE)

library(argparser, quietly=TRUE)

# Create a parser object
parser <- ArgumentParser(description = 'Convert VCF format')
argp = arg_parser("convert vcf format")
argp = add_argument(argp, "--ref", default="", help="Reference genome to use. Must be a valid installed BSgenome package")
argp = add_argument(argp, "--input", help="GRIDSS VCF")
argp = add_argument(argp, "--output", help="Output VCF")
argp = add_argument(argp, "--n_jobs", type="integer", default=-1, help="Number of parallel jobs")
argv = parse_args(argp)


# Define the path to your VCF file
# gs://cromwell-backend-ultima-data-307918/cromwell-execution/SVPipeline/9aaff528-f4e7-439e-b5b5-47ee747e2515/call-GermlineLinkVariants/NA24385_linked.vcf.bgz
# vcf_file <-"/Users/mayalevy/Downloads/gridss/NA24385_linked.vcf.bgz"
#vcf_file <- "/Users/mayalevy/Downloads/gridss/401882-CL10366-Z0082-CTCTGCTGTGCAATGAT_chr1_linked_orig.vcf.bgz"
#vcf_file <- "/Users/mayalevy/Downloads/gridss/diploidSV.vcf.gz"
# gsutil cp modified_vcf_wgs.vcf.bgz gs://ultimagen-users-data/maya/deepvariant/gridss/

# Read the VCF file
vcf <- readVcf(argv$input, "")

# Remove all SVTYPE values from the INFO field
info(vcf)$SVTYPE <- NULL

# Extract breakpoint ranges
vcf_bp <- breakpointRanges(vcf, nominalPosition=FALSE, suffix="_bp", inferMissingBreakends = FALSE)

# Determine the simple event types for each breakpoint
vcf_bp$simpleEvent <- simpleEventType(vcf_bp)

# Initialize a simpleEvent vector with NA to match the length of the original VCF
simpleEvent <- rep(NA, length(vcf))
end_positions <- rep(NA, length(vcf))
svlens <- rep(NA, length(vcf))

# Find the matching indices between the original VCF and vcf_bp
matching_indices <- match(names(vcf), names(vcf_bp))

# Assign the simpleEvent information to the matching positions
simpleEvent[!is.na(matching_indices)] <- vcf_bp$simpleEvent[na.omit(matching_indices)]
svlens[!is.na(matching_indices)] <- vcf_bp$svLen[na.omit(matching_indices)]
end_positions[!is.na(matching_indices)] <- start(vcf)[!is.na(matching_indices)] + vcf_bp$svLen[na.omit(matching_indices)]

# Create a new header line for the END field
end_header <- DataFrame(
  Number = "1",
  Type = "Integer",
  Description = "End position of the variant",
  row.names = "END"
)

# Create a new header line for the SVLEN field
svlen_header <- DataFrame(
  Number = "1",
  Type = "Integer",
  Description = "Length of variant",
  row.names = "SVLEN"
)

# Add the new header line to the VCF header
vcf_header <- header(vcf)
info(vcf_header) <- rbind(info(vcf_header), end_header)
info(vcf_header) <- rbind(info(vcf_header), svlen_header)
header(vcf) <- vcf_header

# Add the simpleEvent info to the original VCF object
info(vcf)$SVTYPE <- simpleEvent
info(vcf)$END <- end_positions
info(vcf)$SVLEN <- svlens

# remove variants with svLen smaller than 30
vcf = vcf[which(is.na(svlens) | abs(svlens) >= 50)]

### Remove MATE variant - keep only one variant  only in case of DEL or INS

# Identify indices of deletions (DEL) in the VCF
del_ins_indices <- which(!is.na(info(vcf)$SVTYPE) & ((info(vcf)$SVTYPE == "DEL") | ((info(vcf)$SVTYPE == "INS"))))

# Extract MATEID for deletions
mate_ids <- info(vcf)$MATEID[del_ins_indices]

# Initialize a logical vector to keep track of variants to remove
remove_indices <- rep(FALSE, length(vcf))

# Create a set to keep track of processed MATEIDs
processed_mate_ids <- character()
count=0
# Iterate over deletion indices to remove only one of each mate pair
for (i in seq_along(del_ins_indices)) {
  current_index <- del_ins_indices[i]
  mate_id <- info(vcf)$MATEID[[current_index]]
  var_id <- names(vcf)[current_index]
  
  if (length(mate_id) > 0 && !is.na(mate_id) && var_id %in% processed_mate_ids) {
    # If the mate ID has already been processed, mark the current variant for removal
    remove_indices[current_index] <- TRUE
  } else {
    # Mark the mate ID as processed
    processed_mate_ids <- c(processed_mate_ids, mate_id)
    info(vcf)$MATEID[[current_index]]  <- NA_character_
  }
}

# Subset the VCF to remove the marked variants
vcf <- vcf[!remove_indices]

# Change SVTYPE from INV to BND for all INV variants
inv_indices <- which(!is.na(info(vcf)$SVTYPE) & info(vcf)$SVTYPE == "INV")
info(vcf)$SVTYPE[inv_indices] <- "BND"

# Collect start and end positions for DEL variants
del_indices <- which(!is.na(info(vcf)$SVTYPE) & info(vcf)$SVTYPE == "DEL")
del_lengths <- abs(info(vcf)$SVLEN[del_indices])

# Separate DEL variants based on length
short_del_indices <- del_indices[del_lengths <= 329]
short_del_lengths <- del_lengths[del_lengths <= 329]
short_del_starts <- start(vcf)[short_del_indices]
short_del_ends <- short_del_starts + short_del_lengths - 1

# Fetch sequences for short DEL variants in one call
short_del_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, names = seqnames(vcf)[short_del_indices], start = short_del_starts, end = short_del_ends)

# Collect indices for INS variants
ins_indices <- which(!is.na(info(vcf)$SVTYPE) & info(vcf)$SVTYPE == "INS")

# Function to process each variant
process_variant <- function(i, vcf, short_del_indices, short_del_seqs, short_del_lengths, ins_indices) {
  result <- list()
  
  if (i %in% short_del_indices) {
    del_idx <- which(short_del_indices == i)
    deletion_length <- short_del_lengths[del_idx]
    
    # Update the REF and ALT fields for short DEL variants
    full_ref_seq <- short_del_seqs[del_idx]
    result$alt <- as.character(ref(vcf)[i])
    result$ref <- DNAString(as.character(full_ref_seq))
    result$end <- start(vcf)[i] + length(result$ref) - 1
    
  } else if (i %in% del_indices && !(i %in% short_del_indices)) {
    # Handle long DEL variants
    result$alt <- CharacterList("<DEL>")
    result$ref <- DNAString(as.character(ref(vcf)[i]))
    result$end <- start(vcf)[i] + length(result$ref) - 1
  } else if (i %in% ins_indices) {
    # Update the REF and ALT fields for INS variants
    mateid_field <- info(vcf)$MATEID[i]
    if (!is.null(mateid_field) && length(mateid_field) > 0 && !is.na(mateid_field[[1]])) {
      result$alt <- CharacterList("<INS>")
    } else {
      result$alt <- gsub("\\]chrX:[0-9]+\\]", "", as.character(alt(vcf)[i]))
    }
    result$ref <- DNAString(as.character(ref(vcf)[i]))
    result$end <- start(vcf)[i] + length(result$ref) - 1
  } else if (is.null(info(vcf)$SVTYPE[i])) {
    alt_field <- as.character(alt(vcf)[i])
    mateid_field <- info(vcf)$MATEID[i]
    if (!is.null(mateid_field) && length(mateid_field[[1]]) != 0 && (grepl(".*\\[.*\\].*", alt_field) || grepl(".*\\].*\\[.*", alt_field))) {
      # Handle breakends BND
      result$svtype <- "BND"
      alt_length <- nchar(alt_field)
      result$svlen <- alt_length
      result$end <- start(vcf)[i] + alt_length
    }
  }
  
  return(result)
}

# Number of cores to use
if (argv$n_jobs == -1) {
  num_cores <- detectCores() - 1
} else {
  num_cores <- argv$n_jobs
}

# Process variants in parallel with progress bar using pblapply
results <- pblapply(seq_along(vcf), process_variant, vcf = vcf, short_del_indices = short_del_indices, short_del_seqs = short_del_seqs, short_del_lengths = short_del_lengths, ins_indices = ins_indices, cl = num_cores)

# Update the VCF object with results

# Assuming results from pblapply are ready

# Initialize the necessary fields
alt_updates <- vector("list", length(vcf))
ref_updates <- vector("list", length(vcf))
end_updates <- rep(NA, length(vcf))
svtype_updates <- rep(NA_character_, length(vcf))
svlen_updates <- rep(NA, length(vcf))

# Extract updates from results
for (i in seq_along(results)) {
  if (!is.null(results[[i]]$alt)) {
    alt_updates[[i]] <- results[[i]]$alt
  }
  if (!is.null(results[[i]]$ref)) {
    ref_updates[[i]] <- results[[i]]$ref
  }
  if (!is.null(results[[i]]$end)) {
    end_updates[i] <- results[[i]]$end
  }
  if (!is.null(results[[i]]$svtype)) {
    svtype_updates[i] <- results[[i]]$svtype
  }
  if (!is.null(results[[i]]$svlen)) {
    svlen_updates[i] <- results[[i]]$svlen
  }
}


ref_updates <- lapply(ref_updates, function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return(DNAString(""))
  } else {
    return(as(x, "DNAString"))
  }
})

# Apply updates to the VCF object
alt(vcf) <- as(alt_updates, "CharacterList")
ref(vcf) <- DNAStringSet(unlist(ref_updates))
info(vcf)$END <- end_updates
info(vcf)$SVTYPE <- svtype_updates
info(vcf)$SVLEN <- svlen_updates
# for (i in seq_along(vcf)) {
#   if (!is.null(results[[i]]$alt)) {
#     alt(vcf)[i] <- results[[i]]$alt
#   }
#   if (!is.null(results[[i]]$ref)) {
#     ref(vcf)[i] <- results[[i]]$ref
#   }
#   if (!is.null(results[[i]]$end)) {
#     info(vcf)$END[i] <- results[[i]]$end
#   }
#   if (!is.null(results[[i]]$svtype)) {
#     info(vcf)$SVTYPE[i] <- results[[i]]$svtype
#   }
#   if (!is.null(results[[i]]$svlen)) {
#     info(vcf)$SVLEN[i] <- results[[i]]$svlen
#   }
# }





### Modify ALT and REF fields for deletions
# for (i in seq_along(vcf)) {
#   if (!is.null(info(vcf)$SVTYPE[i]) && !is.na(info(vcf)$SVTYPE[i]) && info(vcf)$SVTYPE[i] == "DEL") {
# 
#       # Calculate the sequence length of the deletion
#       deletion_length <- abs(info(vcf)$SVLEN[i])
# 
# 
#       # Update the REF field
#       if (deletion_length <= 329) {
#         # If the deletion length is not larger than 329 bases, use the full sequence
#         full_ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,
#                                names = seqnames(vcf)[i],
#                                start = start(vcf)[i],
#                                end = start(vcf)[i] + deletion_length - 1)
#         alt(vcf)[i] <- ref(vcf)[i]
#         ref(vcf)[i] <- DNAStringSet(full_ref_seq)
#       } else {
#         alt(vcf)[i] <- CharacterList("<DEL>")
#       }
#       info(vcf)$END[i] = start(vcf)[i] + length(ref(vcf)[i]) -1
#   } else if (!is.null(info(vcf)$SVTYPE[i]) && !is.na(info(vcf)$SVTYPE[i]) && info(vcf)$SVTYPE[i] == "INS") {
#     ### Modify ALT and REF fields for insertions
#     # Set ALT to REF seq
#     mateid_field = info(vcf)$MATEID[i]
#     if(!is.null(mateid_field) && length(mateid_field[[1]]) != 0) {
#       alt(vcf)[i] <- CharacterList("<INS>")
#     } else {
#       alt(vcf)[i] <- gsub("\\]chrX:[0-9]+\\]", "", alt(vcf)[i])
#     }
#     
#     
#     info(vcf)$END[i] = start(vcf)[i] + length(ref(vcf)[i]) -1
#   } else if(is.null(info(vcf)$SVTYPE[i])){
#     alt_field <- as.character(alt(vcf)[i])
#     mateid_field = info(vcf)$MATEID[i]
#     if(!is.null(mateid_field) && length(mateid_field[[1]]) != 0 && (grepl(".*\\[.*\\].*", alt_field) || grepl(".*\\].*\\[.*", alt_field))) {
#       ### Handle breakends BND
#       info(vcf)$SVTYPE[i] = "BND"
#       # Get the length of the string
#       alt_length = nchar(alt_field)
#       info(vcf)$SVLEN[i] = alt_length
#       info(vcf)$END[i] = start(vcf)[i] + alt_length
#     }
#   }
# }




# Define the path to save the modified VCF file as compressed VCF
#modified_vcf_file <- "/Users/mayalevy/Downloads/gridss/modified_vcf.vcf.gz"

# Write the modified VCF file as compressed VCF
writeVcf(vcf, filename = argv$output, index = TRUE)


parser = argparse.ArgumentParser(description='Convert vcf format.')
parser.add_argument('--input', required=True, help='The input vcf file')
parser.add_argument('--output', required=True, help='The output vcf file')
parser.add_argument('--reference', required=True, help='The reference genome FASTA file')
parser.add_argument("--n_jobs", help="n_jobs of parallel on contigs", type=int, default=-1)
args = parser.parse_args()