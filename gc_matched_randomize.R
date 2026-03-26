#' GC-matched Randomization for Genomic Regions
#'
#' This function generates a set of random genomic regions that match the GC content 
#' distribution of a target set of regions. It allows for excluding blacklist regions 
#' and ensuring no overlap between randomized regions.
#' 
#' contact: 24111510029@m.fudan.edu.cn  2026.03
#'
#' @param A GRanges object. The original target regions to be randomized.
#' @param genome BSgenome object. The reference genome for GC calculation and sampling.
#' @param mask GRanges object (optional). Regions to exclude from sampling (e.g., blacklists).
#' @param non.overlapping Logical. If TRUE, randomized regions will not overlap original regions A.
#' @param per.chromosome Logical. If TRUE, randomizes regions within the same chromosome as original.
#' @param random_overlap Logical. If FALSE, ensures generated random regions do not overlap each other. Default is TRUE.
#' @param gc.tol Numeric. Initial GC tolerance (absolute difference). Default is 0.05.
#' @param max.iter Integer. Maximum iterations per region before widening GC tolerance. Default is 500.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A sorted GRanges object of randomized regions.
#'
#' @export
#' @import GenomicRanges
#' @import Biostrings
#' @import IRanges

library(GenomicRanges)
library(Biostrings)

calc_gc <- function(gr, genome) {
    seqs <- getSeq(genome, gr)
    gc_matrix <- letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE)
    return(as.numeric(rowSums(gc_matrix)))
}

.sample_one_candidate <- function(avail_wide, int_width, chr) {
    frag_weights <- width(avail_wide) - int_width + 1
    chosen_frag  <- sample(length(avail_wide), 1, prob = frag_weights / sum(frag_weights))
    frag         <- avail_wide[chosen_frag]
    max_start    <- start(frag) + width(frag) - int_width
    rnd_start    <- sample(start(frag):max_start, 1)
    
    GRanges(seqnames = chr, ranges = IRanges(rnd_start, rnd_start + int_width - 1))
}

gc_matched_randomize <- function(A,
                                 genome,
                                 mask            = NULL,
                                 non.overlapping = FALSE,
                                 per.chromosome  = FALSE,
                                 random_overlap  = TRUE,
                                 gc.tol          = 0.05,
                                 max.iter        = 500,
                                 verbose         = FALSE) {
    
    # 1. Prepare available pool
    all_chr_lengths <- seqlengths(genome)
    all_chr_lengths <- all_chr_lengths[!is.na(all_chr_lengths)]
    
    chrs_to_build <- if (per.chromosome) {
        unique(as.character(seqnames(A)))
    } else {
        names(all_chr_lengths)
    }
    
    available <- lapply(chrs_to_build, function(chr) {
        chr_len  <- all_chr_lengths[chr]
        avail_gr <- GRanges(seqnames = chr, ranges = IRanges(1, chr_len))
        excl_list <- list()
        if (non.overlapping) {
            a_on_chr <- A[seqnames(A) == chr]
            if (length(a_on_chr) > 0) excl_list <- c(excl_list, list(a_on_chr))
        }
        if (!is.null(mask)) {
            m_on_chr <- mask[seqnames(mask) == chr]
            if (length(m_on_chr) > 0) excl_list <- c(excl_list, list(m_on_chr))
        }
        if (length(excl_list) > 0) {
            excl_merged <- do.call(c, excl_list)
            avail_gr    <- setdiff(avail_gr, excl_merged)
        }
        return(avail_gr)
    })
    names(available) <- chrs_to_build
    
    if (!per.chromosome) {
        chr_avail_len <- sapply(available, function(gr) sum(as.numeric(width(gr))))
        chr_avail_len <- chr_avail_len[chr_avail_len > 0]
    }
    
    gc_orig <- calc_gc(A, genome)
    result_list <- vector("list", length(A))
    failed_idx  <- integer(0)
    current_available <- available
    
    # 2. Iterative Sampling
    for (i in seq_along(A)) {
        int_width   <- width(A[i])
        target_gc   <- gc_orig[i]
        tol_current <- gc.tol
        matched     <- FALSE
        
        for (iter in seq_len(max.iter)) {
            if (per.chromosome) {
                target_chr <- as.character(seqnames(A[i]))
                avail_chr  <- current_available[[target_chr]]
            } else {
                if (!random_overlap) {
                    chr_avail_len_curr <- sapply(current_available, function(gr) sum(as.numeric(width(gr))))
                    valid_chrs <- names(chr_avail_len_curr[chr_avail_len_curr >= int_width])
                    if (length(valid_chrs) == 0) {
                        target_chr <- sample(names(chr_avail_len), 1) # fallback
                        avail_chr  <- available[[target_chr]]
                    } else {
                        target_chr <- sample(valid_chrs, 1, prob = chr_avail_len_curr[valid_chrs])
                        avail_chr  <- current_available[[target_chr]]
                    }
                } else {
                    target_chr <- sample(names(chr_avail_len), 1, prob = chr_avail_len)
                    avail_chr  <- available[[target_chr]]
                }
            }
            
            avail_wide <- avail_chr[width(avail_chr) >= int_width]
            if (length(avail_wide) == 0) {
                if (iter == max.iter && verbose) message("No space on ", target_chr)
                next
            }
            
            candidate <- .sample_one_candidate(avail_wide, int_width, target_chr)
            cand_gc <- calc_gc(candidate, genome)
            
            if (abs(cand_gc - target_gc) <= tol_current) {
                result_list[[i]] <- candidate
                matched <- TRUE
                break
            }
            
            if (iter %% 100 == 0) tol_current <- tol_current + 0.01
        }
        
        if (!matched) {
            failed_idx <- c(failed_idx, i)
            result_list[[i]] <- candidate
        }
        
        if (!random_overlap && matched) {
            target_chr <- as.character(seqnames(result_list[[i]]))
            current_available[[target_chr]] <- setdiff(current_available[[target_chr]], result_list[[i]])
        }
    }
    
    valid     <- result_list[!sapply(result_list, is.null)]
    result_gr <- do.call(c, valid)
    seqlevels(result_gr, pruning.mode = "coarse") <- seqlevels(genome)
    seqlengths(result_gr) <- seqlengths(genome)
    
    return(sort(result_gr))
}
