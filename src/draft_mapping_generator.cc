#include "draft_mapping_generator.h"

#include <vector>

#include "alignment.h"

namespace chromap {

void DraftMappingGenerator::GenerateDraftMappings(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  mapping_metadata.SetMinNumErrors(error_threshold_ + 1);
  mapping_metadata.SetNumBestMappings(0);
  mapping_metadata.SetSecondMinNumErrors(error_threshold_ + 1);
  mapping_metadata.SetNumSecondBestMappings(0);

  // Directly obtain the non-split mapping in ideal case and return without
  // running actual verification.
  const bool is_mapping_generated =
      GenerateNonSplitDraftMappingSupportedByAllMinimizers(
          read_batch, read_index, reference, mapping_metadata);
  if (is_mapping_generated) {
    return;
  }

  // Use more sophicated approach to obtain the mapping.
  // Sort the candidates by their count in descending order.
  // TODO: check if this sorting is necessary.
  mapping_metadata.SortCandidates();

  // For split alignments, SIMD cannot be used.
  if (split_alignment_) {
    GenerateDraftMappingsOnOneStrand(kPositive, read_index, read_batch,
                                     reference, mapping_metadata);

    GenerateDraftMappingsOnOneStrand(kNegative, read_index, read_batch,
                                     reference, mapping_metadata);
    return;
  }

  // For non-split alignments, use SIMD when possible.
  if (mapping_metadata.GetNumPositiveCandidates() < (size_t)num_vpu_lanes_) {
    GenerateDraftMappingsOnOneStrand(kPositive, read_index, read_batch,
                                     reference, mapping_metadata);
  } else {
    GenerateDraftMappingsOnOneStrandUsingSIMD(kPositive, read_index, read_batch,
                                              reference, mapping_metadata);
  }

  if (mapping_metadata.GetNumNegativeCandidates() < (size_t)num_vpu_lanes_) {
    GenerateDraftMappingsOnOneStrand(kNegative, read_index, read_batch,
                                     reference, mapping_metadata);
  } else {
    GenerateDraftMappingsOnOneStrandUsingSIMD(kNegative, read_index, read_batch,
                                              reference, mapping_metadata);
  }
}

bool DraftMappingGenerator::IsValidCandidate(uint32_t rid, uint32_t position,
                                             uint32_t read_length,
                                             const SequenceBatch &reference) {
  const uint32_t reference_length = reference.GetSequenceLengthAt(rid);

  if (position < (uint32_t)error_threshold_ || position >= reference_length ||
      position + read_length + (uint32_t)error_threshold_ >= reference_length) {
    return false;
  }

  return true;
}

bool DraftMappingGenerator::
    GenerateNonSplitDraftMappingSupportedByAllMinimizers(
        const SequenceBatch &read_batch, uint32_t read_index,
        const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  if (split_alignment_) {
    return false;
  }

  const bool has_one_candidate = (mapping_metadata.GetNumCandidates() == 1);

  if (!has_one_candidate) {
    return false;
  }

  const std::vector<Candidate> &positive_candidates =
      mapping_metadata.positive_candidates_;
  const std::vector<Candidate> &negative_candidates =
      mapping_metadata.negative_candidates_;

  std::vector<DraftMapping> &positive_mappings =
      mapping_metadata.positive_mappings_;
  std::vector<DraftMapping> &negative_mappings =
      mapping_metadata.negative_mappings_;

  uint32_t num_all_minimizer_candidates = 0;
  uint32_t all_minimizer_candidate_index = 0;
  Strand all_minimizer_candidate_strand = kPositive;

  for (uint32_t i = 0; i < positive_candidates.size(); ++i) {
    if (positive_candidates[i].count == mapping_metadata.GetNumMinimizers()) {
      all_minimizer_candidate_index = i;
      ++num_all_minimizer_candidates;
    }
  }

  for (uint32_t i = 0; i < negative_candidates.size(); ++i) {
    if (negative_candidates[i].count == mapping_metadata.GetNumMinimizers()) {
      all_minimizer_candidate_index = i;
      all_minimizer_candidate_strand = kNegative;
      ++num_all_minimizer_candidates;
    }
  }

  if (num_all_minimizer_candidates != 1) {
    return false;
  }

  mapping_metadata.SetMinNumErrors(0);
  mapping_metadata.SetNumBestMappings(1);
  mapping_metadata.SetNumSecondBestMappings(0);

  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::vector<Candidate> &candidates =
      all_minimizer_candidate_strand == kPositive ? positive_candidates
                                                  : negative_candidates;

  const uint32_t rid =
      candidates[all_minimizer_candidate_index].GetReferenceSequenceIndex();

  uint32_t position = 0;

  if (all_minimizer_candidate_strand == kPositive) {
    position = positive_candidates[all_minimizer_candidate_index]
                   .GetReferenceSequencePosition();
  } else {
    position = negative_candidates[all_minimizer_candidate_index]
                   .GetReferenceSequencePosition() -
               read_length + 1;
  }

  const bool is_valid_candidate =
      IsValidCandidate(rid, position, read_length, reference);
  if (is_valid_candidate) {
    if (all_minimizer_candidate_strand == kPositive) {
      positive_mappings.emplace_back(
          0, positive_candidates[all_minimizer_candidate_index].position +
                 read_length - 1);
    } else {
      negative_mappings.emplace_back(
          0, negative_candidates[all_minimizer_candidate_index].position);
    }
    return true;
  }

  return false;
}

void DraftMappingGenerator::GenerateDraftMappingsOnOneStrandUsingSIMD(
    const Strand candidate_strand, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &reference,
    MappingMetadata &mapping_metadata) {
  const char *read = read_batch.GetSequenceAt(read_index);
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  const std::vector<Candidate> &candidates =
      candidate_strand == kPositive ? mapping_metadata.positive_candidates_
                                    : mapping_metadata.negative_candidates_;
  std::vector<DraftMapping> &mappings =
      candidate_strand == kPositive ? mapping_metadata.positive_mappings_
                                    : mapping_metadata.negative_mappings_;
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  Candidate valid_candidates[num_vpu_lanes_];
  const char *valid_candidate_starts[num_vpu_lanes_];
  uint32_t valid_candidate_index = 0;
  size_t candidate_index = 0;
  uint32_t candidate_count_threshold = 0;

  while (candidate_index < candidates.size()) {
    if (candidates[candidate_index].count < candidate_count_threshold) {
      break;
    }

    uint32_t rid = candidates[candidate_index].GetReferenceSequenceIndex();
    uint32_t position =
        candidates[candidate_index].GetReferenceSequencePosition();

    if (candidate_strand == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      ++candidate_index;
      continue;
    }

    valid_candidates[valid_candidate_index] = candidates[candidate_index];
    valid_candidate_starts[valid_candidate_index] =
        reference.GetSequenceAt(rid) + position - error_threshold_;
    ++valid_candidate_index;
    ++candidate_index;

    if (valid_candidate_index < (uint32_t)num_vpu_lanes_) {
      continue;
    }

    if (num_vpu_lanes_ == 8) {
      int16_t mapping_edit_distances[num_vpu_lanes_];
      int16_t mapping_end_positions[num_vpu_lanes_];
      for (int li = 0; li < num_vpu_lanes_; ++li) {
        mapping_end_positions[li] = read_length - 1;
      }
      if (candidate_strand == kPositive) {
        BandedAlign8PatternsToText(error_threshold_, valid_candidate_starts,
                                   read, read_length, mapping_edit_distances,
                                   mapping_end_positions);
      } else {
        BandedAlign8PatternsToText(
            error_threshold_, valid_candidate_starts, negative_read.data(),
            read_length, mapping_edit_distances, mapping_end_positions);
      }
      for (int mi = 0; mi < num_vpu_lanes_; ++mi) {
        if (mapping_edit_distances[mi] <= error_threshold_) {
          if (mapping_edit_distances[mi] < min_num_errors) {
            second_min_num_errors = min_num_errors;
            num_second_best_mappings = num_best_mappings;
            min_num_errors = mapping_edit_distances[mi];
            num_best_mappings = 1;
          } else if (mapping_edit_distances[mi] == min_num_errors) {
            num_best_mappings++;
          } else if (mapping_edit_distances[mi] == second_min_num_errors) {
            num_second_best_mappings++;
          } else if (mapping_edit_distances[mi] < second_min_num_errors) {
            num_second_best_mappings = 1;
            second_min_num_errors = mapping_edit_distances[mi];
          }
          if (candidate_strand == kPositive) {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position -
                                      error_threshold_ +
                                      mapping_end_positions[mi]);
          } else {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position - read_length +
                                      1 - error_threshold_ +
                                      mapping_end_positions[mi]);
          }
        } else {
          candidate_count_threshold = valid_candidates[mi].count;
        }
      }
    } else if (num_vpu_lanes_ == 4) {
      int32_t mapping_edit_distances[num_vpu_lanes_];
      int32_t mapping_end_positions[num_vpu_lanes_];
      for (int li = 0; li < num_vpu_lanes_; ++li) {
        mapping_end_positions[li] = read_length - 1;
      }
      if (candidate_strand == kPositive) {
        BandedAlign4PatternsToText(error_threshold_, valid_candidate_starts,
                                   read, read_length, mapping_edit_distances,
                                   mapping_end_positions);
      } else {
        BandedAlign4PatternsToText(
            error_threshold_, valid_candidate_starts, negative_read.data(),
            read_length, mapping_edit_distances, mapping_end_positions);
      }
      for (int mi = 0; mi < num_vpu_lanes_; ++mi) {
        if (mapping_edit_distances[mi] <= error_threshold_) {
          if (mapping_edit_distances[mi] < min_num_errors) {
            second_min_num_errors = min_num_errors;
            num_second_best_mappings = num_best_mappings;
            min_num_errors = mapping_edit_distances[mi];
            num_best_mappings = 1;
          } else if (mapping_edit_distances[mi] == min_num_errors) {
            num_best_mappings++;
          } else if (mapping_edit_distances[mi] == second_min_num_errors) {
            num_second_best_mappings++;
          } else if (mapping_edit_distances[mi] < second_min_num_errors) {
            num_second_best_mappings = 1;
            second_min_num_errors = mapping_edit_distances[mi];
          }
          if (candidate_strand == kPositive) {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position -
                                      error_threshold_ +
                                      mapping_end_positions[mi]);
          } else {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position - read_length +
                                      1 - error_threshold_ +
                                      mapping_end_positions[mi]);
          }
        } else {
          candidate_count_threshold = valid_candidates[mi].count;
        }
      }
    }

    valid_candidate_index = 0;
  }

  for (uint32_t ci = 0; ci < valid_candidate_index; ++ci) {
    uint32_t rid = valid_candidates[ci].GetReferenceSequenceIndex();
    uint32_t position = valid_candidates[ci].GetReferenceSequencePosition();
    if (candidate_strand == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      continue;
    }

    int mapping_end_position;
    int num_errors;
    if (candidate_strand == kPositive) {
      num_errors = BandedAlignPatternToText(
          error_threshold_,
          reference.GetSequenceAt(rid) + position - error_threshold_, read,
          read_length, &mapping_end_position);
    } else {
      num_errors = BandedAlignPatternToText(
          error_threshold_,
          reference.GetSequenceAt(rid) + position - error_threshold_,
          negative_read.data(), read_length, &mapping_end_position);
    }
    if (num_errors <= error_threshold_) {
      if (num_errors < min_num_errors) {
        second_min_num_errors = min_num_errors;
        num_second_best_mappings = num_best_mappings;
        min_num_errors = num_errors;
        num_best_mappings = 1;
      } else if (num_errors == min_num_errors) {
        num_best_mappings++;
      } else if (num_errors == second_min_num_errors) {
        num_second_best_mappings++;
      } else if (num_errors < second_min_num_errors) {
        num_second_best_mappings = 1;
        second_min_num_errors = num_errors;
      }
      if (candidate_strand == kPositive) {
        mappings.emplace_back(num_errors, valid_candidates[ci].position -
                                              error_threshold_ +
                                              mapping_end_position);
      } else {
        mappings.emplace_back(num_errors,
                              valid_candidates[ci].position - read_length + 1 -
                                  error_threshold_ + mapping_end_position);
      }
    }
  }
}

void DraftMappingGenerator::GenerateDraftMappingsOnOneStrand(
    const Strand candidate_strand, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &reference,
    MappingMetadata &mapping_metadata) {
  const char *read = read_batch.GetSequenceAt(read_index);
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  const std::vector<Candidate> &candidates =
      candidate_strand == kPositive ? mapping_metadata.positive_candidates_
                                    : mapping_metadata.negative_candidates_;
  std::vector<DraftMapping> &mappings =
      candidate_strand == kPositive ? mapping_metadata.positive_mappings_
                                    : mapping_metadata.negative_mappings_;
  std::vector<int> &split_sites = candidate_strand == kPositive
                                      ? mapping_metadata.positive_split_sites_
                                      : mapping_metadata.negative_split_sites_;
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  uint32_t candidate_count_threshold = 0;

  for (uint32_t ci = 0; ci < candidates.size(); ++ci) {
    if (candidates[ci].count < candidate_count_threshold) {
      break;
    }

    uint32_t rid = candidates[ci].GetReferenceSequenceIndex();
    uint32_t position = candidates[ci].GetReferenceSequencePosition();
    if (candidate_strand == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      continue;
    }

    int mapping_end_position = read_length;
    int gap_beginning = 0;
    int num_errors = 0;
    const int allow_gap_beginning_ = 20;
    const int mapping_length_threshold = 30;
    int allow_gap_beginning = allow_gap_beginning_ - error_threshold_;
    int actual_num_errors = 0;
    int read_mapping_length = 0;
    int best_mapping_longest_match = 0;
    int longest_match = 0;

    if (split_alignment_) {
      if (candidate_strand == kPositive) {
        num_errors = BandedAlignPatternToTextWithDropOff(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position, &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOff(
              error_threshold_,
              reference.GetSequenceAt(rid) + position - error_threshold_ +
                  allow_gap_beginning,
              read + allow_gap_beginning, read_length - allow_gap_beginning,
              &mapping_end_position, &read_mapping_length);
          if (num_errors > error_threshold_ || mapping_end_position < 0) {
            num_errors = backup_num_errors;
            mapping_end_position = backup_mapping_end_position;
            read_mapping_length = backup_read_mapping_length;
          } else {
            gap_beginning = allow_gap_beginning;
            // Realign the mapping end position as it is the alignment from the
            // whole read.
            mapping_end_position += gap_beginning;
            // I use this adjustment since "position" is based on the whole
            // read, and it will be more consistent with no gap beginning case.
            read_mapping_length += gap_beginning;
          }
        }

        // In stitched-read mode, if the 5'-anchored alignment is too short or fails,
        // try a 3'-anchored alignment as a rescue. This handles cases where the
        // candidate represents the 3' end of the read.
        if (split_alignment_ && stitched_read_mode_ &&
            (read_mapping_length < mapping_length_threshold ||
             num_errors > error_threshold_ || mapping_end_position < 0)) {
          int mapping_end_position3p = 0;
          int read_mapping_length3p = 0;
          int num_errors3p = BandedAlignPatternToTextWithDropOffFrom3End(
              error_threshold_,
              reference.GetSequenceAt(rid) + position - error_threshold_,
              read, read_length,
              &mapping_end_position3p, &read_mapping_length3p);

          double error_rate3p = (read_mapping_length3p > 0)
                                   ? static_cast<double>(num_errors3p) /
                                         static_cast<double>(read_mapping_length3p)
                                   : 1.0;
          const double max_error_rate3p = 0.08;

          if (mapping_end_position3p >= 0 &&
              read_mapping_length3p >= mapping_length_threshold &&
              read_mapping_length3p < (int)read_length &&
              error_rate3p <= max_error_rate3p) {
            // The 3'-anchored alignment succeeded! Use it.
            num_errors = num_errors3p;
            mapping_end_position = mapping_end_position3p;
            read_mapping_length = read_mapping_length3p;
            // For positive strand, 3'-anchored means we're aligning from the 3' end
            // The unaligned portion is at the 5' end
            gap_beginning = read_length - read_mapping_length3p;
          }
        }
      } else {
        num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position,
            &read_mapping_length);
        
#ifdef CHROMAP_DEBUG
        // Diagnose why chr1:194541610 alignment failed
        {
          const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
          if (std::string(dbg_name) ==
                  "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
              rid == 0 && position >= 194541600 && position <= 194541650) {
            std::cerr << "DEBUG NEG ALIGN INITIAL: read=" << dbg_name
                      << " strand=- rid=" << rid << " pos=" << position << "\n"
                      << "  After BandedAlignPatternToTextWithDropOffFrom3End:\n"
                      << "    mapping_end_position=" << mapping_end_position
                      << " read_mapping_length=" << read_mapping_length
                      << " num_errors=" << num_errors
                      << " read_length=" << read_length
                      << " error_threshold=" << error_threshold_ << "\n"
                      << "  Interpretation: Aligned " << read_mapping_length
                      << "bp from 3' end with " << num_errors << " errors.\n";
            if (read_mapping_length < 30) {
              std::cerr << "  PROBLEM: Very short alignment (" << read_mapping_length
                        << "bp). This suggests:\n"
                        << "    1. The read doesn't match well from the 3' end\n"
                        << "    2. This might be a 5' core that needs forward alignment\n"
                        << "    3. Or the candidate position is incorrect\n";
            }
          }
        }
#endif
        
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
              error_threshold_,
              reference.GetSequenceAt(rid) + position - error_threshold_,
              negative_read.data(), read_length - allow_gap_beginning,
              &mapping_end_position, &read_mapping_length);
          if (num_errors > error_threshold_ || mapping_end_position < 0) {
            num_errors = backup_num_errors;
            mapping_end_position = backup_mapping_end_position;
            read_mapping_length = backup_read_mapping_length;
          } else {
            gap_beginning = allow_gap_beginning;
            mapping_end_position += gap_beginning;
            read_mapping_length += gap_beginning;
          }
        }
        
        // In stitched-read mode, if the 3'-anchored alignment (from reverse complement 3' end)
        // is too short or has too many errors, try a 5'-anchored alignment (from reverse complement 5' end)
        // as a rescue. This is the dual-end anchoring strategy for negative-strand candidates.
        if (split_alignment_ && stitched_read_mode_ &&
            (read_mapping_length < mapping_length_threshold ||
             num_errors > error_threshold_ || mapping_end_position < 0)) {
          int mapping_end_position5p = 0;
          int read_mapping_length5p = 0;
          int num_errors5p = BandedAlignPatternToTextWithDropOff(
              error_threshold_,
              reference.GetSequenceAt(rid) + position - error_threshold_,
              negative_read.data(), read_length,
              &mapping_end_position5p, &read_mapping_length5p);

#ifdef CHROMAP_DEBUG
          {
            const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
            if (std::string(dbg_name) ==
                    "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
                rid == 0 && position >= 194541600 && position <= 194541650) {
              std::cerr << "DEBUG NEG 5P RESCUE SEARCH: read=" << dbg_name
                        << " strand=- rid=" << rid << " pos=" << position << "\n"
                        << "  3'-anchored result: read_len=" << read_mapping_length
                        << " num_errors=" << num_errors
                        << " mapping_end=" << mapping_end_position << "\n"
                        << "  5'-anchored attempt: read_len=" << read_mapping_length5p
                        << " num_errors=" << num_errors5p
                        << " mapping_end=" << mapping_end_position5p << "\n";
            }
          }
#endif

          // For stitched Hi-C, treat the 5'-anchored alignment as a "hint":
          // accept long prefixes even if their absolute num_errors exceeds the
          // global error_threshold_, as long as the error *rate* is small.
          double error_rate5p = (read_mapping_length5p > 0)
                                   ? static_cast<double>(num_errors5p) /
                                         static_cast<double>(read_mapping_length5p)
                                   : 1.0;
          const double max_error_rate5p = 0.08;  // allow up to 8% errors for prefix

          if (mapping_end_position5p >= 0 &&
              read_mapping_length5p >= mapping_length_threshold &&
              read_mapping_length5p < (int)read_length &&
              error_rate5p <= max_error_rate5p) {
            // The 5'-anchored alignment succeeded! Use it.
            num_errors = num_errors5p;
            mapping_end_position = mapping_end_position5p;
            read_mapping_length = read_mapping_length5p;

            // CRITICAL FIX: For negative strand, 5'-anchored alignment on RC
            // means we're aligning from the 5' end of the original read.
            // Set gap_beginning to indicate this is a 5' core.
            // The unaligned portion is at the 3' end of the original read.
            gap_beginning = read_length - read_mapping_length5p;

#ifdef CHROMAP_DEBUG
            {
              const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
              if (std::string(dbg_name) ==
                      "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
                  rid == 0 && position >= 194541600 && position <= 194541650) {
                std::cerr << "DEBUG NEG 5P RESCUE OK: read=" << dbg_name
                          << " rid=" << rid << " pos=" << position
                          << " num_errors5p=" << num_errors5p
                          << " read_mapping_length5p=" << read_mapping_length5p
                          << " error_rate5p=" << error_rate5p
                          << " gap_beginning=" << gap_beginning << " (5' core indicator)\n";
              }
            }
#endif
          }
        }
      }

#ifdef CHROMAP_DEBUG
      // Detailed debug for the problematic stitched read at the first-site
      // candidate near chr7:19858730 (rid=56, pos~19858680 in this index).
      {
        const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
        if (std::string(dbg_name) ==
                "LH00708:218:22WYCCLT4:6:1102:18911:1532" &&
            candidate_strand == kPositive && rid == 56 &&
            position == 19858680) {
          std::cerr << "DEBUG SPLIT DETAIL: read=" << dbg_name
                    << " strand=+"
                    << " rid=" << rid
                    << " pos=" << position
                    << " mapping_end_position=" << mapping_end_position
                    << " num_errors=" << num_errors
                    << " gap_beginning=" << gap_beginning
                    << " read_mapping_length=" << read_mapping_length
                    << " error_threshold=" << error_threshold_ << "\n";
        }
      }
#endif

      // Debug for chr1:194541610 investigation
#ifdef CHROMAP_DEBUG
      {
        const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
        if (std::string(dbg_name) ==
                "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
            rid == 0 && position >= 194541600 && position <= 194541650 &&
            candidate_strand == kNegative) {
          std::cerr << "DEBUG ALIGN_LEN CALC: read=" << dbg_name
                    << " strand=- rid=" << rid << " pos=" << position << "\n"
                    << "  BEFORE align_len calc:\n"
                    << "    mapping_end_position=" << mapping_end_position
                    << " num_errors=" << num_errors
                    << " gap_beginning=" << gap_beginning
                    << " read_mapping_length=" << read_mapping_length
                    << " error_threshold=" << error_threshold_
                    << " read_length=" << read_length << "\n";
        }
      }
#endif

      int align_len = mapping_end_position + 1 - error_threshold_ - num_errors -
                      gap_beginning;

#ifdef CHROMAP_DEBUG
      {
        const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
        if (std::string(dbg_name) ==
                "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
            rid == 0 && position >= 194541600 && position <= 194541650 &&
            candidate_strand == kNegative) {
          std::cerr << "  AFTER align_len calc:\n"
                    << "    align_len=" << align_len
                    << " (formula: " << mapping_end_position << " + 1 - "
                    << error_threshold_ << " - " << num_errors << " - "
                    << gap_beginning << ")\n"
                    << "    threshold=" << mapping_length_threshold
                    << " -> " << (align_len >= mapping_length_threshold ? "KEEP" : "DISCARD") << "\n";
        }
      }
#endif

      // In stitched-read mode, if the 5'-anchored split alignment is too short
      // or has too many errors, try a 3'-anchored alignment to discover the
      // ligation junction and then realign from that split site on the 5' side.
      if (split_alignment_ && stitched_read_mode_ &&
          candidate_strand == kPositive &&
          (align_len < mapping_length_threshold ||
           num_errors > error_threshold_ || mapping_end_position < 0)) {
        int mapping_end_position3 = 0;
        int read_mapping_length3 = 0;
        int num_errors3 = BandedAlignPatternToTextWithDropOffFrom3End(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position3, &read_mapping_length3);

#ifdef CHROMAP_DEBUG
        {
          const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
          if (std::string(dbg_name) ==
                  "LH00708:218:22WYCCLT4:6:1102:18911:1532" &&
              rid == 56 && position == 19858680) {
            std::cerr << "DEBUG 3P SEARCH: read=" << dbg_name
                      << " rid=" << rid << " pos=" << position
                      << " num_errors3=" << num_errors3
                      << " mapping_end_position3=" << mapping_end_position3
                      << " read_mapping_length3=" << read_mapping_length3
                      << " align_len5p=" << align_len
                      << " num_errors5p=" << num_errors << "\n";
          }
        }
#endif

        // For stitched Hi-C, treat the 3'-anchored alignment as a "hint":
        // accept long suffixes even if their absolute num_errors exceeds the
        // global error_threshold_, as long as the error *rate* is small.
        double error_rate3 = (read_mapping_length3 > 0)
                                 ? static_cast<double>(num_errors3) /
                                       static_cast<double>(read_mapping_length3)
                                 : 1.0;
        const double max_error_rate3 = 0.08;  // allow up to 8% errors for suffix

        if (mapping_end_position3 >= 0 &&
            read_mapping_length3 >= mapping_length_threshold &&
            read_mapping_length3 < (int)read_length &&
            error_rate3 <= max_error_rate3) {
          // Derive the split site on the read from the 3'-anchored suffix
          // length, e.g., L_read - L_suffix.
          int discovered_split = read_length - read_mapping_length3;

          if (discovered_split > 0 && discovered_split < (int)read_length) {
            // Adaptive "walk-back" safety buffer: try to avoid the noisy
            // junction region that the 3'-anchored alignment may have
            // over-extended into. Use ~5% of read length, capped at 15 bp.
            const int max_safety_bp = 15;
            int safety_buffer =
                std::min(max_safety_bp, static_cast<int>(read_length * 0.05));

            // Ensure we do not trim away almost the entire suffix.
            if (read_length - discovered_split - safety_buffer <
                mapping_length_threshold) {
              safety_buffer = 0;
            }

            auto try_rescue_from =
                [&](int extra_shift, int &out_num_errors,
                    int &out_mapping_end_position,
                    int &out_read_mapping_length) -> bool {
              int mapping_end_position_resc = 0;
              int read_mapping_length_resc = 0;
              int num_errors_resc = BandedAlignPatternToTextWithDropOff(
                  error_threshold_,
                  reference.GetSequenceAt(rid) + position - error_threshold_ +
                      discovered_split + extra_shift,
                  read + discovered_split + extra_shift,
                  read_length - discovered_split - extra_shift,
                  &mapping_end_position_resc, &read_mapping_length_resc);

              double error_rate_resc =
                  (read_mapping_length_resc > 0)
                      ? static_cast<double>(num_errors_resc) /
                            static_cast<double>(read_mapping_length_resc)
                      : 1.0;
              bool resc_strict_ok = (num_errors_resc <= error_threshold_);
              bool resc_long_ok = (read_mapping_length_resc > 50 &&
                                   error_rate_resc <= 0.08);

              if ((resc_strict_ok || resc_long_ok) &&
                  mapping_end_position_resc >= 0) {
                out_num_errors = num_errors_resc;
                out_mapping_end_position = mapping_end_position_resc;
                out_read_mapping_length = read_mapping_length_resc;
                return true;
              }
              return false;
            };

            int chosen_shift = -1;
            int chosen_num_errors = 0;
            int chosen_mapping_end_position = 0;
            int chosen_read_mapping_length = 0;

            // First, try from discovered_split + safety_buffer (conservative).
            if (try_rescue_from(safety_buffer, chosen_num_errors,
                                chosen_mapping_end_position,
                                chosen_read_mapping_length)) {
              chosen_shift = safety_buffer;
            } else {
              // Fallback: try from discovered_split itself (no buffer) in case
              // the 3' alignment did not actually overrun into noisy bases.
              if (try_rescue_from(0, chosen_num_errors,
                                  chosen_mapping_end_position,
                                  chosen_read_mapping_length)) {
                chosen_shift = 0;
              }
            }

            if (chosen_shift >= 0) {
#ifdef CHROMAP_DEBUG
              {
                const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
                if (std::string(dbg_name) ==
                        "LH00708:218:22WYCCLT4:6:1102:18911:1532" &&
                    rid == 56 && position == 19858680) {
                  std::cerr << "DEBUG 3P RESCUE OK: read=" << dbg_name
                            << " rid=" << rid << " pos=" << position
                            << " discovered_split=" << discovered_split
                            << " safety_buffer=" << safety_buffer
                            << " chosen_shift=" << chosen_shift
                            << " num_errors_resc=" << chosen_num_errors
                            << " mapping_end_position_resc="
                            << chosen_mapping_end_position
                            << " read_mapping_length_resc="
                            << chosen_read_mapping_length << "\n";
                }
              }
#endif
              // Successful rescue: update alignment statistics to reflect a
              // 5'-anchored alignment starting at the adjusted split site.
              gap_beginning = discovered_split + chosen_shift;
              mapping_end_position =
                  chosen_mapping_end_position + gap_beginning;
              read_mapping_length =
                  chosen_read_mapping_length + gap_beginning;
              num_errors = chosen_num_errors;
              align_len = mapping_end_position + 1 - error_threshold_ -
                          num_errors - gap_beginning;
            }
          }
        }
      }

      if (align_len >= mapping_length_threshold) {
        actual_num_errors = num_errors;
        num_errors = -(mapping_end_position - error_threshold_ - num_errors -
                       gap_beginning);

#ifdef CHROMAP_DEBUG
        // Investigate why chr6 cores have short alignments
        {
          const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
          if (std::string(dbg_name) ==
                  "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
              ((rid == 55 && (position == 102168919 || position == 84164305)) ||
               (rid == 0 && position >= 194541600 && position <= 194541650))) {
            std::cerr << "DEBUG KEPT ALIGNMENT: read=" << dbg_name
                      << " strand=" << (candidate_strand == kPositive ? '+' : '-')
                      << " rid=" << rid << " pos=" << position << "\n"
                      << "  align_len=" << align_len
                      << " read_length=" << read_length
                      << " gap_beginning=" << gap_beginning
                      << " read_mapping_length=" << read_mapping_length
                      << " mapping_end_position=" << mapping_end_position
                      << " num_errors=" << actual_num_errors
                      << " encoded_num_errors=" << num_errors << "\n";
            if (align_len < (int)read_length * 0.5) {
              std::cerr << "  WARNING: Short alignment! Only " << align_len
                        << "bp of " << read_length << "bp aligned ("
                        << (align_len * 100 / read_length) << "%)\n";
            }
          }
        }
#endif

        if (candidates.size() > 200) {
          if (candidate_strand == kPositive) {
            longest_match = GetLongestMatchLength(
                reference.GetSequenceAt(rid) + position, read, read_length);
          } else {
            longest_match =
                GetLongestMatchLength(reference.GetSequenceAt(rid) + position,
                                      negative_read.data(), read_length);
          }
        }
      } else {
#ifdef CHROMAP_DEBUG
        // Compact, throttled logging of discarded split segments.
        static thread_local int last_read_index = -1;
        static thread_local int discard_count_for_read = 0;
        if ((int)read_index != last_read_index) {
          // New read: reset counter and print a header.
          last_read_index = read_index;
          discard_count_for_read = 0;
          std::cerr << "DEBUG DISCARD SPLIT SUMMARY: read_index=" << read_index
                    << " name=" << read_batch.GetSequenceNameAt(read_index)
                    << " (logging up to 20 discarded splits)\n";
        }
        if (discard_count_for_read < 20) {
          std::cerr << "  split#" << discard_count_for_read
                    << " strand=" << (candidate_strand == kPositive ? '+' : '-')
                    << " rid=" << rid
                    << " pos=" << position
                    << " align_len=" << align_len
                    << " < min_len=" << mapping_length_threshold;
          
          // Detailed investigation for chr1:194541610 discard
          const char *dbg_name = read_batch.GetSequenceNameAt(read_index);
          if (std::string(dbg_name) ==
                  "LH00708:218:22WYCCLT4:6:1102:20578:1532" &&
              rid == 0 && position >= 194541600 && position <= 194541650) {
            std::cerr << "\n    *** INVESTIGATING chr1:194541610 DISCARD ***\n"
                      << "      mapping_end_position=" << mapping_end_position
                      << " num_errors=" << num_errors
                      << " gap_beginning=" << gap_beginning
                      << " read_mapping_length=" << read_mapping_length
                      << " error_threshold=" << error_threshold_
                      << " read_length=" << read_length
                      << "\n      align_len formula: " << mapping_end_position << " + 1 - "
                      << error_threshold_ << " - " << num_errors << " - " << gap_beginning
                      << " = " << align_len << "\n";
          } else {
            std::cerr << "\n";
          }
        }
        ++discard_count_for_read;
#endif
        // In both normal and stitched modes, discard alignments shorter than
        // mapping_length_threshold. This avoids rescuing very short fragments
        // that are highly ambiguous and can dominate candidate ranking.
        num_errors = error_threshold_ + 1;
        actual_num_errors = error_threshold_ + 1;
      }
    } else {
      if (candidate_strand == kPositive) {
        num_errors = BandedAlignPatternToText(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position);
      } else {
        num_errors = BandedAlignPatternToText(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position);
      }
    }

    if (num_errors <= error_threshold_) {
      if (num_errors < min_num_errors) {
        second_min_num_errors = min_num_errors;
        num_second_best_mappings = num_best_mappings;
        min_num_errors = num_errors;
        num_best_mappings = 1;
        if (split_alignment_) {
          if (candidates.size() > 50) {
            candidate_count_threshold = candidates[ci].count;
          } else {
            candidate_count_threshold = candidates[ci].count / 2;
          }
          if (second_min_num_errors < min_num_errors + error_threshold_ / 2 &&
              best_mapping_longest_match > longest_match &&
              candidates.size() > 200) {
            second_min_num_errors = min_num_errors;
          }
        }
        best_mapping_longest_match = longest_match;
      } else if (num_errors == min_num_errors) {
        num_best_mappings++;
      } else if (num_errors == second_min_num_errors) {
        num_second_best_mappings++;
      } else if (num_errors < second_min_num_errors) {
        num_second_best_mappings = 1;
        second_min_num_errors = num_errors;
      }

      if (candidate_strand == kPositive) {
        mappings.emplace_back(
            num_errors,
            candidates[ci].position - error_threshold_ + mapping_end_position);
      } else {
        if (split_alignment_ && mapping_output_format_ != MAPPINGFORMAT_SAM) {
          // TODO: this if condition is suspicious. Check this later.
          mappings.emplace_back(num_errors,
                                candidates[ci].position - gap_beginning);
        } else {
          // Need to minus gap_beginning because mapping_end_position is
          // adjusted by it, but read_length is not.
          // printf("%d %d %d\n", candidates[ci].position, mapping_end_position,
          // gap_beginning);
          mappings.emplace_back(num_errors,
                                candidates[ci].position - read_length + 1 -
                                    error_threshold_ + mapping_end_position);
        }
      }

      if (split_alignment_) {
        split_sites.emplace_back(((actual_num_errors & 0xff) << 24) |
                                 ((gap_beginning & 0xff) << 16) |
                                 (read_mapping_length & 0xffff));
      }
    }
  }
}

}  // namespace chromap
