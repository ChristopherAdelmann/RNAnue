#include "Detect.hpp"

using namespace seqan3::literals;

// WARNING: This split read calling is not working as intended and is currently not handling
// multiple alignments correctly.

Detect::Detect(po::variables_map params) : params(params) {
    readsCount = 0;
    alignedcount = 0;
    splitscount = 0;
    msplitscount = 0;
    nsurvivedcount = 0;
}
//
void Detect::iterate(std::string matched, std::string splitsPath, std::string multSplitsPath) {
    seqan3::sam_file_input alignmentsIn{matched, sam_field_ids{}};

    std::vector<size_t> ref_lengths{};
    std::ranges::transform(alignmentsIn.header().ref_id_info, std::back_inserter(ref_lengths),
                           [](auto const &info) { return std::get<0>(info); });

    std::deque<std::string> &ref_ids = alignmentsIn.header().ref_ids();

    // output file for splits
    seqan3::sam_file_output splitsOut{splitsPath, ref_ids, ref_lengths, sam_field_ids{}};

    // output file for multi splits
    seqan3::sam_file_output multiSplitsOut{multSplitsPath, ref_ids, ref_lengths, sam_field_ids{}};

    std::string currentRecordId = "";

    using RecordType = typename decltype(alignmentsIn)::record_type;
    std::vector<RecordType> currentRecordList;

    for (auto &record : alignmentsIn) {
        if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) {
            continue;
        }
        if (record.sequence().size() <= params["minlen"].as<std::size_t>()) {
            continue;
        }

        // This line depends on the SAM file being sorted by read ID
        std::string recordId = record.id();
        bool isNewRecord = !currentRecordId.empty() && (currentRecordId != recordId);

        if (isNewRecord) {
            process(currentRecordList, splitsOut, multiSplitsOut);
            currentRecordList.clear();
        }

        currentRecordList.push_back(std::move(record));
        currentRecordId = recordId;

        readsCount++;
        if (readsCount % 100000 == 0) {
            Logger::log(LogLevel::INFO, "Processed ", readsCount, " reads");
        }
    }

    if (!currentRecordList.empty()) {
        readsCount++;
        process(currentRecordList, splitsOut, multiSplitsOut);
    }

    Logger::log(LogLevel::INFO, "Processed ", readsCount, " reads");
}

void Detect::process(auto &splitrecords, auto &splitsfile, auto &multsplitsfile) {
    Splits splits;

    seqan3::cigar::operation cigarOp;  // operation of cigar string
    uint32_t cigarOpSize;

    uint32_t startPosRead;
    uint32_t endPosRead;  // absolute position in read/alignment (e.g., 1 to end)
    uint32_t startPosSplit;
    uint32_t endPosSplit;  // position in split read (e.g., XX:i to XY:i / 14 to 20)

    seqan3::sam_tag_dictionary tags;

    int segmentNr = 0;  // intialize active segment
    int splitID = 0;    //

    std::string qname = "";
    std::vector<SamRecord> curated;

    // create object of putative splits
    std::map<int, std::vector<SamRecord>> putative;
    std::vector<std::pair<SamRecord, SamRecord>> splitSegments;

    size_t process_count = 0;

    for (const auto &it : splitrecords) {
        // extract information
        // TODO Handle remaining tags
        qname = it.id();  // QNAME
        // seqan3::sam_flag flag = it->flag();  // SAMFLAG
        //  std::optional<int32_t> refID = it->reference_id();
        std::optional<int32_t> refOffset = it.reference_position();
        // std::optional<uint8_t> qual = it->mapping_quality();  // MAPQ

        // CIGAR string
        std::vector<seqan3::cigar> cigar{it.cigar_sequence()};
        std::vector<seqan3::cigar> cigarSplit{};  // individual cigar for each split
        const auto &seq = it.sequence();          // SEQ

        // TODO Handle remaining tags
        // extract the tags (information about split reads)
        tags = it.tags();
        // auto xhtag = tags.get<"XH"_tag>();
        auto xjtag = tags.get<"XJ"_tag>();
        auto xxtag = tags.get<"XX"_tag>();
        // auto xytag = tags.get<"XY"_tag>();

        if (xjtag == 2) {      // splits in SAMfile need to consists of at least 2 segments
            startPosRead = 1;  // intialise start & end of reads
            endPosRead = 0;

            startPosSplit = xxtag;  // determine the start position within split
            endPosSplit = xxtag - 1;
            // check the cigar string for splits
            for (auto &cig : cigar) {
                // determine size and operator of cigar element
                cigarOpSize = get<uint32_t>(cig);
                cigarOp = get<seqan3::cigar::operation>(cig);

                //
                if (cigarOp == 'N'_cigar_operation) {
                    auto subSeq = seq | seqan3::views::slice(startPosRead - 1, endPosRead);

                    // add properties to tags - lightweight
                    seqan3::sam_tag_dictionary ntags;
                    ntags.get<"XX"_tag>() = startPosSplit;
                    ntags.get<"XY"_tag>() = endPosSplit;
                    ntags["XN"_tag] = splitID;

                    filterSegments(it, refOffset, cigarSplit, subSeq, ntags, curated);

                    // settings for prepare for next split
                    startPosSplit = endPosSplit + 1;
                    startPosRead = endPosRead + 1;
                    cigarSplit = {};  // new split - new CIGAR
                    refOffset.value() +=
                        cigarOpSize + endPosRead + 1;  // adjust leftmost mapping position

                    // change for next iteration
                    segmentNr++;  // counts as segment of split
                } else {
                    // deletion does not account for length in split
                    seqan3::cigar cigarElement{cigarOpSize, cigarOp};
                    // seqan3::debug_stream << "cigarElement: " << cigarElement << std::endl;

                    // exclude soft clipping from alignment
                    if (cigarOp == 'S'_cigar_operation &&
                        params["exclclipping"].as<std::bitset<1>>() == 1) {
                        if (cigarSplit.size() == 0) {
                            startPosRead += cigarOpSize;
                            endPosRead += cigarOpSize;
                            startPosSplit += cigarOpSize;
                            endPosSplit += cigarOpSize;
                        }
                    } else {
                        if (cigarOp != 'D'_cigar_operation) {
                            endPosRead += cigarOpSize;
                            endPosSplit += cigarOpSize;
                        }
                        cigarSplit.push_back(cigarElement);
                    }
                }
            }

            auto subSeq = seq | seqan3::views::slice(startPosRead - 1, endPosRead);

            seqan3::sam_tag_dictionary ntags;
            ntags.get<"XX"_tag>() = startPosSplit;
            ntags.get<"XY"_tag>() = endPosSplit;
            ntags["XN"_tag] = splitID;

            filterSegments(it, refOffset, cigarSplit, subSeq, ntags, curated);
            segmentNr++;

            if (segmentNr == xjtag) {
                std::map<int, std::vector<SamRecord>>::iterator itPutative = putative.begin();
                putative.insert(itPutative, std::make_pair(splitID, curated));
                curated.clear();
                segmentNr = 0;
                ++splitID;
            }
        }
    }

    if (putative.size() > 0) {  // split reads found
        std::vector<SamRecord> p1;
        std::vector<SamRecord> p2;

        seqan3::sam_tag_dictionary p1Tags;
        seqan3::sam_tag_dictionary p2Tags;

        // value for complementarity and hybridization energy
        std::pair<double, double> filters;

        std::string p1Qname;
        std::string p2Qname;

        for (auto &element : putative) {
            std::vector<SamRecord> splits = element.second;
            if (splits.size() > 1) {  // splits consists at least of 2 segments
                for (unsigned i = 0; i < splits.size(); ++i) {
                    for (unsigned j = i + 1; j < splits.size(); ++j) {
                        p1Tags = splits[i].tags();
                        p2Tags = splits[j].tags();

                        auto p1Start = p1Tags.get<"XX"_tag>();
                        auto p1End = p1Tags.get<"XY"_tag>();
                        auto p2Start = p2Tags.get<"XX"_tag>();
                        auto p2End = p2Tags.get<"XY"_tag>();

                        // prevent an overlap between reads position -> same segment of read
                        if ((p2Start >= p1Start && p2Start <= p1End) ||
                            (p1Start >= p2Start && p1Start <= p2End)) {
                            continue;
                        }

                        if (true)  // splitSegments.empty()
                        {
                            TracebackResult initCmpl;
                            std::optional<HybridizationResult> initHyb;
                            {
                                // helper::Timer timer;
                                complementaritySeqAn(splits[i].sequence(), splits[j].sequence());
                                initCmpl =
                                    complementarity(splits[i].sequence(), splits[j].sequence());
                                process_count++;
                            }
                            { initHyb = hybridize(splits[i].sequence(), splits[j].sequence()); }

                            // split read survived complementarity & sitelenratio cutoff
                            if (!initCmpl.a.empty() && initHyb.has_value()) {
                                double &hybEnergy = initHyb.value().energy;
                                if (hybEnergy <= params["nrgmax"].as<double>()) {
                                    filters = std::make_pair(initCmpl.cmpl, hybEnergy);
                                    addComplementarityToSamRecord(splits[i], splits[j], initCmpl);
                                    addHybEnergyToSamRecord(splits[i], splits[j], initHyb.value());
                                    splitSegments.push_back(std::make_pair(splits[i], splits[j]));
                                }
                            }
                        } else {
                            complementaritySeqAn(splits[i].sequence(), splits[j].sequence());

                            TracebackResult cmpl =
                                complementarity(splits[i].sequence(), splits[j].sequence());

                            seqan3::debug_stream << "Sequences ELSE: " << splits[i].sequence()
                                                 << " " << splits[j].sequence() << "\n";
                            std::cout << "cmpl: " << cmpl.score << std::endl;
                            std::optional<HybridizationResult> hyb =
                                hybridize(splits[i].sequence(), splits[j].sequence());

                            if (!cmpl.a.empty() &&
                                hyb.has_value()) {  // check that there is data for
                                                    // complementarity / passed cutoffs
                                double hybEnergy = hyb.value().energy;
                                if (cmpl.cmpl > filters.first) {  // complementarity is higher
                                    splitSegments.clear();
                                    filters = std::make_pair(cmpl.cmpl, hybEnergy);
                                    addComplementarityToSamRecord(splits[i], splits[j], cmpl);
                                    addHybEnergyToSamRecord(splits[i], splits[j], hyb.value());
                                    splitSegments.push_back(std::make_pair(splits[i], splits[j]));
                                } else {
                                    if (cmpl.cmpl == filters.first) {
                                        if (hybEnergy > filters.second) {
                                            splitSegments.clear();
                                            filters = std::make_pair(cmpl.cmpl, hybEnergy);
                                            addComplementarityToSamRecord(splits[i], splits[j],
                                                                          cmpl);
                                            addHybEnergyToSamRecord(splits[i], splits[j],
                                                                    hyb.value());
                                            splitSegments.push_back(
                                                std::make_pair(splits[i], splits[j]));
                                        }
                                    } else {
                                        if (hybEnergy ==
                                            filters.second) {  // same cmpl & hyb -> ambigious read!
                                            addComplementarityToSamRecord(splits[i], splits[j],
                                                                          cmpl);
                                            addHybEnergyToSamRecord(splits[i], splits[j],
                                                                    hyb.value());
                                            splitSegments.push_back(
                                                std::make_pair(splits[i], splits[j]));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "Processed " << process_count << " split reads" << std::endl;

    {
        // TODO Multisplits file not correct, segments size can be higher due to multiple alignments
        // without splits! write to file
        if (splitSegments.size() == 1) {
            writeSamFile(splitsfile, splitSegments);
            splitscount++;
        } else {
            if (splitSegments.size() > 1) {
                writeSamFile(multsplitsfile, splitSegments);
                msplitscount++;
            }
        }
        std::vector<std::pair<SamRecord, SamRecord>>().swap(splitSegments);
    }
}

// write splits back to file
void Detect::writeSamFile(auto &samfile, std::vector<std::pair<SamRecord, SamRecord>> &splits) {
    for (auto &[first, second] : splits) {
        auto &[id1, flag1, ref_id1, ref_offset1, mapq1, cigar1, seq1, tags1] = first;
        samfile.emplace_back(id1, flag1, ref_id1, ref_offset1, mapq1.value(), cigar1, seq1, tags1);

        auto &[id2, flag2, ref_id2, ref_offset2, mapq2, cigar2, seq2, tags2] = second;
        samfile.emplace_back(id2, flag2, ref_id2, ref_offset2, mapq2.value(), cigar2, seq2, tags2);
    }
}

void Detect::filterSegments(const auto &splitrecord, std::optional<int32_t> &refOffset,
                            std::vector<seqan3::cigar> &cigar, std::span<const seqan3::dna5> seq,
                            seqan3::sam_tag_dictionary &tags, std::vector<SamRecord> &curated) {
    SamRecord segment{};  // create new SamRecord

    seqan3::sam_flag flag{0};
    if (static_cast<bool>(splitrecord.flag() & seqan3::sam_flag::on_reverse_strand)) {
        flag = seqan3::sam_flag{seqan3::sam_flag::on_reverse_strand};
    }

    segment.id() = splitrecord.id();
    segment.flag() = flag;
    segment.reference_id() = splitrecord.reference_id();
    segment.reference_position() = refOffset;
    segment.mapping_quality() = splitrecord.mapping_quality();
    segment.cigar_sequence() = cigar;
    segment.sequence().assign(seq.begin(), seq.end());
    segment.tags() = tags;

    if (seq.size() <= (std::size_t)params["minfraglen"].as<int>()) {
        return;
    }

    // std::cout << "length " << seq.size() << std::endl;
    curated.push_back(segment);
}

void Detect::addFilterToSamRecord(SamRecord &rec, std::pair<float, float> filters) {
    rec.tags()["XC"_tag] = filters.first;
    rec.tags()["XE"_tag] = filters.second;
}

void Detect::addComplementarityToSamRecord(SamRecord &rec1, SamRecord &rec2, TracebackResult &res) {
    // number of matches
    rec1.tags()["XM"_tag] = res.matches;
    rec2.tags()["XM"_tag] = res.matches;

    // length of alignment
    rec1.tags()["XL"_tag] = res.length;
    rec2.tags()["XL"_tag] = res.length;

    // complementarity
    rec1.tags()["XC"_tag] = static_cast<float>(res.cmpl);
    rec2.tags()["XC"_tag] = static_cast<float>(res.cmpl);

    // sitelenratio
    rec1.tags()["XR"_tag] = static_cast<float>(res.ratio);
    rec2.tags()["XR"_tag] = static_cast<float>(res.ratio);

    // alignments
    rec1.tags()["XA"_tag] = res.a;
    rec2.tags()["XA"_tag] = res.b;

    // score
    rec1.tags()["XS"_tag] = res.score;
    rec2.tags()["XS"_tag] = res.score;
}

//
void Detect::addHybEnergyToSamRecord(SamRecord &rec1, SamRecord &rec2, HybridizationResult &hyb) {
    rec1.tags()["XE"_tag] = static_cast<float>(hyb.energy);
    rec2.tags()["XE"_tag] = static_cast<float>(hyb.energy);

    if (hyb.crosslinkingResult.has_value()) {
        CrosslinkingResult &crosslinkingResult = hyb.crosslinkingResult.value();
        rec1.tags()["XK"_tag] = static_cast<float>(crosslinkingResult.normCrosslinkingScore);
        rec2.tags()["XK"_tag] = static_cast<float>(crosslinkingResult.normCrosslinkingScore);

        rec1.tags()["XP"_tag] = crosslinkingResult.prefferedCrosslinkingScore;
        rec2.tags()["XP"_tag] = crosslinkingResult.prefferedCrosslinkingScore;

        rec1.tags()["XO"_tag] = crosslinkingResult.nonPrefferedCrosslinkingScore;
        rec2.tags()["XO"_tag] = crosslinkingResult.nonPrefferedCrosslinkingScore;

        rec1.tags()["XW"_tag] = crosslinkingResult.wobbleCrosslinkingScore;
        rec2.tags()["XW"_tag] = crosslinkingResult.wobbleCrosslinkingScore;
    }
}

//
TracebackResult Detect::complementarity(std::vector<seqan3::dna5> &seq1,
                                        std::vector<seqan3::dna5> &seq2) {
    std::string seq1_str;
    std::string seq2_str;

    if (seq1.size() == 0 || seq2.size() == 0) {
        return {"", "", 0, 0, 0, -1.0, 0.0};
    }

    // rna1 (reverse)
    for (unsigned z = seq1.size(); z-- > 0;) {
        seq1_str += seq1[z].to_char();
    }

    // rna
    for (unsigned y = 0; y < seq2.size(); ++y) {
        seq2_str += seq2[y].to_char();
    }

    ScoringMatrix matrix = createScoringMatrix(seq1_str.c_str(), seq2_str.c_str(), 1, -1, 2);
    std::vector<TracebackResult> res = traceback(matrix, seq1_str.c_str(), seq2_str.c_str());

    int idx = -1;
    int cmpl = -1;
    // filter (multiple) alignments
    if (res.size() > 0) {
        for (unsigned i = 0; i < res.size(); ++i) {
            // check if sitelenratio & complementarity exceed cutoffs
            if (params["sitelenratio"].as<double>() <= res[i].ratio &&
                params["cmplmin"].as<double>() <= res[i].cmpl) {
                if (res[i].cmpl > cmpl) {
                    idx = i;
                }
            }
        }
    }

    freeMatrix(&matrix);
    if (idx != -1) {
        return res[idx];
    } else {
        return {"", "", 0, 0, 0, 0.0, 0.0};
    }
}

void Detect::complementaritySeqAn(seqan3::dna5_vector &seq1, seqan3::dna5_vector &seq2) {
    seqan3::nucleotide_scoring_scheme scheme;
    scheme.set_simple_scheme(seqan3::match_score{1}, seqan3::mismatch_score{-1});

    scheme.score('A'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'A'_dna5) = 1;
    scheme.score('G'_dna5, 'C'_dna5) = 1;
    scheme.score('C'_dna5, 'G'_dna5) = 1;
    scheme.score('G'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'G'_dna5) = 1;

    scheme.score('T'_dna5, 'T'_dna5) = -1;
    scheme.score('A'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'G'_dna5) = -1;

    scheme.score('A'_dna5, 'G'_dna5) = -1;
    scheme.score('A'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'A'_dna5) = -1;

    scheme.score('T'_dna5, 'C'_dna5) = -1;
    scheme.score('C'_dna5, 'T'_dna5) = -1;

    const auto config = seqan3::align_cfg::method_local{} |
                        seqan3::align_cfg::scoring_scheme{scheme} |
                        seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-2},
                                                           seqan3::align_cfg::extension_score{-2}} |
                        seqan3::align_cfg::output_alignment{} | seqan3::align_cfg::output_score{} |
                        seqan3::align_cfg::detail::debug{};
    ;

    const auto &seq1Reverse = seq1 | std::views::reverse;

    seqan3::debug_stream << "Seq1: " << seq1 << " Seq2: " << seq2 << std::endl;

    using trace_matrix_t =
        seqan3::detail::two_dimensional_matrix<std::optional<seqan3::detail::trace_directions>>;
    using score_matrix_t = seqan3::detail::row_wise_matrix<std::optional<int>>;
    for (auto const &result : seqan3::align_pairwise(std::tie(seq1Reverse, seq2), config)) {
        seqan3::debug_stream << "Original result: " << result << std::endl;

        score_matrix_t score_matrix{result.score_matrix()};
        trace_matrix_t trace_matrix{result.trace_matrix()};

        const std::optional<int> score = result.score();

        // Find all positions of the maximum score.

        auto it = std::find(score_matrix.begin(), score_matrix.end(), score);
        while (it != score_matrix.end()) {
            size_t row = std::distance(score_matrix.begin(), it) / score_matrix.cols();
            size_t col = std::distance(score_matrix.begin(), it) % score_matrix.cols();

            seqan3::detail::matrix_coordinate index{};
            index.row = row;
            index.col = col;
            const auto result =
                Detect::make_result(std::tie(seq2, seq1Reverse), 0, score, index, trace_matrix);

            seqan3::debug_stream << "Result: " << result.score.value()
                                 << " END FIRST: " << result.end_positions.first
                                 << " END SECOND: " << result.end_positions.second
                                 << " Begin first: " << result.begin_positions.first
                                 << " Begin second: " << result.begin_positions.second << std::endl;

            it = std::find(std::next(it), score_matrix.end(), score);
        }
    }
}

auto Detect::trace_path(
    seqan3::detail::matrix_coordinate const &trace_begin,
    seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> &complete_matrix,
    const size_t row_count, const size_t column_count) {
    using matrix_t = seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions>;
    using matrix_iter_t = std::ranges::iterator_t<matrix_t const>;
    using trace_iterator_t = seqan3::detail::trace_iterator<matrix_iter_t>;
    using path_t = std::ranges::subrange<trace_iterator_t, std::default_sentinel_t>;

    if (trace_begin.row >= row_count || trace_begin.col >= column_count)
        throw std::invalid_argument{
            "The given coordinate exceeds the matrix in vertical or horizontal direction."};

    seqan3::detail::matrix_offset const offset{trace_begin};

    return path_t{trace_iterator_t{complete_matrix.begin() + offset}, std::default_sentinel};
}

template <typename sequence_pair_t, typename index_t, typename score_t,
          typename matrix_coordinate_t>
al_res Detect::make_result(
    [[maybe_unused]] sequence_pair_t &&sequence_pair, [[maybe_unused]] index_t &&id,
    [[maybe_unused]] score_t score, [[maybe_unused]] matrix_coordinate_t end_positions,
    [[maybe_unused]] seqan3::detail::two_dimensional_matrix<
        std::optional<seqan3::detail::trace_directions>> const &alignment_matrix) {
    using std::get;

    al_res result{};

    result.sequence1_id = id;

    result.sequence2_id = id;

    // Cast alignment matrix to
    // seqan3::detail::two_dimensional_matrixseqan3::detail::trace_directions>

    std::vector<seqan3::detail::trace_directions> trace_directions;

    for (auto const &direction : alignment_matrix) {
        trace_directions.push_back(direction.value_or(seqan3::detail::trace_directions::none));
    }

    seqan3::detail::number_rows rows_n{alignment_matrix.rows()};
    seqan3::detail::number_cols cols_n{alignment_matrix.cols()};

    seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> trace_matrix(
        rows_n, cols_n, trace_directions);

    result.score = std::move(score);
    result.end_positions.first = end_positions.col;
    result.end_positions.second = end_positions.row;

    seqan3::detail::aligned_sequence_builder builder{get<0>(sequence_pair), get<1>(sequence_pair)};
    auto aligned_sequence_result = builder(
        trace_path(end_positions, trace_matrix, alignment_matrix.rows(), alignment_matrix.cols()));

    result.begin_positions.first = aligned_sequence_result.first_sequence_slice_positions.first;
    result.begin_positions.second = aligned_sequence_result.second_sequence_slice_positions.first;

    seqan3::debug_stream << aligned_sequence_result.alignment << std::endl;

    return result;
}

std::optional<HybridizationResult> Detect::hybridize(std::span<const seqan3::dna5> seq1,
                                                     std::span<const seqan3::dna5> seq2) {
    auto toString = [](auto &seq) {
        return (seq | seqan3::views::to_char | seqan3::ranges::to<std::string>());
    };

    std::string interactionSeq = toString(seq1) + "&" + toString(seq2);

    vrna_fold_compound_t *vc =
        vrna_fold_compound(interactionSeq.c_str(), NULL, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
    char *structure = new char[interactionSeq.size() + 1];
    float mfe = vrna_cofold(interactionSeq.c_str(), structure);

    auto secondaryStructure = std::string(vrna_cut_point_insert(structure, seq1.size() + 1)) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    std::optional<CrosslinkingResult> crosslinkingResult =
        findCrosslinkingSites(seq1, seq2, secondaryStructure);

    delete[] structure;
    vrna_fold_compound_free(vc);

    return HybridizationResult{mfe, crosslinkingResult};
}

std::optional<InteractionWindow> Detect::getContinuosNucleotideWindows(
    std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
    NucleotidePositionsWindow positionsPair) {
    std::pair<uint16_t, uint16_t> forwardPair =
        std::make_pair(positionsPair.first.first, positionsPair.second.first);
    std::pair<uint16_t, uint16_t> reversePair =
        std::make_pair(positionsPair.first.second, positionsPair.second.second);

    // Check that both pairs are from a continous region
    if (forwardPair.first + 1 != forwardPair.second ||
        reversePair.first - 1 != reversePair.second) {
        return std::nullopt;
    }

    size_t seq1Length = seq1.size();

    bool forwardPairSplit = forwardPair.first < seq1Length && forwardPair.second >= seq1Length;
    bool reversePairSplit = reversePair.first >= seq1Length && reversePair.second < seq1Length;

    if (forwardPairSplit || reversePairSplit) {
        return std::nullopt;
    }

    bool firstPairInSeq1 = forwardPair.second < seq1Length;
    bool secondPairInSeq1 = reversePair.first < seq1Length;

    // Both pairs are in the first sequence
    if (firstPairInSeq1 && secondPairInSeq1) {
        return InteractionWindow{
            seqan3::dna5_vector{seq1[forwardPair.first], seq1[forwardPair.second]},
            seqan3::dna5_vector{seq1[reversePair.first], seq1[reversePair.second]}, forwardPair,
            reversePair, false};
    }

    // Both pairs are in the second sequence
    if (!firstPairInSeq1 && !secondPairInSeq1) {
        return InteractionWindow{seqan3::dna5_vector{seq2[forwardPair.first - seq1Length - 1],
                                                     seq2[forwardPair.second - seq1Length - 1]},
                                 seqan3::dna5_vector{seq2[reversePair.first - seq1Length - 1],
                                                     seq2[reversePair.second - seq1Length - 1]},
                                 forwardPair, reversePair, false};
    }

    // The first pair is in the first sequence and the second pair is in the second sequence
    if (firstPairInSeq1 && !secondPairInSeq1) {
        return InteractionWindow{
            seqan3::dna5_vector{seq1[forwardPair.first], seq1[forwardPair.second]},
            seqan3::dna5_vector{seq2[reversePair.first - seq1Length - 1],
                                seq2[reversePair.second - seq1Length - 1]},
            forwardPair, reversePair, true};
    }

    // Due to sorting of base pairs first pair in second sequence and second pair in first sequence
    // is not possible
    return std::nullopt;
}

std::optional<CrosslinkingResult> Detect::findCrosslinkingSites(
    std::span<const seqan3::dna5> seq1, std::span<const seqan3::dna5> seq2,
    std::vector<seqan3::dot_bracket3> &dotbracket) {
    if (seq1.empty() || seq2.empty() || dotbracket.empty()) {
        Logger::log(LogLevel::WARNING, "Empty input sequences or dot-bracket vector!");
        return std::nullopt;
    }

    std::vector<size_t> openPos;

    std::vector<NucleotidePairPositions> interactionPositions;
    interactionPositions.reserve(dotbracket.size() / 2);

    for (size_t i = 0; i < dotbracket.size(); ++i) {
        if (dotbracket[i].is_pair_open()) {
            openPos.push_back(i);
        } else if (dotbracket[i].is_pair_close()) {
            int open = openPos.back();
            openPos.pop_back();
            interactionPositions.emplace_back(open, i);
        }
    }

    // Sort pairing sites by first position
    std::sort(interactionPositions.begin(), interactionPositions.end(),
              [](const NucleotidePairPositions &a, const NucleotidePairPositions &b) {
                  return a.first < b.first;
              });

    if (!openPos.empty()) {
        Logger::log(LogLevel::WARNING, "Unpaired bases in dot-bracket notation! (" +
                                           std::to_string(openPos.size()) + ")");
        return std::nullopt;
    }

    if (interactionPositions.size() == 0) {
        Logger::log(LogLevel::WARNING, "No interaction sites found!");
        return std::nullopt;
    }

    std::vector<NucleotidePositionsWindow> crosslinkingSites;
    crosslinkingSites.reserve(interactionPositions.size() - 1);

    size_t intraCrosslinkingScore = 0;
    size_t interCrosslinkingScore = 0;

    int preferredCrosslinkingCount = 0;
    int nonPreferredCrosslinkingCount = 0;
    int wobbleCrosslinkingCount = 0;

    // Iterate over all pairs of interaction sites until the second last pair (last pair has no
    // following pair for window)
    for (size_t i = 0; i < interactionPositions.size() - 1; ++i) {
        const NucleotidePositionsWindow window =
            std::make_pair(interactionPositions[i], interactionPositions[i + 1]);

        std::optional<InteractionWindow> interactionWindow =
            getContinuosNucleotideWindows(seq1, seq2, window);

        if (interactionWindow) {
            InteractionWindow &v_interactionWindow = interactionWindow.value();

            const NucleotideWindowPair nucleotidePairs =
                std::make_pair(v_interactionWindow.forwardWindowNucleotides,
                               v_interactionWindow.reverseWindowNucleotides);

            const auto it = crosslinkingScoringScheme.find(nucleotidePairs);
            if (it != crosslinkingScoringScheme.end()) {
                const size_t score = it->second;
                if (v_interactionWindow.isInterFragment) {
                    interCrosslinkingScore += score;
                    if (score == 3) {
                        preferredCrosslinkingCount++;
                    } else if (score == 2) {
                        nonPreferredCrosslinkingCount++;
                    } else if (score == 1) {
                        wobbleCrosslinkingCount++;
                    }
                } else {
                    intraCrosslinkingScore += score;
                }
                crosslinkingSites.emplace_back(window);
            }
        }
    }

    size_t minSeqLength = std::min(seq1.size(), seq2.size());

    if (minSeqLength == 0) {
        Logger::log(LogLevel::WARNING, "Division by zero!");
        return std::nullopt;
    }

    double normCrosslinkingScore = static_cast<double>(interCrosslinkingScore) / minSeqLength;
    return CrosslinkingResult{normCrosslinkingScore, preferredCrosslinkingCount,
                              nonPreferredCrosslinkingCount, wobbleCrosslinkingCount};
}

// calculate rows in SAMfile (exclusing header)
int Detect::countSamEntries(std::string file, std::string command) {
    std::string entries = "awk '!/^@/ {print $0}' " + file + " " + command + " | wc -l";
    const char *call = entries.c_str();
    std::array<char, 128> buffer;

    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(call, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return std::stoi(result);
}

// adds suffix to filename (optional:
std::string Detect::addSuffix(std::string _file, std::string _suffix,
                              std::vector<std::string> _keywords) {
    int keyPos, tmpKeyPos = -1;    // buffer the positions of the keywords
    int dotPos = _file.find(".");  // determine position of the dot
    keyPos = dotPos;
    if (!_keywords.empty()) {
        for (unsigned i = 0; i < _keywords.size(); ++i) {
            tmpKeyPos = _file.find(_keywords[i]);
            if (tmpKeyPos != -1) {  // key could be found
                keyPos = tmpKeyPos;
            }
        }
    }
    std::string newFile = _file.substr(0, keyPos);
    newFile += _suffix;
    newFile += _file.substr(dotPos, _file.size());

    return newFile;
}

void Detect::createDir(fs::path path) {
    if (fs::exists(path)) {
        fs::remove_all(path);
    }
    fs::create_directory(path);
}
//
std::vector<std::vector<fs::path>> Detect::splitInputFile(std::string matched, std::string splits,
                                                          int entries) {
    fs::path matchedFilePath = fs::path(matched);
    fs::path splitsFilePath = fs::path(splits);

    // create tmp folder (for inputs)
    fs::path tmpIn = matchedFilePath.parent_path() / "tmpIn";
    fs::path tmpSplit = splitsFilePath.parent_path() / "tmpSplit";
    fs::path tmpMsplit = splitsFilePath.parent_path() / "tmpMsplit";

    createDir(tmpIn);
    createDir(tmpSplit);
    createDir(tmpMsplit);

    // fs::path tmpInPath = tmpIn / fs::path(matched).filename();
    // std::string tmpInNoExt = tmpInPath.string().substr(0,tmpInPath.size()-4);

    std::string tmpInPath = (tmpIn / "tmp").string();

    std::string call = "awk '/^@/ {hdr = hdr $0 ORS; next;}";
    call +=
        "( (++numLines) % " + std::to_string(entries) + " ) == 1 { if($0 == prev ) {--numLines}";
    call += "else {close(out); out = \"" + tmpInPath +
            "_\" (++numBlocks) \".sam\"; printf \"%s\", hdr > out; numLines = 1;}}";
    call += "{print > out; prev = $0}' " + matched;

    system(call.c_str());

    // list all files in tmpIn directory
    std::vector<fs::path> subfilesInput;  //
    copy(fs::directory_iterator(tmpIn), fs::directory_iterator(), back_inserter(subfilesInput));
    sort(subfilesInput.begin(), subfilesInput.end());  // sort the content

    fs::path splitsPath;
    fs::path msplitsPath;
    std::vector<fs::path> subfilesSplits;
    std::vector<fs::path> subfilesMsplits;

    int counter = 1;
    for (auto &file : subfilesInput) {
        //     std::cout << tmpSplit / file.filename() << std::endl;
        //        fs::path splitsPath = tmpSplit / file.stem();
        //       fs::path msplitsPath = tmpMsplit / file.stem();
        subfilesSplits.push_back(tmpSplit / file.filename());
        subfilesMsplits.push_back(tmpMsplit / file.filename());
        //        subfilesSplits.push_back(fs::path(splitsPath.string().substr(0,splitsPath.size()-2)+"_splits_"+std::to_string(counter)+".sam"));
        //       subfilesMsplits.push_back(fs::path(msplitsPath.string().substr(0,msplitsPath.size()-2)+"_msplits_"+std::to_string(counter)+".sam"));
        counter++;
    }

    std::vector<std::vector<fs::path>> subfiles;
    subfiles.push_back(subfilesInput);
    subfiles.push_back(subfilesSplits);
    subfiles.push_back(subfilesMsplits);

    return subfiles;
}

// start
void Detect::start(pt::ptree sample) {
    // reset counts
    readsCount = 0;
    alignedcount = 0;
    splitscount = 0;
    msplitscount = 0;
    nsurvivedcount = 0;

    // input
    pt::ptree input = sample.get_child("input");
    std::string matched = input.get<std::string>("matched");
    fs::path matchedPath = fs::path(matched);

    int threads = params["threads"].as<int>();
    int entries = std::round(helper::countSamEntries(matchedPath) / threads) + threads;

    // output
    pt::ptree output = sample.get_child("output");
    std::string splits = output.get<std::string>("splits");
    std::string multsplits = output.get<std::string>("multsplits");

    //
    std::vector<std::vector<fs::path>> subfiles = splitInputFile(matched, splits, entries);

#pragma omp parallel for num_threads(params["threads"].as<int>())
    for (unsigned i = 0; i < subfiles[0].size(); ++i) {
        iterate(subfiles[0][i].string(), subfiles[1][i].string(), subfiles[2][i].string());
    }

    // merge results
    std::string callSplits;
    std::string callMsplits;
    if (subfiles[0].size() > 0) {
        for (unsigned i = 0; i < subfiles[0].size(); ++i) {
            if (i == 0) {
                callSplits = "awk '/^@/ {print}' ";
                callSplits += subfiles[1][i].string() + " > " + splits;
                callMsplits = "awk '/^@/ {print}' ";
                callMsplits += subfiles[2][i].string() + " > " + multsplits;
                system(callSplits.c_str());
                system(callMsplits.c_str());
            }
            callSplits = "awk '!/^@/ {print}' ";
            callSplits += subfiles[1][i].string() + " >> " + splits;
            callMsplits = "awk '!/^@/ {print}' ";
            callMsplits += subfiles[2][i].string() + " >> " + multsplits;
            system(callSplits.c_str());
            system(callMsplits.c_str());
        }
        // remove tmp folders
        fs::remove_all(subfiles[0][0].parent_path());
        fs::remove_all(subfiles[1][0].parent_path());
        fs::remove_all(subfiles[2][0].parent_path());
    }

    // generate stats
    if (params["stats"].as<std::bitset<1>>() == 1) {
        fs::path statsfile = fs::path(output.get<std::string>("stats"));
        fs::path splitsPath{splits};
        fs::path multiSplitsPath{multsplits};
        std::fstream fs;
        fs.open(statsfile.string(), std::fstream::app);
        fs << fs::path(matched).stem().string() << "\t";
        fs << helper::countUniqueSamEntries(splitsPath) << "\t";
        fs << helper::countUniqueSamEntries(multiSplitsPath) << "\t";
        fs << "\n";
        fs.close();
    }
}
